#include <sys/types.h>
#include <sys/stat.h>
#include <sys/mman.h>
#include <unistd.h>
#include <inttypes.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <ctype.h>
#include "2bit.h"

uint64_t twobitTell(TwoBit *tb);

/*
    Read nmemb elements, each of size sz from the current file offset 
    into data. Return the number of elements read. On error, the return
    value is either 0 or less than nmemb
*/
size_t twobitRead(void *data, size_t sz, size_t nmemb, TwoBit *tb) {
    if(tb->data) {
        if(memcpy(data, tb->data + tb->offset, nmemb * sz) == NULL) return 0;
        tb->offset += nmemb * sz;
        return nmemb;
    } else {
        return fread(data, sz, nmemb, tb->fp);
    }
}

/*
    Seek to a specific position, which is essentially trivial for memmaped stuff

    Returns: 0 on success, -1 on error
*/
int twobitSeek(TwoBit *tb, uint64_t offset) {
    if(offset >= tb->sz) return -1;
    if(tb->data) {
        tb->offset = offset;
        return 0;
    } else {
        return fseek(tb->fp, (long) offset, SEEK_SET);
    }
}

/*
    Like ftell, but generalized to handle memmaped files

    Returns the offset
*/
uint64_t twobitTell(TwoBit *tb) {
    if(tb->data) return tb->offset;
    return (uint64_t) ftell(tb->fp);
}

/*
    Given a byte containing 4 bases, return the character representation of the offset'th base
*/
char byte2base(uint8_t byte, int offset) {
    int rev = 3 - offset;
    uint8_t mask = 3 << (2 * rev);
    int foo = (mask & byte) >> (2 * rev);
    char bases[4] = "TCAG";
    return bases[foo];
}

/*
    Replace Ts (or whatever else is being used) with N as appropriate
*/
void NMask(char *seq, TwoBit *tb, uint32_t tid, uint32_t start, uint32_t end) {
    uint32_t i, width, pos = 0;
    uint32_t blockStart, blockEnd;

    for(i=0; i<tb->idx->nBlockCount[tid]; i++) {
        blockStart = tb->idx->nBlockStart[tid][i];
        blockEnd = blockStart + tb->idx->nBlockSizes[tid][i];
        if(blockEnd <= start) continue;
        if(blockStart >= end) break;
        if(blockStart < start) {
            blockEnd = (blockEnd < end) ? blockEnd : end;
            pos = 0;
            width = blockEnd - start;
        } else {
            blockEnd = (blockEnd < end) ? blockEnd : end;
            pos = blockStart - start;
            width = blockEnd - blockStart;
        }
        width += pos;
        for(; pos < width; pos++) seq[pos] = 'N';
    }
}

/*
    Replace uppercase with lower-case letters, if required
*/
void softMask(char *seq, TwoBit *tb, uint32_t tid, uint32_t start, uint32_t end) {
    uint32_t i, width, pos = 0;
    uint32_t blockStart, blockEnd;

    if(!tb->idx->maskBlockStart) return;

    for(i=0; i<tb->idx->maskBlockCount[tid]; i++) {
        blockStart = tb->idx->maskBlockStart[tid][i];
        blockEnd = blockStart + tb->idx->maskBlockSizes[tid][i];
        if(blockEnd <= start) continue;
        if(blockStart >= end) break;
        if(blockStart < start) {
            blockEnd = (blockEnd < end) ? blockEnd : end;
            pos = 0;
            width = blockEnd - start;
        } else {
            blockEnd = (blockEnd < end) ? blockEnd : end;
            pos = blockStart - start;
            width = blockEnd - blockStart;
        }
        width += pos;
        for(; pos < width; pos++) {
            if(seq[pos] != 'N') seq[pos] = tolower(seq[pos]);
        }
    }
}

/*
    This is the worker function for twobitSequence, which mostly does error checking
*/
char *constructSequence(TwoBit *tb, uint32_t tid, uint32_t start, uint32_t end) {
    uint32_t sz = end - start + 1, pos = 0;
    uint32_t blockStart, offset;
    char *seq = malloc(sz * sizeof(char)), byte;
    if(!seq) return NULL;

    //There are 4 bases/byte
    blockStart = start/4;
    offset = start % 4;

    if(twobitSeek(tb, tb->idx->offset[tid] + blockStart) != 0) goto error;
    while(pos < sz - 1) {
        if(twobitRead(&byte, 1, 1, tb) != 1) goto error;

        for(; offset<4; offset++) {
            seq[pos++] = byte2base(byte, offset);
            if(pos >= sz - 1) break;
        }
        offset = 0;
    }

    //Null terminate the output
    seq[sz - 1] = '\0';

    //N-mask everything
    NMask(seq, tb, tid, start, end);

    //Soft-mask if requested
    softMask(seq, tb, tid, start, end);

    return seq;

error:
    if(seq) free(seq);
    return NULL;
}

/*
    Given a chromosome, name, and optional range, return the corresponding sequence.

    The start and end or 0-based half-open, so end-start is the number of bases.
    If both start and end are 0, then the whole chromosome is used.

    On error (e.g., a missing chromosome), NULL is returned.
*/
char *twobitSequence(TwoBit *tb, char *chrom, uint32_t start, uint32_t end) {
    uint32_t i, tid=0;

    //Get the chromosome ID
    for(i=0; i<tb->hdr->nChroms; i++) {
        if(strcmp(tb->cl->chrom[i], chrom) == 0) {
            tid = i;
            break;
        }
    }
    if(tid == 0 && strcmp(tb->cl->chrom[i], chrom) != 0) return NULL;

    //Get the start/end if not specified
    if(start == end && end == 0) {
        end = tb->idx->size[tid];
    }

    //Sanity check the bounds
    if(end > tb->idx->size[tid]) return NULL;
    if(start >= end) return NULL;

    return constructSequence(tb, tid, start, end);
}

/*
    Given a tid and a position, set the various mask variables to an appropriate block of Ns.

     * If maskIdx is not -1, these are set to the first overlapping block (or maskIdx is set to the number of N blocks).
     * If maskIdx is not -1 then it's incremented and maskStart/maskEnd set appropriately.

    If the returned interval doesn't overlap the start/end range, then both values will be -1.
*/
void getMask(TwoBit *tb, uint32_t tid, uint32_t start, uint32_t end, uint32_t *maskIdx, uint32_t *maskStart, uint32_t *maskEnd) {
    if(*maskIdx == (uint32_t) -1) {
        for((*maskIdx)=0; (*maskIdx)<tb->idx->nBlockCount[tid]; (*maskIdx)++) {
           *maskStart = tb->idx->nBlockStart[tid][*maskIdx];
           *maskEnd = (*maskStart) + tb->idx->nBlockSizes[tid][*maskIdx];
           if(*maskEnd < start) continue;
           if(*maskEnd >= start) break;
        }
    } else if(*maskIdx >= tb->idx->nBlockCount[tid]) {
        *maskStart = (uint32_t) -1;
        *maskEnd = (uint32_t) -1;
    } else {
        *maskIdx += 1;
        if(*maskIdx >= tb->idx->nBlockCount[tid]) {
            *maskStart = (uint32_t) -1;
            *maskEnd = (uint32_t) -1;
        } else {
            *maskStart = tb->idx->nBlockStart[tid][*maskIdx];
            *maskEnd = (*maskStart) + tb->idx->nBlockSizes[tid][*maskIdx];
        }
    }

    //maskStart = maskEnd = -1 if no overlap
    if(*maskIdx >= tb->idx->nBlockCount[tid] || *maskStart >= end) {
        *maskStart = (uint32_t) -1;
        *maskEnd = (uint32_t) -1;
    }
}

/*
    Determine whether "start + pos" overlaps an N block. If so, update pos, return 1, and update the mask stuff (as well as update blockStart and offset). If there's no overlap, return 0.

    Returning a value of 1 indicates that pos has been updated (i.e., there was an overlap). A value of 0 indicates otherwise.

    A return value of -1 indicates an error (in twobitSeek).
*/
int cycle(TwoBit *tb, uint32_t tid, uint32_t start, uint32_t end, uint32_t *pos, uint32_t *maskIdx, uint32_t *maskStart, uint32_t *maskEnd, uint32_t *blockStart, uint32_t *offset) {
    if(*maskIdx < tb->idx->nBlockCount[tid]) {
        if(end < *maskStart) {
            getMask(tb, tid, start, end, maskIdx, maskStart, maskEnd);
            if(*maskIdx < tb->idx->nBlockCount[tid]) return 0;
        }
        if((start + *pos >= *maskStart) && (start + *pos < *maskEnd)) {
            *pos += (*maskEnd) - (start + *pos);
            *blockStart = (start + *pos)/4;
            *offset = (start + *pos) % 4;
            getMask(tb, tid, start, end, maskIdx, maskStart, maskEnd);
            if(twobitSeek(tb, tb->idx->offset[tid] + *blockStart) != 0) return -1;
            return 1;
        }
    }
    return 0;
}

/*
    Increment the counter for the appropriate base
*/
void increment(char base, uint32_t *A, uint32_t *C, uint32_t *T, uint32_t *G) {
    switch(base) {
    case 'T':
        *T += 1;
        break;
    case 'C':
        *C += 1;
        break;
    case 'A':
        *A += 1;
        break;
    case 'G':
        *G += 1;
        break;
    }
}

void *twobitBasesWorker(TwoBit *tb, uint32_t tid, uint32_t start, uint32_t end, int fraction) {
    void *out;
    uint32_t sz = end - start, pos = 0;
    uint32_t A = 0, C = 0, T = 0, G = 0, len = end - start;
    uint32_t maskIdx = -1, maskStart = -1, maskEnd = -1;
    uint32_t blockStart, offset;
    uint8_t byte, base;
    int rv = 0;

    if(fraction) {
        out = malloc(4 * sizeof(double));
    } else {
        out = malloc(4 * sizeof(uint32_t));
    }
    if(!out) return NULL;

    getMask(tb, tid, start, end, &maskIdx, &maskStart, &maskEnd);

    //There are 4 bases/byte
    blockStart = start/4;
    offset = start % 4;

    if(twobitSeek(tb, tb->idx->offset[tid] + blockStart) != 0) goto error;
    while(pos < sz) {
        //Iterate over the 4 (or less) bases in each byte
        if(twobitRead(&byte, 1, 1, tb) != 1) goto error;
        for(; offset<4; offset++) {
            rv = cycle(tb, tid, start, end, &pos, &maskIdx, &maskStart, &maskEnd, &blockStart, &offset);
            if(rv == -1) goto error;
            if(rv == 1) break;
            base = byte2base(byte, offset);
            increment(base, &A, &C, &T, &G);
            if(++pos >= sz) break;
        }
        if(rv != 1) offset = 0;
    }

    if(fraction) {
        ((double*) out)[0] = ((double) A)/((double) len);
        ((double*) out)[1] = ((double) C)/((double) len);
        ((double*) out)[2] = ((double) T)/((double) len);
        ((double*) out)[3] = ((double) G)/((double) len);
    } else {
        ((uint32_t*) out)[0] = A;
        ((uint32_t*) out)[1] = C;
        ((uint32_t*) out)[2] = T;
        ((uint32_t*) out)[3] = G;
    }

    return out;

error:
    free(out);
    return NULL;
}

void *twobitBases(TwoBit *tb, char *chrom, uint32_t start, uint32_t end, int fraction) {
    uint32_t tid = 0, i;

    //Get the chromosome ID
    for(i=0; i<tb->hdr->nChroms; i++) {
        if(strcmp(tb->cl->chrom[i], chrom) == 0) {
            tid = i;
            break;
        }
    }

    if(tid == 0 && strcmp(tb->cl->chrom[i], chrom) != 0) return NULL;

    //Get the start/end if not specified
    if(start == end && end == 0) {
        end = tb->idx->size[tid];
    }

    //Sanity check the bounds
    if(end > tb->idx->size[tid]) return NULL;
    if(start >= end) return NULL;

    return twobitBasesWorker(tb, tid, start, end, fraction);
}

/*
    Given a chromosome, chrom, return it's length. 0 is used if the chromosome isn't present.
*/
uint32_t twobitChromLen(TwoBit *tb, char *chrom) {
    uint32_t i;
    for(i=0; i<tb->hdr->nChroms; i++) {
        if(strcmp(tb->cl->chrom[i], chrom) == 0) return tb->idx->size[i];
    }
    return 0;
}

/*
    Fill in tb->idx.

    Note that the masked stuff will only be stored if storeMasked == 1, since it uses gobs of memory otherwise.
    On error, tb->idx is left as NULL.
*/
void twobitIndexRead(TwoBit *tb, int storeMasked) {
    uint32_t i, data[2];
    TwoBitMaskedIdx *idx = calloc(1, sizeof(TwoBitMaskedIdx));

    //Allocation and error checking
    if(!idx) return;
    idx->size = malloc(tb->hdr->nChroms * sizeof(uint32_t));
    idx->nBlockCount = calloc(tb->hdr->nChroms, sizeof(uint32_t));
    idx->nBlockStart = calloc(tb->hdr->nChroms, sizeof(uint32_t*));
    idx->nBlockSizes = calloc(tb->hdr->nChroms, sizeof(uint32_t*));
    if(!idx->size) goto error;
    if(!idx->nBlockCount) goto error;
    if(!idx->nBlockStart) goto error;
    if(!idx->nBlockSizes) goto error;
    idx->maskBlockCount = calloc(tb->hdr->nChroms, sizeof(uint32_t));
    if(!idx->maskBlockCount) goto error;
    if(storeMasked) {
        idx->maskBlockStart = calloc(tb->hdr->nChroms, sizeof(uint32_t*));
        idx->maskBlockSizes = calloc(tb->hdr->nChroms, sizeof(uint32_t*));
        if(!idx->maskBlockStart) goto error;
        if(!idx->maskBlockSizes) goto error;
    }
    idx->offset = malloc(tb->hdr->nChroms * sizeof(uint64_t));
    if(!idx->offset) goto error;

    //Read in each chromosome/contig
    for(i=0; i<tb->hdr->nChroms; i++) {
        if(twobitSeek(tb, tb->cl->offset[i]) != 0) goto error;
        if(twobitRead(data, sizeof(uint32_t), 2, tb) != 2) goto error;
        idx->size[i] = data[0];
        idx->nBlockCount[i] = data[1];

        //Allocate the nBlock starts/sizes and fill them in
        idx->nBlockStart[i] = malloc(idx->nBlockCount[i] * sizeof(uint32_t));
        idx->nBlockSizes[i] = malloc(idx->nBlockCount[i] * sizeof(uint32_t));
        if(!idx->nBlockStart[i]) goto error;
        if(!idx->nBlockSizes[i]) goto error;
        if(twobitRead(idx->nBlockStart[i], sizeof(uint32_t), idx->nBlockCount[i], tb) != idx->nBlockCount[i]) goto error;
        if(twobitRead(idx->nBlockSizes[i], sizeof(uint32_t), idx->nBlockCount[i], tb) != idx->nBlockCount[i]) goto error;

        //Get the masked block information
        if(twobitRead(idx->maskBlockCount + i, sizeof(uint32_t), 1, tb) != 1) goto error;

        //Allocate the maskBlock starts/sizes and fill them in
        if(storeMasked) {
            idx->maskBlockStart[i] = malloc(idx->maskBlockCount[i] * sizeof(uint32_t));
            idx->maskBlockSizes[i] = malloc(idx->maskBlockCount[i] * sizeof(uint32_t));
            if(!idx->maskBlockStart[i]) goto error;
            if(!idx->maskBlockSizes[i]) goto error;
            if(twobitRead(idx->maskBlockStart[i], sizeof(uint32_t), idx->maskBlockCount[i], tb) != idx->maskBlockCount[i]) goto error;
            if(twobitRead(idx->maskBlockSizes[i], sizeof(uint32_t), idx->maskBlockCount[i], tb) != idx->maskBlockCount[i]) goto error;
        } else {
            if(twobitSeek(tb, twobitTell(tb) + 8 * idx->maskBlockCount[i]) != 0) goto error;
        }

        //Reserved
        if(twobitRead(data, sizeof(uint32_t), 1, tb) != 1) goto error;

        idx->offset[i] = twobitTell(tb);
    }

    tb->idx = idx;
    return;

error:
    if(idx) {
        if(idx->size) free(idx->size);

        if(idx->nBlockCount) free(idx->nBlockCount);
        if(idx->nBlockStart) {
            for(i=0; i<tb->hdr->nChroms; i++) {
                if(idx->nBlockStart[i]) free(idx->nBlockStart[i]);
            }
            free(idx->nBlockStart[i]);
        }
        if(idx->nBlockSizes) {
            for(i=0; i<tb->hdr->nChroms; i++) {
                if(idx->nBlockSizes[i]) free(idx->nBlockSizes[i]);
            }
            free(idx->nBlockSizes[i]);
        }

        if(idx->maskBlockCount) free(idx->maskBlockCount);
        if(idx->maskBlockStart) {
            for(i=0; i<tb->hdr->nChroms; i++) {
                if(idx->maskBlockStart[i]) free(idx->maskBlockStart[i]);
            }
            free(idx->maskBlockStart[i]);
        }
        if(idx->maskBlockSizes) {
            for(i=0; i<tb->hdr->nChroms; i++) {
                if(idx->maskBlockSizes[i]) free(idx->maskBlockSizes[i]);
            }
            free(idx->maskBlockSizes[i]);
        }

        if(idx->offset) free(idx->offset);

        free(idx);
    }
}

void twobitIndexDestroy(TwoBit *tb) {
    uint32_t i;

    if(tb->idx) {
        if(tb->idx->size) free(tb->idx->size);

        if(tb->idx->nBlockCount) free(tb->idx->nBlockCount);
        if(tb->idx->nBlockStart) {
            for(i=0; i<tb->hdr->nChroms; i++) {
                if(tb->idx->nBlockStart[i]) free(tb->idx->nBlockStart[i]);
            }
            free(tb->idx->nBlockStart);
        }
        if(tb->idx->nBlockSizes) {
            for(i=0; i<tb->hdr->nChroms; i++) {
                if(tb->idx->nBlockSizes[i]) free(tb->idx->nBlockSizes[i]);
            }
            free(tb->idx->nBlockSizes);
        }

        if(tb->idx->maskBlockCount) free(tb->idx->maskBlockCount);
        if(tb->idx->maskBlockStart) {
            for(i=0; i<tb->hdr->nChroms; i++) {
                if(tb->idx->maskBlockStart[i]) free(tb->idx->maskBlockStart[i]);
            }
            free(tb->idx->maskBlockStart);
        }
        if(tb->idx->maskBlockSizes) {
            for(i=0; i<tb->hdr->nChroms; i++) {
                if(tb->idx->maskBlockSizes[i]) free(tb->idx->maskBlockSizes[i]);
            }
            free(tb->idx->maskBlockSizes);
        }

        if(tb->idx->offset) free(tb->idx->offset);

        free(tb->idx);
    }
}

void twobitChromListRead(TwoBit *tb) {
    uint32_t i;
    uint8_t byte;
    char *str = NULL;
    TwoBitCL *cl = calloc(1, sizeof(TwoBitCL));

    //Allocate cl and do error checking
    if(!cl) goto error;
    cl->chrom = calloc(tb->hdr->nChroms, sizeof(char*));
    cl->offset = malloc(sizeof(uint32_t) * tb->hdr->nChroms);
    if(!cl->chrom) goto error;
    if(!cl->offset) goto error;

    for(i=0; i<tb->hdr->nChroms; i++) {
        //Get the string size (not null terminated!)
        if(twobitRead(&byte, 1, 1, tb) != 1) goto error;

        //Read in the string
        str = calloc(1 + byte, sizeof(char));
        if(!str) goto error;
        if(twobitRead(str, 1, byte, tb) != byte) goto error;
        cl->chrom[i] = str;
        str = NULL;

        //Read in the size
        if(twobitRead(cl->offset + i, sizeof(uint32_t), 1, tb) != 1) goto error;
    }

    tb->cl = cl;
    return;

error:
    if(str) free(str);
    if(cl) {
        if(cl->offset) free(cl->offset);
        if(cl->chrom) {
            for(i=0; i<tb->hdr->nChroms; i++) {
                if(cl->chrom[i]) free(cl->chrom[i]);
            }
            free(cl->chrom);
        }
        free(cl);
    }
}

void twobitChromListDestroy(TwoBit *tb) {
    uint32_t i;

    if(tb->cl) {
        if(tb->cl->offset) free(tb->cl->offset);
        if(tb->cl->chrom) {
            for(i=0; i<tb->hdr->nChroms; i++) {
                if(tb->cl->chrom[i]) free(tb->cl->chrom[i]);
            }
            free(tb->cl->chrom);
        }
        free(tb->cl);
    }
}

void twobitHdrRead(TwoBit *tb) {
    //Read the first 16 bytes
    uint32_t data[4];
    TwoBitHeader *hdr = calloc(1, sizeof(TwoBitHeader));

    if(!hdr) return;

    if(twobitRead(data, 4, 4, tb) != 4) goto error;

    //Magic
    hdr->magic = data[0];
    if(hdr->magic != 0x1A412743) {
        fprintf(stderr, "[twobitHdrRead] Received an invalid file magic number (0x%"PRIx32")!\n", hdr->magic);
        goto error;
    }

    //Version
    hdr->version = data[1];
    if(hdr->version != 0) {
        fprintf(stderr, "[twobitHdrRead] The file version is %"PRIu32" while only version 0 is defined!\n", hdr->version);
        goto error;
    }

    //Sequence Count
    hdr->nChroms = data[2];
    if(hdr->nChroms == 0) {
        fprintf(stderr, "[twobitHdrRead] There are apparently no chromosomes/contigs in this file!\n");
        goto error;
    }

    tb->hdr = hdr;
    return;

error:
    if(hdr) free(hdr);
}

void twobitHdrDestroy(TwoBit *tb) {
    if(tb->hdr) free(tb->hdr);
}

void twobitClose(TwoBit *tb) {
    if(tb) {
        if(tb->fp) fclose(tb->fp);
        if(tb->data) munmap(tb->data, tb->sz);
        twobitChromListDestroy(tb);
        twobitIndexDestroy(tb);
        //N.B., this needs to be called last
        twobitHdrDestroy(tb);
        free(tb);
    }
}

TwoBit* twobitOpen(char *fname, int storeMasked) {
    int fd;
    struct stat fs;
    TwoBit *tb = calloc(1, sizeof(TwoBit));
    if(!tb) return NULL;

    tb->fp = fopen(fname, "rb");
    if(!tb->fp) goto error;

    //Try to memory map the whole thing, since these aren't terribly large
    //Since we might be multithreading this in python, use shared memory
    fd = fileno(tb->fp);
    if(fstat(fd, &fs) == 0) {
        tb->sz = (uint64_t) fs.st_size;
        tb->data = mmap(NULL, fs.st_size, PROT_READ, MAP_SHARED, fd, 0);
        if(tb->data) {
            if(madvise(tb->data, fs.st_size, MADV_RANDOM) != 0) {
                munmap(tb->data, fs.st_size);
                tb->data = NULL;
            }
        }
    }

    //Attempt to read in the fixed header
    twobitHdrRead(tb);
    if(!tb->hdr) goto error;

    //Read in the chromosome list
    twobitChromListRead(tb);
    if(!tb->cl) goto error;

    //Read in the mask index
    twobitIndexRead(tb, storeMasked);
    if(!tb->idx) goto error;

    return tb;

error:
    twobitClose(tb);
    return NULL;
}
