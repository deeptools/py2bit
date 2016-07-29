#include <inttypes.h>
#include <stdio.h>

#ifdef __cplusplus
extern "C" {
#endif

typedef struct {
    uint32_t magic; /**<Holds the magic number, should be 0x1A412743 */
    uint32_t version; /**<File version, should be 0 */
    uint32_t nChroms; /**<Number of chromosomes/contigs */
} TwoBitHeader;

typedef struct {
    char **chrom; /**<A list of null terminated chromosomes */
    uint32_t *offset; /**<The file offset for the beginning of each chromosome */
} TwoBitCL;

typedef struct {
    uint32_t *size; /**<The size of a given chromosome/contig */
    uint32_t *nBlockCount; /**<The number of blocks of Ns in a given chromosome/contig */
    uint32_t **nBlockStart; /**<For each chromosome/contig, the list (size nBlockCount) of start positions of the block of Ns */
    uint32_t **nBlockSizes; /**<The size of each block specified above */
    uint32_t *maskBlockCount; /**<The number of blocks of masked sequence in a given chromosome/contig */
    uint32_t **maskBlockStart; /**<For each chromosome/contig, the list (size maskBlockCount) of start positions of the masked sequence blocks */
    uint32_t **maskBlockSizes; /**<The size of each block specified above */
    uint64_t *offset; /**<The offset to the packed 2-bit sequence */
} TwoBitMaskedIdx;

typedef struct {
    FILE *fp;    //The file pointer for the opened file
    uint64_t sz; //File size in bytes (needed for munmap)
    uint64_t offset; //If the file is memory mapped, then this is the current file offset (otherwise ignored)
    void *data;  //The memory mapped file, if it exists.
    TwoBitHeader *hdr; //File header
    TwoBitCL *cl; //Chromosome list with sizes
    TwoBitMaskedIdx *idx; //Index of masked blocks
} TwoBit;

//Open/close functions
TwoBit* twobitOpen(char *fname, int storeMasked);
void twobitClose(TwoBit *tb);

//Return the length of a given chromosome/contig or 0 if not present
uint32_t twobitChromLen(TwoBit *tb, char *chrom);

//Return the sequence of the given range (or a whole chromosome if start=end=0
char *twobitSequence(TwoBit *tb, char *chrom, uint32_t start, uint32_t end);

//We can try to assume that N is always a T and see if that's correct...
double *twobitFrequency(TwoBit *tb, char *chrom, uint32_t start, uint32_t end);

#ifdef __cplusplus
}
#endif

