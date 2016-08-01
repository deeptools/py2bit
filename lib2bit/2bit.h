#include <inttypes.h>
#include <stdio.h>

/*! \mainpage libBigWig
 *
 * \section Introduction
 *
 * lib2bit is a C-based library for accessing [2bit files](https://genome.ucsc.edu/FAQ/FAQformat.html#format7). At the moment, only reading 2bit files is supported (there are no plans to change this, though if someone wants to submit a pull request...). Though it's unlikely to matter, 
 *
 * The motivation for this project is due to needing fast access to 2bit files in [deepTools](https://github.com/fidelram/deepTools). Originally, we were using bx-python for this, which had the benefit of being easy to install and pretty quick. However, that wasn't compatible with python3, so we switched to [twobitreader](https://github.com/benjschiller/twobitreader). While doing everything we needed and working under both python2 and python3, it turns out that it has terrible performance (up to 1000x slow down in `computeGCBias`). Since we'd like to have our cake and eat it too, I began wrote a C library for convenient 2bit access and then [a python wrapper](https://github.com/dpryan79/py2bit) around it to work in python2 and 3.
 *
 * \section Installation
 *
 * 2bit files are very simple and there are no dependencies. Simply typing `make` should suffice for compilation. To install into a specific path (the default is `/usr/local`):
 * 
 *     make install prefix=/some/where/else
 *
 * `lib2bit.so` and `lib2bit.a` will then be in `/some/where/else/lib` and `2bit.h` in `/some/where/else/include`.
 *
 * \section Example
 *
 * See the `test/` directory for an example of using the library.
 */

/*! \file 2bit.h
 *
 * These are all functions and structures exported in lib2bit. There are a few things that could be more efficiently implemented, but at the moment theverything is "fast enough".
 */

#ifdef __cplusplus
extern "C" {
#endif

/*!
 * @brief This structure holds the fixed-sized file header (16 bytes, of which 4 are blank). The version should always be 0. In theory, the endianness of the magic number can change (indicating that everything in the file should be swapped). As I've never actually seen this occur in the wild I've not bothered implementing it, though it'd be simple enough to do so.
 */
typedef struct {
    uint32_t magic; /**<Holds the magic number, should be 0x1A412743 */
    uint32_t version; /**<File version, should be 0 */
    uint32_t nChroms; /**<Number of chromosomes/contigs */
} TwoBitHeader;

/*!
 * @brief This structure holds the chromosome names and the offset to the on-disk beginning of their sequences
 */
typedef struct {
    char **chrom; /**<A list of null terminated chromosomes */
    uint32_t *offset; /**<The file offset for the beginning of each chromosome */
} TwoBitCL;

/*!
 * @brief This structure holds the number, location and size of the hard (N) and soft (lower case) masked blocks.
 *
 * Note that this isn't a great data structure for random access, particularly for the soft-masked blocks. In practice, soft-masking is typically ignored and file access is less random and more blocky. Nonetheless, if performance is not acceptable then this is the structure to change.
 */
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

/*!
 * @brief This is the main structure for holding a 2bit file
 *
 * Note that currently the 2bit file is mmap()ed prior to reading and that this isn't optional.
 */
typedef struct {
    FILE *fp;    /**<The file pointer for the opened file */
    uint64_t sz; /**<File size in bytes (needed for munmap) */
    uint64_t offset; /**<If the file is memory mapped, then this is the current file offset (otherwise ignored) */
    void *data;  /**<The memory mapped file, if it exists. */
    TwoBitHeader *hdr; /**<File header */
    TwoBitCL *cl; /**<Chromosome list with sizes */
    TwoBitMaskedIdx *idx; /**<Index of masked blocks */
} TwoBit;

/*!
 * @brief Opens a local 2bit file
 *
 * @param fname The name of the 2bit file.
 * @param storeMasked Whether soft-masking information should be stored. If this is 1 then soft-masking information will be stored and the `twobitSequence()` function will return lower case letters in soft-masked regions. Note that this has a considerable performance and memory impact.
 * @return A pointer to a TwoBit object.
 * @note The file is memory mapped.
 */
TwoBit* twobitOpen(char *fname, int storeMasked);

/*!
 * @brief Closes a 2bit file and free memory.
 */
void twobitClose(TwoBit *tb);

/*!
 * @brief Returns the length of a given chromosome.
 * 
 * @param tb A pointer to a TwoBit object.
 * @param chrom The chromosome name.
 * @return The chromosome length as a uint32_t. Note that if the chromosome/contig isn't present in the file that 0 is returned.
 */
uint32_t twobitChromLen(TwoBit *tb, char *chrom);

/*!
 * @brief Returns the sequence of a chromosome/contig or range of it.
 *
 * @param tb A pointer to a TwoBit object.
 * @param chrom The chromosome name.
 * @param start The starting position in 0-based coordinates.
 * @param end The end position in 1-based coordinates.
 * @return The sequence or NULL on error. If both start and end are 0 then the sequence for the entire chromosome/contig is returned.
 * @note The result MUST be `free()`d. Care is taken to return reasonable sequences when illegal regions are requested. If the end value is beyond the possible end of the chromosome then it is modified according.
 */
char *twobitSequence(TwoBit *tb, char *chrom, uint32_t start, uint32_t end);

/*!
 * @brief Return the number/fraction of A, C, T, and G in a chromosome/region
 * 
 * @param tb A pointer to a TwoBit object.
 * @param chrom The chromosome name.
 * @param start The starting position in 0-based coordinates.
 * @param end The end position in 1-based coordinates.
 * @param fraction Whether to return the values as fractions (1) or integers (0).
 * @return If fraction is not 0, then 4 `double`s with the fraction of bases as A, C, T and G, respectively. If fraction is 1, integer counts are returned as 4 `uint32_t`s in the aforementioned order.
 * @note On error NULL is returned. The result MUST be `free()`d.
 */

void *twobitBases(TwoBit *tb, char *chrom, uint32_t start, uint32_t end, int fraction);

#ifdef __cplusplus
}
#endif
