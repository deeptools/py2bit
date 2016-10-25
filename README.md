[![Build Status](https://travis-ci.org/dpryan79/py2bit.svg?branch=master)](https://travis-ci.org/dpryan79/py2bit)

# py2bit

A python extension, written in C, for quick access to [2bit](https://genome.ucsc.edu/FAQ/FAQformat.html#format7) files. The extension uses [lib2bit](https://github.com/dpryan79/lib2bit) for file access.

Table of Contents
=================

 * [Installation](#installation)
 * [Usage](#usage)
   * [Load the extension](#load-the-extension)
   * [Open a 2bit file](#open-a-2bit-file)
   * [Access the list of chromosomes and their lengths](#access-the-list-of-chromosomes-and-their-lengths)
   * [Print file information](#print-file-information)
   * [Fetch a sequence](#fetch-a-sequence)
   * [Fetch per-base statistics](#fetch-per-base-statistics)
 * [A note on coordinates](#a-note-on-coordinates)

# Installation

You can install the extension directly from github with:

    pip install git+https://github.com/dpryan79/py2bit

# Usage

Basic usage is as follows:

## Load the extension

    >>> import py2bit

## Open a 2bit file

This will work if your working directory is the py2bit source code directory.

    >>> tb = py2bit.open("test/foo.2bit")

Note that if you would like to include information about soft-masked bases, you need to manually specify that:

    >>> tb = py2bit.open("test/foo.2bit", True)

## Access the list of chromosomes and the lengths

`TwoBit` objects contain a dictionary holding the chromosome/contig lengths, which can be accessed with the `chroms()` method.

    >>> tb.chroms()
    {'chr1': 150L, 'chr2': 100L}

You can directly access a particular chromosome by specifying its name.

    >>> tb.chroms('chr1')
    150L

The lengths are stored as a "long" integer type, which is why there's an `L` suffix. If you specify a nonexistent chromosome then nothing is output.

    >>> tb.chroms("foo")
    >>>

## Print file information

The following information about and contained within a 2bit file can be accessed with the `info()` method:

 * file size, in bytes (`file size`)
 * number of chromosomes/contigs (`nChroms`)
 * total sequence length, in bases (`sequence length`)
 * total number of hard-masked (N) bases (`hard-masked length`)
 * total number of soft-masked (lower case) bases(`soft-masked length`).

Note that `soft-masked length` will only be present if `open("file.2bit", True)` is used, since handling soft-masking increases memory requirements and decreases perfomance.

    >>> tb.info()
    {'file size': 161, 'nChroms': 2, 'sequence length': 250, 'hard-masked length': 150, 'soft-masked length': 8}

## Fetch a sequence

The sequence of a full or partial chromosome/contig can be fetched with the `sequence()` method.

    >>> tb.sequence("chr1")
    'NNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNACGTACGTACGTagctagctGATCGATCGTAGCTAGCTAGCTAGCTGATCNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNN'

By default, the whole chromosome/contig is returned. A specific range can also be requested.

    >>> tb.sequence("chr1", 24, 74)
    NNNNNNNNNNNNNNNNNNNNNNNNNNACGTACGTACGTagctagctGATC

The first number is the (0-based) position on the chromosome/contig where the sequence should begin. The second number is the (1-based) position on the chromosome where the sequence should end.

If it was requested during file opening that soft-masking information be stored, then lower case bases may be present. If a nonexistent chromosome/contig is specified then a runtime error occurs.

## Fetch per-base statistics

It's often required to compute the percentage of 1 or more bases in a chromosome. This can be done with the `bases()` method.

    >>> tb.bases("chr1")
    {'A': 0.08, 'C': 0.08, 'T': 0.08666666666666667, 'G': 0.08666666666666667}

This returns a dictionary with bases as keys and the fraction of the sequence composed of them as values. Note that this will not sum to 1 if there are any hard-masked bases (the chromosome is 2/3 `N` in this case). One can also request this information over a particular region.

    >>> tb.bases("chr1", 24, 74)
    {'A': 0.12, 'C': 0.12, 'T': 0.12, 'G': 0.12}

The start and end position are as with the `sequence()` method described above.

If integer counts are preferred, then they can instead be returned.

    >>> tb.bases("chr1", 24, 74, False)
    {'A': 6, 'C': 6, 'T': 6, 'G': 6}

## Close a file

A `TwoBit` object can be closed with the `close()` method.

    >>> tb.close()

# A note on coordinates

0-based half-open coordinates are used by this python module. So to access the value for the first base on `chr1`, one would specify the starting position as `0` and the end position as `1`. Similarly, bases 100 to 115 would have a start of `99` and an end of `115`. This is simply for the sake of consistency with most other bioinformatics packages.
