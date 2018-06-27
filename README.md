# Sparse Project VCF (spVCF)

Project VCF (aka multi-sample VCF) is the prevailing file format for reporting small genetic variants discovered by high-throughput cohort sequencing. It encodes a two-dimensional matrix with variant sites down the rows and study participants across the columns, filled in with all the genotypes and associated QC measures (read depths, genotype likelihoods, etc.). Large cohorts harbor many rare variants, resulting in a sparse genotype matrix with reference-homozygous genotypes in the vast majority of entries. But the dense pVCF format encodes this very inefficiently -- a growing challenge as sequenced cohorts continue growing larger.

In the last few years, a number of new formats and advanced data structures have been developed to address this and other shortcomings of VCF. In view of VCF's considerable inertia and utility, however, here we explore a minimal evolution to sparsely encode the "lowest-hanging fruit" of repetition, while leaving other aspects of the format unchanged. We can only emphasize that spVCF is far less sophisticated than other efforts, a "strawman" aimed at measuring what could be gained with a modest change.

## Build and test

[![Build Status](https://travis-ci.org/mlin/spVCF.svg?branch=master)](https://travis-ci.org/mlin/spVCF)

```
cmake -DCMAKE_BUILD_TYPE=Release . && make
ctest -V
```

## Usage

The `spvcf` executable encodes an existing pVCF to spVCF, or conversely decodes spVCF to pVCF. The input and output streams are uncompressed VCF text, so you may wish to arrange a pipe with [de]compression programs like `bgzip`.

```
$ ./spvcf -h
spvcf [options] [in.vcf|-]
Read from standard input if input filename is empty or -
Options:
  -o,--output out.vcf    Write to out.vcf instead of standard output
  -d,--decode            Decode from the sparse format instead of encoding to
  -q,--quiet             Suppress statistics printed to standard error
  -h,--help              Show this usage message
```

Examples:

```
$ ./spvcf my.vcf > my.spvcf
$ bgzip -dc my.vcf.gz | ./spvcf | bgzip -c > my.spvcf.gz
$ bgzip -dc my.spvcf.gz | ./spvcf -d > my.decoded.vcf
```
