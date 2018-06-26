# Sparse Project VCF (spVCF)

Project VCF (aka multi-sample VCF) is the prevailing format for reporting small genetic variants discovered by high-throughput sequencing of a cohort of individuals. It encodes a two-dimensional matrix with variant sites down the rows and individuals across the columns, filled in with all the genotypes and associated QC measures (read depths, genotype likelihoods, etc.). Large cohorts harbor many rare variants, resulting in a sparse genotype matrix with reference-homozygous genotypes in the vast majority of entries. But the dense pVCF format encodes this very inefficiently -- a growing challenge as sequenced cohorts continue growing larger.

In the last few years, a number of new formats and advanced data structures have been developed to address this and other shortcomings of VCF. In view of VCF's considerable inertia and utility for interoperable exchange, however, here we explore a minimal evolution to sparsely encode the "lowest-hanging fruit" of repetition, while leaving other aspects of the format undisturbed. We can only emphasize that spVCF is at present a "strawman," far less sophisticated than other efforts, aimed at measuring what could be gained with a modest change.

## Build and test

```
cmake -DCMAKE_BUILD_TYPE=Release . && make
ctest -V
```

## Usage

The `spvcf` executable encodes an existing pVCF to spVCF, or conversely decodes spVCF to pVCF. The input and output streams are uncompressed text, so you may wish to arrange it in a pipe with [de]compression programs like `bgzip`.

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
