# Sparse Project VCF (spVCF)

**Maintainer: Mike Lin [@DNAmlin](https://twitter.com/DNAmlin)**

Project VCF (pVCF; aka multi-sample VCF) is the prevailing file format for small genetic variants discovered by cohort sequencing. It encodes a two-dimensional matrix with variant sites down the rows and study participants across the columns, filled in with all the genotypes and associated QC measures (read depths, genotype likelihoods, etc.). Large cohorts harbor many rare variants, implying a sparse genotype matrix composed largely of reference-homozygous or non-called cells. But the dense pVCF format encodes this inefficiently, growing super-linearly with the cohort size.

spVCF is an evolution of VCF that keeps most aspects of its tab-delimited text format, but presents the genotype matrix sparsely, by run-length encoding repetitive information about reference coverage. This is less sophisticated than some other efforts to address VCF's density and other shortcomings, but perhaps more palatable to existing VCF consumers by virtue of simplicity.

Further resources:

* [bioRxiv preprint](https://www.biorxiv.org/content/10.1101/611954v1) for a short manuscript on the approach & tool
* [doc/SPEC.md](https://github.com/mlin/spVCF/blob/master/doc/SPEC.md) has complete details and a worked example
* [doc/compression_results.md](https://github.com/mlin/spVCF/blob/master/doc/compression_results.md) tests spVCF with *N*=50K exomes, observing up to 15X size reduction for bgzip-compressed pVCF, and scaling much more gently with *N*.
* [slide deck](https://docs.google.com/presentation/d/13lzEkdWAVwcsKofhsiYEdl92xMQgx5_dSOSIyZDggfM/edit?usp=sharing) presented at the GA4GH & MPEG-G Genome Compression Workshop, October 2018.

## `spvcf` utility

[![Build Status](https://travis-ci.org/mlin/spVCF.svg?branch=master)](https://travis-ci.org/mlin/spVCF)

This repository has a command-line utility for encoding pVCF to spVCF and vice versa. The [Releases](https://github.com/mlin/spVCF/releases) page has pre-built executables compatible with most Linux x86-64 hosts, which you can download and `chmod +x spvcf`. To build and test it yourself, clone this repository and:

```
cmake . && make
ctest -V
```

The subcommands `spvcf encode` and `spvcf decode` encode existing pVCF to spVCF and vice versa. The input and output streams are uncompressed VCF text, so you may wish to arrange a pipe with [de]compression programs like `bgzip`.

```
spvcf encode [options] [in.vcf|-]
Reads VCF text from standard input if filename is empty or -

Options:
  -o,--output out.spvcf  Write to out.spvcf instead of standard output
  -n,--no-squeeze        Disable lossy QC squeezing transformation (lossless run-encoding only)
  -p,--period P          Ensure checkpoints (full dense rows) at this period or less (default: 1000)
  -t,--threads N         Use multithreaded encoder with this number of worker threads
  -q,--quiet             Suppress statistics printed to standard error
  -h,--help              Show this help message
```

```
spvcf decode [options] [in.spvcf|-]
Reads spVCF text from standard input if filename is empty or -

Options:
  -o,--output out.vcf  Write to out.vcf instead of standard output
  -q,--quiet           Suppress statistics printed to standard error
  -h,--help            Show this help message
```

Examples:

```
$ ./spvcf encode my.vcf > my.spvcf
$ bgzip -dc my.vcf.gz | ./spvcf encode | bgzip -c > my.spvcf.gz
$ bgzip -dc my.spvcf.gz | ./spvcf decode > my.decoded.vcf
```

There's also `spvcf squeeze` to apply the QC squeezing transformation to a pVCF, without the sparse quote-encoding. This produces valid pVCF that's typically much smaller, although not as small as spVCF.

The multithreaded encoder should be used only if the single-threaded version is a proven bottleneck. It's capable of higher throughput in favorable circumstances, but trades off memory usage and copying. The memory usage scales with threads and period.

### Tabix slicing

If the familiar `bgzip` and `tabix -p vcf` utilities are used to block-compress and index a spVCF file, then `spvcf tabix` can take a genomic range slice from it, extracting spVCF which decodes standalone. (The regular `tabix` utility generates the index, but using it to take the slice would yield a broken fragment.) Internally, this entails decoding a small bit of the spVCF, determined by the `spvcf encode --period` option.

```
spvcf tabix [options] in.spvcf.gz chr1:1000-2000 [chr2 ...]
Requires tabix index present e.g. in.spvcf.gz.tbi. Includes all header lines.

Options:
  -o,--output out.spvcf  Write to out.spvcf instead of standard output
  -h,--help              Show this help message
```

Example:

```
$ ./spvcf encode my.vcf | bgzip -c > my.spvcf.gz
$ tabix -p vcf my.spvcf.gz
$ ./spvcf tabix my.spvcf.gz chr21:5143000-5219900 | ./spvcf decode > slice.vcf
```
