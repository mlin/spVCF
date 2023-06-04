# Sparse Project VCF (spVCF)

**Maintainer: Mike Lin [@DNAmlin](https://twitter.com/DNAmlin)**

Project VCF (pVCF; aka multi-sample VCF or population VCF) is the prevailing file format for small genetic variants discovered by cohort sequencing. It encodes a two-dimensional matrix with variant sites down the rows and study participants across the columns, filled in with all the genotypes and associated QC measures (read depths, genotype likelihoods, etc.). Large cohorts harbor many rare variants, implying a sparse genotype matrix composed largely of reference-identical or non-called cells. But the dense pVCF format encodes this inefficiently, growing super-linearly with the cohort size.

spVCF is an evolution of VCF that keeps most aspects of its tab-delimited text format, but presents the genotype matrix sparsely, by selectively reducing QC measure entropy and run-length encoding repetitive information about reference coverage. This is less sophisticated than some other efforts to address VCF's density and other shortcomings, but perhaps more palatable to existing VCF consumers by virtue of simplicity.

Further resources:

* Our [*Bioinformatics* Applications Note](https://doi.org/10.1093/bioinformatics/btaa1004) describing the approach, tool, and example results
* [doc/SPEC.md](https://github.com/mlin/spVCF/blob/master/doc/SPEC.md) has format details and a worked example
* [doc/compression_results.md](https://github.com/mlin/spVCF/blob/master/doc/compression_results.md) tests spVCF with *N*=50K exomes, observing up to 15X size reduction for bgzip-compressed pVCF, and scaling much more gently with *N*.
* [slide deck](https://docs.google.com/presentation/d/13lzEkdWAVwcsKofhsiYEdl92xMQgx5_dSOSIyZDggfM/edit?usp=sharing) presented at the GA4GH & MPEG-G Genome Compression Workshop, October 2018.
* [spVCF files for the resequenced 1000 Genomes Project cohort](https://github.com/mlin/spVCF/blob/master/doc/1000G_NYGC_GATK.md) (*N*=2,504 WGS)

## `spvcf` utility

[![build](https://github.com/mlin/spVCF/actions/workflows/build.yml/badge.svg?branch=main)](https://github.com/mlin/spVCF/actions/workflows/build.yml)

This repository has a command-line utility for encoding pVCF to spVCF and vice versa. The [Releases](https://github.com/mlin/spVCF/releases) page has pre-built executables compatible with most Linux x86-64 hosts, which you can download and `chmod +x spvcf`.

To build and test it locally, begin with a C++14 Linux development environment with CMake and [libdeflate](https://github.com/ebiggers/libdeflate). Clone this repository and:

```
cmake . && make
ctest -V
```

The subcommands `spvcf encode` and `spvcf decode` encode existing pVCF to spVCF and vice versa. The input and output streams are uncompressed VCF text, so you usually arrange a pipe with `bgzip`. Examples:

```
$ ./spvcf encode cohort.vcf > cohort.spvcf
$ bgzip -dc cohort.vcf.gz | ./spvcf encode | bgzip -c -@ $(nproc)  > cohort.spvcf.gz
$ bgzip -dc cohort.spvcf.gz | ./spvcf decode > cohort.decoded.vcf
```

Details:

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
  --with-missing-fields  Include trailing FORMAT fields with missing values
  -o,--output out.vcf    Write to out.vcf instead of standard output
  -q,--quiet             Suppress statistics printed to standard error
  -h,--help              Show this help message
```

There's also `spvcf squeeze` to apply the QC squeezing transformation to a pVCF, without the sparse quote-encoding. This produces valid pVCF that's typically much smaller, although not as small as spVCF.

The multithreaded encoder should be used only if the single-threaded version is a proven bottleneck. It's capable of higher throughput in favorable circumstances, but trades off memory usage and copying. The memory usage scales with threads, period, and *N*.

### Tabix slicing

If the familiar `bgzip` and `tabix -p vcf` utilities are used to block-compress and index a spVCF file, then `spvcf tabix` can take a genomic range slice from it, extracting spVCF which decodes standalone. (The regular `tabix` utility generates the index, but using it to take the slice would yield a broken fragment.) Example:

```
$ bgzip -dc cohort.vcf.gz | ./spvcf encode | bgzip -c -@ $(nproc)  > cohort.spvcf.gz
$ tabix -p vcf cohort.spvcf.gz
$ ./spvcf tabix cohort.spvcf.gz chr21:5143000-5219900 > slice.spvcf
$ ./spvcf decode slice.spvcf > slice.vcf
```

## Compatibility

spVCF is frequently used with project VCF files generated by [GATK GenotypeGVCFs](https://gatk.broadinstitute.org/hc/en-us/articles/360037057852-GenotypeGVCFs) and [GLnexus](https://github.com/dnanexus-rnd/GLnexus). Other joint-callers' products should work too, but aren't as routinely tested.

GLnexus now has a `--squeeze` command-line option to generate squeezed project VCF directly, which also speeds it up significantly. For large cohorts this should still be piped into `spvcf encode` to clean it up a little and add the run-encoding.

Squeezed project VCF (decoded from spVCF, or generated by `spvcf squeeze`) keeps the declarations of all FORMAT fields, but in most cells omits all except `GT` and `DP`. The rest aren't just marked missing, but omitted completely. Some downstream tools may be confused by this, even though the [VCF specification](https://samtools.github.io/hts-specs/VCFv4.3.pdf) expressly allows it (*"Trailing fields can be dropped..."*). We're advocating for more recognition of this useful, existing feature.
