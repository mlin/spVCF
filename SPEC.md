# Sparse Project VCF: format specification

**Maintainer: Mike Lin [@DNAmlin](https://twitter.com/DNAmlin)**

Project VCF (pVCF; aka multi-sample VCF) is the prevailing file format for small genetic variants discovered by cohort sequencing. It encodes a two-dimensional matrix with variant sites down the rows and study participants across the columns, filled in with all the genotypes and associated QC measures (read depths, genotype likelihoods, etc.). Large cohorts harbor many rare variants, implying a sparse genotype matrix composed largely of reference-homozygous or non-called cells. The dense pVCF format encodes this very inefficiently. See the [VCF specification](http://samtools.github.io/hts-specs/VCFv4.3.pdf) for full details of this format.

[Sparse Project VCF (spVCF)](https://github.com/mlin/spVCF) is a simple scheme to encode the pVCF matrix sparsely, by run-length encoding repetition which arises along both dimensions from pVCF production using tools like [GATK GenotypeGVCFs](https://software.broadinstitute.org/gatk/documentation/tooldocs/3.8-0/org_broadinstitute_gatk_tools_walkers_variantutils_GenotypeGVCFs.php) and [GLnexus](https://github.com/dnanexus-rnd/GLnexus). This is advantageous because the columnar repetition of reference coverage information isn't easily visible to generic compression algorithms. The encoding includes a checkpointing feature to enable random access within a block-compressed spVCF file, using familiar tools like `bgzip` and `tabix`.

In addition to this lossless encoding, spVCF suggests a convention for discarding QC measures in cells where they probably won't be useful, which markedly reduces data volume and increases the compressibility of what remains.

### Lossless encoding

### Internal checkpoints

### QC entropy reduction or "squeezing"
