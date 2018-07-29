# Sparse Project VCF: format specification

**Maintainer: Mike Lin [@DNAmlin](https://twitter.com/DNAmlin)**

Project VCF (pVCF; aka multi-sample VCF) is the prevailing file format for small genetic variants discovered by cohort sequencing. It encodes a two-dimensional matrix with variant sites down the rows and study participants across the columns, filled in with all the genotypes and associated QC measures (read depths, genotype likelihoods, etc.). Large cohorts harbor many rare variants, implying a sparse genotype matrix composed largely of reference-homozygous or non-called cells. The dense pVCF format encodes this very inefficiently. See the [VCF specification](http://samtools.github.io/hts-specs/VCFv4.3.pdf) for full details of this format.

[Sparse Project VCF (spVCF)](https://github.com/mlin/spVCF) is a simple scheme to encode the pVCF matrix sparsely, by run-length encoding repetition which arises along both dimensions from pVCF production using tools like [GATK GenotypeGVCFs](https://software.broadinstitute.org/gatk/documentation/tooldocs/3.8-0/org_broadinstitute_gatk_tools_walkers_variantutils_GenotypeGVCFs.php) and [GLnexus](https://github.com/dnanexus-rnd/GLnexus). This is advantageous because the columnar repetition of reference coverage information isn't easily visible to generic compression algorithms. The encoding includes a checkpointing feature to enable random access within a block-compressed spVCF file, using familiar tools like `bgzip` and `tabix`.

In addition to this lossless encoding, spVCF suggests a convention for discarding QC measures in cells where they probably won't be useful, which markedly reduces data volume and increases the compressibility of what remains.

### Lossless encoding

spVCF adopts identically from pVCF the tab-delimited text scheme with header, and the first nine columns providing all variant-level details. The sparse encoding affects the genotype matrix `V[i,j]`, *i* indexing variant sites and *j* indexing the *N* samples, written across tab-delimited columns ten through 9+*N* of the pVCF text file. Each entry `V[i,j]` is a colon-delimited text string including the genotype and various QC measures (DP, AD, PL, ...).

In the spVCF encoding, entries are replaced with a double-quotation mark `"` if they're identical to the entry *above*: `S[i,j] <- " if i>0 and V[i,j] == V[i-1,j]; V[i,j] otherwise`. Here 'identical' includes all QC measures exactly; such repetition is actually common in pVCF produced by merging gVCF files or other intermediates that convey reference coverage in lengthy bands.

Then, within each row of `S`, consecutive runs double-quotation marks are abbreviated with a text integer, so for example a horizontal run of 42 quotes is written `"42`, appropriately tab-delimited from adjacent entries. The result is a ragged matrix.

### Checkpoints

### QC entropy reduction or "squeezing"
