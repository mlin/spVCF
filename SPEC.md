# Sparse Project VCF: format specification

**Maintainer: Mike Lin [@DNAmlin](https://twitter.com/DNAmlin)**

Project VCF (pVCF; aka multi-sample VCF) is the prevailing file format for small genetic variants discovered by cohort sequencing. It encodes a two-dimensional matrix with variant sites down the rows and study participants across the columns, filled in with all the genotypes and associated QC measures (read depths, genotype likelihoods, etc.). Large cohorts harbor many rare variants, implying a sparse genotype matrix composed largely of reference-homozygous or non-called cells. The dense pVCF format encodes this very inefficiently. See the [VCF specification](http://samtools.github.io/hts-specs/VCFv4.3.pdf) for full details of this format.

[Sparse Project VCF (spVCF)](https://github.com/mlin/spVCF) is a simple scheme to encode the pVCF matrix sparsely, by run-length encoding repetition which arises along both dimensions from pVCF production using tools like [GATK GenotypeGVCFs](https://software.broadinstitute.org/gatk/documentation/tooldocs/3.8-0/org_broadinstitute_gatk_tools_walkers_variantutils_GenotypeGVCFs.php) and [GLnexus](https://github.com/dnanexus-rnd/GLnexus). This is advantageous because the columnar repetition of reference coverage information isn't easily visible to generic compression algorithms. The encoding includes a checkpointing feature to enable random access within a block-compressed spVCF file, using familiar tools like `bgzip` and `tabix`.

In addition to this lossless encoding, spVCF suggests a convention for discarding QC measures in cells where they probably won't be useful, which markedly reduces data volume and increases the compressibility of what remains.

### Sparse encoding

spVCF adopts from pVCF the tab-delimited text format with header, and the first nine columns providing all variant-level details. The sparse encoding concerns the genotype matrix `V[i,j]`, *i* indexing variant sites and *j* indexing the *N* samples, written across tab-delimited columns ten through 9+*N* of the pVCF text file. Each cell `V[i,j]` is a colon-delimited text string including the genotype and various QC measures (DP, AD, PL, etc.).

In the spVCF encoding, cells are replaced with a double-quotation mark `"` if they're identical to the cell *above*: 

```
S[i,j] :=   "    if i>0 and V[i,j] == V[i-1,j],
          V[i,j] otherwise.
```

Here 'identical' covers all QC measures exactly. (Such repetition is common in pVCF produced by merging gVCF files or other intermediates summarizing reference coverage in lengthy bands.)

Then, within each row of `S`, consecutive runs of quotation marks are abbreviated with a text integer, so for example a horizontal run of 42 quotes is written `"42`, tab-delimited from adjacent cells. The result is a ragged tab-delimited matrix.

Worked example:

### Decoding and checkpoints

spVCF is decoded back to pVCF by, first, repeating the first line, identical to the original by construction. On subsequent lines, the decoder copies out explicit cells and upon encountering a quotation mark or an encoded run thereof, repeats the last-emitted cell from the respective column(s).

Decoding a given line of spVCF generally requires contextual state from previous lines, potentially back to the beginning of the file. To expedite random access within a spVCF file, the encoder should also generate periodic *checkpoints*, which are simply pVCF lines copied verbatim without any run-encoding. Subsequent spVCF lines can be decoded by looking back no further than the last checkpoint.

To facilitate finding the last checkpoint, the encoder must prepend an INFO field to the eighth column of each non-checkpoint line, `spVCF_checkpointPOS=12345`, giving the VCF `POS` of the last checkpoint line. The decoder must remove this extra field from the output pVCF. A spVCF line is a checkpoint if and only if it lacks this `spVCF_checkpointPOS` field first in its INFO column.

The first line for each reference contig (chromosome) must be a checkpoint, naturally including the first line in the file.

(Using `POS` for the pointer, rather than line numbers, is convenient for reusing tabix for random access within a block-compressed spVCF file.)

### QC entropy reduction or "squeezing"

Lastly, spVCF suggests the following convention to remove typically-unneeded detail from the matrix, and increase the compressibility of what remains, prior to the sparse encoding discussed above.

In any cell with QC measures indicating zero non-reference reads (typically `AD=d,0` for some *d*, but this depends on how the pVCF-generating pipeline expresses non-reference read depth), keep only `GT` and `DP` and omit any other fields. Also, round `DP` down to a power of two (0, 1, 2, 4, 8, 16, ...).

This "squeezing" requires the encoder to reorder the colon-delimited fields in each cell so that `GT` and `DP` precede any other fields. This makes it valid for a subset of cells to omit the other fields completely, as permitted in VCF.
