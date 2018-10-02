# Sparse Project VCF: format specification

**Maintainer: Mike Lin [@DNAmlin](https://twitter.com/DNAmlin)**

Project VCF (pVCF; aka multi-sample VCF) is the prevailing file format for small genetic variants discovered by cohort sequencing. It encodes a two-dimensional matrix with variant sites down the rows and study participants across the columns, filled in with all the genotypes and associated QC measures (read depths, genotype likelihoods, etc.). Large cohorts harbor many rare variants, implying a sparse genotype matrix composed largely of reference-homozygous or non-called cells. The dense pVCF format encodes this very inefficiently. See the [VCF specification](http://samtools.github.io/hts-specs/VCFv4.3.pdf) for full details of this format.

[Sparse Project VCF (spVCF)](https://github.com/mlin/spVCF) is a simple scheme to encode the pVCF matrix sparsely, by keeping most aspects of the VCF format while run-length encoding repetitive information about reference coverage. The encoding includes a checkpointing feature to facilitate random access within a block-compressed spVCF file, using familiar tools like `bgzip` and `tabix`. Lastly, spVCF suggests an optional convention to strip typically-unneeded QC details from the matrix.

### Sparse encoding

spVCF adopts from pVCF the tab-delimited text format with header, and the first nine columns providing all variant-level details. The sparse encoding concerns the genotype matrix `V[i,j]`, *i* indexing variant sites and *j* indexing the *N* samples, written across tab-delimited columns ten through 9+*N* of the pVCF text file. Each cell `V[i,j]` is a colon-delimited text string including the genotype and various QC measures (DP, AD, PL, etc.).

In the spVCF encoding, cells are first replaced with a double-quotation mark `"` if they're (i) identical to the cell *above*, and (ii) their GT field is reference-identical or non-called:

```
S[i,j] :=   "    if i>0 and V[i,j] == V[i-1,j] and V[i,j]["GT"] in ["0/0","0|0","./.",".|."],
          V[i,j] otherwise.
```

Here 'identical' covers all QC measures exactly. Such exact repetition is common in pVCF produced using tools like [GATK GenotypeGVCFs](https://software.broadinstitute.org/gatk/documentation/tooldocs/3.8-0/org_broadinstitute_gatk_tools_walkers_variantutils_GenotypeGVCFs.php) and [GLnexus](https://github.com/dnanexus-rnd/GLnexus), which merge gVCF or similar files summarizing reference coverage in lengthy bands.

For clarity, the list of "quotable" GTs enumerated above shows diploid genotypes only. In general, quotable GTs are those whose constituent allele calls are either all reference (0), or all non-called (.).

Second, within each row of `S`, consecutive runs of quotation marks are abbreviated with a text integer, so for example a horizontal run of 42 quotes is written `"42` and tab-delimited from adjacent cells. The result is a ragged, tab-delimited matrix.

**Worked example**

```
#CHROM  POS ID REF ALT ... FORMAT       Alice                           Bob                      Carol
22     1000  . A   G   ... GT:DP:AD:PL  0/0:35:35,0:0,117,402           0/0:29:29,0:0,109,387    0/0:22:22,0:0,63,188
22     1012  . CT  C   ... GT:DP:AD:PL  0/0:35:35,0:0,117,402           0/0:31:31,0:0,117,396    0/1:28:17,11:74,0,188
22     1018  . G   A   ... GT:DP:AD:PL  0/0:35:35,0:0,117,402           0/0:31:31,0:0,117,396    1/1:27:0,27:312,87,0
22     1074  . T   C,G ... GT:DP:AD:PL  0/0:33:33,0,0:0,48,62,52,71,94  ./.:0:0,0:.,.,.,.,.,.    1/2:42:4,20,18:93,83,76,87,0,77
```

encodes to

```
#CHROM  POS ID REF ALT ... FORMAT       Alice                           Bob                      Carol
22     1000  . A   G   ... GT:DP:AD:PL  0/0:35:35,0:0,117,402           0/0:29:29,0:0,109,387    0/0:22:22,0:0,63,188
22     1012  . CT  C   ... GT:DP:AD:PL           "                      0/0:31:31,0:0,117,396    0/1:28:17,11:74,0,188
22     1018  . G   A   ... GT:DP:AD:PL           "2                                              1/1:27:0,27:312,87,0
22     1074  . T   C,G ... GT:DP:AD:PL  0/0:33:33,0,0:0,48,62,52,71,94  ./.:0:0,0:.,.,.,.,.,.    1/2:42:4,20,18:93,83,76,87,0,77
```

Here some VCF features have been omitted for brevity, and for clarity the columns have been aligned artificially (for example, in the spVCF there would be only one tab delimiting `"2` and Carol's `1/1`).

Notice that a site with multiple alternate alleles usually breaks the columnar runs of repetition. We'll revisit this below.

### Decoding and checkpoints

spVCF is decoded back to pVCF by, first, copying out the header and the first line, which are identical to the original. On subsequent lines, the decoder copies out explicit cells and, upon encountering a quotation mark or an encoded run thereof, repeats the last-emitted cell from the respective column(s).

Decoding a given line of spVCF requires context from previous lines, potentially back to the beginning of the file. To enable random access within a spVCF file, the encoder should generate periodic *checkpoints*, which are simply pVCF lines copied verbatim without any run-encoding. Subsequent spVCF lines can be decoded by looking back no farther than the last checkpoint.

To facilitate finding the last checkpoint, the encoder prepends an INFO field to the eighth column of each non-checkpoint line, `spVCF_checkpointPOS=12345`, giving the VCF `POS` of the last checkpoint line. The decoder must remove this extra field from the output pVCF. A spVCF line is a checkpoint if and only if it lacks this `spVCF_checkpointPOS` field first in its INFO column. The first line for each reference contig (chromosome) must be a checkpoint, naturally including the first line of the file. 

With checkpoints, it's possible to reuse the familiar `bgzip` and `tabix` utilities with spVCF files. Compression and indexing use the original utilities as-is, while random access (genomic range slicing) requires specialized logic to construct self-contained spVCF from the whole original, locating a checkpoint and decoding from it as needed. The decoder seeking a checkpoint must accommodate the possibility that multiple VCF lines could share `POS` with the desired checkpoint.

### Optional: QC entropy reduction or "squeezing"

Lastly, spVCF suggests the following convention to remove typically-unneeded detail from the matrix, and increase the compressibility of what remains, prior to the sparse encoding. In any cell with QC measures indicating zero non-reference reads (typically `AD=d,0` for some *d*, but this depends on how the pVCF-generating pipeline expresses non-reference read depth), report only `GT` and `DP` and omit any other fields. Also, round `DP` down to a power of two (0, 1, 2, 4, 8, 16, ...).

This "squeezing" requires the encoder to reorder the colon-delimited fields in each cell so that `GT` and `DP` precede any other fields. Then it's valid for a subset of cells to omit remaining fields completely, as permitted by VCF. The FORMAT specification in column 9 of each line must reflect this reordering. Notice that not all reference-identical genotype calls are necessarily squeezed, namely if the QC data indicate even one non-reference read.

The optional squeezing transformation can be applied to any pVCF, usually to great benefit, whether or not the spVCF sparse encoding is also used.

Revisiting the worked example above reveals another benefit of squeezing, that it extends repetitive runs through sites with multiple alternate alleles.

```
#CHROM  POS ID REF ALT ... FORMAT       Alice                  Bob                      Carol
22     1000  . A   G   ... GT:DP:AD:PL  0/0:32                 0/0:16                   0/0:16
22     1012  . CT  C   ... GT:DP:AD:PL            "2                                    0/1:28:17,11:74,0,188
22     1018  . G   A   ... GT:DP:AD:PL            "2                                    1/1:27:0,27:312,87,0
22     1074  . T   C,G ... GT:DP:AD:PL            "            ./.:0                    1/2:42:4,20,18:93,83,76,87,0,77

```

