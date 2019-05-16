# spVCF for the resequenced 1000 Genomes Project

### Updated: 15 May 2019

The 2,504 individuals from the 1000 Genomes Project phase 3 cohort were [recently resequenced](https://twitter.com/notSoJunkDNA/status/1125401248348495873) at the New York Genome Center (to high depth on modern instruments). Accompanying this awesome data drop were "working" versions of [project VCF](http://ftp.1000genomes.ebi.ac.uk/vol1/ftp/data_collections/1000G_2504_high_coverage/working/20190425_NYGC_GATK/) files [joint-called using GATK](https://software.broadinstitute.org/gatk/best-practices/workflow?id=11145), with total size of roughly 1,250GB. We encoded and recompressed these using [spVCF](https://www.biorxiv.org/content/10.1101/611954v1) to produce files about 1/5 the size with (we contend) minimal loss of useful information.

* [chr21 spvcf.gz on figshare](https://figshare.com/articles/1000G_2504_high_coverage_20190425_NYGC_GATK_chr21_spvcf_gz/8132261) (3.7GB compared to 18GB original pvcf.gz)

**Help wanted:** please [get in touch](https://twitter.com/DNAmlin) if you could host & distribute the full dataset on our behalf (about 250GB). Right now we have it stored in a couple places which understandably are reluctant to provide free-for-all public downloads (the 80% reduced size notwithstanding!).

We're really grateful to all involved in the sponsorship & production of the new sequencing results, a tremendous resource for global collaboration amongst genome geeks everywhere!

### Example usage

The easiest way to start using the spVCF file is with [our command-line tool](https://github.com/mlin/spVCF/releases) to decode it on-the-fly for existing tools that consume project VCF files. Here are some example invocations, supposing you've just downloaded `1000G_2504_high_coverage.20190425_NYGC_GATK.chr21.spvcf.gz`.

```bash
# Fetch & enable spvcf utility
wget https://github.com/mlin/spVCF/releases/download/v1.0.0/spvcf
chmod +x spvcf

# For brevity below:
mv 1000G_2504_high_coverage.20190425_NYGC_GATK.chr21.spvcf.gz chr21.spvcf.gz

# Generate tabix index of the spVCF file
tabix -p vcf chr21.spvcf.gz

# Slice & decode a range for use with plain bcftools
./spvcf tabix chr21.spvcf.gz chr21:7000000-8000000 | bgzip -c > chr21.slice.spvcf.gz
gzip -dc chr21.slice.spvcf.gz | ./spvcf decode | bgzip -c > chr21.slice.vcf.gz
bcftools stats chr21.slice.vcf.gz

# Repeat with a much better streaming approach
bcftools stats <(./spvcf tabix chr21.spvcf.gz chr21:7000000-8000000 | ./spvcf decode)
```

These invocations report transition/transversion ratio of 1.17 which is super low, suggesting that downstream analyses need to look at the filters carefully. Updated and additional variant call sets for these genomes are no doubt coming from numerous groups, and we hope spVCF might facilitate their exchange.

### Afterthought

spVCF reduces project VCF file size in two ways, first through selective removal & rounding of QC measures ("squeezing"), and second by run-length encoding the remainder. For this dataset, nearly all of the 80% size reduction arises from the lossy squeezing step; restated, the `vcf.gz` decoded from a given `spvcf.gz` is only slightly larger (about 10%). The run encoding is ineffective with "only" *N*=2,504 genomes because the variant sites are spaced far enough apart that there are not many long identical runs (in the way spVCF sees them) to encode. At *N*=50K the run-encoding halves the size after squeezing, and quarters it at *N*=100K. So, while this dataset didn't exhibit spVCF to its full potential, we hope it provides a useful familiarization exercise for things to come.

### Appendix: workflow

The following [WDL](http://openwdl.org/) workflow generated the spVCF files. The task downloads and encodes the pVCF for one chromosome. The workflow scatters this task across the chromosomes.

```wdl
version 1.0

task spvcf_encode_1000G_chrom {
  input {
    String chrom
  }
  command <<<
    set -eux -o pipefail
    apt-get update && apt-get -y install aria2 pigz tabix
    aria2c https://github.com/mlin/spVCF/releases/download/v1.0.0/spvcf
    chmod +x spvcf
    aria2c -s 4 -j 4 -x 4 --file-allocation=none "http://ftp.1000genomes.ebi.ac.uk/vol1/ftp/data_collections/1000G_2504_high_coverage/working/20190425_NYGC_GATK/CCDG_13607_B01_GRM_WGS_2019-02-19_~{chrom}.recalibrated_variants.vcf.gz"
    pigz -dc *.vcf.gz | ./spvcf encode | bgzip -c > "1000G_2504_high_coverage.20190425_NYGC_GATK.~{chrom}.spvcf.gz"
  >>>
  output {
    File spvcf_gz = glob("*.spvcf.gz")[0]
  }
  runtime {
    docker: "ubuntu:18.04"
    preemptible: 2
    cpu: 2
    disks: "local-disk 160 HDD"
  }
}

workflow spvcf_encode_1000G {
  input {
    Array[String] chroms = [
         "chr1", "chr2", "chr3", "chr4",
         "chr5", "chr6", "chr7", "chr8",
         "chr9","chr10","chr11","chr12",
        "chr13","chr14","chr15","chr16",
        "chr17","chr18","chr19","chr20",
        "chr21","chr22", "chrX", "chrY",
        "others"
    ]
  }
  scatter (chrom in chroms) {
    call spvcf_encode_1000G_chrom as encode1 {
      input:
        chrom = chrom
    }
  }
  output {
    Array[File] spvcf_gz = encode1.spvcf_gz
  }
}
```
