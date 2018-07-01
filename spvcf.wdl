task spvcf {
    File pvcf_gz
    Boolean squeeze = false
    String release = "v0.2.3"

    parameter_meta {
        pvcf_gz: "stream"
    }

    command {
        set -ex -o pipefail

        apt-get update -qq && apt-get install -y -qq pigz wget
        wget -nv https://github.com/mlin/spVCF/releases/download/${release}/spvcf
        wget -nv https://github.com/dnanexus-rnd/GLnexus/raw/master/cli/dxapplet/resources/usr/local/bin/bgzip
        chmod +x spvcf bgzip

        nm=$(basename "${pvcf_gz}" .vcf.gz)
        mkdir out
        pigz -dc "${pvcf_gz}" | ./spvcf ${if squeeze then '-S' else ''} | ./bgzip -@ $(nproc) > out/$nm.spvcf.gz
    }

    runtime {
        docker: "ubuntu:18.04"
    }

    output {
        File spvcf_gz = glob("out/*.gz")[0]
    }
}
