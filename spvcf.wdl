task spvcf {
    File in_gz
    Boolean squeeze = false
    Boolean decode = false
    String release = "v0.2.3"

    parameter_meta {
        in_gz: "stream"
    }

    command {
        set -ex -o pipefail

        apt-get update -qq && apt-get install -y -qq pigz wget
        wget -nv https://github.com/mlin/spVCF/releases/download/${release}/spvcf
        wget -nv https://github.com/dnanexus-rnd/GLnexus/raw/master/cli/dxapplet/resources/usr/local/bin/bgzip
        chmod +x spvcf bgzip

        nm=$(basename "${in_gz}" .vcf.gz)
        nm=$(basename "$nm" .spvcf.gz)
        nm="$nm.${if decode then 'vcf.gz' else 'spvcf.gz'}"
        mkdir out
        pigz -dc "${in_gz}" | ./spvcf ${if squeeze then '-S' else ''} ${if decode then '-d' else ''} | ./bgzip -@ $(nproc) > out/$nm
    }

    runtime {
        docker: "ubuntu:18.04"
    }

    output {
        File out_gz = glob("out/*.gz")[0]
    }
}
