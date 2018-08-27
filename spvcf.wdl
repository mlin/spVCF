task spvcf {
    File in_gz
    Boolean squeeze = false
    Boolean decode = false
    Boolean multithread = false
    String release = "v0.3.0"

    parameter_meta {
        in_gz: "stream"
    }

    command {
        set -ex -o pipefail

        apt-get update -qq && apt-get install -y -qq pigz wget
        wget -nv https://github.com/mlin/spVCF/releases/download/${release}/spvcf
        wget -nv https://github.com/dnanexus-rnd/GLnexus/raw/master/cli/dxapplet/resources/usr/local/bin/bgzip
        chmod +x spvcf bgzip

        threads_arg=""
        if [ "${multithread}" == "true" ]; then
            threads_arg="--threads $(nproc)"
        fi

        nm=$(basename "${in_gz}" .vcf.gz)
        nm=$(basename "$nm" .spvcf.gz)
        nm="$nm.${if decode then 'vcf.gz' else 'spvcf.gz'}"
        mkdir out
        pigz -dc "${in_gz}" | ./spvcf ${if decode then 'decode' else 'encode'} ${if squeeze then '-S' else ''} $threads_arg | ./bgzip -@ $(nproc) > out/$nm
    }

    runtime {
        docker: "ubuntu:18.04"
    }

    output {
        File out_gz = glob("out/*.gz")[0]
    }
}
