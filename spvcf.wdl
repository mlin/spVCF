task spvcf_encode {
    File vcf_gz
    Boolean multithread = false
    String release = "v0.7.0"

    parameter_meta {
        vcf_gz: "stream"
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

        nm=$(basename "${vcf_gz}" .vcf.gz)
        nm="$nm.spvcf.gz"
        mkdir out
        pigz -dc "${vcf_gz}" | ./spvcf encode $threads_arg | ./bgzip -@ $(nproc) > "out/$nm"
    }

    runtime {
        docker: "ubuntu:18.04"
    }

    output {
        File spvcf_gz = glob("out/*.gz")[0]
    }
}

task spvcf_decode {
    File spvcf_gz
    String release = "v0.7.0"

    parameter_meta {
        spvcf_gz: "stream"
    }

    command {
        set -ex -o pipefail

        apt-get update -qq && apt-get install -y -qq pigz wget
        wget -nv https://github.com/mlin/spVCF/releases/download/${release}/spvcf
        wget -nv https://github.com/dnanexus-rnd/GLnexus/raw/master/cli/dxapplet/resources/usr/local/bin/bgzip
        chmod +x spvcf bgzip

        nm=$(basename "${spvcf_gz}" .spvcf.gz)
        nm="$nm.vcf.gz"
        mkdir out
        pigz -dc "${spvcf_gz}" | ./spvcf decode | ./bgzip -@ $(nproc) > "out/$nm"
    }

    runtime {
        docker: "ubuntu:18.04"
    }

    output {
        File vcf_gz = glob("out/*.gz")[0]
    }
}

task spvcf_squeeze {
    File vcf_gz
    Boolean multithread = false
    String release = "v0.7.0"

    parameter_meta {
        vcf_gz: "stream"
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

        nm=$(basename "${vcf_gz}" .vcf.gz)
        nm="$nm.squeeze.vcf.gz"
        mkdir out
        pigz -dc "${vcf_gz}" | ./spvcf squeeze $threads_arg | ./bgzip -@ $(nproc) > "out/$nm"
    }

    runtime {
        docker: "ubuntu:18.04"
    }

    output {
        File squeeze_vcf_gz = glob("out/*.gz")[0]
    }
}
