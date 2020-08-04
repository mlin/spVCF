version 1.0

import "spvcf.wdl" as tasks

workflow test_spvcf {
    input {
        File vcf_gz # pVCF
    }

    # spVCF-encode the pVCF
    call tasks.spvcf_encode {
        input:
            vcf_gz = vcf_gz
    }

    # also squeeze the pVCF (without run-encoding)
    call tasks.spvcf_squeeze {
        input:
            vcf_gz = vcf_gz
    }

    # decode the spVCF back to squeezed pVCF
    call tasks.spvcf_decode {
        input:
            spvcf_gz = spvcf_encode.spvcf_gz
    }

    # verify decoded pVCF is identical to the squeezed pVCF
    # (modulo arbitrary differences in compression block framing)
    call verify_identical_gz_content {
        input:
            gz1 = spvcf_squeeze.squeeze_vcf_gz,
            gz2 = spvcf_decode.vcf_gz
    }

    output {
        File spvcf_gz = spvcf_encode.spvcf_gz
        File squeeze_vcf_gz = spvcf_squeeze.squeeze_vcf_gz
        File decoded_vcf_gz = spvcf_decode.vcf_gz
    }
}

task verify_identical_gz_content {
    input {
        File gz1
        File gz2
    }

    command <<<
        set -euxo pipefail
        apt-get -qq update && apt-get install -y tabix
        cmp --silent <(bgzip -dc "~{gz1}") <(bgzip -dc "~{gz2}")
    >>>

    runtime {
        docker: "ubuntu:20.04"
        cpu: 4
        memory: "4 GB"
        disks: "local-disk ~{ceil(size(gz2,'GB')+size(gz1,'GiB'))+4} SSD"
    }
}
