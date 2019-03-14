import "spvcf.wdl" as tasks

workflow test_spvcf {
    File vcf_gz # pVCF

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
    File gz1
    File gz2

    command {
        set -e -o pipefail
        cmp --silent <(gzip -dc "${gz1}") <(gzip -dc "${gz2}")
    }
}
