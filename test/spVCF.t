#!/bin/bash
set -o pipefail

HERE="$(dirname $0)"
BASH_TAP_ROOT="$HERE"
source "$HERE/bash-tap-bootstrap"

EXE=$HERE/../spvcf
D=/tmp/spVCFTests
rm -rf $D
mkdir -p $D

plan tests 7

pigz -dc "$HERE/data/small.vcf.gz" > $D/small.vcf
"$EXE" -o $D/small.ssvcf $D/small.vcf
is "$?" "0" "filename I/O"
is "$(cat $D/small.ssvcf | wc -c)" "36842025" "filename I/O output size"

pigz -dc "$HERE/data/small.vcf.gz" | "$EXE" -q > $D/small.ssvcf
is "$?" "0" "piped I/O"
is "$(cat $D/small.ssvcf | wc -c)" "36842025" "piped I/O output size"

"$EXE" -d -o $D/small.roundtrip.vcf $D/small.ssvcf
is "$?" "0" "decode"
is "$(cat $D/small.roundtrip.vcf | wc -c)" "54007969" "roundtrip decode"

is "$(cat $D/small.vcf | grep -v ^# | sha256sum)" \
   "$(cat $D/small.roundtrip.vcf | grep -v ^# | sha256sum)" \
   "roundtrip fidelity"

rm -rf $D
