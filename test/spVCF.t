#!/bin/bash
set -o pipefail

HERE="$(dirname $0)"
BASH_TAP_ROOT="$HERE"
source "$HERE/bash-tap-bootstrap"

EXE=$HERE/../spvcf
D=/tmp/spVCFTests
rm -rf $D
mkdir -p $D

plan tests 11

pigz -dc "$HERE/data/small.vcf.gz" > $D/small.vcf
"$EXE" -o $D/small.spvcf $D/small.vcf
is "$?" "0" "filename I/O"
is "$(cat $D/small.spvcf | wc -c)" "36878360" "filename I/O output size"

pigz -dc "$HERE/data/small.vcf.gz" | "$EXE" -q > $D/small.spvcf
is "$?" "0" "piped I/O"
is "$(cat $D/small.spvcf | wc -c)" "36878360" "piped I/O output size"

"$EXE" -d -o $D/small.roundtrip.vcf $D/small.spvcf
is "$?" "0" "decode"
is "$(cat $D/small.roundtrip.vcf | wc -c)" "54007969" "roundtrip decode"

is "$(cat $D/small.vcf | grep -v ^# | sha256sum)" \
   "$(cat $D/small.roundtrip.vcf | grep -v ^# | sha256sum)" \
   "roundtrip fidelity"

"$EXE" -S -p 500 -o $D/small.squeezed.spvcf $D/small.vcf
is "$?" "0" "squeeze"
is "$(cat $D/small.squeezed.spvcf | wc -c)" "18647725" "squeezed output size"

"$EXE" -d -q -o $D/small.squeezed.roundtrip.vcf $D/small.squeezed.spvcf
is "$?" "0" "squeezed roundtrip decode"
is "$(cat $D/small.vcf | grep -v ^# | sed -r 's/(\t[^:]+):[^\t]+/\1/g' | sha256sum)" \
   "$(cat $D/small.squeezed.roundtrip.vcf | grep -v ^# | sed -r 's/(\t[^:]+):[^\t]+/\1/g' | sha256sum)" \
   "squeezed roundtrip GT fidelity"

#rm -rf $D
