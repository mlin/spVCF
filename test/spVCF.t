#!/bin/bash
set -o pipefail

HERE="$(dirname $0)"
BASH_TAP_ROOT="$HERE"
source "$HERE/bash-tap-bootstrap"

EXE=$HERE/../spvcf
D=/tmp/spVCFTests
rm -rf $D
mkdir -p $D

plan tests 25

pigz -dc "$HERE/data/small.vcf.gz" > $D/small.vcf
"$EXE" encode --no-squeeze -o $D/small.spvcf $D/small.vcf
is "$?" "0" "filename I/O"
is "$(cat $D/small.spvcf | wc -c)" "37095284" "filename I/O output size"

pigz -dc "$HERE/data/small.vcf.gz" | "$EXE" encode -n -q > $D/small.spvcf
is "$?" "0" "piped I/O"
is "$(cat $D/small.spvcf | wc -c)" "37095284" "piped I/O output size"

"$EXE" decode -o $D/small.roundtrip.vcf $D/small.spvcf
is "$?" "0" "decode"
is "$(cat $D/small.roundtrip.vcf | wc -c)" "54007969" "roundtrip decode"

is "$(cat $D/small.vcf | grep -v ^# | sha256sum)" \
   "$(cat $D/small.roundtrip.vcf | grep -v ^# | sha256sum)" \
   "roundtrip fidelity"

is "$(egrep -o "spVCF_checkpointPOS=[0-9]+" $D/small.spvcf | uniq | cut -f2 -d = | tr '\n' ' ')" \
   "5030088 5142716 5232967 5252665 5273811 " \
   "checkpoint positions"

"$EXE" encode -p 500 -o $D/small.squeezed.spvcf $D/small.vcf
is "$?" "0" "squeeze"
is "$(cat $D/small.squeezed.spvcf | wc -c)" "17553478" "squeezed output size"

"$EXE" decode -q -o $D/small.squeezed.roundtrip.vcf $D/small.squeezed.spvcf
is "$?" "0" "squeezed roundtrip decode"
is "$(cat $D/small.vcf | grep -v ^# | sed -r 's/(\t[^:]+):[^\t]+/\1/g' | sha256sum)" \
   "$(cat $D/small.squeezed.roundtrip.vcf | grep -v ^# | sed -r 's/(\t[^:]+):[^\t]+/\1/g' | sha256sum)" \
   "squeezed roundtrip GT fidelity"

"$EXE" squeeze -q -o $D/small.squeezed_only.vcf $D/small.vcf
is "$?" "0" "squeeze (only)"
is "$(cat $D/small.squeezed_only.vcf | grep -v ^# | sha256sum)" \
   "$(cat $D/small.squeezed.roundtrip.vcf | grep -v ^# | sha256sum)" \
   "squeeze (only) fidelity"

is "$(egrep -o "spVCF_checkpointPOS=[0-9]+" $D/small.squeezed.spvcf | uniq | cut -f2 -d = | tr '\n' ' ')" \
   "5030088 5085728 5142746 5225415 5232998 5243839 5252707 5264483 5273919 " \
   "squeezed checkpoint positions"

bgzip -c $D/small.squeezed.spvcf > $D/small.squeezed.spvcf.gz
tabix $D/small.squeezed.spvcf.gz
"$EXE" tabix -o $D/small.squeezed.slice.spvcf $D/small.squeezed.spvcf.gz chr21:5143000-5226000
is "$?" "0" "tabix slice"

is "$(egrep -o "spVCF_checkpointPOS=[0-9]+" $D/small.squeezed.slice.spvcf | uniq -c | tr -d ' ' | tr '\n' ' ')" \
   "497spVCF_checkpointPOS=5143363 28spVCF_checkpointPOS=5225415 " \
   "slice checkpoint"

"$EXE" decode $D/small.squeezed.slice.spvcf > $D/small.squeezed.slice.vcf
is "$?" "0" "decode tabix slice"

bgzip -c $D/small.squeezed.roundtrip.vcf > $D/small.squeezed.roundtrip.vcf.gz
tabix $D/small.squeezed.roundtrip.vcf.gz
tabix $D/small.squeezed.roundtrip.vcf.gz chr21:5143000-5226000 > $D/small.squeezed.roundtrip.slice.vcf
is "$(cat $D/small.squeezed.slice.vcf | grep -v ^# | sha256sum)" \
   "$(cat $D/small.squeezed.roundtrip.slice.vcf | grep -v ^# | sha256sum)" \
   "slice fidelity"

"$EXE" tabix -o $D/small.squeezed.slice_chr21.spvcf $D/small.squeezed.spvcf.gz chr21
is "$(cat $D/small.squeezed.slice_chr21.spvcf | sha256sum)" \
   "$(cat $D/small.squeezed.spvcf | sha256sum)" \
   "chromosome slice"

pigz -dc "$HERE/data/small.vcf.gz" | "$EXE" encode -n -t $(nproc) - > $D/small.mt.spvcf
is "$?" "0" "multithreaded encode"
is "$(cat $D/small.mt.spvcf | wc -c)" "37091691" "multithreaded output size"

time "$EXE" decode -o $D/small.mt.roundtrip.vcf $D/small.mt.spvcf
is "$?" "0" "decode from multithreaded"

is "$(cat $D/small.vcf | grep -v ^# | sha256sum)" \
   "$(cat $D/small.mt.roundtrip.vcf | grep -v ^# | sha256sum)" \
   "multithreaded roundtrip fidelity"

is "$(egrep -o "spVCF_checkpointPOS=[0-9]+" $D/small.mt.spvcf | uniq | cut -f2 -d = | tr '\n' ' ')" \
   "5030088 5139057 5232476 5252288 5273680 " \
   "multithreaded checkpoint positions"

rm -rf $D
