#!/bin/bash

# req 
#  olapBed

bed=$1;
cpg=$2;
win=$3;     # space separated window size list in quotes: "0 10 100"
outtmp=$4;  #out_fn template: <outtmp>_d[dist].cnt
fmt=2;      # printing number of significant aftr 

tmp=$(mktemp $0.XXXXXX.tmp)
for d in $win; do
  echo "$d" >&2
  awk '{OFS="\t";$2-=($2<D)?0:D;$3+=D;print}' D=$d $cpg > $tmp 
  olapBed -s $bed $tmp 2>$0.$bed.$d.err| awk '{printf("%s\t%."F"f\n",$4,$NF)}' F=$fmt >${outtmp}_d${d}.cnt
done;
rm $tmp

ls -l ${dip/.bed}_d${d}.cnt  >&2

