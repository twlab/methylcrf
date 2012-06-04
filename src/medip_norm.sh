#!/bin/bash

# req 
#  olapBed

dip=$1;   # dipfn
cpg=$2;   # all cpg file

# sanity check
for i in $dip $cpg; do 
  if [ ! -e $i ]; then echo "$i not exist!" >&2; exit 1; fi;
done;


# calc reads over every CpG
dipreadbed=cpg_${dip/.bed/}_read.bed
olapBed -s $dip $cpg|awk '{OFS="\t";$5=$NF;NF=6;print}' > $dipreadbed
ls -l $dipreadbed >&2

# normalize to 75th percentile
p75=$(expr 25 \* $(awk '($5!=0){i++}END{print i}' $dipreadbed) / 100  )
p75cnt=$(awk '($5!=0){print $5}' $dipreadbed |sort -rn|head -$p75|tail -n 1)
echo "p75:$p75 P75cnt:$p75cnt" >&2

## ~8Gb
dipcpgbed=${dip/.bed/}_$cpg
olapBed -c $cpg $dip | awk '{OFS="\t";$5=$NF;NF=5;print}' > $dipcpgbed

# split count amongst cpgs and normalize (mutiplicatively) so 75th percentile will be 10
# NOTE: MeDIP has chrM, cpg doesnt
awk '($5){$5=(10/P)*(1/$5)}{OFS="\t";;print}' P=$p75cnt $dipcpgbed
 



