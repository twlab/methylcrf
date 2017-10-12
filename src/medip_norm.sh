#!/bin/bash
#------------------------------------------------------------------------#
# Copyright 2012                                                         #
# Author: stevens _at_ cse.wustl.edu                                     #
#                                                                        #
# This program is free software: you can redistribute it and/or modify   #
# it under the terms of the GNU General Public License as published by   #
# the Free Software Foundation, either version 3 of the License, or      #
# (at your option) any later version.                                    #
#                                                                        #
# This program is distributed in the hope that it will be useful,        #
# but WITHOUT ANY WARRANTY; without even the implied warranty of         #
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the          #
# GNU General Public License for more details.                           #
#                                                                        #
# You should have received a copy of the GNU General Public License      #
# along with this program.  If not, see <http://www.gnu.org/licenses/>.  #
#------------------------------------------------------------------------#

# req 
#  mapBed

dip=$1;   # dipfn
cpg=$2;   # all cpg file

# sanity check
for i in $dip $cpg; do 
  if [ ! -e $i ]; then echo "$i not exist!" >&2; exit 1; fi;
done;

# sort both dip/mre and cpg bed file, so I can use mapbed instead of olapbed.
sort -k1,1V -k2,2n -o $dip $dip
sort -k1,1V -k2,2n -o $cpg $cpg

# calc reads over every CpG
dipreadbed=cpg_${dip/.bed/}_read.bed
mapBed -b $dip -a $cpg -o sum | awk '{OFS="\t";$4=".";$5=$NF;NF=6;print}' > $dipreadbed
# olapBed -s $dip $cpg|awk '{OFS="\t";$4=".";$5=$NF;NF=6;print}' > $dipreadbed
ls -l $dipreadbed >&2

# normalize to 75th percentile
p75=$(expr 25 \* $(awk '($5!=0){i++}END{print i}' $dipreadbed) / 100  )
p75cnt=$(awk '($5!=0){print $5}' $dipreadbed |sort -rn|head -$p75|tail -n 1)
echo "p75:$p75 P75cnt:$p75cnt" >&2

## ~8Gb
dipcpgbed=${dip/.bed/}_$(basename $cpg)
mapBed -b $cpg -a $dip -o count | awk '{OFS="\t";$5=$NF;NF=5;print}' > $dipcpgbed
# olapBed -c $cpg $dip | awk '{OFS="\t";$5=$NF;NF=5;print}' > $dipcpgbed

# split count amongst cpgs and normalize (mutiplicatively) so 75th percentile will be 10
# NOTE: MeDIP has chrM, cpg doesnt
awk 'BEGIN{OFS="\t"}($5){$5=(10/P)*(1/$5)}{print}' P=$p75cnt $dipcpgbed
# awk '($5){$5=(10/P)*(1/$5)}{OFS="\t";;print}' P=$p75cnt $dipcpgbed  # origninal, seems wrong?
 


