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

bed=$1;
cpg=$2;
win=$3;     # space separated window size list in quotes: "0 10 100"
outtmp=$4;  #out_fn template: <outtmp>_d[dist].cnt
fmt=2;      # printing number of significant aftr 
sort -k1,1V -k2,2n -o $bed $bed
sort -k1,1V -k2,2n -o $cpg $cpg

tmp=$(mktemp $0.XXXXXX.tmp)
for d in $win; do
  echo "$d" >&2
  awk '{OFS="\t";$2-=($2<D)?$2:D;$3+=D;print}' D=$d $cpg > $tmp  # it should be this. The original code is $2-=($2<D)?0:D. Alternatively, use bedtools slop?
  # awk '{OFS="\t";$2-=($2<D)?0:D;$3+=D;print}' D=$d $cpg > $tmp
  # sort -k1,1V -k2,2n -o $tmp $tmp
  out=${outtmp}_d${d}.cnt
  mapBed -b $bed -a $tmp -o sum | awk '{printf("%s\t%."F"f\n",$4,$NF)}' F=$fmt >$out
  # olapBed -s $bed $tmp | awk '{printf("%s\t%."F"f\n",$4,$NF)}' F=$fmt >$out
#  ls -l $out  >&2
done;
rm $tmp


