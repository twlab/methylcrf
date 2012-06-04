#!/bin/bash

# req 
#  olapBed

tmp=$(mktemp $0.XXXXXX.tmp)
cat $1/* |cut -f1,2,3|sort -u >$tmp
olapBed -c $tmp $2 |awk '{OFS="\t";print $4,!!$7}' 
rm $tmp
