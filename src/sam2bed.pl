#!/usr/bin/perl  -an
#------------------------------------------------------------------------#
# Copyright 2012                                                         #
# Author: stevens _at_ cse.wustl.edu (from Ting Wang original)           #
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

use strict;
use autodie;

# Name says it all...
#  filters out reads that are:
#    un-mapped
#    reads that align to chromo's not all numbers,X,Y,M,MT or starting with chr
#    with -q, reads below minq thresh
#  chrMT is converted to chrM


use vars qw(%unmapped $contigcnt $reads $hqreads $mappedreads %a);
BEGIN{ 
  while ($ARGV[0] =~ /^-/) {
    my $a=shift @ARGV;
    if ($a eq '-r') {$a{$a}=1          }
    else            {$a{$a}=shift @ARGV} 
  }
  # -c <number of bases to remove from begining of read>
  # -r  # report stats to stderr
  # -q <min mapq valu to accept>
}

## find call
#bam2mre_recheck.sh <max call to check> <bam> 
#(echo ~mstevens/data/genome/hg19/local
#
#
## get call 
#a=($(grep $b MRE.callcheck.121103.sum |sed -e 's/.MREcall[^ ]* / /'))
#
## combine and handle
# bam2bed.pl -q <minq> -c <call> <bam> 2>${b}_sam.err |\
# MRE_handler.pl - $fdir/MRE_${enz}enz_$sz.bed <libnm> |\
# sort -Vk1,3 > out.bed


# hdr line
next if (/^@/);

$reads++;
if (defined $a{'-q'} && $F[4] < $a{'-q'}) {next}
$hqreads++;

# doc says if set, everything is unreliable
# http://genome.sph.umich.edu/wiki/SAM
# MRE_handler sam_2_bed check chr name and whether there >= 10 fields instead
if ((2**2)&$F[1]) { $unmapped{"$F[1] $F[2]"}++ if ($a{-r});next }

my $chr = $F[2];
if ($chr =~ m/^(\d+|[XYM]|MT)$/) { $chr ="chr$chr"}
elsif ($chr =~ /^GL\d+/)         { $contigcnt++; next}           # dont know what else to do
elsif ($chr !~ /^chr/)           { warn "whatchrom:$_ \n"; next} # dont know what else to do

#match legacy
$chr =~ s/chrMT/chrM/;

my $strand=(2**4)&$F[1]? "-":"+";

my ($bgn,$end);
if ($strand eq '-') {
  $bgn = $F[3] -$a{-c}; 
  $end = $F[3] +length($F[9]) -1; 
} else {
  $bgn = $F[3]; 
  $end = $F[3] +length($F[9]) -1 +$a{-c}; 
}

#                     0-based         use call as score
print join("\t",$chr,$bgn-1,$end,$F[9],$a{-c},$strand),"\n";

END{ 
  if ($a{-r}) {
    print STDERR "reads         : $reads\n";
    print STDERR "hqreads       : $hqreads\n";
    print STDERR "unmapped:\n"; 
    print STDERR "$_: $unmapped{$_}\n" for (keys %unmapped);
    print STDERR "$contigcnt reads mapped to contigs\n";
  } 
}
