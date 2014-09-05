#!/usr/bin/env perl
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
use List::Util qw(sum reduce);our ($a, $b); 
use List::MoreUtils qw(uniq);
use Data::Dumper;

my $usage = '
MRE_handler.pl [options] <expr bed file or - > <virtual digest file> <prefix for additional output files>

[options]
 -vdat             : virtial digest file is a dumped struct (ie, not a bed file)
 -vdump            : process virtual digest file and dump to stdout
expr bed file      : should have the call issue already fixed
virtual digest file: file of all the virtual enzyme fragments combined. should include the whole 
                     recognition site on both sides. [should change this to include just the cut 
                     fragment].  and have the enzyme name as the id


MRE script, does the following:

Main output: bed file: for each CpG site, generate a score (read per site per million mapped per enzyme)

Additional Files:
1) expr.filtered.bed: possible fragments,  ie those matching the fragment file


2) expr.report: generate additional statitics:
   No of total reads;
   No of reads passed quality filter (.bed file)
   No of reads passed MRE filter (.filtered.bed file)
   No of solo reads (no paired reads from expected fragments)
   Distribution of size of all fragments with both ends mapped
   How unbalancing are the reads: No of reads mapped to one end (more) vs mapped to the other end (fewer)

NOTE: BStUI gives information about 2 CpGs.  These are both counted.  CGCG and GCGC can overlap, so if a read could belong to either, each read is fractionally counted.

';


my $vdat;
while ($ARGV[0] =~ /(-[^ ])/) {
  my $arg = shift @ARGV;
  if ($arg eq '-vdat') {
    $vdat="load"
  } elsif ($arg eq '-vdump') {
    $vdat="dump";
  }
}

die $usage unless (@ARGV > 2);
my ( $exp_f, $vdigest, $name, $printidx) = @ARGV;

# enzyme cut sit2s
#my %enz =( 'CCGG',1, 'GCGG',1, 'GCGC',1, 'ACGT',1, 'CGCG',2,);
my $total;
my $unique;
my $q          = 0;
my $sam_f      = $name.".sam";
my $filtered_f = $name.".filtered.bed";
my $cg_f       = $name."_MRE.bed";
my $report_f   = $name.".report";
my $normconst  = 1000000;

open ( REPORT, ">$report_f" );

# hash that keys in paired ends $PET{chr|start|strand} = chr|start|strand where start is where the read maps to

#PE      # hash for paired RE sites
#cnt     # count for reads at end of a fragment
#CpG     # count for CpG by CCGG reads (may include reads of two directions)
#reads	 # number of reads from CCGG site, sum of %CCGG_cnt
my (%cnt,%CpG,%reads);


warn "build_PE_hash\n";
my ($PE,$enzcnt);
if ($vdat eq 'load') {
  my $VAR1 = do { 
    no strict 'vars';
    eval qx(cat $vdigest);
  };
 ($PE,$enzcnt) = @$VAR1;
1;
} else {
  ($PE,$enzcnt) =build_PE_hash($vdigest);
  if ($vdat eq 'dump') {
  print Dumper([$PE,$enzcnt]);
  exit;
  }
}

warn "filter reads by mre sites\n";
# readin bed, write xxx_full.index.bed, gens CGCG_CpG;
filter_reads_by_mre_sites( $exp_f, $filtered_f );

warn "calculate CpG score\n";
calculate_CpG_score( $cg_f );

warn "fragment stats\n";
fragment_stats();

# site: refers to the begining of the read you would expect (0-based)
# PE would be the read you expect on the - strand  (1-based)
sub build_PE_hash {
  open ( IN, $_[0] );
  my (%rtn,%enzcnt);

  my (@line, $site, $pe_site);
  while ( <IN> ) {
    @line   = split /\t/;
    my $enz = uc($line[3]); $enz = 'GCGG' if ($enz eq 'CCGC');
    $enzcnt{$enz}++;

    $site    = join "|",$line[0],$line[1],"+";
    $pe_site = join "|",$line[0],$line[2],"-"; # 0-base it?

    # 120822 choose the smallest..(equivalent to assuming complete digestion
    # NOTE: PE site is always - (even if a gcgg site in -),
    #       so, can safely rmv old pe_site (if the smaller frag exists it will be in the list
    #if (!exists $rtn{$enz}{ $site } || (split '\|',$rtn{$enz}{ $site })[1] > $line[2]) { $rtn{$enz}{ $site }  =$pe_site}
    if (exists $rtn{$enz}{ $site }) {
      my $pe =$rtn{$enz}{ $site };
      if ((split '\|',$pe)[1] > $line[2]) { undef $rtn{$enz}{$pe}; }
      else                                { $pe_site = $pe}
    }
    $rtn{$enz}{ $site }    =$pe_site;
    $rtn{$enz}{ $pe_site } =$site;

#  1+ -> 7-
#  7- -> 1+
#  [1     7]
# A[CGTA CG]TA CG T 
# T[GCAT GC]AT GC A
    #
#       [5      11]
# A CGTA[CG TA CG]T 
# T GCAT[GC AT GC]A
    #
#              9
# A CGTA CG TA[CG T 
# T GCAT GC AT[GC A

  }
  close IN;

  return (\%rtn,\%enzcnt);
} 

#NOTE: could have olap between CGCG and GCGC
#      so should only count fractional read and site for each enzyme
#     ie, you don't know which, so spread across all -should fix by indexing MRE's libs
sub filter_reads_by_mre_sites {
  my ( $in_f, $out_f ) = @_;

  open ( IN, $in_f );
  open ( OUT, ">$out_f" );

  my @fnnm = ('Unknown'); push @fnnm, (keys %$PE) if ($printidx);

  my %fh = map {
    my $fn=$name."_".$_."_fullindex.txt";
    open my $fh, '>',$fn;
    ($_,$fh)  
  }  @fnnm;

  while ( <IN> ) {
    chomp;
    my @line = split /\t/;

    my ($site,$CpG);
    if ( $line[5] eq '+' ) {
      $site = $line[0].'|'.$line[1].'|+';
      $CpG  = $line[0].'|'.$line[1];
    } elsif ( $line[5] eq "-" ) {
      $site = $line[0].'|'.$line[2].'|-';
      $CpG  = $line[0].'|'.($line[2]-2); # converts to 0-based at bgn of CG
    } else { die "dont know strand [$line[5]] line:$.  [$_]"}

   chomp; print join("\t",$_,$site,$CpG),"\n";
    my $line = join("\t", @line[0..3,5])."\n";

    my @cuts = grep { exists $PE->{$_}{$site} } keys %$PE;
    if  (@cuts) {
      print OUT $_,"\n";

      # add to each enzyme stats and files too
      for my $enz (@cuts) {
        $cnt{$enz}{ $site } += 1/@cuts;
        $reads{$enz}        += 1/@cuts;
        $CpG{$enz}{ $CpG }++;  # this is averaged later in calculated CpG score
        print {$fh{$enz}} $line if (exists $fh{$enz}); 

        #MS: also report the 1st (or last if "-") CG (it must have been !mC too)
        if ($enz eq 'CGCG') {
          my $psn = $line[5] eq '+' ? $line[1]-2 : $line[2];
          $CpG{$enz}{ join "|", $line[0], $psn }++;
        }
      }
    } else{ print {$fh{'Unknown'}}  $line;$reads{'unknown'}++}
  }

  my $sum = reduce{ $a+$reads{$b} } (0,keys %$PE);

  print REPORT "mre filtered reads:\t$sum\n";
  for (keys %$PE) { print REPORT "\t$_ reads:\t", $reads{$_}, "\n"; }
  print REPORT "\tUnknown reads: $reads{'unknown'}\n";

  close IN; close OUT;  
}

sub calculate_CpG_score {
  my ( $out_f ) = @_;

  my %mil= map{($_, $reads{$_}/$normconst) } keys %$PE;

  my @CpG=sort uniq ( map {keys %{$CpG{$_}}} keys %$PE );

  print REPORT "Sampled CpG sites:\t", scalar @CpG, "\n";

  open ( OUT, ">$out_f" );
  # need to take avg for any CpG sampled by multiple enzymes
  for ( my $i=0; $i<=$#CpG; $i++ ) {

    my @rcme = map { $CpG{$_}{ $CpG[$i] } / $mil{$_} } 
                (grep {exists $CpG{$_}{$CpG[$i]}} keys %$PE);

    my $rcme= sum @rcme; $rcme /= scalar @rcme;

    my @line = split /\||\n/, $CpG[$i];
    printf OUT "%s\t%.4f\n",join("\t", $line[0],$line[1], $line[1]+2,"."), $rcme;
  }
  close OUT;
} 

sub fragment_stats {
  my $high_end   = 0;
  my $low_end    = 0;
  my $solo_site  = 0;
  my $solo_reads = 0;
  my $pair_site  = 0;
  my @frags;
  my %paired; # keep a record of ends that got paired

  for my $enz (keys %$PE) {
    for my $key ( keys %{$cnt{$enz}} ) {
      my $pe = $PE->{$enz}{$key};

      if ( defined( $paired{$key} ) ) {next}

      if ( defined( $cnt{$enz}{$pe} ) ) {
        # a valid fragment
        $pair_site++;
        my @left  = split /\|/, $key;
        my @right = split /\|/, $pe;
        push @frags, abs($right[1]-$left[1]);
        $paired{ $pe } = 1;
        if ( $cnt{$enz}{$key} >= $cnt{$enz}{$pe} ) {
          $high_end += $cnt{$enz}{$key};
          $low_end  += $cnt{$enz}{$pe};
        } else {
          $high_end += $cnt{$enz}{$pe};
          $low_end  += $cnt{$enz}{$key};
        }
      } else {
        $solo_reads += $cnt{$enz}{$key};
        $solo_site++;
      }
    }
  }

  print REPORT "solo ends:\t$solo_site\n";
  print REPORT "    reads on solo ends:\t$solo_reads\n";
  print REPORT "fragments:\t$pair_site\n";
  print REPORT "    reads on higher end of fragments:\t$high_end\n";
  print REPORT "    reads on lower end of fragments:\t$low_end\n";
  print REPORT "fragment size distribution:\n";

  histograph( \@frags, 20, 40, 400 );

}

sub histograph {
  my ( $data_r, $win, $min, $max ) = @_;

  my @histo;
  my @scale;
  my $total = 0;
  my $histo_min = 0;
  my $histo_max = 0;

  for ( my $i=$min, my $j=0; $i<=$max; $i+=$win, $j++ ) {
    $scale[$j][0] = $i;
    $scale[$j][1] = $i+$win;
    $histo[$j] = 0;
  }
  $scale[$#scale][1] = $max;

  for ( my $i=0; $i<=$#$data_r; $i++ ) {
    $total++;
    if ( $data_r->[$i] <= $min )     { $histo_min++;
    } elsif ( $data_r->[$i] > $max ) { $histo_max++;
    } else {
      for ( my $j=0; $j<=$#scale; $j++ ) {
        if ( $data_r->[$i]>$scale[$j][0] && $data_r->[$i]<=$scale[$j][1] ) {
          $histo[$j]++;
          last;
        }
      }
    }
  }

  printf REPORT "%s\t%10s\t%10s\t|\n", "Scale", "Count", "Percent";
  printf REPORT "%s\t%10d\t%10.2f\t|", "<=$min", $histo_min, $histo_min/$total*100;
  print_bar ( int($histo_min/$total*100) );
  for ( my $i=0; $i<=$#histo; $i++ ) {
    printf REPORT "%s\t%10d\t%10.2f\t|", $scale[$i][1], $histo[$i], $histo[$i]/$total*100;
    print_bar ( int($histo[$i]/$total*100) );
  }
  printf REPORT "%s\t%10d\t%10.2f\t|", ">$max", $histo_max, $histo_max/$total*100;
  print_bar ( int($histo_max/$total*100) );
}
  
sub print_bar {
  my ( $cnt ) = @_;
  for ( my $i=0; $i<$cnt; $i++ ) { print REPORT "*"; }
  print REPORT "\n";
}
  
# Legacy Files:
#  if ( $database eq "hg18" || $database eq "hg19" || $database eq "mm9" ) {
#    $chromSize_f = "/home/comp/twlab/twang/twlab-shared/genomes/${database}/${database}_chrom_sizes";
#    $MRE_f       = "/home/comp/twlab/twang/twlab-shared/genomes/${database}/MRE/FiveMRE_frags.bed";
#  } elsif ( $database eq "pig" || $database eq "dog" ) {
#    $chromSize_f = "/home/comp/twlab/mxie/${database}/${database}_chrom_sizes";
#    $MRE_f       = "/home/comp/twlab/mxie/${database}/MRE/FiveMRE_frags.bed";
#  } elsif ( $database eq "rat" ) {
#    $chromSize_f = "/home/comp/twlab/mxie/Rat/seq/Rat_chrom_sizes";
#    $MRE_f       = "/home/comp/twlab/mxie/Rat/MRE/FiveMRE_frags.bed";
#  } elsif ( $database eq "zebrafish" ) {
#    $chromSize_f = "/home/comp/twlab/mxie/Zebrafish/danRer7/danRer7_chrom_sizes";
#    $MRE_f       = "/home/comp/twlab/mxie/Zebrafish/danRer7/MRE/FiveMRE_frags.bed";
#  }

