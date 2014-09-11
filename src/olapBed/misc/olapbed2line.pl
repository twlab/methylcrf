#!/usr/bin/perl -n


# olapBed returns one line per matching qry line.  This script will
#  make one line for each matching db line, ie, qry<tab>dbline
# USAGE: olapBed .....   | olapbed2line.pl
#        -d delimeter between types of data olapBed adds, [dflt: <tab>]

INIT {
  $dlm="\t";
  if ($ARGV[0] eq '-d') {shift; $dlm = shift;}
}

chomp;
my ($dbline,@F) = split $dlm;

my @v; for (@F) { push @v, [split ",", $_] }

for my $mdx (0..$#{$v[0]} ) {
 print join($dlm, $dbline,( map{$v[$_][$mdx]} (0..$#v))),"\n";
}


