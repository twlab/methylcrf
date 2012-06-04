#!/usr/bin/env perl

use strict;
# NOTES #
# should have used BS midpt for training too, so I could avoid caring here
#########

# hiested from http://www.perlmonks.org/bare/?node_id=178170
# NOTE: IO::FILE allows handling fh's in either OO-mode or procedural
#        --maybe switch to filehandle --it's on the system, so maybe it's std
package IO::File::AutoChomp;
use base 'IO::File';
sub getline  { my $v = readline($_[0]) || die "readline failed: $!";chomp $v;return $v}
sub getlines { my ($self) = @_; return map { chomp; $_ } $self->SUPER::getlines(); }



package main;

my $bin     = "./crfsgd";
my @windist = (0,10,100,1000,10000);
my (
  $mdldir,  # has one model file per crf
  $gdatdir, # has specifies-specific files: 1 cpg.bin for each crf, cpg.bed, gdtat.tbl
  $dipfn,   # MeDIP-seq file
  $mrefn,   # MRE-seq file
  $fragdir, # dir with MRE enzyme frag files -will detect if it's instead a filename and use that
  $gapsz,   # size of gap
  $eid,     # prefix for output files
  $startfrom # control which processes run
) = @ARGV;

my ($crffn,$cutfn,$cpgfn,$gdatfn,$fragfn) = 
   ("$mdldir/crf.list","$mdldir/cut.list","$gdatdir/cpg.bed","$gdatdir/gdata.tbl",$fragdir);

if (-d $fragdir) { $fragfn  = ${eid}."_".qx(basename $fragdir); chomp $fragfn; }

# startfrom:
# 0 all (format DIP/MRE bed files)
# 1 DIP/MRE avgwin (25Gb)
# 2 make sampled MRE fragment file
# 3 make tables  (1.5hr)
# 4 predict
# 5 combine


## 0: Get and Format DIP/MRE data ##

# make correct bed files from handler output (dipfn.bed,mrefn.bed)
if ($startfrom <=0) {
  print STDERR scalar localtime(), " format $dipfn and $mrefn\n";
  format_DIPMRE($dipfn,$mrefn,$cpgfn) 
}

## 1: generate windowed files ({e}_DIP_d{d}.cnt, {e}_MRE_d{d}.cnt)
if ($startfrom <=1) {
  print STDERR scalar localtime(), " make avg windows for $dipfn and $mrefn\n";
  # makes e_DIPMRE.cnt
  gen_avgwindows($eid,"$dipfn.norm.bed","$mrefn.bed",$cpgfn, "0 10 100 1000 10000") ;
}

## 2: make MRE frag file
if ($startfrom <=2) {
  print STDERR scalar localtime(), " make MRE frag in $fragdir\n";
  make_fragfn($fragdir,$cpgfn,$fragfn);
}
 


## stuff needed for steps (3,4,5) ##

# crf list
my @crf = map {chomp;$_} qx(cat $crffn|grep -v "^#"); 

# genomic cpg.bed (step 2,3,4)
my $CPG = IO::File::AutoChomp->new($cpgfn,'r') or die "can't open $cpgfn: $!";

# files have cpgid [0|1] -telling weather CpG is in crf (this is a bad idea, i need to redo it)
# step (2,4)
my @CLS =map { 
    my $fn = "$gdatdir/${_}_cpg.bin";
    IO::File::AutoChomp->new($fn, 'r') or die "can't open $fn: $!";
  } @crf;

# CRF table fn's (step 2,3,4)
my @tblfn = map{"${eid}_${_}.tbl"} (@crf);

# get cut list hash (step 2,4[midpt])
my (%cut,$id);
for (qx(cat $cutfn)) {
  chomp;
  if (/[a-zA-Z]/) { $id=$_                 }
  else            { push @{$cut{$id}},$_+0 }
}



## 2: Make Tables ##



if ($startfrom <=3) {
  print STDERR scalar localtime(), " make table\n";
 # local $| = 1;
  # get common filehandles
  my $EFN  = IO::File::AutoChomp->new("${eid}_DIPMRE.cnt", 'r') or die "can't open DIPMRE.cnt $_: $!";
  my $FRAG = IO::File::AutoChomp->new($fragfn,'r') or die "can't open $fragfn: $!";
  my $GDAT = IO::File::AutoChomp->new($gdatfn,'r') or die "can't open $gdatfn: $!";
  my @gdat_hdr = split /\t/,$GDAT->getline();

  # per CRF fh's (these should close as they go out of scope)
  my @OUT = map { IO::File->new($_, 'w') or die "can't open $_: $!"; } @tblfn ;

  my (@pchr,@plocn);
  while (my @cpg = split /\t/, <$CPG>) {

    my @dipmre = (split /\t/,  $EFN->getline());chomp $dipmre[-1];
    my $frag   = (split /\t/, $FRAG->getline())[1];
    my @gdat   = (split /\t/, $GDAT->getline()); chomp @gdat[-1];

    # each crf
    for my $cidx (0..$#crf) {
      # print this cpg or not?
      my @cls = split /\t/, $CLS[$cidx]->getline();
      if ($cls[0] ne $cpg[3]) {die "cls cpg != cpg ".join(" ", @cpg)."  ".join(" ", @cls)}
      if (!$cls[1]) {next}

      # this crf's fh 
      my $OUT   = $OUT[$cidx];
      my $crfnm = $crf[$cidx];

      # add newline for gaps
      if ( $pchr[$cidx] && ($pchr[$cidx] ne $cpg[0] || ($cpg[1] -$plocn[$cidx] >$gapsz)) ) 
      {$OUT->say("");}

      # MRE and MeDIP cols
      $OUT->print( ($dipmre[$_]           ==-1 ? -1 : binon( $cut{ "${crfnm}__H1ES_DIP_d$windist[$_]" },$dipmre[$_]          )),"\t" ) for (0..$#windist);
      $OUT->print( ($dipmre[ $_+@windist ]==-1 ? -1 : binon( $cut{ "${crfnm}__H1ES_MRE_d$windist[$_]" },$dipmre[$_+@windist] )),"\t" ) for (0..$#windist);

      # mre frag 
      die "no frag at cpg $cpg[3]" unless defined $frag;
      $OUT->print( $frag, "\t");

      # all other cols (no autochomp) 
      #$gdat =~ s/.*?\t//; #(cut off chrx.x)
      #$OUT->print($gdat);
      $OUT->print( ($gdat[$_]==-1 ? -1 : binon( $cut{ "${crfnm}__$gdat_hdr[$_]" },$gdat[$_] )),"\t" ) for (1..$#gdat_hdr-1);
      $OUT->print( ($gdat[-1]==-1 ? -1 : binon( $cut{ "${crfnm}__$gdat_hdr[-1]" },$gdat[-1] )),"\n" );

      ($pchr[$cidx],$plocn[$cidx]) = ($cpg[0],$cpg[2])
    }
  }
  # add empty last line (maybe didn't work..)
    for (@OUT) { $_->say("");}
}


## 3: Predict ##
if ($startfrom <= 4) {
  print STDERR scalar localtime(), " predict\n";
  predict($bin,\@tblfn,\@crf,$mdldir, $eid) 
}
  

## 4: Combine ##
# midpt map
my $midpt = midpt_map(\@crf,\%cut,1);

# reset fh posn
$CPG->seek(0,0);
for (@CLS) { $_->seek(0,0)}

#prints to stdout
print STDERR scalar localtime(), " making ensemble prediction\n";
combine(\@CLS,\@tblfn,$midpt, $CPG);
print STDERR scalar localtime(), " fin\n";






sub format_DIPMRE {
  my ($dfn,$mfn,$cfn) = @_;
  qx(awk '(c != \$1){i=0}{OFS="\t";print \$1,\$2,\$3,"dip."\$1"."(++i),1;c=\$1}' $dfn |sort -Vk1,3 > $dfn.bed) == 0 
    or die "$0: make dip.bed " . ($? >> 8);
  qx(awk '(c != \$1){i=0}{OFS="\t";print \$1,\$2,\$3,"mre."\$1"."(++i),\$4;c=\$1}' $mfn |sort -Vk1,3 > $mfn.bed) == 0 
    or die "$0: make mre.bed " . ($? >> 8);

# do medip norm
  qx(medip_norm.sh $dfn.bed $cfn  >$dfn.norm.bed 2>$dfn.norm.err) == 0 
    or die "$0: medip_norm " . ($? >> 8);

 qx(ls -l $dfn.bed $mfn.bed $dfn.norm.bed >&2);
}

sub gen_avgwindows {
  my ($e,$dfn,$mfn,$cfn,$win) = @_;
  qx(bed2avgwinbin.sh $dfn $cfn "$win" ${e}_DIP 2>${e}_DIP_avgwin.err) == 0
    or die "$0: bed2avgwingin(DIP) " . ($? >> 8);
  qx(bed2avgwinbin.sh $mfn $cfn "$win" ${e}_MRE 2>${e}_MRE_avgwin.err) == 0
    or die "$0: bed2avgwingin(MRE) " . ($? >> 8);

  my $dip = join " ", map{"<\(cut -f2 ${e}_DIP_d${_}.cnt)"} (split /\W/,$win);
  my $mre = join " ", map{"<\(cut -f2 ${e}_MRE_d${_}.cnt)"} (split /\W/,$win);
  my $out = "${e}_DIPMRE.cnt"; 
  system("bash","-c", "paste $dip $mre > $out") == 0 or die "$0: gen_avgwin(paste) " . ($? >> 8);
   qx(ls -l $out >&2);
  # del cnt files??
}

sub make_fragfn {
  my ($dir,$cpg,$outfn) = @_;
  qx(dir2frag.sh $dir $cpg > $outfn) == 0 or die "$0: medip_norm " . ($? >> 8);
  qx(ls -l $outfn >&2;)
}
sub predict {
  my ($bin, $tfn,$crf,$mdldir,$e) = @_;
  for (0..$#$tfn) {
    my $mdl="$mdldir/$crf->[$_].mdl.gz";
    my $out="$tfn->[$_].out";
    print STDERR "predict: [$crf->[$_]] [$mdl]\n";
    qx($bin -t $mdl $tfn->[$_] 2>$out.err |awk '(NF){print \$NF}' >$out) == 0
    or die "$0: crfsgd " . ($? >> 8);

    qx(ls -l $out >&2);
}
}
sub combine {
  my ($CLS,$tblfn,$midmap,$CPG) = @_;

# prediction fh's
  my @pred= map { 
    IO::File::AutoChomp->new("${_}.out", 'r') or die "can't open $_.out: $!";
    } @$tblfn;

    while (my @l = split /\t/, <$CPG>) {
      $#l=4; # remove anything after score
      my ($v,$c);
      for my $cidx (0..$#$tblfn) {
        if ( (split /\t/, $CLS->[$cidx]->getline())[1] ) {
          if (defined (my $p = $pred[$cidx]->getline() )){
              $v += $midmap->[$cidx]{ $p+0 };$c++;
            } else {die "$tblfn->[$cidx] has no val  for $l[3]!"}
#          my $s =  $pred[$cidx]->getline() +0 ;$c++; my $m =$midmap->[$cidx]{$s}; $v+=$m;
            #print "$tblfn->[$cidx] $s ($m) tot: $v\n" if $l[3] eq 'chr1.1388';
          }
        }
      die "@l has no predictions!" if !$c;
      $l[4] = $v/$c;
      print join("\t",@l),"\n";
   #  exit if $l[3] eq 'chr1.1388';
    }
 
# sanity check all pred's are eof 
  for (@$CLS) { die "$_ not finished" if !$_->eof()}
}


sub binon {
  my ($cut,$v, $nobin) = @_;
  if ($v eq $nobin) {return $v};
  #NOTE: v can be < cut[0] if binning was done on small samples
  my $idx=1; while ($idx <= $#$cut && $v >= $cut->[$idx]) {$idx++}
  return  $cut->[$idx-1];
}


sub midpt_map{
  my ($crf, $cut, $max) = @_;
  my @rtn;
  for my $cidx (0..$#$crf) {
    my $k=$crf->[$cidx]."__H1ES_BS";
    for ( 0 .. $#{$cut->{$k}}-1 ) {
      my $v = $cut->{$k}[$_]; 
      $rtn[$cidx]{$v} = ( $v + $cut->{$k}[$_+1] )/2;
    }
    $rtn[$cidx]{ $cut->{$k}[-1] } =  ( $cut->{$k}[-1] + $max )/2;
  }
  return \@rtn;
}

# Not really sure the best way to do this:
# 1) loop through crfs w/ cpg's inside
# 2) loop through cpg w/ crf's inside
sub make_crf {
  my ($CPG,$CLASS,$gapsz,$crfnm,$windist,$cut,$EFN,$FRAG,$GDAT,$gdat_hdr,$OUT) = @_;
  # local $| = 1;

  my $c;
  my ($pchr,$plocn);
  while (my @cpg = split /\t/, <$CPG>) {
    $c++;
    # MRE and MeDIP cols
    #my @mre = map{ ( split /\t/, readline($efn->{'MRE'}[$_]) )[1]+0 } (0..$#windist); #force numerical contxt
    #my @dip = map{ (split /\t/, $efn->{'DIP'}[$_]->getline())[1] } (0..$#windist);
    #my $frag = (split /\t/, $FRAG->getline())[1];

    my $dipmre = $EFN->getline();
    my $frag   = $FRAG->getline();
    my $gdat   = $GDAT->getline();

    # print this cpg or not?
    #if (!(split /\t/,$CLASS->getline())[1]) {next}
    my @cls = (split /\t/,$CLASS->getline());
    if ($cls[0] ne $cpg[3]) {die "cls cpg != cpg ".join(" ", @cpg)."  ".join(" ", @cls)}
    if (!$cls[1]){next}

    # add newline for gaps
    if ( $pchr && ($pchr ne $cpg[0] || ($cpg[1] -$plocn >$gapsz)) ) 
    {$OUT->say("");}

    # MRE and MeDIP cols
    my @dipmre = split /\t/, $dipmre;chomp $dipmre;
    $OUT->print( ($dipmre[$_]           ==-1 ? -1 : binon( $cut->{ "${crfnm}__H1ES_DIP_d$windist[$_]" },$dipmre[$_]          )),"\t" ) for (0..$#windist);
    $OUT->print( ($dipmre[ $_+@windist ]==-1 ? -1 : binon( $cut->{ "${crfnm}__H1ES_MRE_d$windist[$_]" },$dipmre[$_+@windist] )),"\t" ) for (0..$#windist);

    # mre frag 
    $OUT->print( (split /\t/, $frag)[1], "\t");

    # all other cols (no autochomp) 
    #$gdat =~ s/.*?\t//; #(cut off chrx.x)
    #$OUT->print($gdat);
    my @gdat = split /\t/,$gdat;
    $OUT->print( ($gdat[$_]==-1 ? -1 : binon( $cut->{ "${crfnm}__$gdat_hdr->[$_]" },$gdat[$_] )),"\t" ) for (1..$#$gdat_hdr-1);
    $OUT->print( ($gdat[-1]==-1 ? -1 : binon( $cut->{ "${crfnm}__$gdat_hdr->[-1]" },$gdat[-1] )),"\n" );

    ($pchr,$plocn) = ($cpg[0],$cpg[2]);
  }
  # add empty last line
  $OUT->say("");
}





## once through crf's
#if ($startfrom < 3) {
#
#  # get filehandles
## DIP/MRE
#  my $EFN  = IO::File::AutoChomp->new("${eid}_DIPMRE.cnt", 'r') or die "can't open DIPMRE.cnt $_: $!";
#  my $FRAG = IO::File::AutoChomp->new($fragfn,'r') or die "can't open $fragfn: $!";
#  my $GDAT = IO::File::AutoChomp->new($gdatfn,'r') or die "can't open $gdatfn: $!";
#
#  # each crf
#  for my $cidx (0..$#crf) {
#    print STDERR "make $crf[$cidx] table\n"; 
#
#    # reset global fh's
#    #for my $k (keys %efn) { $_->seek(0,0) for @{$efn{$k}}; }
#    $CPG->seek(0,0);
#    $EFN->seek(0,0);  
#    $FRAG->seek(0,0);
#    $GDAT->seek(0,0);my @gdat_hdr = split /\t/,$GDAT->getline();
# 
#    # get crf-specific fh's (these should close as they go out of scope)
#    my $CLASS    = $CLS[$cidx];
#    my $OUT      = IO::File->new($tblfn[$cidx], 'w') or die "can't open $_: $!";
#
#    make_crf($CPG,$CLASS,$gapsz,$crf[$cidx],\@windist,\%cut,$EFN,$FRAG,$GDAT,\@gdat_hdr,$OUT);
#    $OUT->close();
# exit;
#  }
#}

