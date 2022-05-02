#! /usr/bin/perl -w
# filter variants whose ref alleles overlap in VCF
# ============================================================

@preLine = ();
$preChr = "";
$preEnd = 0;
$preStart = 0;

use List::Util qw(max);

while (<>) {
  
  chomp $_;
  
  if ($_ =~ m/^\#/) {
    print $_, "\n"; next;
  }
  
  @line = split /\t/, $_;
  $thisChr = $line[0];
  $thisStart = $line[1] - 1;
  $thisEnd = $thisStart + length($line[3]);
  
  if (($thisChr eq $preChr && $thisStart >= $preEnd) || $thisChr ne $preChr) {
    
    if ($#preLine == 0) {
      @line = split /\t/, $preLine[0];
      $line[2] = $line[0]."_".$line[1]."_".$line[3];
      print join("\t", @line), "\n";
    } else {
      print STDERR join("\n", @preLine), "\n";
    }
    
    @preLine = ($_);
    $preChr = $thisChr;
    $preStart = $thisStart;
    $preEnd = $thisEnd;
    
  } elsif ($thisChr eq $preChr && $thisStart < $preEnd) {
    
    push(@preLine, $_);
    $preChr = $thisChr;
    $preEnd = max(($preEnd, $thisEnd));
    
  }
  
}

if ($#preLine == 0) {
  @line = split /\t/, $preLine[0];
  $line[2] = $line[0]."_".$line[1]."_".$line[3];
  print join("\t", @line), "\n";
}
