#! /usr/bin/perl -w
# find allele information
# 1. whether to flip
# 2. whether to replace reference
# ============================================================

%flip = ("AC" => "GT",
         "AG" => "CT",
         "CT" => "AG",
         "GT" => "AC");

while (<>) {
  
  chomp $_;
  @line = split /\t/, $_;
  
  # sort and join alleles
  $chipAlleles = join("", sort(@line[2..3]));
  $genomeAlleles = join("", sort(@line[7..8]));
  
  if (!defined($flip{$chipAlleles}) || !defined($flip{$genomeAlleles})) { next; }
  
  if ($chipAlleles eq $genomeAlleles) {
    print join("\t", @line), "\tpass\n";
  } elsif ($flip{$chipAlleles} eq $genomeAlleles) {
    print join("\t", @line), "\tflip\n";
  } else {
    next; # not flippable
  }
  
}
