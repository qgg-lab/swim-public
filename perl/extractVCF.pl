#! /usr/bin/perl -w
# extract variants from VCF according to a variant list
# basically an intersection operation
# ============================================================

($file1, $varList) = @ARGV;

open VAR, "<$varList";
%var = ();

while (<VAR>) {
  
  chomp $_;
  if ($_ =~ /^\#/) {
    if ($_ =~ /^\#\#/) {
      print $_, "\n"; next;
    } else {
      next;
    }
  }
  
  @line = split /\t/, $_;
  $var{$line[2]} = $_;
  
}

close VAR;

open VCF, "<$file1";

while (<VCF>) {
  
  chomp $_;
  if ($_ =~ /^\#/) {
    if ($_ =~ /^\#\#/) {
      next;
    } else {
      print $_, "\n"; next;
    }
  }
  
  @line = split /\t/, $_;
  $id = $line[0]."_".$line[1]."_".$line[3];
  if (defined($var{$id})) {
    print $var{$id}, "\t";
    print join("\t", @line[8..$#line]), "\n";
  }
  
}

close VCF;
