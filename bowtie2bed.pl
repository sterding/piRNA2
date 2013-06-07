#!/usr/bin/perl -w
#
########################################
#
# File Name:
#   bowtie2bed.pl
# 
# Description:
#   Convert bowtie concise output to bed.
# 
# Usage:
#   $0 <bowtie-out> <out.bed>
# 
# Author:
#   Xi Wang, wang-xi05@mails.thu.edu.cn
# 
# Date:
#   Wed Oct 14 15:38:30 CST 2009
#
########################################

use strict;
my $usage = "$0 <bowtie-out> <out.bed>\n";
my $infile = shift || die $usage;
my $outfile = shift || die $usage;
open(IN, $infile) || die "Can't open $infile for reading!\n";
open(OUT, ">$outfile") || die "Can't open $outfile for writing!\n";

=head
bowtie output examples:

GA-EAS46_1_209DH:5:330:347:125   -       chrVII  1013392 GGACGAAAAGTTTGAATGGTTCCAAGTGTGGCAAGC    GAEF@AUFP>NIOOQKOBESDOLU`NLQSQO\```N        0       26:A>G
GA-EAS46_1_209DH:5:330:85:6      -       chrVII  591540  AAAGACCAAACTACCTTTTGTGCGGACTTGGTGTCG    FGDAB@AHEDAGAENLFQWODEQRJEQQKMH\H^RI        0
GA-EAS46_1_209DH:5:330:166:14    -       chrV    280087  GAAATACCAATGATCTTACAACCGGCGGCTTTTCCG    DAEIG?GFKH?BFMCGKFHB>DK?B@LMDQZHIP^A        0       34:C>A
GA-EAS46_1_209DH:5:330:106:862   +       chrIII  52072   CGCGATTTCTAGGAGATCTGCGTATTTAGACGCCGT    PVL_>V]XKR@NI>GCWHJKGCEGMLGAF?J?=B>G        0       29:C>A

=cut

my @col;
my ($chr, $start, $end, $name, $strand);

while(<IN>)
{
  chomp;
  @col = split '\t', $_;
  $chr = $col[2];
  my $l = length($col[4]);
  $start = $col[3];
  $end = $start + $l;
  $strand = $col[1];
  if(@col < 8)
  {
    $name = "U0";
  }
  else
  {
    my @tmp = split ',', $col[7];
    $name = "U".@tmp;
  }
  print OUT "$chr\t$start\t$end\t$name\t0\t$strand\n";
}

close IN;
close OUT;
