#!/usr/bin/perl

# input: cDNA sequence in fasta format
# output: tab format for each sequence, with ORF position
# usage: cat xxx.fa | perl getORF.pl > xxx.fa.orf.tab

use Bio::SeqIO;
my $in=Bio::SeqIO->new(-fh => \*STDIN, -format=>'fasta');

my $minORFLength=12;

while($seq = $in->next_seq())
{
    my $s=$seq->seq();

    # positions of all start/stop codons (0-based)
    my @STARTcodon; my @STOPcodon;
    while ($s=~/ATG/gi) {push @STARTcodon, pos($s)-3;}
    while ($s=~/TA[AG]|TGA/gi) {push @STOPcodon, pos($s);}

    # sort
    @STARTcodon = sort {$a<=>$b} @STARTcodon;
    @STOPcodon = sort {$a<=>$b} @STOPcodon;

    $found=0;

    # defining ORF for paired start/stop codons
    my @ORFs;
    my $maxlength=0; my $maxlength_start=0;
    foreach $start (@STARTcodon)
    {
        foreach $stop (@STOPcodon)
        {
            if(($stop-$start)>=$minORFLength && ($stop-$start)%3==0){
                #print $seq->id, "\t", $start, "\t", $stop-$start, "\n";
                #if(($stop-$start)>$maxlength) {$maxlength = ($stop-$start); $maxlength_start=$start;}
                push @ORFs, [$start, $stop-$start];
                $found=1;
                last;
            }
        }
    }

    if($found==0){
        print $seq->id, "\t", -1, "\t", -1, "\n";
    }
    else{
        # sort by length
        @ORFs=sort {$b->[1] <=> $a->[1]} @ORFs;
        for (@ORFs){
            print $seq->id, "\t", $_->[0], "\t", $_->[1], "\n";
        }
    }
}
