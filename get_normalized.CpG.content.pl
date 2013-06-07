#!/usr/bin/perl

# input: chr/start/end/strand like GTF file for gencode gene (e.g. download from http://www.gencodegenes.org/releases/3c.html)
# output: tab format for each transcript, with CpG content
# usage: awk 'BEGIN{OFS="\t"}{if($3=="transcript")print $1,$4,$5,$7,$10,$12,$14,$16,$18,$20,$22,$24,$26}' ../data/gencode.v3c.annotation.GRCh37.gtf | sed 's/;//g;s/"//g;' | perl get_normalized.CpG.content.pl > gencode_v3c_hg19_transcripts.gtf.cpg.tab

use lib "/home/dongx/src/ensembl67/ensembl/";
use Bio::EnsEMBL::Registry;

my $registry = 'Bio::EnsEMBL::Registry';

$registry->load_registry_from_db(
    -host => 'ensembldb.ensembl.org',
    -user => 'anonymous',
    -version => 67  # change into the correponding version of ENSEMBL for each synchronized Gencode. See: http://genome.ucsc.edu/cgi-bin/hgTrackUi?hgsid=186973795&c=chr1&g=wgEncodeGencode
);
my $slice_adaptor = $registry->get_adaptor( 'Mouse', 'Core', 'Slice' );

while($line = <STDIN>){
    chomp($line);
    @t=split("\t", $line);
    # chr	start	end	strand gene_id   transcript_id gene_type   gene_status gene_name   transcript_type transcript_status   transcript_name level
    $chr=$t[0];
    if($chr eq 'chr'){
        print join("\t", @t, "normalized.CpG.content"), "\n";
        next;
    }
    $tss = ($t[3] eq "-")?$t[2]:$t[1];
    $start=$tss-1500;
    $end=$tss+1500;

    $chr =~ s/^chr//g;
    $chr="Mt" if($chr eq 'M');

    my $s = $slice_adaptor->fetch_by_region( 'chromosome', $chr, $start, $end, ($t[3] eq '+')?1:-1)->seq;

    my @n1 = ($s=~/CG/g); my $n1=@n1;
    my $n2 = length($s);
    my @n3 = (($s=~/C/g), ($s=~/G/g)); my $n3=@n3;

    print join("\t", @t, sprintf("%.3f", 4*$n2*$n1/($n3**2))), "\n";
    # Normalized CpG fraction was computed as (observed CpG)/(expected CpG), where expected CpG was calculated as (GC content/2)^2.
}
