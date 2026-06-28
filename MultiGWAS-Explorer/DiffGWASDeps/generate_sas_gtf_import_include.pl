#!/usr/bin/env perl
use strict;
use warnings;

use Getopt::Long qw(GetOptions);

my $dataset = 'gtf_hg38';
my $remote_basename = '';

GetOptions(
    'dataset=s'         => \$dataset,
    'remote-basename=s' => \$remote_basename,
) or die usage();

die "--remote-basename is required\n" unless length $remote_basename;

my $text = <<'SAS';
filename gtfdata zip "~/__REMOTE_BASENAME__" gzip;

data __DATASET__;
  infile gtfdata dlm='09'x dsd firstobs=2 truncover lrecl=1048576;
  length chr_text chr_raw seqname source $64 feature type $32 score $32 strand $4 frame $8 ensembl $64;
  length gene gene_name gene_id transcript_name transcript_id gene_type transcript_type exon_id exon_number genesymbol $128;
  input
    chr
    chr_text :$64.
    chr_raw :$64.
    seqname :$64.
    source :$64.
    feature :$32.
    type :$32.
    start
    end
    st
    en
    bp1
    bp2
    txStart
    txEnd
    score :$32.
    strand :$4.
    frame :$8.
    gene :$128.
    gene_name :$128.
    gene_id :$128.
    transcript_name :$128.
    transcript_id :$128.
    gene_type :$128.
    transcript_type :$128.
    exon_id :$128.
    exon_number :$128.
    ensembl :$64.
    genesymbol :$128.
    protein_coding
    original_protein_coding
  ;
run;
SAS

$text =~ s/__REMOTE_BASENAME__/$remote_basename/g;
$text =~ s/__DATASET__/$dataset/g;
print $text;

sub usage {
    return <<"USAGE";
Usage:
  perl generate_sas_gtf_import_include.pl --remote-basename subset.tsv.gz [--dataset gtf_hg38]
USAGE
}
