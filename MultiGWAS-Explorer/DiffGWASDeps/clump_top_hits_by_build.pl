#!/usr/bin/env perl
use strict;
use warnings;
use Getopt::Long qw(GetOptions);
use Text::CSV;
use File::Basename qw(dirname);
use File::Path qw(make_path);
use File::Spec;
use File::Temp qw(tempdir);

my %opt = (build => 'hg38', populations => 'EUR EAS', population_rule => 'ANY',
    r2_threshold => 0.1, window_kb => 1000, max_leads => 0);
GetOptions(
    'candidates=s' => \$opt{candidates}, 'build=s' => \$opt{build},
    'panel-dir=s' => \$opt{panel_dir}, 'pfile=s' => \$opt{pfile},
    'plink2=s' => \$opt{plink2}, 'populations=s' => \$opt{populations},
    'population-rule=s' => \$opt{population_rule},
    'r2-threshold=f' => \$opt{r2_threshold}, 'window-kb=f' => \$opt{window_kb},
    'max-leads=i' => \$opt{max_leads},
    'output-leads=s' => \$opt{output_leads},
    'output-audit=s' => \$opt{output_audit},
    'output-cache=s' => \$opt{output_cache},
    'output-status=s' => \$opt{output_status},
) or die "Invalid local top-hit clumping option\n";
for my $key (qw(candidates plink2 output_leads output_audit output_cache output_status)) {
    die "Missing --$key\n" unless length($opt{$key} // '');
}
die "--build must be hg19 or hg38\n" unless $opt{build} =~ /^(?:hg19|hg38)$/;
die "--panel-dir is required for hg38\n"
    if $opt{build} eq 'hg38' && !length($opt{panel_dir} // '');
die "--pfile is required for hg19\n"
    if $opt{build} eq 'hg19' && !length($opt{pfile} // '');
die "Invalid r2, window, or max-leads\n"
    unless $opt{r2_threshold} >= 0 && $opt{r2_threshold} <= 1
        && $opt{window_kb} > 0 && $opt{max_leads} >= 0;
make_path(dirname($opt{output_leads})) unless -d dirname($opt{output_leads});
my $deps = dirname(__FILE__);
my $clumper = File::Spec->catfile($deps, 'clump_top_hits_with_plink2_1kg.pl');
my @common = ('--plink2', $opt{plink2}, '--populations', $opt{populations},
    '--population-rule', $opt{population_rule}, '--signal-column', 'focus_signal',
    '--signal-threshold', 1, '--r2-threshold', $opt{r2_threshold},
    '--window-kb', $opt{window_kb});

if ($opt{build} eq 'hg19') {
    run($^X, $clumper, '--candidates', $opt{candidates},
        '--pfile', $opt{pfile}, '--max-leads', $opt{max_leads}, @common,
        '--output-leads', $opt{output_leads}, '--output-audit', $opt{output_audit},
        '--output-cache', $opt{output_cache}, '--output-status', $opt{output_status});
    exit 0;
}

my $reader = Text::CSV->new({ binary => 1, auto_diag => 2 });
open my $input, '<:raw', $opt{candidates} or die "Open $opt{candidates}: $!\n";
my $header = $reader->getline($input) or die "Empty candidate CSV\n";
my %idx = map { uc($header->[$_]) => $_ } 0 .. $#$header;
for my $column (qw(SNP CHR BP FOCUS_SIGNAL)) {
    die "Candidate CSV lacks $column\n" unless exists $idx{$column};
}
my (%by_chr, @candidate_rows);
while (my $row = $reader->getline($input)) {
    my $chr = $row->[$idx{CHR}] // '';
    $chr =~ s/^chr//i;
    $chr = 'X' if $chr eq '23';
    die "No hg38 1KG chromosome panel for candidate chromosome $chr\n"
        unless $chr eq 'X' || ($chr =~ /^\d+$/ && $chr >= 1 && $chr <= 22);
    push @{ $by_chr{$chr} }, $row;
    push @candidate_rows, $row;
}
close $input;
die "No top-hit candidates\n" unless @candidate_rows;
my @ranked = sort {
    ($a->[$idx{FOCUS_SIGNAL}] + 0) <=> ($b->[$idx{FOCUS_SIGNAL}] + 0)
        || ($a->[$idx{CHR}] + 0) <=> ($b->[$idx{CHR}] + 0)
        || ($a->[$idx{BP}] + 0) <=> ($b->[$idx{BP}] + 0)
} @candidate_rows;
my %candidate_rank = map { uc($ranked[$_][$idx{SNP}]) => $_ + 1 } 0 .. $#ranked;
my $tmp = tempdir('top_hit_ld_hg38.XXXXXX',
    DIR => dirname($opt{output_leads}), CLEANUP => 1);
my (@lead_rows, $lead_header, @audit_rows, @cache_rows, @status_rows);
my ($audit_header, $cache_header, $status_header);
for my $chr (sort { $a eq 'X' ? 1 : $b eq 'X' ? -1 : $a <=> $b } keys %by_chr) {
    my $prefix = File::Spec->catfile($opt{panel_dir}, "chr${chr}_hg38");
    run($^X, File::Spec->catfile($deps, 'prepare_plink2_1kg_hg38_chr.pl'),
        '--chr', $chr, '--output-dir', $opt{panel_dir}, '--plink2', $opt{plink2});
    my $part = File::Spec->catfile($tmp, "chr$chr");
    my $part_candidates = "$part.candidates.csv";
    open my $out, '>:raw', $part_candidates or die "Write $part_candidates: $!\n";
    my $writer = Text::CSV->new({ binary => 1, eol => "\n" });
    $writer->print($out, $header);
    $writer->print($out, $_) for @{ $by_chr{$chr} };
    close $out or die "Close $part_candidates: $!\n";
    print "[local LD] Clumping " . scalar(@{ $by_chr{$chr} })
        . " candidate(s) on hg38 chr$chr with EUR/EAS 1KG genotypes.\n";
    run($^X, $clumper, '--candidates', $part_candidates,
        '--pfile', $prefix, '--max-leads', 0, @common,
        '--output-leads', "$part.leads.csv", '--output-audit', "$part.audit.tsv",
        '--output-cache', "$part.cache.tsv", '--output-status', "$part.status.tsv");
    my ($h, $rows) = read_csv("$part.leads.csv");
    $lead_header ||= $h;
    die "Inconsistent LD lead columns across chromosomes\n"
        unless join("\t", @$h) eq join("\t", @$lead_header);
    push @lead_rows, @$rows;
    append_tsv("$part.audit.tsv", \$audit_header, \@audit_rows);
    append_tsv("$part.cache.tsv", \$cache_header, \@cache_rows);
    append_tsv("$part.status.tsv", \$status_header, \@status_rows);
}
my %lead_idx = map { uc($lead_header->[$_]) => $_ } 0 .. $#$lead_header;
@lead_rows = sort {
    ($a->[$lead_idx{FOCUS_SIGNAL}] + 0) <=> ($b->[$lead_idx{FOCUS_SIGNAL}] + 0)
        || ($a->[$lead_idx{CHR}] + 0) <=> ($b->[$lead_idx{CHR}] + 0)
        || ($a->[$lead_idx{BP}] + 0) <=> ($b->[$lead_idx{BP}] + 0)
} @lead_rows;
splice @lead_rows, $opt{max_leads} if $opt{max_leads} && @lead_rows > $opt{max_leads};
my %lead_rank = map { uc($lead_rows[$_][$lead_idx{SNP}]) => $_ + 1 } 0 .. $#lead_rows;
for my $row (@lead_rows) {
    $row->[$lead_idx{LD_LEAD_RANK}] = $lead_rank{uc($row->[$lead_idx{SNP}])};
}
write_csv($opt{output_leads}, $lead_header, \@lead_rows);
my @audit_header_fields = split /\t/, $audit_header, -1;
my %audit_idx = map { $audit_header_fields[$_] => $_ } 0 .. $#audit_header_fields;
for my $line (@audit_rows) {
    my @fields = split /\t/, $line, -1;
    my $candidate = uc($fields[$audit_idx{candidate_snp}] // '');
    my $lead = uc($fields[$audit_idx{lead_snp}] // '');
    $fields[$audit_idx{candidate_rank}] = $candidate_rank{$candidate} // '';
    if (length($lead) && !$lead_rank{$lead}) {
        $fields[$audit_idx{selection_action}] = 'NOT_REACHED';
        $fields[$audit_idx{lead_snp}] = '';
        $fields[$audit_idx{lead_rank}] = '';
    } elsif (length $lead) {
        $fields[$audit_idx{lead_rank}] = $lead_rank{$lead};
    }
    $line = join("\t", @fields);
}
write_tsv($opt{output_audit}, $audit_header, \@audit_rows);
write_tsv($opt{output_cache}, $cache_header, \@cache_rows);
write_tsv($opt{output_status}, $status_header, \@status_rows);
print "[local LD] Selected " . scalar(@lead_rows)
    . " lead(s) using build-matched PLINK2/1KG hg38 panels.\n";

sub read_csv {
    my ($path) = @_;
    open my $fh, '<:raw', $path or die "Read $path: $!\n";
    my $csv = Text::CSV->new({ binary => 1, auto_diag => 2 });
    my $h = $csv->getline($fh) or die "Empty CSV $path\n";
    my @rows;
    push @rows, $_ while $_ = $csv->getline($fh);
    close $fh;
    return ($h, \@rows);
}
sub append_tsv {
    my ($path, $header_ref, $rows) = @_;
    open my $fh, '<:raw', $path or die "Read $path: $!\n";
    my $h = <$fh> // die "Empty TSV $path\n";
    chomp $h;
    $h =~ s/\r$//;
    $$header_ref ||= $h;
    die "Inconsistent LD TSV columns across chromosomes\n" unless $$header_ref eq $h;
    while (my $line = <$fh>) {
        chomp $line;
        $line =~ s/\r$//;
        push @$rows, $line if length $line;
    }
    close $fh;
}
sub write_csv {
    my ($path, $h, $rows) = @_;
    open my $fh, '>:raw', $path or die "Write $path: $!\n";
    my $csv = Text::CSV->new({ binary => 1, eol => "\n" });
    $csv->print($fh, $h);
    $csv->print($fh, $_) for @$rows;
    close $fh or die "Close $path: $!\n";
}
sub write_tsv {
    my ($path, $h, $rows) = @_;
    open my $fh, '>:raw', $path or die "Write $path: $!\n";
    print {$fh} "$h\n", map { "$_\n" } @$rows;
    close $fh or die "Close $path: $!\n";
}
sub run {
    my (@cmd) = @_;
    system(@cmd) == 0 or die "Command failed (exit $?): @cmd\n";
}
