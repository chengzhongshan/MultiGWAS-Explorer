#!/usr/bin/env perl
use strict;
use warnings;
use Getopt::Long qw(GetOptions);
use File::Spec;
use Digest::MD5 qw(md5_hex);
use JSON::PP qw(decode_json);

my %opt = (mode => 'heatmap', population => 'EUR', window_bp => '1e6');
GetOptions(
    'targets-csv=s'   => \$opt{targets_csv},
    'runner-config=s' => \$opt{runner_config},
    'output-html=s'   => \$opt{output_html},
    'single-runner=s' => \$opt{single_runner},
    'window-bp=s'     => \$opt{window_bp},
    'population=s'    => \$opt{population},
    'mode=s'          => \$opt{mode},
) or die "Invalid signed-LD dispatch option\n";
for my $required (qw(targets_csv runner_config output_html single_runner)) {
    die "Missing --$required\n" unless length($opt{$required} // '');
}
die "Signed-LD mode must be heatmap or both\n"
    unless $opt{mode} =~ /^(?:heatmap|both)$/i;
die "Missing selected-hit CSV: $opt{targets_csv}\n" unless -s $opt{targets_csv};
die "Missing SAS runner: $opt{single_runner}\n" unless -f $opt{single_runner};

open my $fh, '<', $opt{targets_csv} or die "Open $opt{targets_csv}: $!\n";
my $header = <$fh> // die "Empty selected-hit CSV: $opt{targets_csv}\n";
chomp $header;
$header =~ s/\r$//;
my @header = split /,/, $header, -1;
my %idx = map { $header[$_] => $_ } 0 .. $#header;
die "Selected-hit CSV needs SNP, CHR, and BP columns\n"
    unless exists $idx{SNP} && exists $idx{CHR} && exists $idx{BP};
my @targets;
my %seen;
while (my $line = <$fh>) {
    chomp $line;
    $line =~ s/\r$//;
    next unless length $line;
    my @fields = split /,/, $line, -1;
    my $snp = $fields[$idx{SNP}] // '';
    $snp =~ s/^"|"$//g;
    die "Unsafe selected SNP: $snp\n" unless $snp =~ /^[A-Za-z0-9_.:-]+$/;
    next if $seen{lc $snp}++;
    push @targets, $snp;
}
close $fh;
die "Selected-hit CSV contains no SNPs: $opt{targets_csv}\n" unless @targets;
print "[signed-LD] Rendering " . scalar(@targets)
    . " listed loci separately; the list is not proof of LD-independent leads.\n";

my ($volume, $dir, $base) = File::Spec->splitpath($opt{output_html});
die "Output HTML must end in .html: $opt{output_html}\n" unless $base =~ /\.html$/i;
my $output_stem = $base;
$output_stem =~ s/\.html$//i;
open my $config_fh, '<', $opt{runner_config}
    or die "Cannot read $opt{runner_config}: $!\n";
local $/;
my $runner_config = decode_json(<$config_fh>);
close $config_fh;
die "Runner configuration must be a JSON object\n"
    unless ref($runner_config) eq 'HASH';
my %render_config = map { $_ => $runner_config->{$_} }
    grep { /^(?:GTF_|LOCAL_|DISPLAY_GWAS|GROUP_TRACKS$|PAIR_DEFS$|DATA_GZ$|SOURCE_LONG_GZ$|SOURCE_MODE$|EXTRACTOR_CONFIG_JSON$|REFERENCE_BUILD$)/ }
    keys %$runner_config;
my $render_request = JSON::PP->new->canonical(1)->encode({
    version => 1, config => \%render_config,
    data => file_signature($runner_config->{DATA_GZ}),
    source_long => file_signature($runner_config->{SOURCE_LONG_GZ}),
    single_runner => file_signature($opt{single_runner}),
    dispatcher => file_signature(__FILE__),
    window_bp => $opt{window_bp}, population => uc($opt{population}),
    mode => lc($opt{mode}),
});
my @links;
for my $snp (@targets) {
    (my $safe_snp = $snp) =~ s/[^A-Za-z0-9._-]/_/g;
    my $target_html = @targets == 1 ? $base : "${output_stem}_${safe_snp}.html";
    my $target_csv = "${output_stem}_${safe_snp}_top_hit.csv";
    my $target_path = File::Spec->catpath($volume, $dir, $target_html);
    (my $target_png = $target_path) =~ s/\.html$/.png/i;
    my $request_file = "$target_path.request.md5";
    my $request_key = md5_hex(join("\0", $render_request, $snp));
    if (-s $target_path && -s $target_png && -s $request_file) {
        open my $request_fh, '<', $request_file or die "Cannot read $request_file: $!\n";
        my $saved = <$request_fh> // '';
        close $request_fh;
        $saved =~ s/\s+\z//;
        if ($saved eq $request_key) {
            print "[signed-LD] Reusing completed plot for $snp\n";
            push @links, [$snp, $target_html];
            next;
        }
    }
    local %ENV = %ENV;
    delete @ENV{qw(DATA_GZ REMOTE_DATA_BASENAME)};
    $ENV{RUNNER_CONFIG_JSON} = $opt{runner_config};
    $ENV{TARGET_SNP} = $snp;
    $ENV{LOCAL_WINDOW_BP} = $opt{window_bp};
    $ENV{OUTPUT_HTML_BASENAME} = $target_html;
    $ENV{SINGLE_SNP_ALLOW_GENERIC_OUTPUT_BASENAME} = 1;
    $ENV{SINGLE_SNP_TOP_HITS_CSV_BASENAME} = $target_csv;
    $ENV{GTF_LABEL_SNPS} = $snp;
    $ENV{GTF_LD_DISPLAY_MODE} = lc $opt{mode};
    $ENV{GTF_LD_R2_CACHE} = '';
    $ENV{GTF_LD_REFERENCE_SNP} = $snp;
    $ENV{GTF_LD_HEATMAP_LEGEND_TITLE} =
        "Signed LD r2 to $snp (" . uc($opt{population}) . ", 1000G Phase 3 / PLINK2)";
    $ENV{OPEN_RESULT} = 0 if @targets > 1;
    $ENV{CLEAN_ODA_MACROS} = 0;
    # The parent automation may retain a shared genome-wide upload.  These
    # target-specific inputs must be removed after each successful locus.
    $ENV{KEEP_REMOTE_PLOT_DATA} = 0;
    $ENV{CLEAN_ODA_INPUT} = 1;
    print "[signed-LD] Rendering $snp with its own PLINK2 LD reference\n";
    system('/bin/bash', $opt{single_runner}) == 0
        or die "Signed-LD SAS runner failed for $snp: $?\n";
    die "Signed-LD SAS runner did not create $target_path\n" unless -s $target_path;
    die "Signed-LD SAS runner did not create $target_png\n" unless -s $target_png;
    open my $request_fh, '>', $request_file or die "Cannot write $request_file: $!\n";
    print {$request_fh} "$request_key\n";
    close $request_fh or die "Cannot close $request_file: $!\n";
    push @links, [$snp, $target_html];
}

sub file_signature {
    my ($path) = @_;
    return { path => ($path // ''), size => 0, mtime => 0 }
        unless defined($path) && length($path) && -f $path;
    my @stat = stat($path);
    return { path => $path, size => $stat[7], mtime => $stat[9] };
}

if (@links > 1) {
    open my $out, '>', $opt{output_html} or die "Write $opt{output_html}: $!\n";
    print {$out} '<!DOCTYPE html><html lang="en"><head><meta charset="utf-8"><title>Signed-LD local GTF plots</title></head><body>', "\n";
    print {$out} '<h1>Signed-LD local GTF plots</h1><ul>', "\n";
    print {$out} '<li><a href="', $_->[1], '">', $_->[0], '</a></li>', "\n" for @links;
    print {$out} '</ul></body></html>', "\n";
    close $out or die "Close $opt{output_html}: $!\n";
    print "[signed-LD] Wrote plot index: $opt{output_html}\n";
}
