#!/usr/bin/env perl
use strict;
use warnings;
use Getopt::Long qw(GetOptions);
use File::Spec;
use Fcntl qw(:flock);
use Digest::MD5 qw(md5_hex);
use Digest::SHA qw(sha256_hex);
use JSON::PP qw(decode_json);

my %opt = (mode => 'heatmap', population => 'EUR', window_bp => '1e6');
GetOptions(
    'targets-csv=s'   => \$opt{targets_csv},
    'target-snps=s'   => \$opt{target_snps},
    'runner-config=s' => \$opt{runner_config},
    'output-html=s'   => \$opt{output_html},
    'single-runner=s' => \$opt{single_runner},
    'window-bp=s'     => \$opt{window_bp},
    'population=s'    => \$opt{population},
    'mode=s'          => \$opt{mode},
) or die "Invalid signed-LD dispatch option\n";
for my $required (qw(runner_config output_html single_runner)) {
    die "Missing --$required\n" unless length($opt{$required} // '');
}
die "Specify either --targets-csv or --target-snps\n"
    if !!length($opt{targets_csv} // '') == !!length($opt{target_snps} // '');
die "Signed-LD mode must be heatmap or both\n"
    unless $opt{mode} =~ /^(?:heatmap|both)$/i;
die "Missing selected-hit CSV: $opt{targets_csv}\n"
    if length($opt{targets_csv} // '') && !-s $opt{targets_csv};
die "Missing SAS runner: $opt{single_runner}\n" unless -f $opt{single_runner};

my @targets;
my %seen;
if (length($opt{targets_csv} // '')) {
    open my $fh, '<', $opt{targets_csv} or die "Open $opt{targets_csv}: $!\n";
    my $header = <$fh> // die "Empty selected-hit CSV: $opt{targets_csv}\n";
    chomp $header;
    $header =~ s/\r$//;
    my @header = split /,/, $header, -1;
    my %idx = map { $header[$_] => $_ } 0 .. $#header;
    die "Selected-hit CSV needs SNP, CHR, and BP columns\n"
        unless exists $idx{SNP} && exists $idx{CHR} && exists $idx{BP};
    while (my $line = <$fh>) {
        chomp $line;
        $line =~ s/\r$//;
        next unless length $line;
        my @fields = split /,/, $line, -1;
        my $snp = $fields[$idx{SNP}] // '';
        $snp =~ s/^"|"$//g;
        add_target($snp);
    }
    close $fh;
} else {
    add_target($_) for split /,/, $opt{target_snps};
}
sub add_target {
    my ($snp) = @_;
    $snp =~ s/^\s+|\s+$//g;
    die "Unsafe selected SNP: $snp\n" unless $snp =~ /^[A-Za-z0-9_.:-]+$/;
    return if $seen{lc $snp}++;
    push @targets, $snp;
}
die "No selected SNPs\n" unless @targets;
print "[signed-LD] Rendering " . scalar(@targets)
    . " listed loci separately; the list is not proof of LD-independent leads.\n";

my ($volume, $dir, $base) = File::Spec->splitpath($opt{output_html});
die "Output HTML must end in .html: $opt{output_html}\n" unless $base =~ /\.html$/i;
my $output_stem = $base;
$output_stem =~ s/\.html$//i;
my $lock_path = File::Spec->catpath($volume, $dir, "${output_stem}.dispatch.lock");
open my $lock_fh, '>>', $lock_path or die "Open $lock_path: $!\n";
flock($lock_fh, LOCK_EX | LOCK_NB)
    or die "Another signed-LD dispatcher is already using $opt{output_html}; wait for it to finish before resuming\n";
open my $config_fh, '<', $opt{runner_config}
    or die "Cannot read $opt{runner_config}: $!\n";
local $/;
my $runner_config = decode_json(<$config_fh>);
close $config_fh;
die "Runner configuration must be a JSON object\n"
    unless ref($runner_config) eq 'HASH';
my %render_config = map { $_ => $runner_config->{$_} }
    grep { $_ ne 'GTF_LD_R2_CACHE_BY_SNP'
        && /^(?:GTF_|LOCAL_|DISPLAY_GWAS|GROUP_TRACKS$|PAIR_DEFS$|DATA_GZ$|SOURCE_LONG_GZ$|SOURCE_MODE$|EXTRACTOR_CONFIG_JSON$|REFERENCE_BUILD$)/ }
    keys %$runner_config;
my $deps_dir = File::Spec->catpath((File::Spec->splitpath(__FILE__))[0,1], '');
my @plot_sources = qw(
    run_sas_oda_single_snp_with_gtf.sas
    SNP_Local_Manhattan_With_GTF.sas
    Lattice_gscatter_over_bed_track.sas
    map_grp_assoc2gene4covidsexgwas.sas
    Multgscatter_with_gene_exons.sas
    adj_grpnum4close_gene_bed_regs.sas
    extract_single_snp_wide_diff_gwas.pl
    gnuplot/prepare_indexed_merged_locus.pl
    extract_gencode_gtf_subset.pl
    resolve_plink2_local_ld.pl
    augment_gwas_with_ld_r2.pl
    generate_sas_wide_import_include.pl
    generate_sas_gtf_import_include.pl
);
my %plot_source_hashes = map {
    my $path = File::Spec->catfile($deps_dir, $_);
    open my $source, '<:raw', $path or die "Cannot read plot source $path: $!\n";
    local $/;
    my $digest = sha256_hex(<$source>);
    close $source;
    $_ => $digest;
} @plot_sources;
my $render_request = JSON::PP->new->canonical(1)->encode({
    version => 2, config => \%render_config,
    data => file_signature($runner_config->{DATA_GZ}),
    source_long => file_signature($runner_config->{SOURCE_LONG_GZ}),
    plot_sources => \%plot_source_hashes,
    window_bp => $opt{window_bp}, population => uc($opt{population}),
    mode => lc($opt{mode}),
});
my @links;
my $shared_macros_ready = 0;
my $progress_path = File::Spec->catpath($volume, $dir, "${output_stem}.progress.json");
my $partial_path = File::Spec->catpath($volume, $dir, "${output_stem}.partial.html");
my %status = map { $_ => 'pending' } @targets;
my $request_id = md5_hex($render_request);
my $failure = '';
my %status_write_warning;
write_progress();
for my $snp (@targets) {
    (my $safe_snp = $snp) =~ s/[^A-Za-z0-9._-]/_/g;
    my $target_html = @targets == 1 ? $base : "${output_stem}_${safe_snp}.html";
    my $target_csv = "${output_stem}_${safe_snp}_top_hit.csv";
    my $target_path = File::Spec->catpath($volume, $dir, $target_html);
    (my $target_png = $target_path) =~ s/\.html$/.png/i;
    my $request_file = "$target_path.request.md5";
    my $ld_cache = $runner_config->{GTF_LD_R2_CACHE_BY_SNP}{lc $snp}
        // $runner_config->{GTF_LD_R2_CACHE_BY_SNP}{$snp}
        // '';
    my $request_key = md5_hex(join("\0", $render_request, $snp,
        JSON::PP->new->canonical(1)->encode(file_signature($ld_cache))));
    if (-s $target_path && -s $target_png && -s $request_file) {
        open my $request_fh, '<', $request_file or die "Cannot read $request_file: $!\n";
        my $saved = <$request_fh> // '';
        close $request_fh;
        $saved =~ s/\s+\z//;
        if ($saved eq $request_key) {
            print "[signed-LD] Reusing completed plot for $snp\n";
            push @links, [$snp, $target_html];
            $status{$snp} = 'complete';
            write_progress();
            next;
        }
    }
    # A failed rerender must never leave an old matching sidecar behind.
    unlink $request_file if -e $request_file;
    $status{$snp} = 'running';
    write_progress();
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
    $ENV{GTF_LD_R2_CACHE} = $ld_cache;
    $ENV{GTF_LD_REFERENCE_SNP} = $snp;
    $ENV{GTF_LD_HEATMAP_LEGEND_TITLE} =
        "Signed LD r2 to $snp (" . uc($opt{population}) . ", 1000G Phase 3 / PLINK2)";
    $ENV{OPEN_RESULT} = 0 if @targets > 1;
    $ENV{CLEAN_ODA_MACROS} = 0;
    # The first successful locus uploads and checks the shared macro files.
    # Later loci keep using those remote files while uploading only their
    # distinct GWAS/GTF inputs. Cached plots do not establish remote readiness.
    $ENV{SINGLE_SNP_REUSE_SHARED_MACROS} = $shared_macros_ready ? 1 : 0;
    $ENV{INCLUDE_PREFLIGHT_ENABLED} = 0 if $shared_macros_ready;
    $ENV{SAS_ODA_REUSE_VERIFIED_MACRO_BOOTSTRAP_HELPER} = $shared_macros_ready ? 1 : 0;
    # The parent automation may retain a shared genome-wide upload.  These
    # target-specific inputs must be removed after each successful locus.
    $ENV{KEEP_REMOTE_PLOT_DATA} = 0;
    $ENV{CLEAN_ODA_INPUT} = 1;
    print "[signed-LD] Rendering $snp with its own PLINK2 LD reference\n";
    my $rc = system('/bin/bash', $opt{single_runner});
    if ($rc != 0 || !-s $target_path || !-s $target_png) {
        my $exit_code = $rc == -1 ? 1 : ($rc >> 8 || 1);
        $failure = "Signed-LD SAS runner failed for $snp (exit $exit_code)";
        $failure .= '; missing HTML or PNG' if $rc == 0;
        $status{$snp} = 'failed';
        write_progress();
        print STDERR "[signed-LD] $failure. Re-run the same pipeline command to resume remaining loci. Progress: $progress_path\n";
        exit $exit_code;
    }
    $shared_macros_ready = 1;
    atomic_write($request_file, "$request_key\n");
    push @links, [$snp, $target_html];
    $status{$snp} = 'complete';
    write_progress();
}

sub file_signature {
    my ($path) = @_;
    my $stat_path = $path // '';
    if ($^O eq 'cygwin' && $stat_path =~ m{^([A-Za-z]):[/\\](.*)$}) {
        $stat_path = '/cygdrive/' . lc($1) . '/' . $2;
        $stat_path =~ s{\\}{/}g;
    }
    return { path => ($path // ''), size => 0, mtime => 0 }
        unless defined($path) && length($path) && -f $stat_path;
    my @stat = stat($stat_path);
    return { path => $path, size => $stat[7], mtime => $stat[9] };
}

if (@links > 1) {
    atomic_write($opt{output_html}, plot_index_html(\@links));
    print "[signed-LD] Wrote plot index: $opt{output_html}\n";
}
unlink $partial_path if -e $partial_path;
write_progress();

sub write_progress {
    my $complete = scalar grep { $status{$_} eq 'complete' } @targets;
    my $report = {
        version => 1, request_id => $request_id, output_html => $opt{output_html},
        total => scalar(@targets), complete => $complete,
        remaining => scalar(@targets) - $complete,
        failure => $failure,
        loci => [ map { { snp => $_, status => $status{$_} } } @targets ],
    };
    write_status_file($progress_path, JSON::PP->new->canonical(1)->pretty(1)->encode($report));
    if ($complete < @targets && @links) {
        write_status_file($partial_path, plot_index_html(\@links));
    }
    # The full index is a completion signal; keep it only when every locus
    # belongs to the current request and has verified output files.
    unlink $opt{output_html} if @targets > 1 && $complete < @targets && -e $opt{output_html};
}

sub write_status_file {
    my ($path, $contents) = @_;
    if (eval { atomic_write($path, $contents); 1 }) {
        delete $status_write_warning{$path};
        return 1;
    }
    my $error = $@ || 'unknown write error';
    warn "[signed-LD] Could not refresh optional status file $path: $error"
        unless $status_write_warning{$path}++;
    return 0;
}

sub plot_index_html {
    my ($links) = @_;
    my $html = '<!DOCTYPE html><html lang="en"><head><meta charset="utf-8"><title>Signed-LD local GTF plots</title></head><body>' . "\n";
    $html .= '<h1>Signed-LD local GTF plots</h1><ul>' . "\n";
    $html .= '<li><a href="' . $_->[1] . '">' . $_->[0] . '</a></li>' . "\n" for @$links;
    return $html . '</ul></body></html>' . "\n";
}

sub atomic_write {
    my ($path, $contents) = @_;
    my $tmp = "$path.tmp.$$";
    open my $out, '>', $tmp or die "Write $tmp: $!\n";
    print {$out} $contents or die "Write $tmp: $!\n";
    close $out or die "Close $tmp: $!\n";
    my $error = '';
    for my $attempt (1 .. 15) {
        return if rename $tmp, $path;
        $error = "$!";
        last unless $^O =~ /^(?:cygwin|MSWin32)$/i
            && ($!{EBUSY} || $!{EACCES} || $!{EPERM});
        select undef, undef, undef, 0.2 if $attempt < 15;
    }
    unlink $tmp or warn "Cannot remove failed status write $tmp: $!\n";
    die "Rename $tmp to $path: $error\n";
}
