#!/usr/bin/env perl
use strict;
use warnings;

use Cwd qw(abs_path);
use File::Basename qw(dirname);
use File::Path qw(make_path);
use File::Spec;
use File::Temp qw(tempdir);
use Getopt::Long qw(GetOptions);
use Text::CSV;

my ($candidates, $pfile, $bfile, $plink2, $populations, $output_cache,
    $output_status, $work_prefix, $snp_column, $signal_column);
my $signal_threshold = 1;
my $r2_threshold = 0.1;
my $window_kb = 1000;
my $threads = 4;
my $keep_work = 0;

$plink2 = $ENV{PLINK2} || 'plink2';
$populations = 'EUR,EAS';
$snp_column = 'SNP';
$signal_column = '';

GetOptions(
    'candidates=s'       => \$candidates,
    'pfile=s'            => \$pfile,
    'bfile=s'            => \$bfile,
    'plink2=s'           => \$plink2,
    'populations=s'      => \$populations,
    'snp-column=s'       => \$snp_column,
    'signal-column=s'    => \$signal_column,
    'signal-threshold=f' => \$signal_threshold,
    'r2-threshold=f'     => \$r2_threshold,
    'window-kb=f'        => \$window_kb,
    'threads=i'          => \$threads,
    'output-cache=s'     => \$output_cache,
    'output-status=s'    => \$output_status,
    'work-prefix=s'      => \$work_prefix,
    'keep-work!'         => \$keep_work,
) or die usage();

die usage() unless defined($candidates) && (defined($pfile) xor defined($bfile))
    && defined($output_cache) && defined($output_status);
die "--r2-threshold must be between 0 and 1\n"
    unless $r2_threshold >= 0 && $r2_threshold <= 1;
die "--window-kb must be positive\n" unless $window_kb > 0;
die "--threads must be positive\n" unless $threads > 0;
die "Candidate file does not exist: $candidates\n" unless -s $candidates;

my $input_prefix = defined($pfile) ? $pfile : $bfile;
my @required = defined($pfile)
    ? ("$input_prefix.pgen", (-e "$input_prefix.pvar.zst" ? "$input_prefix.pvar.zst" : "$input_prefix.pvar"), "$input_prefix.psam")
    : ("$input_prefix.bed", "$input_prefix.bim", "$input_prefix.fam");
for my $path (@required) {
    die "PLINK reference component is unavailable: $path\n" unless -s $path;
}

my @pops = grep { length } map {
    my $x = uc($_);
    $x =~ s/^\s+|\s+$//g;
    $x;
} split /[\s,]+/, $populations;
die "At least one population is required\n" unless @pops;

my $tmp_obj;
my $work_dir;
if (defined($work_prefix) && length($work_prefix)) {
    $work_dir = dirname(File::Spec->rel2abs($work_prefix));
    make_path($work_dir) unless -d $work_dir;
} else {
    $tmp_obj = tempdir('plink2_1kg_top_hits.XXXXXX', TMPDIR => 1, CLEANUP => !$keep_work);
    $work_prefix = File::Spec->catfile($tmp_obj, 'top_hits');
    $work_dir = $tmp_obj;
}

my ($candidate_ids, $candidate_order) = load_candidates();
die "No candidate variants passed the requested filters\n" unless @$candidate_order;
my $id_file = "$work_prefix.ids.txt";
open my $id_fh, '>:raw', $id_file or die "Cannot write $id_file: $!\n";
print {$id_fh} "$_\n" for @$candidate_order;
close $id_fh;

my $subset = "$work_prefix.reference_subset";
my @extract_cmd = plink_base($input_prefix);
push @extract_cmd, '--extract', native_path($id_file), '--snps-only', 'just-acgt',
    '--max-alleles', '2', '--make-pgen', '--threads', $threads,
    '--out', native_path($subset);
run_cmd(@extract_cmd);

my %reference_present = read_pvar_ids("$subset.pvar");
open my $cache_fh, '>:raw', $output_cache or die "Cannot write $output_cache: $!\n";
open my $status_fh, '>:raw', $output_status or die "Cannot write $output_status: $!\n";
print {$cache_fh} join("\t", qw(query_snp proxy_snp ld_population proxy_r2 source)), "\n";
print {$status_fh} join("\t", qw(variant ld_population reference_status alt_count obs_count maf source)), "\n";

my ($estimable_total, $edge_total) = (0, 0);
for my $pop (@pops) {
    my $tag = lc($pop);
    my $freq_prefix = "$work_prefix.$tag.freq";
    my @freq_cmd = (plink_executable(), '--pfile', native_path($subset),
        '--keep-if', 'SuperPop', '==', $pop, '--freq', 'counts',
        '--threads', $threads, '--out', native_path($freq_prefix));
    run_cmd(@freq_cmd);
    my %frequency = read_acount("$freq_prefix.acount");

    my %estimable;
    for my $id (@$candidate_order) {
        my ($status, $alt, $obs, $maf) = ('ABSENT_REFERENCE', '', '', '');
        if ($reference_present{$id} && exists $frequency{$id}) {
            ($alt, $obs) = @{$frequency{$id}};
            if ($obs > 0 && $alt > 0 && $alt < $obs) {
                $status = 'POLYMORPHIC';
                $maf = $alt / $obs;
                $maf = 1 - $maf if $maf > 0.5;
                $estimable{$id} = 1;
                ++$estimable_total;
                print {$cache_fh} join("\t", uc($id), uc($id), $pop, 1, 'PLINK2_1KG_DIRECT'), "\n";
                ++$edge_total;
            } elsif ($obs > 0) {
                $status = 'MONOMORPHIC';
                $maf = 0;
            } else {
                $status = 'NO_GENOTYPE_OBSERVATIONS';
            }
        }
        print {$status_fh} join("\t", uc($id), $pop, $status, $alt, $obs,
            (length($maf) ? sprintf('%.8g', $maf) : ''), 'PLINK2_1KG_DIRECT'), "\n";
    }

    next unless keys(%estimable) > 1;
    my $ld_prefix = "$work_prefix.$tag.ld";
    my @ld_cmd = (plink_executable(), '--pfile', native_path($subset),
        '--keep-if', 'SuperPop', '==', $pop, '--mac', 1,
        '--r2-phased', 'allow-ambiguous-allele', 'cols=chrom,pos,id',
        '--ld-window-kb', $window_kb, '--ld-window-r2', $r2_threshold,
        '--threads', $threads, '--out', native_path($ld_prefix));
    run_cmd(@ld_cmd);
    my $report = -s "$ld_prefix.vcor" ? "$ld_prefix.vcor"
               : -s "$ld_prefix.vcor.zst" ? "$ld_prefix.vcor.zst" : '';
    die "PLINK2 did not produce an LD report for $pop\n" unless length $report;
    $edge_total += append_vcor_edges($report, $pop, \%estimable, $cache_fh);
}

close $cache_fh;
close $status_fh;
print "LD_SOURCE\tPLINK2_1KG_DIRECT\n";
print "CANDIDATES\t", scalar(@$candidate_order), "\n";
print "POPULATIONS\t", join(',', @pops), "\n";
print "ESTIMABLE_VARIANT_POPULATION_ROWS\t$estimable_total\n";
print "DIRECTED_CACHE_EDGES\t$edge_total\n";
print "OUTPUT_CACHE\t$output_cache\nOUTPUT_STATUS\t$output_status\n";

sub load_candidates {
    open my $fh, '<:raw', $candidates or die "Cannot read $candidates: $!\n";
    my $first = <$fh>;
    die "Empty candidate file: $candidates\n" unless defined $first;
    seek($fh, 0, 0) or die "Cannot rewind $candidates: $!\n";
    my $sep = index($first, "\t") >= 0 ? "\t" : ',';
    my $csv = Text::CSV->new({ binary => 1, auto_diag => 2, sep_char => $sep });
    my $header = $csv->getline($fh) or die "Cannot read candidate header\n";
    my %idx = map { uc(trim($header->[$_])) => $_ } 0 .. $#$header;
    die "Candidate file lacks $snp_column\n" unless exists $idx{uc $snp_column};
    if (length $signal_column) {
        die "Candidate file lacks $signal_column\n" unless exists $idx{uc $signal_column};
    }
    my (%ids, @order);
    while (my $row = $csv->getline($fh)) {
        my $id = trim($row->[$idx{uc $snp_column}]);
        next unless length $id;
        if (length $signal_column) {
            my $value = numeric($row->[$idx{uc $signal_column}]);
            next unless defined($value) && $value > 0 && $value <= $signal_threshold;
        }
        my $key = uc($id);
        next if $ids{$key}++;
        push @order, $id;
    }
    close $fh;
    return (\%ids, \@order);
}

sub read_pvar_ids {
    my ($path) = @_;
    open my $fh, '<:raw', $path or die "Cannot read $path: $!\n";
    my %ids;
    while (my $line = <$fh>) {
        next if $line =~ /^##/;
        $line =~ s/[\r\n]+\z//;
        my @f = split /\t/, $line, -1;
        next if $f[0] =~ /^#/;
        $ids{$f[2]} = 1 if defined($f[2]) && length($f[2]);
    }
    close $fh;
    return %ids;
}

sub read_acount {
    my ($path) = @_;
    open my $fh, '<:raw', $path or die "Cannot read $path: $!\n";
    my $header = <$fh> // die "Empty allele-count report: $path\n";
    $header =~ s/[\r\n]+\z//;
    my @h = split /\t/, $header, -1;
    s/^#// for @h;
    my %idx = map { uc($h[$_]) => $_ } 0 .. $#h;
    for my $required (qw(ID ALT_CTS OBS_CT)) {
        die "$path lacks $required\n" unless exists $idx{$required};
    }
    my %result;
    while (my $line = <$fh>) {
        $line =~ s/[\r\n]+\z//;
        my @f = split /\t/, $line, -1;
        my $alt = numeric($f[$idx{ALT_CTS}]);
        my $obs = numeric($f[$idx{OBS_CT}]);
        next unless defined($alt) && defined($obs);
        $result{$f[$idx{ID}]} = [$alt, $obs];
    }
    close $fh;
    return %result;
}

sub append_vcor_edges {
    my ($path, $pop, $estimable, $out) = @_;
    my $fh;
    if ($path =~ /\.zst\z/) {
        open $fh, '-|', ($ENV{ZSTD} || 'zstd'), '-dc', $path
            or die "Cannot decompress $path: $!\n";
    } else {
        open $fh, '<:raw', $path or die "Cannot read $path: $!\n";
    }
    my $header = <$fh> // die "Empty LD report: $path\n";
    $header =~ s/[\r\n]+\z//;
    my @h = split /\t/, $header, -1;
    s/^#// for @h;
    my %idx = map { uc($h[$_]) => $_ } 0 .. $#h;
    my $r2_name = exists($idx{PHASED_R2}) ? 'PHASED_R2' : 'R2';
    for my $required (qw(ID_A ID_B), $r2_name) {
        die "$path lacks $required\n" unless exists $idx{$required};
    }
    my $count = 0;
    while (my $line = <$fh>) {
        $line =~ s/[\r\n]+\z//;
        my @f = split /\t/, $line, -1;
        my ($a, $b) = @f[$idx{ID_A}, $idx{ID_B}];
        my $r2 = numeric($f[$idx{$r2_name}]);
        next unless defined($r2) && $estimable->{$a} && $estimable->{$b};
        next if $a eq $b;
        print {$out} join("\t", uc($a), uc($b), $pop, $r2, 'PLINK2_1KG_DIRECT'), "\n";
        print {$out} join("\t", uc($b), uc($a), $pop, $r2, 'PLINK2_1KG_DIRECT'), "\n";
        $count += 2;
    }
    close $fh;
    return $count;
}

sub plink_base {
    my ($prefix) = @_;
    my @cmd = (plink_executable());
    if (defined $pfile) {
        push @cmd, '--pfile', native_path($prefix);
        push @cmd, 'vzs' if -e "$prefix.pvar.zst";
    } else {
        push @cmd, '--bfile', native_path($prefix);
    }
    return @cmd;
}

sub plink_executable {
    return $plink2 if -f $plink2;
    return $plink2;
}

sub run_cmd {
    my (@cmd) = @_;
    system(@cmd);
    die "PLINK2 command failed (exit $?): @cmd\n" if $? != 0;
}

sub numeric {
    my ($value) = @_;
    return undef unless defined $value;
    $value = trim($value);
    return undef unless $value =~ /\A[+-]?(?:\d+(?:\.\d*)?|\.\d+)(?:[Ee][+-]?\d+)?\z/;
    return 0 + $value;
}

sub trim {
    my ($value) = @_;
    $value //= '';
    $value =~ s/^\s+|\s+$//g;
    return $value;
}

sub native_path {
    my ($path) = @_;
    my $abs = abs_path($path) || File::Spec->rel2abs($path);
    if ($abs =~ m{^/mnt/([A-Za-z])/(.*)$}) {
        my ($drive, $rest) = (uc($1), $2);
        $rest =~ s{/}{\\}g;
        return "$drive:\\$rest";
    }
    if ($^O =~ /cygwin/i && $abs =~ m{^/}) {
        if (open my $cygpath, '-|', 'cygpath', '-w', $abs) {
            my $win = <$cygpath> // '';
            close $cygpath;
            $win =~ s/[\r\n]+\z//;
            return $win if length $win;
        }
    }
    return $abs;
}

sub usage {
    return <<'USAGE';
Usage: build_plink2_1kg_ld_cache.pl --candidates FILE
  (--pfile PREFIX | --bfile PREFIX) --output-cache FILE --output-status FILE
  [--plink2 EXE] [--populations EUR,EAS] [--snp-column SNP]
  [--signal-column P --signal-threshold 5e-8]
  [--r2-threshold 0.1] [--window-kb 1000] [--threads 4]

Builds a direct, population-specific 1000 Genomes LD cache.  Only biallelic
ACGT variants are retained.  The status report explicitly distinguishes
polymorphic, monomorphic, and absent variants; only polymorphic variants get
self-edges and are considered LD-estimable.
USAGE
}
