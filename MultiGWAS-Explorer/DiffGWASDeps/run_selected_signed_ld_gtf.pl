#!/usr/bin/env perl
use strict;
use warnings;
use Getopt::Long qw(GetOptions);
use File::Spec;
use Fcntl qw(:flock);
use Digest::MD5 qw(md5_hex);
use Digest::SHA qw(sha256_hex);
use JSON::PP qw(decode_json);
use IO::Uncompress::Gunzip qw($GunzipError);

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
my %target_location;
my %target_p;
my $explicit_targets = length($opt{target_snps} // '') ? 1 : 0;
if (length($opt{targets_csv} // '')) {
    open my $fh, '<', $opt{targets_csv} or die "Open $opt{targets_csv}: $!\n";
    my $header = <$fh> // die "Empty selected-hit CSV: $opt{targets_csv}\n";
    chomp $header;
    $header =~ s/\r$//;
    my $sep = index($header, "\t") >= 0 ? "\t" : ',';
    my @header = split /\Q$sep\E/, $header, -1;
    my %idx = map { $header[$_] => $_ } 0 .. $#header;
    my %idx_lc = map { lc($header[$_]) => $_ } 0 .. $#header;
    my ($p_column) = grep { exists $idx_lc{$_} }
        qw(common_assoc_p focus_signal p top_p min_p p_value);
    die "Selected-hit CSV needs SNP, CHR, and BP columns\n"
        unless exists $idx{SNP} && exists $idx{CHR} && exists $idx{BP};
    while (my $line = <$fh>) {
        chomp $line;
        $line =~ s/\r$//;
        next unless length $line;
        my @fields = split /\Q$sep\E/, $line, -1;
        my $snp = $fields[$idx{SNP}] // '';
        $snp =~ s/^"|"$//g;
        add_target($snp, $fields[$idx{CHR}], $fields[$idx{BP}],
            defined($p_column) ? $fields[$idx_lc{$p_column}] : undef);
    }
    close $fh;
} else {
    add_target($_) for split /,/, $opt{target_snps};
}
sub add_target {
    my ($snp, $chr, $bp, $p) = @_;
    $snp =~ s/^\s+|\s+$//g;
    die "Unsafe selected SNP: $snp\n" unless $snp =~ /^[A-Za-z0-9_.:-]+$/;
    if (defined $p) {
        $p =~ s/^"|"$//g;
        if ($p =~ /^\s*(?:\d+(?:\.\d*)?|\.\d+)(?:[Ee][+-]?\d+)?\s*$/
            && $p >= 0 && $p <= 1) {
            $target_p{lc $snp} = 0 + $p
                if !exists($target_p{lc $snp}) || $p < $target_p{lc $snp};
        }
    }
    return if $seen{lc $snp}++;
    push @targets, $snp;
    if (defined($chr) && defined($bp) && $bp =~ /^\d+(?:\.\d+)?$/) {
        $chr =~ s/^"|"$//g;
        $chr =~ s/^chr//i;
        $chr = 23 if uc($chr) eq 'X';
        $target_location{lc $snp} = { chr => $chr, bp => 0 + $bp } if length $chr;
    }
}
die "No selected SNPs\n" unless @targets;

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
my $runner_config = do {
    local $/;
    decode_json(<$config_fh>);
};
close $config_fh;
die "Runner configuration must be a JSON object\n"
    unless ref($runner_config) eq 'HASH';
my $coordinate_cache = "$opt{output_html}.target_coordinates.tsv";
load_target_locations_from_tables([$coordinate_cache], \@targets, \%target_location);
resolve_target_locations_from_prior_plot_csvs($opt{output_html}, \@targets, \%target_location);
resolve_target_locations_from_ld_caches($runner_config, \@targets, \%target_location);
resolve_missing_target_locations($runner_config, \@targets, \%target_location);
write_target_location_cache($coordinate_cache, \@targets, \%target_location);
my $requested_ld_reference = $explicit_targets
    ? ($runner_config->{GTF_LD_REFERENCE_SNP} || $runner_config->{LOCAL_LD_REFERENCE_SNP} || '') : '';
my @loci = build_overlapping_loci(\@targets, \%target_location,
    $opt{window_bp}, \%target_p, $explicit_targets, $requested_ld_reference);
my @locus_ids = map { $_->{id} } @loci;
my %locus_by_id = map { $_->{id} => $_ } @loci;
my $merged_target_count = scalar(@targets) - scalar(@loci);
print "[signed-LD] Grouped " . scalar(@targets) . " target SNP(s) into "
    . scalar(@loci) . " non-overlapping genomic locus/loci using a +/-$opt{window_bp}-bp window.\n";
print "[signed-LD] Merged $merged_target_count nearby target SNP(s), "
    . ($explicit_targets ? 'retaining every requested SNP as a plot label.'
                         : 'labeling only the smallest-P lead in each region.') . "\n"
    if $merged_target_count;
my %render_config = map { $_ => $runner_config->{$_} }
    grep { $_ ne 'GTF_LD_R2_CACHE_BY_SNP'
        && /^(?:GTF_|LOCAL_|DISPLAY_GWAS|GROUP_TRACKS$|PAIR_DEFS$|DATA_GZ$|SOURCE_LONG_GZ$|SOURCE_MODE$|EXTRACTOR_CONFIG_JSON$|REFERENCE_BUILD$)/ }
    keys %$runner_config;
my $deps_dir = File::Spec->catpath((File::Spec->splitpath(__FILE__))[0,1], '');
my $perl_planner = File::Spec->catfile($deps_dir, 'plan_local_gtf_labels.pl');
my $perl_planner_available = eval { require $perl_planner; 1 };
warn "[signed-LD] Perl label planner unavailable ($@); SAS ODA will plan labels.\n"
    unless $perl_planner_available;
my @plot_sources = qw(
    plan_local_gtf_labels.pl
    Plan_Local_GTF_Target_Labels.sas
    run_sas_oda_single_snp_with_gtf.sas
    SNP_Local_Manhattan_With_GTF.sas
    Lattice_gscatter_over_bed_track.sas
    rank4grps.sas
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
    my $digest;
    if (!-f $path && $_ eq 'plan_local_gtf_labels.pl') {
        $digest = sha256_hex('SAS fallback: Perl planner unavailable');
    } else {
        open my $source, '<:raw', $path or die "Cannot read plot source $path: $!\n";
        local $/;
        $digest = sha256_hex(<$source>);
        close $source;
    }
    $_ => $digest;
} @plot_sources;
my $render_request = JSON::PP->new->canonical(1)->encode({
    version => 4, config => \%render_config,
    data => file_signature($runner_config->{DATA_GZ}),
    source_long => file_signature($runner_config->{SOURCE_LONG_GZ}),
    plot_sources => \%plot_source_hashes,
    window_bp => $opt{window_bp}, population => uc($opt{population}),
    mode => lc($opt{mode}),
    loci => [ map {
        +{ reference => $_->{reference}, targets => $_->{targets},
           chr => $_->{chr}, bp => $_->{bp}, plot_window_bp => $_->{plot_window_bp} }
    } @loci ],
});
my @links;
my $shared_macros_ready = 0;
my $progress_path = File::Spec->catpath($volume, $dir, "${output_stem}.progress.json");
my $partial_path = File::Spec->catpath($volume, $dir, "${output_stem}.partial.html");
my %status = map { $_ => 'pending' } @locus_ids;
my $request_id = md5_hex($render_request);
my $failure = '';
my %status_write_warning;
write_progress();
for my $locus (@loci) {
    my $snp = $locus->{reference};
    my @label_snps = $explicit_targets ? @{ $locus->{targets} } : ($snp);
    my $label_text = join(',', @label_snps);
    my $plot_window_bp = $locus->{plot_window_bp};
    my $label_plan;
    if ($perl_planner_available && lc($ENV{GTF_LABEL_PLAN_BACKEND} || '') ne 'sas') {
        $label_plan = eval { LocalGTFLabelPlanner::plan(
            targets => [map { +{ snp => $_,
                bp => exists($target_location{lc $_})
                    ? $target_location{lc $_}{bp} : undef } } @label_snps],
            reference_bp => $locus->{bp}, window_bp => $plot_window_bp,
            design_width => $runner_config->{GTF_DESIGN_WIDTH} || 950,
            design_height => $runner_config->{GTF_DESIGN_HEIGHT} || 1000,
            font_size => $runner_config->{GTF_LABEL_FONT_SIZE} || 10,
            layout => $runner_config->{GTF_LABEL_LAYOUT} || 'auto',
        ) };
        warn "[signed-LD] Perl label planning failed ($@); SAS ODA will plan labels.\n"
            if $@;
    }
    $label_plan ||= {
        layout => $runner_config->{GTF_LABEL_LAYOUT} || 'auto',
        positions => '', headroom_frac => '', center_offset => '',
        font_size => $runner_config->{GTF_LABEL_FONT_SIZE} || 10,
        reason => 'sas_fallback',
    };
    my $locus_id = $locus->{id};
    (my $safe_snp = $snp) =~ s/[^A-Za-z0-9._-]/_/g;
    my $target_html = @loci == 1 ? $base : "${output_stem}_${safe_snp}.html";
    my $target_csv = "${output_stem}_${safe_snp}_top_hit.csv";
    my $target_path = File::Spec->catpath($volume, $dir, $target_html);
    (my $target_png = $target_path) =~ s/\.html$/.png/i;
    my $request_file = "$target_path.request.md5";
    my $ld_cache = $runner_config->{GTF_LD_R2_CACHE_BY_SNP}{lc $snp}
        // $runner_config->{GTF_LD_R2_CACHE_BY_SNP}{$snp}
        // '';
    my $request_key = md5_hex(join("\0", $render_request, $snp, $label_text, $plot_window_bp,
        $label_plan->{layout}, $label_plan->{positions}, $label_plan->{font_size},
        $label_plan->{headroom_frac}, $label_plan->{center_offset},
        JSON::PP->new->canonical(1)->encode(file_signature($ld_cache))));
    if (-s $target_path && -s $target_png && -s $request_file) {
        open my $request_fh, '<', $request_file or die "Cannot read $request_file: $!\n";
        my $saved = <$request_fh> // '';
        close $request_fh;
        $saved =~ s/\s+\z//;
        if ($saved eq $request_key) {
            print "[signed-LD] Reusing completed plot for $snp\n";
            push @links, [$label_text, $target_html, locus_region($locus)];
            $status{$locus_id} = 'complete';
            write_progress();
            next;
        }
    }
    # A failed rerender must never leave an old matching sidecar behind.
    unlink $request_file if -e $request_file;
    $status{$locus_id} = 'running';
    write_progress();
    local %ENV = %ENV;
    delete @ENV{qw(DATA_GZ REMOTE_DATA_BASENAME)};
    $ENV{RUNNER_CONFIG_JSON} = $opt{runner_config};
    $ENV{TARGET_SNP} = $snp;
    $ENV{LOCAL_WINDOW_BP} = $plot_window_bp;
    $ENV{OUTPUT_HTML_BASENAME} = $target_html;
    $ENV{SINGLE_SNP_ALLOW_GENERIC_OUTPUT_BASENAME} = 1;
    $ENV{SINGLE_SNP_TOP_HITS_CSV_BASENAME} = $target_csv;
    $ENV{GTF_LABEL_SNPS} = $label_text;
    $ENV{GTF_LABEL_LAYOUT} = $label_plan->{layout};
    $ENV{GTF_LABEL_POSITIONS} = $label_plan->{positions};
    $ENV{GTF_LABEL_FONT_SIZE} = $label_plan->{font_size};
    $ENV{GTF_LABEL_HEADROOM_FRAC} = $label_plan->{headroom_frac};
    $ENV{GTF_LABEL_CENTER_OFFSET} = $label_plan->{center_offset};
    $ENV{GTF_LD_DISPLAY_MODE} = lc $opt{mode};
    $ENV{GTF_LD_R2_CACHE} = $ld_cache;
    $ENV{GTF_LD_REFERENCE_SNP} = $snp;
    $ENV{GTF_LD_HEATMAP_LEGEND_TITLE} =
        "Signed LD r2 x sign(Z) to $snp (" . uc($opt{population}) . ", 1000G Phase 3 / PLINK2)";
    $ENV{OPEN_RESULT} = 0 if @loci > 1;
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
    print "[signed-LD] Rendering locus labels [$label_text] with $label_plan->{layout} layout ($label_plan->{reason}) and PLINK2 LD reference $snp"
        . ($plot_window_bp != 0 + $opt{window_bp} ? " and expanded half-window $plot_window_bp bp" : '') . "\n";
    my $rc = system('/bin/bash', $opt{single_runner});
    if ($rc != 0 || !-s $target_path || !-s $target_png) {
        my $exit_code = $rc == -1 ? 1 : ($rc >> 8 || 1);
        $failure = "Signed-LD SAS runner failed for $snp (exit $exit_code)";
        $failure .= '; missing HTML or PNG' if $rc == 0;
        $status{$locus_id} = 'failed';
        write_progress();
        print STDERR "[signed-LD] $failure. Re-run the same pipeline command to resume remaining loci. Progress: $progress_path\n";
        exit $exit_code;
    }
    $shared_macros_ready = 1;
    atomic_write($request_file, "$request_key\n");
    push @links, [$label_text, $target_html, locus_region($locus)];
    $status{$locus_id} = 'complete';
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

sub resolve_missing_target_locations {
    my ($config, $targets, $locations) = @_;
    my %missing = map { lc($_) => $_ } grep { !exists $locations->{lc $_} } @$targets;
    return unless %missing;
    for my $candidate ($config->{DATA_GZ}, $config->{SOURCE_LONG_GZ}) {
        next unless defined($candidate) && length($candidate);
        my $path = local_path($candidate);
        next unless -s $path;
        my $fh;
        if ($path =~ /\.gz$/i) {
            $fh = IO::Uncompress::Gunzip->new($path);
            next unless $fh;
        } else {
            open $fh, '<', $path or next;
        }
        my $header = <$fh> // '';
        $header =~ s/[\r\n]+\z//;
        my $sep = index($header, "\t") >= 0 ? "\t" : ',';
        my @h = split /\Q$sep\E/, $header, -1;
        my %idx = map { uc($h[$_]) => $_ } 0 .. $#h;
        unless (exists($idx{SNP}) && exists($idx{CHR}) && exists($idx{BP})) {
            close $fh;
            next;
        }
        while (my $line = <$fh>) {
            $line =~ s/[\r\n]+\z//;
            my @f = split /\Q$sep\E/, $line, -1;
            my $key = lc($f[$idx{SNP}] // '');
            next unless exists $missing{$key};
            my ($chr, $bp) = @f[$idx{CHR}, $idx{BP}];
            next unless defined($chr) && defined($bp) && $bp =~ /^\d+(?:\.\d+)?$/;
            $chr =~ s/^chr//i;
            $chr = 23 if uc($chr) eq 'X';
            next unless length $chr;
            $locations->{$key} = { chr => $chr, bp => 0 + $bp };
            delete $missing{$key};
            last unless %missing;
        }
        close $fh;
        last unless %missing;
    }
    warn "[signed-LD] Coordinates were not found for " . join(', ', values %missing)
        . "; those targets will remain separate loci.\n" if %missing;
}

sub resolve_target_locations_from_ld_caches {
    my ($config, $targets, $locations) = @_;
    my %wanted = map { lc($_) => 1 } grep { !exists $locations->{lc $_} } @$targets;
    return unless %wanted;
    my @caches;
    if (ref($config->{GTF_LD_R2_CACHE_BY_SNP}) eq 'HASH') {
        push @caches, values %{ $config->{GTF_LD_R2_CACHE_BY_SNP} };
    }
    push @caches, $config->{GTF_LD_R2_CACHE}, $config->{LOCAL_LD_CACHE_TSV};
    my %seen_cache;
    for my $cache (grep { defined($_) && length($_) && !$seen_cache{$_}++ } @caches) {
        my $path = local_path($cache);
        next unless -s $path;
        open my $fh, '<', $path or next;
        my $header = <$fh> // '';
        $header =~ s/[\r\n]+\z//;
        my @h = split /\t/, $header, -1;
        my %idx = map { lc($h[$_]) => $_ } 0 .. $#h;
        unless (exists($idx{query_snp}) && exists($idx{query_chr}) && exists($idx{query_bp})) {
            close $fh;
            next;
        }
        while (my $line = <$fh>) {
            $line =~ s/[\r\n]+\z//;
            my @f = split /\t/, $line, -1;
            my $key = lc($f[$idx{query_snp}] // '');
            next unless $wanted{$key};
            my ($chr, $bp) = @f[$idx{query_chr}, $idx{query_bp}];
            next unless defined($chr) && length($chr) && defined($bp) && $bp =~ /^\d+(?:\.\d+)?$/;
            $chr =~ s/^chr//i;
            $chr = 23 if uc($chr) eq 'X';
            $locations->{$key} = { chr => $chr, bp => 0 + $bp };
            delete $wanted{$key};
            last unless %wanted;
        }
        close $fh;
        last unless %wanted;
    }
}

sub resolve_target_locations_from_prior_plot_csvs {
    my ($output_html, $targets, $locations) = @_;
    return unless grep { !exists $locations->{lc $_} } @$targets;
    my ($volume, $dir, $base) = File::Spec->splitpath($output_html);
    $base =~ s/\.html$//i;
    my $prefix = File::Spec->catpath($volume, $dir, $base);
    my @tables = glob("${prefix}_*_top_hit.csv");
    load_target_locations_from_tables(\@tables, $targets, $locations);
}

sub load_target_locations_from_tables {
    my ($tables, $targets, $locations) = @_;
    my %wanted = map { lc($_) => 1 } grep { !exists $locations->{lc $_} } @$targets;
    return unless %wanted;
    for my $path (@$tables) {
        next unless defined($path) && -s $path;
        open my $fh, '<', $path or next;
        my $header = <$fh> // '';
        $header =~ s/[\r\n]+\z//;
        my $sep = index($header, "\t") >= 0 ? "\t" : ',';
        my @h = split /\Q$sep\E/, $header, -1;
        my %idx = map { uc($h[$_]) => $_ } 0 .. $#h;
        unless (exists($idx{SNP}) && exists($idx{CHR}) && exists($idx{BP})) {
            close $fh;
            next;
        }
        while (my $line = <$fh>) {
            $line =~ s/[\r\n]+\z//;
            my @f = split /\Q$sep\E/, $line, -1;
            my $key = lc($f[$idx{SNP}] // '');
            next unless $wanted{$key};
            my ($chr, $bp) = @f[$idx{CHR}, $idx{BP}];
            next unless defined($chr) && length($chr) && defined($bp) && $bp =~ /^\d+(?:\.\d+)?$/;
            $chr =~ s/^chr//i;
            $chr = 23 if uc($chr) eq 'X';
            $locations->{$key} = { chr => $chr, bp => 0 + $bp };
            delete $wanted{$key};
            last unless %wanted;
        }
        close $fh;
        last unless %wanted;
    }
}

sub write_target_location_cache {
    my ($path, $targets, $locations) = @_;
    my @known = grep { exists $locations->{lc $_} } @$targets;
    return unless @known;
    my $contents = "SNP\tCHR\tBP\n" . join('', map {
        my $loc = $locations->{lc $_};
        join("\t", $_, $loc->{chr}, $loc->{bp}) . "\n";
    } @known);
    eval { atomic_write($path, $contents); 1 }
        or warn "[signed-LD] Could not persist target-coordinate cache $path: " . ($@ || 'unknown error');
}

sub local_path {
    my ($path) = @_;
    return $path unless $^O eq 'cygwin' && $path =~ m{^([A-Za-z]):[/\\](.*)$};
    my ($drive, $rest) = (lc($1), $2);
    $rest =~ s{\\}{/}g;
    for my $prefix ("/mnt/$drive", "/cygdrive/$drive") {
        my $candidate = "$prefix/$rest";
        return $candidate if -e $candidate;
    }
    return $path;
}

sub build_overlapping_loci {
    my ($targets, $locations, $window_text, $pvalues, $explicit, $reference_override) = @_;
    my $window = 0 + $window_text;
    die "--window-bp must be a positive number\n" unless $window > 0;
    my %input_order = map { lc($targets->[$_]) => $_ } 0 .. $#$targets;
    my @located = map {
        +{ snp => $_, %{ $locations->{lc $_} }, input_order => $input_order{lc $_} }
    } grep { exists $locations->{lc $_} } @$targets;
    @located = sort {
        chromosome_sort_key($a->{chr}) cmp chromosome_sort_key($b->{chr})
            || $a->{bp} <=> $b->{bp} || $a->{input_order} <=> $b->{input_order}
    } @located;
    my @groups;
    for my $row (@located) {
        my $start = $row->{bp} - $window;
        my $end = $row->{bp} + $window;
        if (!@groups || normalize_chr($groups[-1]{chr}) ne normalize_chr($row->{chr})
            || $start > $groups[-1]{window_end}) {
            push @groups, { chr => $row->{chr}, window_end => $end, members => [] };
        }
        push @{ $groups[-1]{members} }, $row;
        $groups[-1]{window_end} = $end if $end > $groups[-1]{window_end};
    }
    for my $snp (@$targets) {
        next if exists $locations->{lc $snp};
        push @groups, { chr => '', members => [ { snp => $snp,
            input_order => $input_order{lc $snp} } ] };
    }
    my @loci;
    for my $group (@groups) {
        my @members = sort { $a->{input_order} <=> $b->{input_order} } @{ $group->{members} };
        my $reference = $members[0];
        if ($explicit && length($reference_override // '')) {
            my ($matching_reference) = grep {
                lc($_->{snp}) eq lc($reference_override)
            } @members;
            $reference = $matching_reference if $matching_reference;
        } elsif (!$explicit) {
            ($reference) = sort {
                (exists($pvalues->{lc $a->{snp}}) ? 0 : 1)
                    <=> (exists($pvalues->{lc $b->{snp}}) ? 0 : 1)
                || ($pvalues->{lc $a->{snp}} // 1) <=> ($pvalues->{lc $b->{snp}} // 1)
                || $a->{input_order} <=> $b->{input_order}
            } @members;
        }
        my $plot_window = $window;
        if (defined $reference->{bp}) {
            for my $member (@members) {
                next unless defined $member->{bp};
                my $required = abs($member->{bp} - $reference->{bp}) + $window;
                $plot_window = $required if $required > $plot_window;
            }
        }
        push @loci, {
            id => lc($reference->{snp}), reference => $reference->{snp},
            targets => [ map { $_->{snp} } @members ], chr => ($reference->{chr} // ''),
            bp => ($reference->{bp} // ''), plot_window_bp => int($plot_window + 0.5),
            first_order => $reference->{input_order},
        };
    }
    return sort { $a->{first_order} <=> $b->{first_order} } @loci;
}

sub normalize_chr {
    my ($chr) = @_;
    $chr = '' unless defined $chr;
    $chr =~ s/^chr//i;
    return uc($chr) eq 'X' ? '23' : uc($chr);
}

sub locus_region {
    my ($locus) = @_;
    my $chr = $locus->{chr} // '';
    return '' unless length $chr;
    my @bp = sort { $a <=> $b } map { $target_location{lc $_}{bp} }
        grep { exists $target_location{lc $_} } @{ $locus->{targets} };
    return '' unless @bp;
    $chr = 'X' if normalize_chr($chr) eq '23';
    return 'chr' . $chr . ':' . $bp[0]
        . (@bp > 1 && $bp[-1] != $bp[0] ? '-' . $bp[-1] : '');
}

sub chromosome_sort_key {
    my ($chr) = @_;
    $chr = normalize_chr($chr);
    return $chr =~ /^\d+$/ ? sprintf('0%04d', $chr) : "1$chr";
}

if (@links > 1) {
    atomic_write($opt{output_html}, plot_index_html(\@links));
    print "[signed-LD] Wrote plot index: $opt{output_html}\n";
}
unlink $partial_path if -e $partial_path;
write_progress();

sub write_progress {
    my $complete = scalar grep { $status{$_} eq 'complete' } @locus_ids;
    my $report = {
        version => 1, request_id => $request_id, output_html => $opt{output_html},
        total => scalar(@loci), complete => $complete,
        remaining => scalar(@loci) - $complete,
        failure => $failure,
        loci => [ map {
            my $locus = $locus_by_id{$_};
            +{ snp => $locus->{reference}, targets => $locus->{targets}, status => $status{$_} }
        } @locus_ids ],
    };
    write_status_file($progress_path, JSON::PP->new->canonical(1)->pretty(1)->encode($report));
    if ($complete < @loci && @links) {
        write_status_file($partial_path, plot_index_html(\@links));
    }
    # The full index is a completion signal; keep it only when every locus
    # belongs to the current request and has verified output files.
    unlink $opt{output_html} if @loci > 1 && $complete < @loci && -e $opt{output_html};
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
    $html .= '<h1>Signed-LD local GTF plots</h1>' . "\n";
    $html .= '<p>Targets in separate genomic regions have separate plots.</p>' . "\n";
    $html .= '<ul>' . "\n";
    for my $link (@$links) {
        my $label = (length($link->[2] // '') ? "$link->[2] | " : '') . $link->[0];
        $html .= '<li><a href="' . html_escape($link->[1]) . '">'
            . html_escape($label) . '</a></li>' . "\n";
    }
    return $html . '</ul></body></html>' . "\n";
}

sub html_escape {
    my ($text) = @_;
    $text //= '';
    $text =~ s/&/&amp;/g;
    $text =~ s/</&lt;/g;
    $text =~ s/>/&gt;/g;
    $text =~ s/"/&quot;/g;
    return $text;
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
