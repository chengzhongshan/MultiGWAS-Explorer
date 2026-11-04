package OneKGPhase3MAF;
use strict;
use warnings;
use Exporter qw(import);
use Digest::MD5 qw(md5_hex);
use File::Basename qw(dirname);
use File::Path qw(make_path);
use File::Spec;
use File::Temp qw(tempdir);
use Cwd qw(abs_path);

our @EXPORT_OK = qw(query_phase3_maf);
my @POPS = qw(AFR AMR EAS EUR SAS);

sub query_phase3_maf {
    my (%arg) = @_;
    my $build = lc($arg{build} || 'hg38');
    die "1KG Phase 3 MAF fallback supports hg38 and hg19, not $build\n"
        unless $build eq 'hg38' || $build eq 'hg19';
    my $plink2 = $arg{plink2} || 'plink2';
    my $panel_dir = $arg{panel_dir} || File::Spec->catdir(dirname(__FILE__), '..', 'cache', 'plink2_1kg_hg38');
    my $pfile = $arg{pfile} || File::Spec->catfile(dirname(__FILE__), '..', 'cache', 'plink2_1kg_phase3', 'all_phase3');
    my $cache_dir = $arg{cache_dir} || File::Spec->catdir(dirname(__FILE__), '..', 'cache', 'top_hit_1kg_maf');
    make_path($cache_dir) unless -d $cache_dir;
    my %by_chr;
    for my $row (@{ $arg{rows} || [] }) {
        my $snp = $row->{SNP} // '';
        my $chr = $row->{CHR} // '';
        $chr =~ s/^chr//i;
        $chr = 'X' if uc($chr) eq 'X' || $chr eq '23';
        next unless $snp =~ /^rs\d+$/i;
        next unless $chr eq 'X' || ($chr =~ /^\d+$/ && $chr >= 1 && $chr <= 22);
        next unless defined($row->{BP}) && $row->{BP} =~ /^\d+$/;
        $by_chr{$chr}{uc $snp}{id} = lc($snp);
        $by_chr{$chr}{uc $snp}{positions}{$row->{BP}} = 1;
    }
    my %result;
    for my $chr (sort { $a eq 'X' ? 1 : $b eq 'X' ? -1 : $a <=> $b } keys %by_chr) {
        my $prefix = $build eq 'hg38'
            ? File::Spec->catfile($panel_dir, "chr${chr}_hg38") : $pfile;
        if ($build eq 'hg38') {
            _run($^X, File::Spec->catfile(dirname(__FILE__), 'prepare_plink2_1kg_hg38_chr.pl'),
                '--chr', $chr, '--output-dir', $panel_dir, '--plink2', $plink2);
        }
        die "Missing 1KG Phase 3 panel $prefix\n"
            unless -s "$prefix.pgen" && (-s "$prefix.pvar.zst" || -s "$prefix.pvar") && -s "$prefix.psam";
        my @ids = sort map { $by_chr{$chr}{$_}{id} } keys %{ $by_chr{$chr} };
        my @reference_stats = map { my @s = stat($_); join(':', $_, @s ? @s[7,9] : ('missing','missing')) }
            "$prefix.pgen", (-s "$prefix.pvar.zst" ? "$prefix.pvar.zst" : "$prefix.pvar"), "$prefix.psam";
        my @position_signature = map {
            $_ . ':' . join(',', sort { $a <=> $b } keys %{ $by_chr{$chr}{$_}{positions} })
        } sort keys %{ $by_chr{$chr} };
        my $key = md5_hex(join("\0", 'position-checked-v1', $build, $chr,
            @position_signature, @reference_stats, $plink2));
        my $cache = File::Spec->catfile($cache_dir, "chr${chr}_${key}.tsv");
        unless (-s $cache) {
            my $tmp = tempdir('phase3_maf.XXXXXX', DIR => $cache_dir, CLEANUP => 1);
            my $ids_file = File::Spec->catfile($tmp, 'ids.txt');
            open my $ids_fh, '>:raw', $ids_file or die "Write $ids_file: $!\n";
            print {$ids_fh} "$_\n" for @ids;
            close $ids_fh;
            my $subset = File::Spec->catfile($tmp, 'subset');
            my @input = ($plink2, '--pfile', _native($prefix));
            push @input, 'vzs' if -s "$prefix.pvar.zst";
            push @input, ('--chr', $chr) if $build eq 'hg19';
            my @subset_cmd = (@input, '--extract', _native($ids_file), '--snps-only', 'just-acgt',
                '--max-alleles', 2, '--make-pgen', '--out', _native($subset));
            my $subset_status = system(@subset_cmd);
            my $has_subset = $subset_status == 0;
            unless ($has_subset) {
                my $log = "$subset.log";
                open my $log_fh, '<', $log or die "Cannot inspect failed 1KG extraction $log: $!\n";
                local $/;
                my $log_text = <$log_fh> // '';
                close $log_fh;
                die "1KG Phase 3 frequency command failed (exit $subset_status): @subset_cmd\n"
                    unless $log_text =~ /No variants remaining after main filters/;
                print STDERR "[1kg-maf] chr$chr: requested IDs absent from biallelic SNP reference; retaining as UNKNOWN\n";
            }
            my %values;
            for my $pop ($has_subset ? @POPS : ()) {
                my $out = File::Spec->catfile($tmp, lc($pop));
                _run($plink2, '--pfile', _native($subset), '--keep-if', 'SuperPop',
                    '==', $pop, '--nonfounders', '--freq', 'counts', '--out', _native($out));
                _read_counts("$out.acount", $pop, \%values);
            }
            my %reference_positions = $has_subset
                ? %{ _read_reference_positions("$subset.pvar") } : ();
            my $pending = "$cache.tmp.$$";
            open my $out, '>:raw', $pending or die "Write $pending: $!\n";
            print {$out} join("\t", 'SNP', @POPS), "\n";
            for my $id (@ids) {
                my $key = uc($id);
                my @observed = @{ $reference_positions{$key} || [] };
                my @expected = keys %{ $by_chr{$chr}{$key}{positions} };
                my $position_ok = @observed == 1 && @expected == 1
                    && $observed[0] == $expected[0];
                print {$out} join("\t", $key, map {
                    $position_ok && defined $values{$key}{$_}
                        ? sprintf('%.12g', $values{$key}{$_}) : ''
                } @POPS), "\n";
            }
            close $out or die "Close $pending: $!\n";
            rename $pending, $cache or die "Install $cache: $!\n";
        }
        open my $in, '<:raw', $cache or die "Read $cache: $!\n";
        <$in>;
        while (my $line = <$in>) {
            $line =~ s/[\r\n]+\z//;
            my ($id, @maf) = split /\t/, $line, -1;
            my %per_pop;
            for my $i (0 .. $#POPS) {
                $per_pop{$POPS[$i]} = 0 + $maf[$i]
                    if defined $maf[$i] && $maf[$i] =~ /^(?:\d+(?:\.\d*)?|\.\d+)(?:e[+-]?\d+)?$/i;
            }
            my @known = sort { $per_pop{$a} <=> $per_pop{$b} } keys %per_pop;
            $result{$id} = {
                maf => @known ? $per_pop{$known[0]} : undef,
                max_maf => @known ? $per_pop{$known[-1]} : undef,
                per_pop => \%per_pop,
                complete => @known == @POPS ? 1 : 0,
            };
        }
        close $in;
    }
    return \%result;
}

sub _read_reference_positions {
    my ($file) = @_;
    open my $fh, '<:raw', $file or die "Read $file: $!\n";
    my %positions;
    while (my $line = <$fh>) {
        next if $line =~ /^#/;
        $line =~ s/[\r\n]+\z//;
        my @f = split /\t/, $line, -1;
        next unless @f >= 3;
        push @{ $positions{uc $f[2]} }, $f[1];
    }
    close $fh;
    return \%positions;
}

sub _read_counts {
    my ($file, $pop, $out) = @_;
    open my $fh, '<:raw', $file or die "Read $file: $!\n";
    my $header = <$fh> // die "Empty PLINK frequency file $file\n";
    $header =~ s/[\r\n]+\z//;
    my @h = split /\t/, $header, -1;
    s/^#// for @h;
    my %idx = map { uc($h[$_]) => $_ } 0 .. $#h;
    die "PLINK frequency file lacks ID, ALT_CTS, or OBS_CT: $file\n"
        unless exists $idx{ID} && exists $idx{ALT_CTS} && exists $idx{OBS_CT};
    while (my $line = <$fh>) {
        $line =~ s/[\r\n]+\z//;
        my @f = split /\t/, $line, -1;
        my $obs = $f[$idx{OBS_CT}] + 0;
        next unless $obs > 0;
        my $af = ($f[$idx{ALT_CTS}] + 0) / $obs;
        my $maf = $af <= 0.5 ? $af : 1 - $af;
        $out->{uc $f[$idx{ID}]}{$pop} = $maf;
    }
    close $fh;
}

sub _native {
    my ($path) = @_;
    my $abs = abs_path($path) || File::Spec->rel2abs($path);
    if ($^O =~ /cygwin/i && $abs =~ m{^/}) {
        open my $fh, '-|', 'cygpath', '-w', $abs or die "cygpath: $!\n";
        my $value = <$fh> // '';
        close $fh;
        $value =~ s/[\r\n]+$//;
        return $value if length $value;
    }
    return $abs;
}

sub _run {
    my @cmd = @_;
    system(@cmd) == 0 or die "1KG Phase 3 frequency command failed (exit $?): @cmd\n";
}
1;
