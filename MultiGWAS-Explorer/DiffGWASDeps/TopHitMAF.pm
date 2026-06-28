package TopHitMAF;

use strict;
use warnings;
use Exporter qw(import);
use IO::Uncompress::Gunzip qw($GunzipError);

our @EXPORT_OK = qw(
  numeric
  format_num
  average_defined
  maf_from_effect_af
  derive_effect_af
  default_population_map
  parse_population_map
  infer_population_codes_for_text
  load_gnomad_lookup
  lookup_gnomad_maf
  unique_values
);

sub numeric {
    my ($x) = @_;
    return undef unless defined $x;
    $x =~ s/^\s+|\s+$//g;
    return undef unless length $x;
    return undef if $x =~ /^(?:NA|NaN|null|\.)$/i;
    return undef unless $x =~ /^[+-]?(?:\d+(?:\.\d*)?|\.\d+)(?:[eE][+-]?\d+)?$/;
    return 0 + $x;
}

sub format_num {
    my ($v) = @_;
    return '' unless defined $v;
    return sprintf('%.12g', $v);
}

sub average_defined {
    my @vals = grep { defined $_ } @_;
    return undef unless @vals;
    my $sum = 0;
    $sum += $_ for @vals;
    return $sum / @vals;
}

sub maf_from_effect_af {
    my ($af) = @_;
    return undef unless defined $af;
    return undef if $af < 0 || $af > 1;
    return $af <= 0.5 ? $af : (1 - $af);
}

sub derive_effect_af {
    my (%args) = @_;
    my $fa = numeric($args{frq_a});
    my $fu = numeric($args{frq_u});
    return average_defined($fa, $fu);
}

sub default_population_map {
    return {
        AFR        => 'AFR',
        AFRICAN    => 'AFR',
        AFRAM      => 'AFR',
        AFRAME     => 'AFR',
        AMR        => 'AMR',
        LAT        => 'AMR',
        LATINO     => 'AMR',
        HISPANIC   => 'AMR',
        HISPANICLATINO => 'AMR',
        ASN        => 'EAS',
        ASIAN      => 'EAS',
        EAS        => 'EAS',
        EASTASIAN  => 'EAS',
        EUR        => 'NFE',
        EUROPEAN   => 'NFE',
        NFE        => 'NFE',
        NONFINNISH => 'NFE',
        FIN        => 'FIN',
        FINNISH    => 'FIN',
        SAS        => 'SAS',
        SOUTHASIAN => 'SAS',
        ASJ        => 'ASJ',
        AJ         => 'ASJ',
        JEWISH     => 'ASJ',
        OTH        => 'OTH',
        OTHER      => 'OTH',
        ALL        => 'AF',
        GLOBAL     => 'AF',
        META       => 'AF',
        MIXED      => 'AF',
        TOTAL      => 'AF',
    };
}

sub parse_population_map {
    my ($raw) = @_;
    my %map = %{ default_population_map() };
    return \%map unless defined $raw;

    if (ref($raw) eq 'HASH') {
        for my $key (keys %{$raw}) {
            next unless defined $key && length $key;
            my $val = $raw->{$key};
            next unless defined $val && length $val;
            $map{ _norm_token($key) } = uc($val);
        }
        return \%map;
    }

    for my $entry (grep { length } split /\s*,\s*/, $raw) {
        my ($key, $val) = split /\s*=\s*/, $entry, 2;
        next unless defined $key && defined $val && length $key && length $val;
        $map{ _norm_token($key) } = uc($val);
    }
    return \%map;
}

sub infer_population_codes_for_text {
    my ($text, $map) = @_;
    $map ||= default_population_map();
    my @tokens = grep { length } map { _norm_token($_) } split /[^A-Za-z0-9]+/, ($text // '');
    my @codes;
    for my $tok (@tokens) {
        next unless exists $map->{$tok};
        push @codes, $map->{$tok};
    }
    return unique_values(@codes);
}

sub unique_values {
    my %seen;
    return grep { defined $_ && length $_ && !$seen{$_}++ } @_;
}

sub load_gnomad_lookup {
    my (%args) = @_;
    my $path = $args{file} || '';
    return undef unless length $path;
    return undef unless -s $path;

    my $fh;
    if ($path =~ /\.gz$/i) {
        $fh = IO::Uncompress::Gunzip->new($path)
          or die "Unable to read gnomAD lookup $path: $GunzipError\n";
    }
    else {
        open($fh, '<', $path) or die "Unable to read gnomAD lookup $path: $!\n";
    }

    my $header = <$fh>;
    die "gnomAD lookup file is empty: $path\n" unless defined $header;
    chomp $header;
    $header =~ s/\r$//;
    my @cols = split /\t/, $header, -1;
    my @norm = map { _norm_token($_) } @cols;
    my %idx = map { $norm[$_] => $_ } 0 .. $#norm;

    my $snp_idx = _first_existing_idx(\%idx, qw(SNP RSID ID VARIANTID VARIANT));
    my $chr_idx = _first_existing_idx(\%idx, qw(CHR CHROM CHROMOSOME));
    my $bp_idx  = _first_existing_idx(\%idx, qw(BP POS POSITION BASEPAIR));

    my (%by_snp, %by_pos);
    while (my $line = <$fh>) {
        chomp $line;
        $line =~ s/\r$//;
        next unless length $line;
        my @f = split /\t/, $line, -1;
        my %row;
        @row{@norm} = @f;
        if (defined $snp_idx) {
            my $snp = uc($f[$snp_idx] // '');
            $by_snp{$snp} = \%row if length $snp;
        }
        if (defined $chr_idx && defined $bp_idx) {
            my $chr = _normalize_chr($f[$chr_idx] // '');
            my $bp = $f[$bp_idx] // '';
            my $key = join(':', $chr, $bp);
            $by_pos{$key} = \%row if length($chr) && length($bp);
        }
    }
    close $fh;

    return {
        file     => $path,
        by_snp   => \%by_snp,
        by_pos   => \%by_pos,
    };
}

sub lookup_gnomad_maf {
    my (%args) = @_;
    my $lookup = $args{lookup} || return {};
    my $record = $args{record} || {};
    my $pop_codes = $args{pop_codes} || [];

    my $row;
    my $snp = uc($record->{SNP} // '');
    if (length $snp && exists $lookup->{by_snp}{$snp}) {
        $row = $lookup->{by_snp}{$snp};
    }
    if (!$row) {
        my $key = join(':', _normalize_chr($record->{CHR}), ($record->{BP} // ''));
        $row = $lookup->{by_pos}{$key} if length($key) > 1 && exists $lookup->{by_pos}{$key};
    }
    return {} unless $row;

    my @codes = @{$pop_codes};
    @codes = ('AF') unless @codes;
    my %pop_maf;
    for my $code (@codes) {
        my $af = _row_af_for_population($row, $code);
        next unless defined $af;
        my $maf = maf_from_effect_af($af);
        next unless defined $maf;
        $pop_maf{$code} = $maf;
    }
    if (!%pop_maf) {
        my $af = _row_af_for_population($row, 'AF');
        my $maf = maf_from_effect_af($af);
        $pop_maf{AF} = $maf if defined $maf;
    }
    return {} unless %pop_maf;

    my @ordered = sort { $pop_maf{$a} <=> $pop_maf{$b} || $a cmp $b } keys %pop_maf;
    my $selected = $ordered[0];
    return {
        maf          => $pop_maf{$selected},
        selected_pop => $selected,
        pop_mafs     => \%pop_maf,
    };
}

sub _first_existing_idx {
    my ($idx, @candidates) = @_;
    for my $cand (@candidates) {
        my $norm = _norm_token($cand);
        return $idx->{$norm} if exists $idx->{$norm};
    }
    return undef;
}

sub _row_af_for_population {
    my ($row, $pop) = @_;
    my $norm = _norm_token($pop);
    my @candidates = unique_values(
        'MAF',
        'AF',
        "AF$norm",
        "${norm}AF",
        "MAF$norm",
        "${norm}MAF",
        "GNOMAD${norm}AF",
        "GNOMAD${norm}MAF",
        "POPMAXAF",
    );
    for my $cand (@candidates) {
        next unless exists $row->{$cand};
        my $val = numeric($row->{$cand});
        return $val if defined $val;
    }
    return undef;
}

sub _normalize_chr {
    my ($chr) = @_;
    $chr //= '';
    $chr =~ s/^\s+|\s+$//g;
    $chr =~ s/^chr//i;
    return 'X' if $chr =~ /^(?:23|X)$/i;
    return $chr;
}

sub _norm_token {
    my ($text) = @_;
    $text //= '';
    $text = uc($text);
    $text =~ s/[^A-Z0-9]+//g;
    return $text;
}

1;
