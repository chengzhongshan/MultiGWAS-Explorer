#!/usr/bin/env perl
use strict;
use warnings;
use DBI;
use Getopt::Long qw(GetOptions);
use Text::CSV;
use Time::HiRes qw(time);

my ($candidates, $cache, $sqlite, $output_leads, $output_audit);
my $snp_column = 'SNP';
my $chr_column = 'CHR';
my $bp_column = 'BP';
my $signal_column = '';
my $signal_columns = '';
my $signal_threshold = 1;
my $populations = 'EUR ASN';
my $r2_threshold = 0.2;
my $population_rule = 'ANY';
my $fallback_distance_bp = 1_000_000;
my $max_leads = 0;
my $rebuild_sqlite = 0;

GetOptions(
    'candidates=s'          => \$candidates,
    'cache=s'               => \$cache,
    'sqlite=s'              => \$sqlite,
    'output-leads=s'        => \$output_leads,
    'output-audit=s'        => \$output_audit,
    'snp-column=s'          => \$snp_column,
    'chr-column=s'          => \$chr_column,
    'bp-column=s'           => \$bp_column,
    'signal-column=s'       => \$signal_column,
    'signal-columns=s'      => \$signal_columns,
    'signal-threshold=f'    => \$signal_threshold,
    'populations=s'         => \$populations,
    'r2-threshold=f'        => \$r2_threshold,
    'population-rule=s'     => \$population_rule,
    'fallback-distance-bp=i'=> \$fallback_distance_bp,
    'max-leads=i'           => \$max_leads,
    'rebuild-sqlite!'       => \$rebuild_sqlite,
) or die usage();

die usage() unless defined $candidates && defined $cache && defined $sqlite
    && defined $output_leads && defined $output_audit;
die "Specify --signal-column or --signal-columns, not both\n"
    if length($signal_column) && length($signal_columns);
die "Specify --signal-column or --signal-columns\n"
    unless length($signal_column) || length($signal_columns);
die "--r2-threshold must be between 0.2 and 1 for the downloadable HaploReg archive\n"
    unless $r2_threshold >= 0.2 && $r2_threshold <= 1;
$population_rule = uc $population_rule;
die "--population-rule must be ANY or ALL\n" unless $population_rule =~ /\A(?:ANY|ALL)\z/;
my @populations = grep { length } map { uc } split /[\s,]+/, $populations;
die "At least one population is required\n" unless @populations;

my $started = time;
my ($rows, $header) = load_candidate_rows();
my $db_started = time;
my $dbh = prepare_cache_database();
my $db_seconds = time - $db_started;

my %pruned;
my @selected;
my %audit;
my $lookup = $dbh->prepare(
    'SELECT proxy_snp, ld_population, proxy_r2 FROM edges '
  . 'WHERE query_snp=? AND proxy_r2>=? AND ld_population IN ('
  . join(',', ('?') x @populations) . ')'
);
my $presence = $dbh->prepare(
    'SELECT DISTINCT ld_population FROM edges WHERE query_snp=? AND proxy_snp=query_snp '
  . 'AND ld_population IN (' . join(',', ('?') x @populations) . ')'
);

ROW:
for my $row (@$rows) {
    next if $pruned{$row->{_SNP_KEY}};
    last if $max_leads > 0 && @selected >= $max_leads;
    my $lead_rank = @selected + 1;
    push @selected, $row;

    $presence->execute($row->{_SNP_KEY}, @populations);
    my %present = map { $_->[0] => 1 } @{$presence->fetchall_arrayref};
    my $query_status = keys(%present) == @populations ? 'OK_LOCAL_CACHE'
                     : keys(%present) ? 'PARTIAL_LOCAL_CACHE'
                     : 'NO_LD_RESPONSE';
    my $method = keys(%present) ? 'LD' : 'DISTANCE_FALLBACK';
    $audit{$row->{_SNP_KEY}} = {
        action => 'SELECTED_LEAD', lead => $row->{_SNP_KEY}, lead_rank => $lead_rank,
        query_status => $query_status, method => $method, population => '', r2 => '',
    };

    if (keys %present) {
        $lookup->execute($row->{_SNP_KEY}, $r2_threshold, @populations);
        my %proxy;
        while (my $edge = $lookup->fetchrow_arrayref) {
            my ($proxy, $population, $r2) = @$edge;
            next if $proxy eq $row->{_SNP_KEY};
            $proxy{$proxy}{populations}{$population} = 1;
            if (!defined $proxy{$proxy}{r2} || $r2 > $proxy{$proxy}{r2}) {
                $proxy{$proxy}{r2} = $r2;
                $proxy{$proxy}{population} = $population;
            }
        }
        for my $proxy_snp (keys %proxy) {
            next if $population_rule eq 'ALL'
                && keys(%{$proxy{$proxy_snp}{populations}}) < @populations;
            next unless exists $row->{_candidate_by_snp}{$proxy_snp};
            my $proxy_row = $row->{_candidate_by_snp}{$proxy_snp};
            next if $proxy_row->{_RANK} <= $row->{_RANK} || $pruned{$proxy_snp};
            $pruned{$proxy_snp} = 1;
            $audit{$proxy_snp} = {
                action => 'PRUNED_LD', lead => $row->{_SNP_KEY}, lead_rank => $lead_rank,
                query_status => '', method => 'LD',
                population => keys(%{$proxy{$proxy_snp}{populations}}) > 1 ? 'MULTI' : $proxy{$proxy_snp}{population},
                r2 => $proxy{$proxy_snp}{r2},
            };
        }
    }
    elsif ($fallback_distance_bp > 0) {
        my $half = $fallback_distance_bp / 2;
        for my $candidate (@$rows) {
            next if $candidate->{_RANK} <= $row->{_RANK} || $pruned{$candidate->{_SNP_KEY}};
            next unless $candidate->{_CHR_KEY} eq $row->{_CHR_KEY};
            next unless abs($candidate->{_BP_NUM} - $row->{_BP_NUM}) <= $half;
            $pruned{$candidate->{_SNP_KEY}} = 1;
            $audit{$candidate->{_SNP_KEY}} = {
                action => 'PRUNED_DISTANCE_FALLBACK', lead => $row->{_SNP_KEY}, lead_rank => $lead_rank,
                query_status => '', method => 'DISTANCE_FALLBACK', population => 'NA', r2 => '',
            };
        }
    }
}

for my $row (@$rows) {
    $audit{$row->{_SNP_KEY}} //= {
        action => 'NOT_REACHED', lead => '', lead_rank => '', query_status => '',
        method => '', population => '', r2 => '',
    };
}

write_leads($header, \@selected);
write_audit($rows, \%audit);
my $elapsed = time - $started;
my $pruned_ld = scalar grep { $audit{$_}{action} eq 'PRUNED_LD' } keys %audit;
my $pruned_distance = scalar grep { $audit{$_}{action} eq 'PRUNED_DISTANCE_FALLBACK' } keys %audit;
print "OUTPUT_LEADS\t$output_leads\nOUTPUT_AUDIT\t$output_audit\n";
print "CANDIDATES\t", scalar(@$rows), "\nSELECTED_LEADS\t", scalar(@selected), "\n";
print "PRUNED_LD\t$pruned_ld\nPRUNED_DISTANCE_FALLBACK\t$pruned_distance\n";
printf "SQLITE_PREP_SECONDS\t%.3f\nELAPSED_SECONDS\t%.3f\n", $db_seconds, $elapsed;

sub load_candidate_rows {
    open my $fh, '<:raw', $candidates or die "Cannot read $candidates: $!\n";
    my $csv = Text::CSV->new({ binary => 1, auto_diag => 2 });
    my $header = $csv->getline($fh) or die "Empty candidate CSV: $candidates\n";
    my %idx = map { uc($header->[$_]) => $_ } 0 .. $#$header;
    for my $required ($snp_column, $chr_column, $bp_column) {
        die "Candidate CSV lacks $required\n" unless exists $idx{uc $required};
    }
    my @signals = length($signal_column) ? ($signal_column) : split /[\s,]+/, $signal_columns;
    for my $signal (@signals) {
        die "Candidate CSV lacks signal column $signal\n" unless exists $idx{uc $signal};
    }
    my (@data, %seen);
    while (my $values = $csv->getline($fh)) {
        my $snp = uc(trim($values->[$idx{uc $snp_column}]));
        next unless $snp =~ /\ARS\d+\z/ || length $snp;
        next if $seen{$snp}++;
        my @p = grep { defined } map { numeric($values->[$idx{uc $_}]) } @signals;
        next unless @p;
        my $signal = $p[0];
        $signal = $_ < $signal ? $_ : $signal for @p;
        next unless $signal > 0 && $signal <= $signal_threshold;
        my $bp = numeric($values->[$idx{uc $bp_column}]);
        next unless defined $bp;
        push @data, {
            _VALUES => $values, _SNP_KEY => $snp,
            _CHR_KEY => uc(trim($values->[$idx{uc $chr_column}])),
            _BP_NUM => $bp, _SIGNAL => $signal,
        };
    }
    close $fh;
    @data = sort { $a->{_SIGNAL} <=> $b->{_SIGNAL} || $a->{_CHR_KEY} cmp $b->{_CHR_KEY}
                   || $a->{_BP_NUM} <=> $b->{_BP_NUM} || $a->{_SNP_KEY} cmp $b->{_SNP_KEY} } @data;
    my %by_snp;
    for my $i (0 .. $#data) {
        $data[$i]{_RANK} = $i + 1;
        $by_snp{$data[$i]{_SNP_KEY}} = $data[$i];
    }
    $_->{_candidate_by_snp} = \%by_snp for @data;
    return (\@data, $header);
}

sub prepare_cache_database {
    my $source_size = -s $cache;
    my $source_mtime = (stat $cache)[9];
    my $dbh = DBI->connect("dbi:SQLite:dbname=$sqlite", '', '', { RaiseError => 1, AutoCommit => 1 });
    $dbh->do('PRAGMA journal_mode=OFF');
    $dbh->do('PRAGMA synchronous=OFF');
    $dbh->do('PRAGMA temp_store=MEMORY');
    my $reuse = 0;
    if (!$rebuild_sqlite) {
        my ($has_meta) = $dbh->selectrow_array("SELECT count(*) FROM sqlite_master WHERE type='table' AND name='cache_meta'");
        if ($has_meta) {
            my $meta = $dbh->selectrow_hashref('SELECT source_size, source_mtime FROM cache_meta LIMIT 1');
            $reuse = $meta && $meta->{source_size} == $source_size && $meta->{source_mtime} == $source_mtime;
        }
    }
    return $dbh if $reuse;

    $dbh->do('DROP TABLE IF EXISTS edges');
    $dbh->do('DROP TABLE IF EXISTS cache_meta');
    $dbh->do('CREATE TABLE edges(query_snp TEXT NOT NULL, proxy_snp TEXT NOT NULL, ld_population TEXT NOT NULL, proxy_r2 REAL NOT NULL)');
    $dbh->do('CREATE TABLE cache_meta(source_size INTEGER, source_mtime INTEGER)');
    open my $fh, '<:raw', $cache or die "Cannot read $cache: $!\n";
    my $header = <$fh>;
    die "Empty LD cache: $cache\n" unless defined $header;
    my $insert = $dbh->prepare('INSERT INTO edges VALUES(?,?,?,?)');
    $dbh->begin_work;
    while (my $line = <$fh>) {
        $line =~ s/[\r\n]+\z//;
        my ($query, $proxy, $population, $r2) = split /\t/, $line, 5;
        next unless defined $r2;
        $insert->execute($query, $proxy, $population, $r2);
    }
    close $fh;
    $dbh->commit;
    $dbh->do('CREATE INDEX edges_query_population ON edges(query_snp,ld_population,proxy_r2)');
    $dbh->do('INSERT INTO cache_meta VALUES(?,?)', undef, $source_size, $source_mtime);
    return $dbh;
}

sub write_leads {
    my ($header, $selected) = @_;
    open my $fh, '>:raw', $output_leads or die "Cannot write $output_leads: $!\n";
    my $csv = Text::CSV->new({ binary => 1, eol => "\n" });
    my %header_index = map { uc($header->[$_]) => $_ } 0 .. $#$header;
    $csv->print($fh, [@$header, qw(INDEPENDENCE_METHOD LD_POPULATIONS LD_POPULATION_RULE LD_QUERY_STATUS LD_LEAD_RANK LD_R2_THRESHOLD LD_FALLBACK_BP)]);
    for my $i (0 .. $#$selected) {
        my $a = $audit{$selected->[$i]{_SNP_KEY}};
        my @values = @{$selected->[$i]{_VALUES}};
        $values[$header_index{FOCUS_SIGNAL}] = $selected->[$i]{_SIGNAL} if exists $header_index{FOCUS_SIGNAL};
        $csv->print($fh, [@values, $a->{method}, join(' ', @populations), $population_rule,
                          $a->{query_status}, $i + 1, $r2_threshold, $fallback_distance_bp]);
    }
    close $fh;
}

sub write_audit {
    my ($rows, $audit) = @_;
    open my $fh, '>:raw', $output_audit or die "Cannot write $output_audit: $!\n";
    print {$fh} join("\t", qw(CHR BP candidate_snp lead_snp selection_action independence_method query_status ld_population candidate_rank signal lead_rank prune_r2)), "\n";
    for my $row (@$rows) {
        my $a = $audit->{$row->{_SNP_KEY}};
        print {$fh} join("\t", $row->{_CHR_KEY}, $row->{_BP_NUM}, $row->{_SNP_KEY}, $a->{lead},
            $a->{action}, $a->{method}, $a->{query_status}, $a->{population}, $row->{_RANK},
            $row->{_SIGNAL}, $a->{lead_rank}, $a->{r2}), "\n";
    }
    close $fh;
}

sub numeric {
    my ($value) = @_;
    return undef unless defined $value;
    $value = trim($value);
    return undef unless $value =~ /\A[+-]?(?:\d+(?:\.\d*)?|\.\d+)(?:[Ee][+-]?\d+)?\z/;
    return 0 + $value;
}
sub trim { my ($v) = @_; $v //= ''; $v =~ s/^\s+|\s+$//g; return $v; }
sub usage {
    return <<'USAGE';
Usage: clump_top_hits_with_haploreg_cache.pl --candidates candidates.csv
  --cache candidate_ld.tsv --sqlite candidate_ld.sqlite
  --signal-column DIFF_P | --signal-columns "P1 P2 ..."
  --signal-threshold 5e-8 --populations "EUR ASN" --r2-threshold 0.2
  --output-leads leads.csv --output-audit audit.tsv
USAGE
}
