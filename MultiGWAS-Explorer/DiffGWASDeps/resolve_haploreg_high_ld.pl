#!/usr/bin/env perl
use strict;
use warnings;

use Fcntl qw(:flock);
use Getopt::Long qw(GetOptions);
use HTTP::Tiny;
use File::Basename qw(dirname);
use File::Path qw(make_path);

my ($query_snps, $population, $min_r2, $web_cache, $output);
my @local_cache;
my $web_fallback = 1;
$population = 'EUR';
$min_r2 = 0.8;

GetOptions(
    'query-snps=s'    => \$query_snps,
    'population=s'    => \$population,
    'min-r2=f'        => \$min_r2,
    'local-cache=s@'  => \@local_cache,
    'web-cache=s'     => \$web_cache,
    'web-fallback!'   => \$web_fallback,
    'output=s'        => \$output,
) or die usage();

die usage() unless defined($query_snps) && length($query_snps);
$population = uc($population // 'EUR');
die "--population must be AFR, AMR, ASN, or EUR\n"
    unless $population =~ /\A(?:AFR|AMR|ASN|EUR)\z/;
die "--min-r2 must be between 0 and 1\n"
    unless $min_r2 > 0 && $min_r2 <= 1;

my @queries;
my %is_query;
for my $snp (split /[,\s]+/, $query_snps) {
    next unless $snp =~ /\Ars\d+\z/i;
    my $key = uc $snp;
    next if $is_query{$key}++;
    push @queries, $key;
}
die "No valid rsIDs were supplied with --query-snps\n" unless @queries;

push @local_cache, $web_cache if defined($web_cache) && length($web_cache);
my (%proxy_for, %source_for);
for my $cache (@local_cache) {
    next unless defined($cache) && -s $cache;
    load_normalized_cache(
        path       => $cache,
        population => $population,
        min_r2     => $min_r2,
        queries    => \%is_query,
        proxy_for  => \%proxy_for,
        source_for => \%source_for,
    );
}

my @missing = grep { !keys %{ $proxy_for{$_} || {} } } @queries;
my @new_rows;
if (@missing && $web_fallback) {
    for my $query (@missing) {
        my @rows = query_haploreg_web(
            query      => $query,
            population => $population,
            min_r2     => $min_r2,
        );
        for my $row (@rows) {
            my ($proxy, $r2, $dprime) = @{$row};
            next if $proxy eq $query;
            $proxy_for{$query}{$proxy} = $r2
                if !exists($proxy_for{$query}{$proxy}) || $r2 > $proxy_for{$query}{$proxy};
            $source_for{$query} = 'HAPLOREG_WEB';
            push @new_rows, [$query, $proxy, $population, $r2, $dprime, $min_r2];
        }
    }
    append_web_cache($web_cache, \@new_rows)
        if @new_rows && defined($web_cache) && length($web_cache);
}

my (%seen, @proxies);
for my $query (@queries) {
    for my $proxy (sort {
        ($proxy_for{$query}{$b} <=> $proxy_for{$query}{$a}) || ($a cmp $b)
    } keys %{ $proxy_for{$query} || {} }) {
        next if $is_query{$proxy};
        next if $seen{$proxy}++;
        push @proxies, $proxy;
    }
}

if (defined($output) && length($output)) {
    my $dir = dirname($output);
    make_path($dir) if length($dir) && !-d $dir;
    open my $out, '>:raw', $output or die "Cannot write $output: $!\n";
    print {$out} join("\t", qw(query_snp proxy_snp ld_population proxy_r2 source)), "\n";
    for my $query (@queries) {
        for my $proxy (sort {
            ($proxy_for{$query}{$b} <=> $proxy_for{$query}{$a}) || ($a cmp $b)
        } keys %{ $proxy_for{$query} || {} }) {
            print {$out} join("\t", $query, $proxy, $population,
                $proxy_for{$query}{$proxy}, ($source_for{$query} || 'LOCAL_CACHE')), "\n";
        }
    }
    close $out or die "Cannot close $output: $!\n";
}

print "LD_SNPS\t", join(',', @proxies), "\n";
print "LD_POPULATION\t$population\n";
print "LD_MIN_R2\t$min_r2\n";
for my $query (@queries) {
    my $count = scalar keys %{ $proxy_for{$query} || {} };
    print "QUERY\t$query\t$count\t", ($source_for{$query} || 'NONE'), "\n";
}

sub load_normalized_cache {
    my (%args) = @_;
    my $path = $args{path};
    if ($path =~ /\.sqlite\z/i) {
        require DBI;
        my $dbh = DBI->connect("dbi:SQLite:dbname=$path", '', '',
            { RaiseError => 1, AutoCommit => 1 });
        my $sth = $dbh->prepare(
            'SELECT proxy_snp, proxy_r2 FROM edges '
          . 'WHERE query_snp=? AND ld_population=? AND proxy_r2>=?'
        );
        for my $query (keys %{ $args{queries} }) {
            $sth->execute($query, $args{population}, $args{min_r2});
            while (my ($proxy, $r2) = $sth->fetchrow_array) {
                $proxy = uc($proxy // '');
                next unless $proxy =~ /\ARS\d+\z/;
                $args{proxy_for}{$query}{$proxy} = 0 + $r2;
                $args{source_for}{$query} = 'LOCAL_SQLITE';
            }
        }
        $dbh->disconnect;
        return;
    }

    open my $fh, '<:raw', $path or die "Cannot read LD cache $path: $!\n";
    my $header = <$fh> // '';
    $header =~ s/[\r\n]+\z//;
    my @cols = split /\t/, $header, -1;
    my %idx = map { lc($cols[$_]) => $_ } 0 .. $#cols;
    unless (exists($idx{query_snp}) && exists($idx{proxy_snp})
        && exists($idx{proxy_r2})) {
        warn "WARNING: Ignoring non-normalized LD cache $path; expected query_snp, proxy_snp, and proxy_r2 columns.\n";
        close $fh;
        return;
    }
    while (my $line = <$fh>) {
        $line =~ s/[\r\n]+\z//;
        my @f = split /\t/, $line, -1;
        my $query = uc($f[$idx{query_snp}] // '');
        next unless $args{queries}{$query};
        if (exists($idx{ld_population})) {
            next unless uc($f[$idx{ld_population}] // '') eq $args{population};
        }
        my $r2 = $f[$idx{proxy_r2}] // '';
        next unless $r2 =~ /\A(?:\d+(?:\.\d*)?|\.\d+)(?:e[+-]?\d+)?\z/i;
        next if $r2 < $args{min_r2};
        my $proxy = uc($f[$idx{proxy_snp}] // '');
        next unless $proxy =~ /\ARS\d+\z/;
        $args{proxy_for}{$query}{$proxy} = 0 + $r2;
        $args{source_for}{$query} = 'LOCAL_TSV';
    }
    close $fh;
}

sub query_haploreg_web {
    my (%args) = @_;
    my $content = join('&',
        'query=' . lc($args{query}),
        'ldThresh=' . $args{min_r2},
        'ldPop=' . $args{population},
        'output=text',
        'submit=submit',
    );
    my $http = HTTP::Tiny->new(timeout => 45, verify_SSL => 1);
    my $res = $http->post(
        'https://pubs.broadinstitute.org/mammals/haploreg/haploreg.php/post',
        {
            headers => {'content-type' => 'application/x-www-form-urlencoded'},
            content => $content,
        },
    );
    unless ($res->{success}) {
        warn "WARNING: HaploReg web fallback failed for $args{query} "
           . "(HTTP " . ($res->{status} // 'unknown') . ").\n";
        return;
    }
    my @lines = split /\r?\n/, ($res->{content} // '');
    my $header = shift @lines // '';
    my @cols = split /\t/, $header, -1;
    my %idx = map { lc($cols[$_]) => $_ } 0 .. $#cols;
    unless (exists($idx{rsid}) && exists($idx{r2})) {
        warn "WARNING: HaploReg response for $args{query} lacks rsID/r2 columns.\n";
        return;
    }
    my @rows;
    for my $line (@lines) {
        next unless length $line;
        my @f = split /\t/, $line, -1;
        my $proxy = uc($f[$idx{rsid}] // '');
        my $r2 = $f[$idx{r2}] // '';
        next unless $proxy =~ /\ARS\d+\z/;
        next unless $r2 =~ /\A(?:\d+(?:\.\d*)?|\.\d+)(?:e[+-]?\d+)?\z/i;
        next if $r2 < $args{min_r2};
        my $dprime = exists($idx{"d'"}) ? ($f[$idx{"d'"}] // '') : '';
        push @rows, [$proxy, 0 + $r2, $dprime];
    }
    return @rows;
}

sub append_web_cache {
    my ($path, $rows) = @_;
    my $dir = dirname($path);
    make_path($dir) if length($dir) && !-d $dir;
    my $exists = -s $path;
    open my $fh, '>>:raw', $path or die "Cannot append HaploReg web cache $path: $!\n";
    flock($fh, LOCK_EX) or die "Cannot lock HaploReg web cache $path: $!\n";
    print {$fh} join("\t", qw(query_snp proxy_snp ld_population proxy_r2 proxy_dprime cache_min_r2)), "\n"
        unless $exists;
    my %written;
    for my $row (@{$rows}) {
        my $key = join("\t", @{$row}[0,1,2]);
        next if $written{$key}++;
        print {$fh} join("\t", @{$row}), "\n";
    }
    close $fh or die "Cannot close HaploReg web cache $path: $!\n";
}

sub usage {
    return <<'USAGE';
Usage:
  perl resolve_haploreg_high_ld.pl --query-snps rs1,rs2 [options]

Options:
  --population EUR        HaploReg population: AFR, AMR, ASN, or EUR (default EUR)
  --min-r2 0.8            High-LD threshold (default 0.8)
  --local-cache FILE      Normalized HaploReg TSV or SQLite edges cache; repeatable
  --web-cache FILE        Persistent normalized cache for successful web fallbacks
  --[no-]web-fallback     Query HaploReg when local caches lack a query (default on)
  --output FILE           Optional resolved query/proxy audit TSV
USAGE
}
