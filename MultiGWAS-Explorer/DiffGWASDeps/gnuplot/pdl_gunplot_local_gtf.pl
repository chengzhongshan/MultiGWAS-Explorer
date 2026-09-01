#!/usr/bin/env perl
use strict;
use warnings;
use Getopt::Long;
use IO::Uncompress::Gunzip qw(gunzip $GunzipError);
use File::Spec;
use File::Basename;

sub usage {
    print <<'USAGE';
Usage: pdl_gunplot_local_gtf.pl --data file.tsv.gz --snp rsID --outdir out [--window 1000000] [--gtf path.gtf]

Generates a local locus plot with gene-track using gnuplot. If --gtf is omitted
gene track will be skipped unless gene columns are present.
USAGE
    exit 1;
}

my %opts;
GetOptions(
    'data=s'  => \$opts{data},
    'snp=s'   => \$opts{snp},
    'outdir=s'=> \$opts{outdir},
    'window=i'=> \$opts{window},
    'gtf=s'   => \$opts{gtf},
) or usage();
usage() unless $opts{data} && $opts{snp};
$opts{outdir} ||= '.';
mkdir $opts{outdir} unless -d $opts{outdir};
$opts{window} ||= 10000000;

my $data = $opts{data};
my $gz = IO::Uncompress::Gunzip->new($data) or die "gunzip failed: $GunzipError";
my $hdr = <$gz>;
chomp $hdr;
my @cols = split /\t/, $hdr;
my %col_idx = map { $cols[$_] => $_ } 0..$#cols;
die "Missing CHR/BP or SNP columns in input" unless exists $col_idx{'CHR'} && exists $col_idx{'BP'} && exists $col_idx{'SNP'};

my ($target_chr, $target_bp);
my @records;
while (my $line = <$gz>) {
    chomp $line;
    next unless $line;
    my @f = split /\t/, $line, -1;
    push @records, \@f;
    if ($f[$col_idx{'SNP'}] && $f[$col_idx{'SNP'}] eq $opts{snp}) {
        $target_chr = $f[$col_idx{'CHR'}];
        $target_bp  = $f[$col_idx{'BP'}] + 0;
    }
}
close $gz;
die "Target SNP $opts{snp} not found in data" unless $target_chr && $target_bp;

my $win = $opts{window};
my $start = $target_bp - $win;
my $end   = $target_bp + $win;

my @locus_rows = grep { $_->[$col_idx{'CHR'}] eq $target_chr && $_->[$col_idx{'BP'}] >= $start && $_->[$col_idx{'BP'}] <= $end } @records;
die "No data in locus window" unless @locus_rows;

# detect P columns similarly to the Manhattan helper (broad match)
my @pcols = grep { $_ !~ /^(CHR|CHROM|BP|POS|SNP|ID)$/i && $_ =~ /p/i } @cols;

my $outpng = File::Spec->catfile($opts{outdir}, "$opts{snp}_local.png");
my $gp = File::Spec->catfile($opts{outdir}, "$opts{snp}_local.gp");
open my $gpfh, '>', $gp or die "open $gp: $!";
print $gpfh "set terminal pngcairo size 1200,800 enhanced font 'Arial,10'\n";
print $gpfh "set output '$outpng'\n";

# prepare data file
my $dat = File::Spec->catfile($opts{outdir}, "$opts{snp}_locus.dat");
open my $dfh, '>', $dat or die "open $dat: $!";
my $maxy = 0;
for my $r (@locus_rows) {
    my ($chr,$bp,@rest) = @$r;
    for my $pc (@pcols) {
        my $pidx = $col_idx{$pc};
        my $p = $r->[$pidx];
        next unless defined $p && $p ne '' && $p > 0;
        my $y = -log($p)/log(10);
        $maxy = $y if $y > $maxy;
        print $dfh "$bp\t$y\t$pc\t$r->[$col_idx{'SNP'}]\n";
    }
}
close $dfh;
close $dfh;

# if gtf provided, parse genes overlapping region
my @gene_objs;
if ($opts{gtf} && -e $opts{gtf}) {
    open my $gfh, '<', $opts{gtf} or warn "open gtf: $!";
    while (<$gfh>) {
        next if /^#/;
        chomp;
        my @g = split /\t/;
        next unless $g[2] eq 'gene';
        my $gchr = $g[0];
        my $gs = $g[3] + 0;
        my $ge = $g[4] + 0;
        next unless $gchr eq $target_chr && $ge >= $start && $gs <= $end;
        my ($name) = ($_ =~ /gene_name "([^"]+)"/);
        $name ||= ($_ =~ /gene_id "([^"]+)"/) || 'gene';
        push @gene_objs, { name=>$name, start=>$gs, end=>$ge };
    }
    close $gfh;
}

# generate gnuplot script: scatter points and gene rectangles as set object
my $yrange_top = int($maxy * 1.1) + 1;
print $gpfh "set xrange [$start:$end]\n";
print $gpfh "set xlabel 'Position (bp)' font ',11'\n";
print $gpfh "set ylabel '-log10(P)' font ',11'\n";
# allow a small negative area for gene boxes below x-axis
print $gpfh "set yrange [-0.6:$yrange_top]\n";
print $gpfh "set format y '%.0f'\n";
print $gpfh "set tics nomirror\n";
print $gpfh "set key off\n";
print $gpfh "plot '$dat' using 1:2 with points pt 7 ps 0.8 lc rgb '#2c7bb6' notitle\n";

my $label_y_base = 0;
if (@gene_objs) {
    my $count = 0;
    for my $g (@gene_objs) {
        my $gname = $g->{name};
        my $gs = $g->{start};
        my $ge = $g->{end};
        my $col = ($count % 2) ? '#cccccc' : '#bbbbbb';
        # draw gene boxes slightly below the x-axis (negative y region)
        print $gpfh "set object rect from $gs, -0.5 to $ge, -0.15 fc rgb '$col' fillstyle solid 0.8 behind\n";
        print $gpfh "set label '$gname' at ( ($gs+$ge)/2 ), -0.55 center font ',9'\n";
        $count++;
    }
}

# annotate SNP label at the bottom
print $gpfh "set label '$opts{snp}' at $target_bp, -0.58 center rotate by 90 font ',10'\n";

close $gpfh;
system('gnuplot', $gp) == 0 or warn "gnuplot failed: $?";
print "Local locus image: $outpng\n";
