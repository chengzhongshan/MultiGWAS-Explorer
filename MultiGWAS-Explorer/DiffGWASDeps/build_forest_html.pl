#!/usr/bin/env perl
use strict;
use warnings;

use File::Basename qw(basename dirname);
use File::Spec;
use Getopt::Long qw(GetOptions);
use MIME::Base64 qw(encode_base64);

my $manifest = '';
my $output = '';
my $title = 'Top-hit forest plots';
my $image_dir = '';
my $top_hits_csv = '';

GetOptions(
    'manifest=s'     => \$manifest,
    'output=s'       => \$output,
    'title=s'        => \$title,
    'image-dir=s'    => \$image_dir,
    'top-hits-csv=s' => \$top_hits_csv,
) or die usage();

die usage() unless length($manifest) && length($output);
die "Forest manifest does not exist or is empty: $manifest\n" unless -s $manifest;
$image_dir = dirname($manifest) unless length $image_dir;

my ($fallback_snp, $fallback_gene) = read_first_top_hit($top_hits_csv);
my @panels = read_manifest($manifest);
die "Forest manifest contains no panel rows: $manifest\n" unless @panels;

for my $panel (@panels) {
    my $image_name = trim($panel->{png_file});
    die "Forest manifest has an empty png_file value\n" unless length $image_name;
    die "Forest manifest png_file must be a basename, not a path: $image_name\n"
        unless basename($image_name) eq $image_name;
    my $image_path = File::Spec->catfile($image_dir, $image_name);
    die "Forest panel image does not exist or is empty: $image_path\n"
        unless -s $image_path;

    my $label = trim($panel->{track_label});
    if (!length($label) && lc(trim($panel->{track_id})) eq 'single_snp') {
        $label = $fallback_snp;
        $label .= " ($fallback_gene)"
            if length($fallback_gene) && uc($fallback_gene) ne 'NA';
    }
    $label = trim($panel->{track_id}) unless length $label;
    $label = 'Forest plot' unless length $label;

    my $mime = $image_name =~ /\.svg$/i ? 'image/svg+xml' : 'image/png';
    $panel->{label} = $label;
    $panel->{data_uri} = "data:$mime;base64," . encode_base64(slurp_binary($image_path), '');
}

my $escaped_title = html_escape($title);
my $html = <<"HTML";
<!doctype html>
<html lang="en"><head><meta charset="utf-8">
<meta name="viewport" content="width=device-width,initial-scale=1">
<title>$escaped_title</title>
<style>
body{margin:0;padding:24px;font-family:Arial,sans-serif;background:#fff;color:#1f2937;}
h1{font-size:24px;line-height:1.25;margin:0 0 20px;}
.grid{display:grid;grid-template-columns:repeat(auto-fit,minmax(360px,1fr));gap:24px;align-items:start;}
.panel{margin:0;padding:14px;border:1px solid #d9dee7;border-radius:10px;background:#fff;box-shadow:0 1px 4px rgba(0,0,0,.06);}
.panel figcaption{font-size:18px;font-weight:700;margin-bottom:10px;text-align:center;}
.panel img{width:100%;height:auto;display:block;}
</style></head><body>
<h1>$escaped_title</h1>
<div class="grid">
HTML

for my $panel (@panels) {
    my $label = html_escape($panel->{label});
    my $uri = $panel->{data_uri};
    $html .= qq{<figure class="panel"><figcaption>$label</figcaption><img src="$uri" alt="$label forest plot"></figure>\n};
}
$html .= "</div></body></html>\n";

my $tmp = "$output.tmp.$$";
open my $out, '>:raw', $tmp or die "Cannot write temporary forest HTML $tmp: $!\n";
print {$out} $html;
close $out or die "Cannot close temporary forest HTML $tmp: $!\n";
unlink $output if -e $output;
rename $tmp, $output or die "Cannot install forest HTML $output: $!\n";
die "Generated forest HTML is empty: $output\n" unless -s $output;

print "Self-contained forest HTML: $output\n";
print "Embedded forest panels: " . scalar(@panels) . "\n";

sub read_manifest {
    my ($path) = @_;
    open my $fh, '<:encoding(UTF-8)', $path
        or die "Cannot read forest manifest $path: $!\n";
    my $header = <$fh>;
    die "Forest manifest has no header: $path\n" unless defined $header;
    chomp $header;
    my @header = map { trim($_) } split /\t/, $header, -1;
    my @rows;
    while (my $line = <$fh>) {
        chomp $line;
        next unless $line =~ /\S/;
        my @fields = split /\t/, $line, -1;
        my %row;
        for my $i (0 .. $#header) {
            $row{$header[$i]} = $fields[$i] // '';
        }
        push @rows, \%row if trim($row{png_file} // '') =~ /\.(?:png|svg)$/i;
    }
    close $fh;
    return @rows;
}

sub read_first_top_hit {
    my ($path) = @_;
    return ('', '') unless defined $path && length($path) && -s $path;
    open my $fh, '<:encoding(UTF-8)', $path
        or die "Cannot read top-hit CSV $path: $!\n";
    my $header_line = <$fh>;
    my $data_line = <$fh>;
    close $fh;
    return ('', '') unless defined $header_line && defined $data_line;
    chomp($header_line, $data_line);
    my @header = parse_csv_line($header_line);
    my @data = parse_csv_line($data_line);
    my %index;
    @index{@header} = (0 .. $#header);
    my $snp = exists($index{SNP}) ? trim($data[$index{SNP}] // '') : '';
    my $gene = exists($index{gene}) ? trim($data[$index{gene}] // '') : '';
    return ($snp, $gene);
}

sub parse_csv_line {
    my ($line) = @_;
    my @fields;
    my $field = '';
    my $quoted = 0;
    my @chars = split //, ($line // '');
    for (my $i = 0; $i < @chars; $i++) {
        my $ch = $chars[$i];
        if ($quoted) {
            if ($ch eq '"') {
                if ($i + 1 < @chars && $chars[$i + 1] eq '"') {
                    $field .= '"';
                    $i++;
                }
                else {
                    $quoted = 0;
                }
            }
            else {
                $field .= $ch;
            }
        }
        elsif ($ch eq '"') {
            $quoted = 1;
        }
        elsif ($ch eq ',') {
            push @fields, $field;
            $field = '';
        }
        else {
            $field .= $ch;
        }
    }
    push @fields, $field;
    return @fields;
}

sub slurp_binary {
    my ($path) = @_;
    open my $fh, '<:raw', $path or die "Cannot read forest image $path: $!\n";
    local $/;
    my $data = <$fh>;
    close $fh;
    return defined($data) ? $data : '';
}

sub trim {
    my ($value) = @_;
    $value = '' unless defined $value;
    $value =~ s/^\s+|\s+$//g;
    return $value;
}

sub html_escape {
    my ($value) = @_;
    $value = '' unless defined $value;
    $value =~ s/&/&amp;/g;
    $value =~ s/</&lt;/g;
    $value =~ s/>/&gt;/g;
    $value =~ s/"/&quot;/g;
    $value =~ s/'/&#39;/g;
    return $value;
}

sub usage {
    return <<'USAGE';
Usage:
  perl build_forest_html.pl --manifest forest.manifest.tsv --output forest.html
       [--title TEXT] [--image-dir DIR] [--top-hits-csv FILE]

Builds a self-contained forest HTML page by embedding every manifest PNG/SVG.
For a blank single-SNP manifest label, SNP and gene are read from the first
row of --top-hits-csv.
USAGE
}
