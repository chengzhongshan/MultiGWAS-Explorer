#!/usr/bin/env perl
use strict;
use warnings;
use Getopt::Long qw(GetOptions);
use File::Path qw(make_path);
use File::Basename qw(basename);
use File::Copy qw(copy);
use JSON::PP;
use Digest::SHA;
my $dir;
GetOptions('output-dir=s'=>\$dir) && $dir or die "Usage: perl $0 --output-dir DIR\n";
my (@sections,@manifest);
my %titles=(manhattan=>'Genome-wide Manhattan',local_manhattan=>'Local Manhattan',local_gtf=>'Local gene tracks',forest=>'Forest plot');
my %order=(manhattan=>0,local_manhattan=>1,local_gtf=>2,forest=>3);
for my $backend (qw(gunplot sas)) {
    my $report="$dir/image_validation_$backend.json";
    next unless -f $report;
    open my $fh,'<',$report or die "$report: $!";
    my $data=decode_json(do {local $/;<$fh>}); close $fh;
    die "Image validation did not pass: $report\n" unless $data->{status} eq 'PASS';
    my $folder="figures/$backend";
    make_path("$dir/$folder");
    my $label=$backend eq 'sas' ? 'SAS ODA' : 'gnuplot';
    my $section=qq{<h2 id="$backend">$label</h2>\n};
    for my $im (sort {($order{$a->{family}}//9)<=>($order{$b->{family}}//9) || $a->{path} cmp $b->{path}} @{$data->{images}}) {
        my $name=basename($im->{path});
        copy($im->{path},"$dir/$folder/$name") or die "Copy $im->{path}: $!\n";
        my $sha=checksum($im->{path});
        die "Copied image differs from source: $name\n" unless checksum("$dir/$folder/$name") eq $sha;
        push @manifest,{backend=>$label,family=>$im->{family},source=>$im->{path},copy=>"$folder/$name",sha256=>$sha};
        my $url=esc("$folder/$name");
        $section.='<figure><figcaption><strong>'.esc($label.' - '.($titles{$im->{family}}//$im->{family})).'</strong><br>'.esc($name).'</figcaption>'
          .qq{<a href="$url"><img loading="lazy" src="$url" alt="}.esc($name).qq{"></a></figure>\n};
    }
    push @sections,$section;
}
die "No image validation reports found in $dir\n" unless @sections;
open my $out,'>',"$dir/results.html" or die $!;
print {$out} '<!doctype html><html lang="en"><meta charset="utf-8"><title>Public GWAS test results</title>'
 .'<style>body{font-family:sans-serif;margin:2rem}img{max-width:100%;height:auto}figure{margin:2rem 0}figcaption{overflow-wrap:anywhere}</style>'
 .'<h1>Public GWAS test results</h1><nav><a href="#sas">SAS ODA figures</a> | <a href="#gunplot">gnuplot figures</a></nav><p>Click a figure to open it at full size. Copied files are verified against their source with SHA-256.</p>'
 .join("\n",@sections).'</html>';
close $out or die $!;
open my $mf,'>',"$dir/results_gallery_manifest.json" or die $!;
print {$mf} JSON::PP->new->canonical->pretty->encode({images=>\@manifest});
close $mf or die $!;
print "Results gallery: $dir/results.html\n";
sub esc {my $s=shift;$s=~s/&/&amp;/g;$s=~s/</&lt;/g;$s=~s/>/&gt;/g;$s=~s/"/&quot;/g;return $s}
sub checksum {open my $fh,'<:raw',$_[0] or die $!;my $s=Digest::SHA->new(256)->addfile($fh)->hexdigest;close $fh;return $s}
