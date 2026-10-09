#!/usr/bin/env perl
use strict;
use warnings;
use Getopt::Long qw(GetOptions);
use File::Path qw(make_path);
use File::Basename qw(basename);
use File::Copy qw(copy);
use JSON::PP;
use Digest::SHA;

my ($dir, $existing_manifest);
GetOptions('output-dir=s'=>\$dir, 'existing-manifest'=>\$existing_manifest)
  && $dir or die "Usage: perl $0 --output-dir DIR [--existing-manifest]\n";
my @manifest;
my %order=(manhattan=>0,local_manhattan=>1,local_gtf=>2,forest=>3);
if ($existing_manifest) {
    my $path="$dir/results_gallery_manifest.json";
    open my $fh,'<',$path or die "$path: $!";
    my $data=decode_json(do {local $/;<$fh>}); close $fh;
    @manifest=@{$data->{images} || []};
    die "No images in $path\n" unless @manifest;
    for my $im (@manifest) {
        my $copy="$dir/$im->{copy}";
        die "Gallery image missing: $copy\n" unless -s $copy;
        die "Gallery checksum mismatch: $copy\n" unless checksum($copy) eq $im->{sha256};
    }
} else {
    for my $backend (qw(gnuplot sas)) {
        my $report="$dir/image_validation_$backend.json";
        next unless -f $report;
        open my $fh,'<',$report or die "$report: $!";
        my $data=decode_json(do {local $/;<$fh>}); close $fh;
        die "Image validation did not pass: $report\n" unless $data->{status} eq 'PASS';
        my $folder="figures/$backend";
        make_path("$dir/$folder");
        my $label=$backend eq 'sas' ? 'SAS ODA' : 'gnuplot';
        for my $im (@{$data->{images}}) {
            my $name=basename($im->{path});
            copy($im->{path},"$dir/$folder/$name") or die "Copy $im->{path}: $!\n";
            my $sha=checksum($im->{path});
            die "Copied image differs from source: $name\n"
              unless checksum("$dir/$folder/$name") eq $sha;
            push @manifest,{backend=>$label,family=>$im->{family},source=>$im->{path},
                            copy=>"$folder/$name",sha256=>$sha};
        }
    }
    die "No image validation reports found in $dir\n" unless @manifest;
    open my $mf,'>',"$dir/results_gallery_manifest.json" or die $!;
    print {$mf} JSON::PP->new->canonical->pretty->encode({images=>\@manifest});
    close $mf or die $!;
}

my @sections;
for my $backend ('SAS ODA','gnuplot') {
    my @images=sort {($order{$a->{family}}//9)<=>($order{$b->{family}}//9)
                     || $a->{copy} cmp $b->{copy}}
               grep {$_->{backend} eq $backend} @manifest;
    next unless @images;
    my $section=qq{<h2 id="}.($backend eq 'SAS ODA' ? 'sas' : 'gnuplot').qq{">$backend</h2>\n};
    for my $im (@images) {
        my $name=basename($im->{copy});
        my $family=$im->{family};
        my $caption=caption($backend,$family,$name);
        my $url=esc($im->{copy});
        $section.='<figure class="'.esc($family).'"><figcaption><strong>'.esc($caption)
          .'</strong></figcaption><a href="'.$url.'"><img loading="lazy" src="'.$url
          .'" alt="'.esc($caption).'"></a></figure>'."\n";
    }
    push @sections,$section;
}
open my $out,'>',"$dir/results.html" or die $!;
print {$out} '<!doctype html><html lang="en"><meta charset="utf-8"><meta name="viewport" content="width=device-width,initial-scale=1"><title>Public schizophrenia GWAS by sex: plots</title>'
 .'<style>body{font:16px/1.5 system-ui,sans-serif;max-width:1280px;margin:0 auto;padding:1rem 1.5rem;color:#202124}nav{margin-bottom:2rem}figure{box-sizing:border-box;margin:1.5rem auto 2.5rem;padding:1rem;border:1px solid #ddd;border-radius:8px;background:#fff}figure.manhattan{max-width:1250px}figure.local_manhattan{max-width:1150px}figure.local_gtf{max-width:1000px}figure.forest{max-width:900px}figure img{display:block;max-width:100%;height:auto;margin:.7rem auto 0}figcaption{font-size:1.05rem}a{color:#1155a0}</style>'
 .'<h1>Public schizophrenia GWAS by sex</h1><nav><a href="#sas">SAS ODA figures</a> | <a href="#gnuplot">gnuplot figures</a></nav><p>Association height is -log10(P). Z-score and signed-LD labels describe point color. Click any figure for its full-resolution PNG.</p>'
 .join("\n",@sections).'</html>';
close $out or die $!;
print "Results gallery: $dir/results.html\n";

sub caption {
    my ($backend,$family,$name)=@_;
    my ($snp)=$name=~/(rs\d+)/;
    my $suffix=$snp ? " ($snp)" : '';
    $suffix .= " (exported panel $1)" if $name =~ /_part(\d+)\.png\z/;
    return "$backend: Genome-wide Manhattan" if $family eq 'manhattan';
    return "$backend: Local Manhattan, ".($backend eq 'gnuplot' ? 'Z-score colors' : 'chromosome colors').$suffix
      if $family eq 'local_manhattan';
    return "$backend: Local Manhattan and gene tracks, signed LD r^2 x sign(Z)".$suffix
      if $family eq 'local_gtf';
    return "$backend: Forest plot".($name =~ /EUR_FEMALE/ ? ' (female)' : $name =~ /EUR_MALE/ ? ' (male)' : ' (combined)')
      if $family eq 'forest';
    return "$backend: $family$suffix";
}
sub esc {my $s=shift;$s=~s/&/&amp;/g;$s=~s/</&lt;/g;$s=~s/>/&gt;/g;$s=~s/"/&quot;/g;return $s}
sub checksum {open my $fh,'<:raw',$_[0] or die $!;my $s=Digest::SHA->new(256)->addfile($fh)->hexdigest;close $fh;return $s}
