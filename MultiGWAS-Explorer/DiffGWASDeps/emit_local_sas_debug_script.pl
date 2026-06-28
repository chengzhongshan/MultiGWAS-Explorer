#!/usr/bin/env perl
use strict;
use warnings;
use Getopt::Long qw(GetOptions);
use File::Basename qw(basename dirname);
use File::Path qw(make_path);

my $mode = '';
my $input = '';
my $output = '';
my $workdir = '';
my $deps_dir = '';
my $data_gz = '';
my $gtf_subset_gz = '';
my $top_hits_csv = '';
my $gtf_macro_upload = '';
my $gtf_local_dataset = 'gtf_hg38';
my $manhattan_macro = '';
my $output_html_basename = '';

GetOptions(
    'mode=s'                => \$mode,
    'input=s'               => \$input,
    'output=s'              => \$output,
    'workdir=s'             => \$workdir,
    'deps-dir=s'            => \$deps_dir,
    'data-gz=s'             => \$data_gz,
    'gtf-subset-gz=s'       => \$gtf_subset_gz,
    'top-hits-csv=s'        => \$top_hits_csv,
    'gtf-macro-upload=s'    => \$gtf_macro_upload,
    'gtf-local-dataset=s'   => \$gtf_local_dataset,
    'manhattan-macro=s'     => \$manhattan_macro,
    'output-html-basename=s'=> \$output_html_basename,
) or die usage();

die "--mode is required\n" unless length $mode;
die "--input is required\n" unless length $input;
die "--output is required\n" unless length $output;
die "--workdir is required\n" unless length $workdir;
die "--deps-dir is required\n" unless length $deps_dir;
die "--data-gz is required\n" unless length $data_gz;

my %valid_mode = map { $_ => 1 } qw(manhattan local_manhattan local_gtf);
die "Unsupported --mode '$mode'\n" unless $valid_mode{$mode};

my $text = slurp($input);

my $local_debug_root = path_to_windows($workdir);
my $local_debug_deps = path_to_windows($deps_dir);
my $local_debug_data_gz = path_to_windows($data_gz);
my $local_debug_output_html = path_to_windows(join_unix($workdir, $output_html_basename || basename($output)));

my $local_debug_top_hits_csv = length($top_hits_csv) ? path_to_windows($top_hits_csv) : '';
my $local_debug_gtf_subset_gz = length($gtf_subset_gz) ? path_to_windows($gtf_subset_gz) : '';
my $local_debug_gtf_macro = length($gtf_macro_upload) ? path_to_windows($gtf_macro_upload) : '';
my $local_debug_manhattan_macro = '';

if ($mode eq 'manhattan' || $mode eq 'local_manhattan') {
    die "--manhattan-macro is required for mode $mode\n" unless length $manhattan_macro;
    my $macro_src = slurp($manhattan_macro);
    $macro_src =~ s/filename\s+gout\s+"~\/&outputfigname\.\.png";/filename gout "&local_debug_output_dir\\&outputfigname..png";/g;
    my $macro_out = join_unix(dirname($output), basename($manhattan_macro, '.sas') . '.local_debug.sas');
    write_text($macro_out, $macro_src);
    $local_debug_manhattan_macro = path_to_windows($macro_out);
}

my $preamble = build_preamble(
    mode                     => $mode,
    local_debug_root         => $local_debug_root,
    local_debug_deps         => $local_debug_deps,
    local_debug_output_html  => $local_debug_output_html,
    local_debug_data_gz      => $local_debug_data_gz,
    local_debug_top_hits_csv => $local_debug_top_hits_csv,
    local_debug_gtf_subset_gz=> $local_debug_gtf_subset_gz,
    local_debug_gtf_macro    => $local_debug_gtf_macro,
    local_debug_manhattan_macro => $local_debug_manhattan_macro,
);

$text = apply_mode_rewrites(
    mode                     => $mode,
    text                     => $text,
    local_debug_top_hits_csv => $local_debug_top_hits_csv,
);

$text = $preamble . "\n" . $text;

write_text($output, $text);
print "$output\n";

sub apply_mode_rewrites {
    my (%args) = @_;
    my $mode = $args{mode};
    my $text = $args{text};
    my $local_debug_top_hits_csv = $args{local_debug_top_hits_csv};

    if ($mode eq 'manhattan') {
        $text =~ s/%include\s+"~\/Manhattan4DiffGWASs_png\.sas";/%include "&local_debug_manhattan_macro";/g;
        $text =~ s/filename\s+mhdata\s+zip\s+"~\/[^"]+"\s+gzip;/filename mhdata zip "&local_debug_wide_data_gz" gzip;/g;
        $text =~ s/file\s+"~\/[^"]+_png\.html"\s+lrecl=32767;/file "&local_debug_output_html" lrecl=32767;/g;
    }
    elsif ($mode eq 'local_manhattan') {
        $text =~ s/%include\s+"~\/get_top_signal_within_dist\.sas";/%include "&local_debug_deps\\get_top_signal_within_dist.sas";/g;
        $text =~ s/%include\s+"~\/[^"]*get_genecode_gtf_data[^"]*\.sas";/%include "&local_debug_gtf_macro";/g;
        $text =~ s/%include\s+"~\/Manhattan4DiffGWASs_png\.sas";/%include "&local_debug_manhattan_macro";/g;
        $text =~ s/filename\s+mhdata\s+zip\s+"~\/[^"]+"\s+gzip;/filename mhdata zip "&local_debug_wide_data_gz" gzip;/g;
        $text =~ s/outfile="~\/&local_top_hits_csv_basename"/outfile="&local_debug_top_hits_csv"/g;
        $text =~ s/file\s+"~\/[^"]+\.html"\s+lrecl=32767;/file "&local_debug_output_html" lrecl=32767;/g;
        $text =~ s/%let\s+gtf_dsd=.*?;/%let gtf_dsd=$gtf_local_dataset;/g;
        $text =~ s/%let\s+fm_libpath=.*?;/%let fm_libpath=;/g;
        $text =~ s/%let\s+gtf_local_dsd=.*?;/%let gtf_local_dsd=$gtf_local_dataset;/g;
    }
    elsif ($mode eq 'local_gtf') {
        $text =~ s/%let\s+gtf_dsd=.*?;/%let gtf_dsd=$gtf_local_dataset;/g;
        $text =~ s/%let\s+fm_libpath=.*?;/%let fm_libpath=;/g;
        $text =~ s/%let\s+gtf_local_dsd=.*?;/%let gtf_local_dsd=$gtf_local_dataset;/g;
        $text =~ s/ods\s+html5\s+file="~\/[^"]+"\s*\n\s*options\(bitmap_mode='inline'\)\s*\n\s*style=HTMLBlue;/ods html5 path="&local_debug_output_dir" (url=none) file="&local_debug_output_name"\n  options(bitmap_mode='inline')\n  style=HTMLBlue;/g;
        $text =~ s/%include\s+"~\/get_top_signal_within_dist\.sas";/%include "&local_debug_deps\\get_top_signal_within_dist.sas";/g;
        $text =~ s/%include\s+"~\/[^"]*get_genecode_gtf_data[^"]*\.sas";/%include "&local_debug_gtf_macro";/g;
        $text =~ s/%include\s+"~\/adj_grpnum4close_gene_bed_regs\.sas";/%include "&local_debug_deps\\adj_grpnum4close_gene_bed_regs.sas";/g;
        $text =~ s/%include\s+"~\/Multgscatter_with_gene_exons\.sas";/%include "&local_debug_deps\\Multgscatter_with_gene_exons.sas";/g;
        $text =~ s/%include\s+"~\/map_grp_assoc2gene4covidsexgwas\.sas";/%include "&local_debug_deps\\map_grp_assoc2gene4covidsexgwas.sas";/g;
        $text =~ s/%include\s+"~\/SNP_Local_Manhattan_With_GTF\.sas";/%include "&local_debug_deps\\SNP_Local_Manhattan_With_GTF.sas";/g;
        $text =~ s/%include\s+"~\/Lattice_gscatter_over_bed_track\.sas";/%include "&local_debug_deps\\Lattice_gscatter_over_bed_track.sas";/g;
        $text =~ s/filename\s+gtfdata\s+zip\s+"~\/[^"]+"\s+gzip;/filename gtfdata zip "&local_debug_gtf_subset_gz" gzip;/g;
        $text =~ s/filename\s+mhdata\s+zip\s+"~\/[^"]+"\s+gzip;/filename mhdata zip "&local_debug_wide_data_gz" gzip;/g;
        $text =~ s/filename\s+_reqcsv\s+"~\/&lth_input_csv";/filename _reqcsv "&local_debug_root\\&lth_input_csv";/g;
        $text =~ s/proc\s+import\s+datafile="~\/&lth_input_csv"/proc import datafile="&local_debug_root\\&lth_input_csv"/g;
        $text =~ s/outfile="~\/&local_top_hits_csv_basename"/outfile="&local_debug_root\\&local_top_hits_csv_basename"/g;
        $text =~ s/file\s+"~\/[^"]+\.html"\s+lrecl=32767;/file "&local_debug_output_html" lrecl=32767;/g;
    }

    if (length $local_debug_top_hits_csv) {
        $text =~ s/%let\s+lth_input_csv=.*?;/%let lth_input_csv=@{[basename_any($local_debug_top_hits_csv)]};/g if $mode eq 'local_gtf';
    }

    return $text;
}

sub build_preamble {
    my (%args) = @_;
    my @lines = (
        '/* Local desktop-SAS debug paths emitted automatically by the pipeline. */',
        qq{%let local_debug_root=$args{local_debug_root};},
        qq{%let local_debug_deps=&local_debug_root\\DiffGWASDeps;},
        qq{%let local_debug_output_dir=&local_debug_root;},
        qq{%let local_debug_output_html=$args{local_debug_output_html};},
    );

    if ($args{mode} eq 'local_gtf') {
        push @lines,
            qq{%let local_debug_output_name=} . basename_any($args{local_debug_output_html}) . q{;},
            qq{%let local_debug_gtf_macro=$args{local_debug_gtf_macro};},
            qq{%let local_debug_gtf_subset_gz=$args{local_debug_gtf_subset_gz};},
            qq{%let local_debug_wide_data_gz=$args{local_debug_data_gz};};
    }
    elsif ($args{mode} eq 'local_manhattan') {
        push @lines,
            qq{%let local_debug_gtf_macro=$args{local_debug_gtf_macro};},
            qq{%let local_debug_manhattan_macro=$args{local_debug_manhattan_macro};},
            qq{%let local_debug_top_hits_csv=$args{local_debug_top_hits_csv};},
            qq{%let local_debug_wide_data_gz=$args{local_debug_data_gz};};
    }
    else {
        push @lines,
            qq{%let local_debug_manhattan_macro=$args{local_debug_manhattan_macro};},
            qq{%let local_debug_wide_data_gz=$args{local_debug_data_gz};};
    }

    return join("\n", @lines) . "\n\n";
}

sub slurp {
    my ($path) = @_;
    open my $fh, '<', $path or die "Cannot read $path: $!\n";
    local $/;
    my $text = <$fh>;
    close $fh;
    return $text;
}

sub write_text {
    my ($path, $text) = @_;
    my $dir = dirname($path);
    make_path($dir) if length($dir) && !-d $dir;
    open my $fh, '>', $path or die "Cannot write $path: $!\n";
    print {$fh} $text;
    close $fh or die "Cannot close $path: $!\n";
}

sub join_unix {
    my ($base, $child) = @_;
    return $child unless defined $base && length $base;
    $base =~ s{[\\/]+$}{};
    $child =~ s{^[\\/]+}{};
    return "$base/$child";
}

sub path_to_windows {
    my ($path) = @_;
    return '' unless defined $path;
    if ($path =~ m{^/mnt/([A-Za-z])/(.*)$}) {
        my ($drive, $rest) = (uc($1), $2);
        $rest =~ s{/}{\\}g;
        return "${drive}:\\${rest}";
    }
    $path =~ s{/}{\\}g;
    return $path;
}

sub basename_any {
    my ($path) = @_;
    return '' unless defined $path;
    $path =~ s{\\}{/}g;
    return basename($path);
}

sub usage {
    return <<"USAGE";
Usage:
  perl emit_local_sas_debug_script.pl --mode manhattan|local_manhattan|local_gtf \\
    --input rendered_oda.sas --output local_debug.sas --workdir /mnt/g/... --deps-dir /mnt/g/.../DiffGWASDeps \\
    --data-gz /mnt/e/... [--gtf-subset-gz ...] [--top-hits-csv ...] [--gtf-macro-upload ...]
    [--gtf-local-dataset gtf_hg38|gtf_hg19|gtf_t2t] [--manhattan-macro ...]
USAGE
}
