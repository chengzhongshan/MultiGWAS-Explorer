#!/usr/bin/env perl
use strict;
use warnings;
use FindBin qw($Bin);
use File::Spec;
use File::Temp qw(tempdir);
use Test::More;

my @wrappers = qw(
  run_sas_oda_local_top_hits_manhattan_download_png.sh
  run_sas_oda_local_top_hits_with_gtf_download_html.sh
);
my $tmp = tempdir('sas_oda_upload_timeout_XXXXXX', TMPDIR => 1, CLEANUP => 1);
my $data = File::Spec->catfile($tmp, 'large subset.tsv.gz');
open my $data_fh, '>', $data or die "Cannot create $data: $!\n";
truncate($data_fh, 118_496_045) or die "Cannot size $data: $!\n";
close $data_fh or die "Cannot close $data: $!\n";

for my $wrapper (@wrappers) {
    my $path = File::Spec->catfile($Bin, $wrapper);
    open my $fh, '<', $path or die "Cannot read $path: $!\n";
    local $/;
    my $source = <$fh>;
    close $fh or die "Cannot close $path: $!\n";

    my ($function) = $source =~ /(^oda_upload_timeout_for_file\(\) \{.*?^\})/ms;
    ok(defined($function), "$wrapper defines the size-aware timeout calculator");
    like(
        $source,
        qr/run_oda_upload_helper\s+"\$\{DATA_GZ\}"\s+\\\s*\n\s*--upload-file\s+"\$\{DATA_GZ\}"/,
        "$wrapper uses the size-aware timeout for its large data upload",
    );
    next unless defined $function;

    my $test_script = File::Spec->catfile($tmp, "$wrapper.test.sh");
    open my $out, '>', $test_script or die "Cannot write $test_script: $!\n";
    print {$out} <<"SHELL";
#!/usr/bin/env bash
set -euo pipefail
$function
ODA_HELPER_TIMEOUT_SECONDS=300
unset ODA_DATA_UPLOAD_TIMEOUT_SECONDS ODA_UPLOAD_TIMEOUT_BASE_SECONDS ODA_UPLOAD_TIMEOUT_SECONDS_PER_MB || true
test "\$(oda_upload_timeout_for_file "\$1")" = 1440
ODA_DATA_UPLOAD_TIMEOUT_SECONDS=900
test "\$(oda_upload_timeout_for_file "\$1")" = 900
SHELL
    close $out or die "Cannot close $test_script: $!\n";
    is(system('/usr/bin/bash', $test_script, $data), 0,
        "$wrapper gives the 118 MB upload enough time and honors an explicit override");
}

done_testing();
