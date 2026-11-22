#!/usr/bin/env perl
use strict;
use warnings;
use FindBin qw($Bin);
use JSON::PP qw(decode_json);
use Digest::SHA qw(sha256_hex);
use Test::More;

sub text {open my $fh,'<:raw',$_[0] or die $!;local$/;my $data=<$fh>;$data=~s/\r\n/\n/g;return $data;}
my $root="$Bin/..";
my $manifest=decode_json(text("$root/docs/SAS_RUNTIME_PORT_MANIFEST.json"));
is($manifest->{upstream_commit},'cd0caa928bd37da3a7720a136d6949c261cbb545','exact reviewed MCP4SAS source revision');
is($manifest->{hash_policy},'SHA-256 after CRLF-to-LF normalization','portable text-hash convention');
for my $file(@{$manifest->{files}}) {
    like($file->{path},qr{^(?:DiffGWASDeps/|tests/|tools/|run_sas_codes_or_script_in_ODA\.pl$|server\.pl$)},'source map stays inside runtime/test paths');
    unlike($file->{path},qr{(?:^|/)\.\.(?:/|$)},'no source-map traversal');
    is(sha256_hex(text("$root/$file->{path}")),$file->{sha256},"recorded port hash $file->{path}");
}
for my $alias(@{$manifest->{compatibility_entrypoints}}) {
    like(text("$root/$alias"),qr/run_sas_codes_or_script_in_ODA\.pl/,'compatibility alias retains canonical dispatch');
}
my $validation=decode_json(text("$root/docs/SAS_RUNTIME_VALIDATION_20261121.json"));
ok($validation->{synthetic_live}{passed},'synthetic SAS integration record is successful');
is($validation->{multi_gwas}{checks},11295,'separate numerical validation count retained');
is(scalar(grep {exists $_->{numerical_equal}} @{$validation->{multi_gwas}{comparison}}),4,'four numerical tables compared');
is(scalar(grep {$_->{pixel_identical}} @{$validation->{multi_gwas}{comparison}}),1,'PNG appearance comparison is separate');
done_testing;
