#!/usr/bin/env perl
use strict;
use warnings;
use FindBin qw($Bin);
use lib "$Bin/../DiffGWASDeps";
use SAS_ODA_Runner;
use Test::More;
my $runner = bless {persistent => 1, session_id => 'offline-test'}, 'SAS_ODA_Runner';
my ($request, $timeout, $calls);
{
    no warnings 'redefine';
    local *SAS_ODA_Runner::_call_persistent_session_server = sub {
        my ($self, $payload, $seconds) = @_;
        ($request, $timeout) = ($payload, $seconds);
        ++$calls;
        return {status => 'ok', value => {uploads => [], downloads => []}};
    };
    for my $archive (0, 1) {
        my $payload = {archive_transfers => $archive, downloads => [
            {remote_path => '~/a', local_path => 'a'}, {remote_path => '~/b', local_path => 'b'}]};
        my $result = $runner->transfer_many($payload);
        is($request->{cmd}, 'bulk_transfer', 'persistent request uses real bulk dispatch');
        is($request->{session_id}, 'offline-test', 'existing session identity is preserved');
        is($request->{payload}, $payload, 'complete payload and archive flag reach the session server');
        is_deeply($result, {uploads => [], downloads => []}, 'successful result is propagated');
        cmp_ok($timeout, '>', 0, 'file-transfer timeout is supplied');
    }
    is($calls, 2, 'one RPC per transfer batch, not one request per file');
}
{
    no warnings 'redefine';
    local *SAS_ODA_Runner::_call_persistent_session_server = sub {
        return {status => 'error', error => 'integrity check failed'};
    };
    like($runner->transfer_many({}), qr/^PYTHON ERROR: integrity check failed$/,
        'bulk failure remains a reported failure');
}
open my $fh, '<:raw', "$Bin/../DiffGWASDeps/sas_oda_session_server.py" or die $!;
my $standalone = do {local $/; <$fh>};
is($standalone, SAS_ODA_Runner::_session_server_source(),
    'checked-in Python helper exactly matches the generated session implementation');
like($standalone, qr/SERVER_API_VERSION = '2026-11-12-persistent-archive'/,
    'new protocol rejects stale helpers that lack bulk dispatch');
like($standalone, qr/elif cmd == 'bulk_transfer':.*?transfer_many\(sess,/s,
    'Python bulk handler operates on the existing SASPy session');
done_testing();
