package SAS_ODA_MCP;
use strict;
use warnings;
use FindBin qw($Bin);
use File::Basename qw(basename dirname);
use File::Spec;
use File::Path qw(make_path remove_tree);
use File::Temp qw(tempdir);
use Cwd qw(abs_path getcwd);
use JSON::PP qw(encode_json decode_json);
use File::Copy qw(copy);
use POSIX qw(WNOHANG);
sub text_result {
    my ($text) = @_;
    return { content => [{ type => 'text', text => ($text // '') }] };
}

sub is_truthy {
    my ($value) = @_;
    return 0 unless defined $value;
    return $value =~ /^(?:1|true|yes|y|on)$/i ? 1 : 0;
}

sub is_falsey {
    my ($value) = @_;
    return 0 unless defined $value;
    return $value =~ /^(?:0|false|no|n|off)$/i ? 1 : 0;
}

sub as_list {
    my ($value) = @_;
    return () unless defined $value;
    return @$value if ref($value) eq 'ARRAY';
    return ($value);
}

sub runner_path {
    my $root = dirname(dirname(abs_path(__FILE__)));
    my $runner = File::Spec->catfile($root, 'run_sas_codes_or_script_in_ODA.pl');
    return $runner if -f $runner;
    die "Cannot find SAS ODA runner under $root\n";
}

sub pid_is_running {
    my ($pid) = @_;
    return 0 unless defined($pid) && $pid =~ /^\d+$/;
    return kill(0, $pid) ? 1 : 0;
}

sub read_text_file {
    my ($path) = @_;
    return '' unless defined($path) && length($path) && -f $path;
    open(my $fh, '<', $path) or return '';
    local $/;
    my $text = <$fh>;
    close $fh;
    return defined($text) ? $text : '';
}

sub write_text_file {
    my ($path, $text) = @_;
    open(my $fh, '>', $path) or die "Cannot write $path: $!\n";
    print {$fh} ($text // '');
    close $fh;
}

sub create_job_directories {
    my $root = getcwd();
    my $result_dir = tempdir('SAS_Output_ODA_XXXXXX', DIR => $root, CLEANUP => 0);
    my $scratch_dir = eval { tempdir('tmpXXXXXX', DIR => $root, CLEANUP => 0) };
    if (!$scratch_dir) {
        my $creation_error = $@ || 'unknown temporary-directory creation error';
        my $errors;
        remove_tree($result_dir, { error => \$errors });
        my $cleanup_detail = $errors && @{$errors} ? ' Could not remove the unused results directory.' : '';
        die "Cannot create worker scratch directory: $creation_error$cleanup_detail\n";
    }
    return ($result_dir, $scratch_dir);
}

sub cleanup_worker_scratch {
    my ($scratch_dir) = @_;
    return 'Worker scratch directory path is missing.'
      unless defined($scratch_dir) && length($scratch_dir);

    my $absolute_dir = File::Spec->rel2abs($scratch_dir);
    my $root = File::Spec->rel2abs(getcwd());
    return "Refusing to remove worker scratch directory outside $root: $absolute_dir"
      unless dirname($absolute_dir) eq $root && basename($absolute_dir) =~ /^tmp[A-Za-z0-9_-]+$/;

    my $errors;
    remove_tree($absolute_dir, { error => \$errors });
    my @messages = $errors && @{$errors}
      ? map { my ($path, $message) = %$_; "$path: $message" } @{$errors}
      : ();
    push @messages, "$absolute_dir still exists" if -e $absolute_dir;
    return @messages
      ? "Could not clean worker scratch directory: " . join('; ', @messages)
      : '';
}

sub find_status_for_pid {
    my ($pid) = @_;
    for my $pid_file (
        glob(File::Spec->catfile('SAS_Output_ODA_*', 'output.html.info.pid')),
        glob(File::Spec->catfile('tmp*', 'output.html.info.pid')),
    ) {
        next unless -f $pid_file;
        my $stored = read_text_file($pid_file);
        chomp $stored;
        next unless defined($stored) && $stored eq "$pid";

        my $out_file = $pid_file;
        $out_file =~ s/\.pid$/.txt/;
        return ($pid_file, $out_file);
    }
    return;
}

sub write_pid_file {
    my ($path, $pid) = @_;
    open(my $fh, '>', $path) or die "Cannot write $path: $!\n";
    print {$fh} $pid;
    close $fh;
}

sub reap_child_if_finished {
    my ($pid) = @_;
    return 0 unless defined($pid) && $pid =~ /^\d+$/;
    my $reaped = waitpid($pid, WNOHANG);
    return ($reaped && $reaped > 0) ? 1 : 0;
}

sub run_sas_oda_tool {
    my ($args) = @_;
    $args ||= {};

    if (defined $args->{pid}) {
        my ($pid_file, $out_file) = find_status_for_pid($args->{pid});
        return text_result("ERROR: PID $args->{pid} not found or already completed.")
          unless $pid_file;

        reap_child_if_finished($args->{pid});
        if (pid_is_running($args->{pid})) {
            return text_result(
                "STATUS: RUNNING (PID $args->{pid})\n"
              . "Output file: $out_file\n"
              . "Ask the AI agent to check status again in a moment."
            );
        }

        my $completion_file = File::Spec->catfile(dirname($pid_file), 'worker.result.json');
        my $completion = eval { decode_json(read_text_file($completion_file)) };
        $completion = {} unless ref($completion) eq 'HASH';
        my $log_path = $completion->{log_path} || $out_file;
        my $content = read_text_file($log_path);
        my $runner_stdout = File::Spec->catfile(dirname($pid_file), 'runner.stdout.txt');
        if ((!defined($content) || $content eq '') && -f $runner_stdout) {
            $log_path = $runner_stdout;
            $content = read_text_file($runner_stdout);
        }
        my $known = defined $completion->{exit_code};
        my $cleanup_error = $completion->{cleanup_error} // '';
        my $worker_error = $completion->{worker_error} // '';
        my $failed = !$known || $completion->{exit_code} != 0
          || length($cleanup_error) || length($worker_error);
        my $label = !$known ? 'UNKNOWN' : ($failed ? 'FAILED' : 'COMPLETE');
        my $exit_detail = $known ? "Exit code: $completion->{exit_code}\n"
            : "Worker stopped without a completion record; success is unconfirmed.\n";
        unlink $completion_file if -f $completion_file;
        unlink "$completion_file.tmp" if -f "$completion_file.tmp";
        unlink $pid_file if -f $pid_file;
        my $response = text_result(
            "STATUS: $label (PID $args->{pid})\n$exit_detail\n"
          . "SAS log for debugging saved to: $log_path\n\n"
          . ($content // '')
          . (length($cleanup_error) ? "\nWARNING: $cleanup_error\n" : '')
          . (length($worker_error) ? "\nERROR: $worker_error\n" : '')
        );
        $response->{isError} = $failed ? JSON::PP::true : JSON::PP::false;
        return $response;
    }

    my ($result_dir, $tmpdir) = create_job_directories();
    my $actual_log = File::Spec->catfile($result_dir, 'output.html.info.txt');
    my $out_file = defined($args->{output_file}) && length($args->{output_file})
        ? File::Spec->rel2abs($args->{output_file}) : $actual_log;
    make_path(dirname($out_file)) unless -d dirname($out_file);
    my $pid_file = File::Spec->catfile($result_dir, 'output.html.info.pid');
    # The CLI interprets --output-prefix as a directory, not a filename stem.
    my $output_prefix = $result_dir;

    my @runner_args = ('--output-prefix', $output_prefix);

    if (is_truthy($args->{kill_saspy_sessions})) {
        push @runner_args, '--kill-saspy-sessions';
    }

    my $persistent = exists($args->{persistent}) ? !is_falsey($args->{persistent}) : 1;
    my $session_id = $args->{session_id} // 'mysession';
    if ($persistent) {
        push @runner_args, '--persistent', '--session-id', $session_id;
    }

    if (defined($args->{run_timeout_seconds}) && $args->{run_timeout_seconds} =~ /^\d+$/) {
        push @runner_args, '--run-timeout-seconds', int($args->{run_timeout_seconds});
    }
    push @runner_args, '--no-run-timeout' if is_truthy($args->{no_run_timeout});

    my $sas_input = $args->{sas_codes_or_file};
    if (defined($sas_input) && length($sas_input)) {
        if (-f $sas_input) {
            push @runner_args, '--file', $sas_input;
        } else {
            my $code_file = File::Spec->catfile($tmpdir, 'input.sas');
            open(my $fh, '>', $code_file) or die "Cannot write $code_file: $!\n";
            print {$fh} $sas_input;
            close $fh;
            push @runner_args, '--file', $code_file;
        }
    }

    push @runner_args, map { ('--upload-file', $_) } as_list($args->{upload_file});
    push @runner_args, '--no-archive-transfers'
      if exists($args->{archive_transfers}) && is_falsey($args->{archive_transfers});
    push @runner_args, '--no-skip-upload-if-same'
      if exists($args->{skip_upload_if_same}) && is_falsey($args->{skip_upload_if_same});
    push @runner_args, map { ('--download-file', $_) } as_list($args->{download_file});
    push @runner_args, map { ('--download-local-path', $_) } as_list($args->{download_local_path});
    push @runner_args, map { ('--delete-file', $_) } as_list($args->{delete_file});
    push @runner_args, map { ('--delete-file-rgx', $_) } as_list($args->{delete_file_rgx});
    push @runner_args, ('--delete-dir', $args->{delete_dir}) if defined($args->{delete_dir}) && length($args->{delete_dir});
    push @runner_args, map { ('--file-info', $_) } as_list($args->{file_info});
    push @runner_args, ('--classify-sas-log', $args->{classify_sas_log})
      if defined($args->{classify_sas_log}) && length($args->{classify_sas_log});
    push @runner_args, '--no-cleanup-empty-output-dir'
      if exists($args->{cleanup_empty_output_dir}) && is_falsey($args->{cleanup_empty_output_dir});
    push @runner_args, ('--dir4listing', $args->{dir4listing}) if defined($args->{dir4listing}) && length($args->{dir4listing});
    push @runner_args, ('--saspy-cfgname', $args->{saspy_cfgname}) if defined($args->{saspy_cfgname}) && length($args->{saspy_cfgname});
    push @runner_args, ('--saspy-cfgfile', $args->{saspy_cfgfile}) if defined($args->{saspy_cfgfile}) && length($args->{saspy_cfgfile});
    push @runner_args, '--check-saspy-connection-only' if is_truthy($args->{check_saspy_connection_only});

    my $pid = fork();
    if (!defined $pid) {
        my $cleanup_error = cleanup_worker_scratch($tmpdir);
        my $errors;
        remove_tree($result_dir, { error => \$errors });
        my $message = "ERROR: Could not fork SAS ODA worker: $!";
        $message .= "\n$cleanup_error" if length $cleanup_error;
        $message .= "\nCould not remove results directory $result_dir." if $errors && @{$errors};
        return text_result($message);
    }

    if ($pid == 0) {
        $ENV{PIPELINE_SAS_ODA_ACCOUNT} = $args->{sas_oda_account}
          if defined($args->{sas_oda_account}) && length($args->{sas_oda_account});
        $ENV{PIPELINE_SAS_ODA_PASSWORD} = $args->{sas_oda_password}
          if defined($args->{sas_oda_password}) && length($args->{sas_oda_password});
        $ENV{PIPELINE_FORCE_SAS_ODA_AUTH_PROMPT} = 1
          if is_truthy($args->{prompt_sas_oda_auth});
        $ENV{SAS_ODA_WORKER_TEMP_DIR} = $tmpdir;
        my $python_temp_dir = $tmpdir;
        if ($^O eq 'cygwin' && $python_temp_dir =~ m{^([A-Za-z]):[\\/](.*)$}) {
            my ($drive, $rest) = (lc($1), $2);
            $rest =~ s{\\}{/}g;
            $python_temp_dir = "/mnt/$drive/$rest";
        }
        $ENV{TMPDIR} = $python_temp_dir;
        $ENV{OPEN_RESULT} //= 0;

        my $stdout_log = File::Spec->catfile($tmpdir, 'runner.stdout.txt');
        open(STDOUT, '>', $stdout_log) or die "Cannot redirect stdout to $stdout_log: $!\n";
        open(STDERR, '>&STDOUT') or die "Cannot redirect stderr: $!\n";

        my $runner = runner_path();
        my $wait_status = system { $^X } $^X, $runner, @runner_args;
        my $exit_code = $wait_status == -1 ? 127
            : ($wait_status & 127) ? 128 + ($wait_status & 127) : $wait_status >> 8;
        close STDOUT;
        close STDERR;
        my $log_path = $actual_log;
        my @worker_errors;
        if (!-f $actual_log && -f $stdout_log) {
            if (!copy($stdout_log, $actual_log)) {
                push @worker_errors, "Cannot preserve worker diagnostics in $actual_log: $!";
                $exit_code = 1 if $exit_code == 0;
            }
        }
        if ($out_file ne $actual_log) {
            if (-f $log_path && !copy($log_path, $out_file)) {
                push @worker_errors, "Cannot copy job log to $out_file: $!";
                $exit_code = 1 if $exit_code == 0;
            }
        }
        my $cleanup_error = cleanup_worker_scratch($tmpdir);
        my $completion_file = File::Spec->catfile($result_dir, 'worker.result.json');
        write_text_file("$completion_file.tmp", encode_json({
            exit_code => $exit_code, log_path => $log_path, cleanup_error => $cleanup_error,
            worker_error => join('; ', @worker_errors),
        }));
        rename "$completion_file.tmp", $completion_file or die "Cannot save completion: $!";
        exit $exit_code;
    }

    write_pid_file($pid_file, $pid);
    return text_result(
        "QUERYING: SAS ODA worker started\n"
      . "PID: $pid\n"
      . "Output file: $out_file\n"
      . "Results directory: $result_dir\n"
      . "Worker scratch is removed automatically when the job finishes.\n"
      . "Ask the AI agent to check status with: {\"pid\": $pid}\n"
      . "For long SAS jobs, check no more than every 30 seconds."
    );
}

1;
