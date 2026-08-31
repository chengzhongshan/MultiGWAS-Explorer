#!/usr/bin/env perl
use strict;
use warnings;
use utf8;
use feature qw(state);
use Getopt::Long qw(GetOptions);
use JSON::PP ();
use Digest::SHA qw(sha256_hex);
use Time::HiRes qw(time sleep);
use POSIX qw(strftime WNOHANG);
use File::Spec;
use File::Path qw(make_path);
use File::Basename qw(dirname);
use File::Find qw(find);
use File::Glob qw(bsd_glob GLOB_NOSORT);
use Cwd qw(abs_path getcwd);
use Text::ParseWords qw(shellwords);
use Text::CSV;
use IO::Socket::INET;
use Encode qw(encode);

# Automated Codex benchmark runner for MultiGWAS-Explorer.
#
# This is the Perl counterpart of the earlier Python runner.  It is designed
# for Cygwin/Linux/macOS Perl and should also work under modern Strawberry Perl
# where fork() emulation is available.
#
# Evaluator-only gold actions are never sent to the AI clients.

my %VOLATILE_JSON_KEYS = map { $_ => 1 } qw(
    timestamp created_at updated_at start_time end_time elapsed elapsed_seconds
    runtime runtime_seconds duration duration_ms pid output_file output_path
    tmpdir temp_dir session_id thread_id
);
my $FAILURE_RE = qr/\b(?:error|failed|failure|not found|missing|invalid|unknown|does not exist|cannot find|unresolved|rejected)\b/i;
my $CLARIFY_RE = qr/\b(?:clarif|contradict|conflict|choose|which backend|cannot both|incompatible)\w*/i;
my $CODEX_UNAVAILABLE_RE = qr/
    (?:you(?:'|\x{2019})ve\s+hit\s+your\s+usage\s+limit)
  | (?:usage\s+limit(?:\s+(?:reached|exceeded))?)
  | (?:purchase\s+more\s+credits)
  | (?:insufficient\s+(?:credits|quota))
  | (?:quota\s+(?:reached|exceeded))
  | (?:rate\s+limit\s+(?:reached|exceeded))
  | (?:too\s+many\s+requests)
  | (?:authentication\s+(?:required|failed))
  | (?:not\s+logged\s+in)
  | (?:unauthorized|forbidden)
/ix;

my $JSON = JSON::PP->new->utf8(0)->allow_nonref(1);
my $JSON_CANON = JSON::PP->new->utf8(0)->canonical(1)->allow_nonref(1);
my $JSON_PRETTY = JSON::PP->new->utf8(0)->canonical(1)->pretty(1)->allow_nonref(1);

sub utcnow {
    return strftime('%Y-%m-%dT%H:%M:%SZ', gmtime(time()));
}

sub slurp {
    my ($path) = @_;
    open my $fh, '<:raw', $path or die "Cannot read $path: $!\n";
    local $/;
    my $x = <$fh>;
    close $fh;
    return $x;
}

sub spew {
    my ($path, $data) = @_;
    my $dir = dirname($path);
    make_path($dir) if $dir && !-d $dir;
    open my $fh, '>:raw', $path or die "Cannot write $path: $!\n";
    $data = encode('UTF-8', $data) if utf8::is_utf8($data);
    print {$fh} $data;
    close $fh;
}

sub load_json {
    my ($path) = @_;
    return $JSON->decode(slurp($path));
}

sub write_json {
    my ($path, $obj) = @_;
    spew($path, $JSON_PRETTY->encode($obj));
}

sub load_jsonl {
    my ($path) = @_;
    open my $fh, '<:encoding(UTF-8)', $path or die "Cannot read $path: $!\n";
    my @rows;
    while (my $line = <$fh>) {
        $line =~ s/^\x{FEFF}//;
        $line =~ s/\s+$//;
        next unless length $line;
        push @rows, $JSON->decode($line);
    }
    close $fh;
    return \@rows;
}

sub sha256_file {
    my ($path) = @_;
    open my $fh, '<:raw', $path or die "Cannot hash $path: $!\n";
    my $sha = Digest::SHA->new(256);
    $sha->addfile($fh);
    close $fh;
    return $sha->hexdigest;
}

sub percentile {
    my ($vals, $p) = @_;
    my @v = sort { $a <=> $b } map { 0 + $_ } @$vals;
    return undef unless @v;
    return $v[0] if @v == 1;
    my $k = ($#v) * $p;
    my $lo = int($k);
    my $hi = ($k == int($k)) ? int($k) : int($k) + 1;
    return $v[$lo] if $lo == $hi;
    return $v[$lo] * ($hi - $k) + $v[$hi] * ($k - $lo);
}

sub median {
    my ($vals) = @_;
    return percentile($vals, 0.5);
}

sub iqr {
    my ($vals) = @_;
    return (percentile($vals, 0.25), percentile($vals, 0.75));
}

sub wilson {
    my ($success, $n, $z) = @_;
    $z //= 1.959963984540054;
    return (undef, undef) unless $n && $n > 0;
    my $phat = $success / $n;
    my $den = 1 + $z*$z/$n;
    my $center = ($phat + $z*$z/(2*$n))/$den;
    my $half = $z * sqrt(($phat*(1-$phat) + $z*$z/(4*$n))/$n)/$den;
    my $lo = $center - $half; $lo = 0 if $lo < 0;
    my $hi = $center + $half; $hi = 1 if $hi > 1;
    return ($lo, $hi);
}

sub fmt_pct {
    my ($success, $n) = @_;
    return 'NA' unless $n;
    my ($lo, $hi) = wilson($success, $n);
    return sprintf('%.1f%% (%.1f-%.1f%%)', 100*$success/$n, 100*$lo, 100*$hi);
}

sub fmt_num {
    my ($x, $digits) = @_;
    $digits //= 3;
    return 'NA' unless defined $x;
    return sprintf("%.${digits}f", $x);
}

sub resolve_placeholders {
    my ($text, $ph) = @_;
    return $text unless defined $text;
    for my $k (keys %$ph) {
        my $needle = '{' . $k . '}';
        my $v = defined $ph->{$k} ? "$ph->{$k}" : '';
        $text =~ s/\Q$needle\E/$v/g;
    }
    return $text;
}

sub resolve_obj {
    my ($obj, $ph) = @_;
    if (!ref $obj) {
        return defined($obj) ? resolve_placeholders("$obj", $ph) : undef;
    }
    if (ref($obj) eq 'ARRAY') {
        return [ map { resolve_obj($_, $ph) } @$obj ];
    }
    if (ref($obj) eq 'HASH') {
        return { map { $_ => resolve_obj($obj->{$_}, $ph) } keys %$obj };
    }
    return $obj;
}

sub unresolved_placeholders {
    my ($obj) = @_;
    my %found;
    my $visit;
    $visit = sub {
        my ($x) = @_;
        if (!ref $x) {
            return unless defined $x;
            while ("$x" =~ /\{([A-Z][A-Z0-9_]*)\}/g) { $found{$1} = 1 }
        } elsif (ref($x) eq 'ARRAY') {
            $visit->($_) for @$x;
        } elsif (ref($x) eq 'HASH') {
            $visit->($_) for values %$x;
        }
    };
    $visit->($obj);
    return [ sort keys %found ];
}

sub shell_quote {
    my ($s) = @_;
    return "''" if !defined($s) || $s eq '';
    return $s if $s =~ /^[A-Za-z0-9_\.\-\/:=,@+]+$/;
    $s =~ s/'/'"'"'/g;
    return "'$s'";
}

sub command_string {
    my ($cmd) = @_;
    return join(' ', map { shell_quote($_) } @$cmd);
}

sub which_cmd {
    my ($cmd) = @_;
    return $cmd if File::Spec->file_name_is_absolute($cmd) && -x $cmd;
    for my $dir (File::Spec->path()) {
        for my $name ($cmd, "$cmd.exe", "$cmd.cmd", "$cmd.bat") {
            my $p = File::Spec->catfile($dir, $name);
            return $p if -f $p && -x _;
        }
    }
    return undef;
}

sub run_capture {
    my (%a) = @_;
    my $cmd = $a{cmd};
    my $cwd = $a{cwd};
    my $timeout = $a{timeout} // 7200;
    my $stdout_path = $a{stdout};
    my $stderr_path = $a{stderr};
    make_path(dirname($stdout_path));
    make_path(dirname($stderr_path));

    my $t0 = time();
    my $pid = fork();
    die "fork failed: $!\n" unless defined $pid;
    if ($pid == 0) {
        chdir $cwd or do { print STDERR "Cannot chdir to $cwd: $!\n"; exit 127 };
        open STDOUT, '>:raw', $stdout_path or exit 127;
        open STDERR, '>:raw', $stderr_path or exit 127;
        { no warnings 'exec'; exec { $cmd->[0] } @$cmd; }
        exit 127;
    }

    my $timed_out = 0;
    my $status;
    while (1) {
        my $w = waitpid($pid, WNOHANG);
        if ($w == $pid) { $status = $?; last }
        if (time() - $t0 > $timeout) {
            $timed_out = 1;
            kill 'TERM', $pid;
            sleep 2;
            my $w2 = waitpid($pid, WNOHANG);
            if ($w2 == 0) { kill 'KILL', $pid; waitpid($pid, 0) }
            $status = $?;
            open my $efh, '>>:encoding(UTF-8)', $stderr_path;
            print {$efh} "\nBENCHMARK_TIMEOUT after $timeout seconds\n";
            close $efh;
            last;
        }
        sleep 0.2;
    }
    my $elapsed = time() - $t0;
    my $exit_code = $timed_out ? 124 : (($status // 0) >> 8);
    return ($exit_code, $elapsed, $timed_out ? JSON::PP::true : JSON::PP::false);
}

sub capture_text {
    my (%a) = @_;
    my $tmp = File::Spec->catfile(File::Spec->tmpdir(), "multigwas_bench_$$." . int(rand(1e9)));
    my $err = "$tmp.err";
    my ($code) = run_capture(cmd=>$a{cmd}, cwd=>$a{cwd}, timeout=>$a{timeout}//60, stdout=>$tmp, stderr=>$err);
    my $txt = (-e $tmp ? slurp($tmp) : '') . (-e $err ? slurp($err) : '');
    unlink $tmp, $err;
    return ($code, $txt);
}

sub git_output {
    my ($repo, @args) = @_;
    my ($code, $txt) = capture_text(cmd=>['git', @args], cwd=>$repo, timeout=>30);
    $txt =~ s/\s+$//;
    return $code == 0 ? $txt : "ERROR: $txt";
}

sub git_diff_fingerprint {
    my ($repo) = @_;
    my ($code1, $diff) = capture_text(cmd=>['git','diff','--binary','HEAD','--'], cwd=>$repo, timeout=>30);
    my ($code2, $names) = capture_text(cmd=>['git','diff','--name-only','HEAD','--'], cwd=>$repo, timeout=>30);
    if ($code1 != 0 || $code2 != 0) {
        return { sha256 => 'ERROR', tracked_changed_files => [] };
    }
    my @names = sort grep { length } map { s!\\!/!g; s/\s+$//r } split /\r?\n/, $names;
    return { sha256 => sha256_hex($diff), tracked_changed_files => \@names };
}

sub command_version {
    my ($cmd, $cwd) = @_;
    return 'NOT_FOUND' unless which_cmd($cmd);
    for my $args (['--version'], ['version']) {
        my ($code, $txt) = capture_text(cmd=>[$cmd, @$args], cwd=>$cwd, timeout=>30);
        if (length $txt) {
            my ($line) = split /\r?\n/, $txt;
            return $line if defined $line && length $line;
        }
    }
    return which_cmd($cmd) // 'NOT_FOUND';
}

sub parse_mcp_url {
    my ($url) = @_;
    if ($url =~ m!^https?://([^/:]+)(?::(\d+))?!i) {
        my $host = $1;
        my $port = defined($2) ? 0+$2 : ($url =~ /^https:/i ? 443 : 80);
        return ($host, $port);
    }
    return ('127.0.0.1', 8080);
}

sub wait_for_port {
    my ($host, $port, $timeout) = @_;
    my $deadline = time() + $timeout;
    while (time() < $deadline) {
        my $sock = IO::Socket::INET->new(PeerAddr=>$host, PeerPort=>$port, Proto=>'tcp', Timeout=>1);
        if ($sock) { close $sock; return 1 }
        sleep 0.5;
    }
    return 0;
}

sub parse_jsonl_file {
    my ($path) = @_;
    return [] unless -e $path;
    open my $fh, '<:encoding(UTF-8)', $path or return [];
    my @ev;
    while (my $line = <$fh>) {
        $line =~ s/\s+$//;
        next unless length $line;
        my $obj = eval { $JSON->decode($line) };
        push @ev, $obj if ref($obj) eq 'HASH';
    }
    close $fh;
    return \@ev;
}

sub recursive_text {
    my ($x) = @_;
    my @parts;
    my $walk;
    $walk = sub {
        my ($v) = @_;
        if (!ref $v) { push @parts, "$v" if defined $v }
        elsif (ref($v) eq 'ARRAY') { $walk->($_) for @$v }
        elsif (ref($v) eq 'HASH') { $walk->($_) for values %$v }
    };
    $walk->($x);
    return join("\n", @parts);
}

sub normalize_bool {
    my ($v) = @_;
    return $v if ref($v) eq 'JSON::PP::Boolean';
    if (!ref $v && defined $v) {
        return JSON::PP::true if "$v" =~ /^(?:1|true|yes|y|on)$/i;
        return JSON::PP::false if "$v" =~ /^(?:0|false|no|n|off)$/i;
    }
    return $v;
}

sub is_bool {
    my ($v) = @_;
    return ref($v) eq 'JSON::PP::Boolean';
}

sub normalize_scalar {
    my ($v) = @_;
    $v = normalize_bool($v);
    if (!ref $v && defined $v) {
        my $s = "$v";
        $s =~ s!\\!/!g;
        $s =~ s/^\s+|\s+$//g;
        return $s;
    }
    return $v;
}

sub split_csvish {
    my ($v) = @_;
    return [] unless defined $v;
    my @xs = ref($v) eq 'ARRAY' ? map { "$_" } @$v : split /,/, "$v";
    @xs = sort grep { length } map { s/^\s+|\s+$//gr } @xs;
    return \@xs;
}

sub numeric_equal {
    my ($a, $b) = @_;
    return 0 unless defined($a) && defined($b);
    return 0 unless "$a" =~ /^[-+]?\d*\.?\d+(?:[eE][-+]?\d+)?$/ && "$b" =~ /^[-+]?\d*\.?\d+(?:[eE][-+]?\d+)?$/;
    my $fa = 0 + $a; my $fb = 0 + $b;
    my $tol = 1e-12 + 1e-12 * (abs($fa) > abs($fb) ? abs($fa) : abs($fb));
    return abs($fa-$fb) <= $tol;
}

sub scalar_equal {
    my ($a, $b) = @_;
    $a = normalize_scalar($a); $b = normalize_scalar($b);
    if (is_bool($a) || is_bool($b)) {
        my $aa = normalize_bool($a); my $bb = normalize_bool($b);
        return (is_bool($aa) && is_bool($bb) && ($aa ? 1:0) == ($bb ? 1:0));
    }
    return 1 if defined($a) && defined($b) && "$a" eq "$b";
    return 1 if numeric_equal($a,$b);
    return lc(defined($a)?"$a":'') eq lc(defined($b)?"$b":'');
}

sub canonical_action {
    my ($entry, $args, $source) = @_;
    $source //= '';
    $args ||= {};
    my %aliases = (
        spec_file=>'spec', get_common_associations=>'common', raw_column_alias_config=>'raw_column_alias_config',
        local_gtf_window_bp=>'local_gtf_window_bp', target_snp_genes=>'target_snp_genes', target_snps=>'target_snps',
        display_gwas=>'display_gwas', spec_out=>'spec_out', gwas_dir=>'gwas_dir', generate_spec_only=>'generate_spec_only',
        preview_spec=>'preview_spec', force=>'force', mode=>'mode'
    );
    my %out;
    for my $k (keys %$args) {
        my $kk = $k; $kk =~ s/-/_/g;
        my $nk = $aliases{$kk} // $kk;
        next if $nk eq 'pid' || $nk eq 'output_file';
        $out{$nk} = normalize_scalar($args->{$k});
    }
    my (@plots,@steps);
    if (exists $out{plots}) { push @plots, @{split_csvish(delete $out{plots})} }
    if (exists $out{step}) {
        for my $s (@{split_csvish(delete $out{step})}) {
            if ($s =~ /^plot_(.+)$/) { push @plots, $1 } else { push @steps, $s }
        }
    }
    for my $key (keys %out) {
        if ($key =~ /^plot_(.+)$/ && is_bool(normalize_bool($out{$key})) && normalize_bool($out{$key})) {
            push @plots, $1; delete $out{$key};
        }
    }
    for my $key (qw(merge_raw sort_long diff_pairs standardize_diff extract_wide_subset)) {
        if (exists $out{$key} && is_bool(normalize_bool($out{$key})) && normalize_bool($out{$key})) {
            push @steps, $key; delete $out{$key};
        }
    }
    if (@plots) { my %u; $out{plots} = [ sort grep { !$u{$_}++ } @plots ] }
    if (@steps) { my %u; $out{steps} = [ sort grep { !$u{$_}++ } @steps ] }
    $out{common} = normalize_bool($out{common}) if exists $out{common};
    for my $key (qw(force generate_spec_only preview_spec no_real keep_workdir)) {
        $out{$key} = normalize_bool($out{$key}) if exists $out{$key};
    }
    return { entry=>$entry, args=>\%out, source=>$source };
}

sub parse_flags {
    my ($tokens) = @_;
    my %args;
    for (my $i=0; $i<@$tokens; $i++) {
        my $tok = $tokens->[$i];
        next unless $tok =~ /^--(.+)/;
        my $kv = $1;
        if ($kv =~ /^([^=]+)=(.*)$/s) {
            (my $k=$1)=~s/-/_/g; my $v=$2; if (exists $args{$k}) { $args{$k}=ref($args{$k}) eq 'ARRAY' ? [@{$args{$k}},$v] : [$args{$k},$v] } else { $args{$k}=$v }
        } else {
            (my $k=$kv)=~s/-/_/g;
            if ($i+1 < @$tokens && $tokens->[$i+1] !~ /^-/) { $args{$k}=$tokens->[++$i] }
            else { $args{$k}=JSON::PP::true }
        }
    }
    return \%args;
}

sub parse_flags_from_command {
    my ($command) = @_;
    my %args;
    while ($command =~ /--([A-Za-z0-9][A-Za-z0-9-]*)(?:=([^\s'";]+)|\s+(?!-)([^\s'";]+))?/g) {
        (my $key=$1)=~s/-/_/g;
        my $value=defined($2)?$2:defined($3)?$3:JSON::PP::true;
        if (exists $args{$key}) {
            $args{$key}=ref($args{$key}) eq 'ARRAY' ? [@{$args{$key}},$value] : [$args{$key},$value];
        } else {
            $args{$key}=$value;
        }
    }
    return \%args;
}

sub extract_actions_from_command {
    my ($command) = @_;
    my @actions;
    for my $seg (split /\s*(?:&&|;)\s*/, $command) {
        next unless $seg =~ /\S/;
        my @tokens = eval { shellwords($seg) };
        @tokens = split /\s+/, $seg if $@ || !@tokens;
        my $joined = join(' ', @tokens);
        my $flags = parse_flags_from_command($seg);
        next if exists $flags->{help} || exists $flags->{version};
        if ($joined =~ /\bperl(?:\.exe)?\s+\S*auto_prepare_and_run_diff_gwas_with_gunplot\.pl\b/) {
            push @actions, canonical_action('gnuplot_pipeline',$flags,'shell');
        } elsif ($joined =~ /\bperl(?:\.exe)?\s+\S*auto_prepare_and_run_diff_gwas\.pl\b/) {
            push @actions, canonical_action('sas_pipeline',$flags,'shell');
        } elsif ($joined =~ /\bbash(?:\.exe)?(?:['"])?\s+\S*check_pipeline_install\.sh\b/) {
            push @actions, canonical_action('install_check',{},'shell');
        } elsif ($joined =~ /\bperl(?:\.exe)?\s+\S*test_top_hit_maf_filter\.pl\b/) {
            push @actions, canonical_action('maf_test',{ no_real=>exists($flags->{no_real})?JSON::PP::true:JSON::PP::false, keep_workdir=>exists($flags->{keep_workdir})?JSON::PP::true:JSON::PP::false },'shell');
        } elsif ($joined =~ /\bperl(?:\.exe)?\s+\S*run_sas_codes_or_script_in_ODA\.pl\b/ && exists $flags->{check_sas_oda_login_only}) {
            push @actions, canonical_action('sas_login_check',{},'shell');
        }
    }
    return \@actions;
}

sub codex_parse {
    my ($events) = @_;
    my (@actions,@tools,@commands,@file_changes,@final_msgs,@tool_texts);
    my $terminal_ok = 0; my $usage = {};
    for my $ev (@$events) {
        $usage = $ev->{usage} || $usage if ($ev->{type}//'') eq 'turn.completed';
        my $item = ref($ev->{item}) eq 'HASH' ? $ev->{item} : next;
        my $it = $item->{type} // '';
        if ($it eq 'agent_message') {
            push @final_msgs, "$item->{text}" if defined $item->{text};
        } elsif ($it eq 'command_execution') {
            my $cmd = "$item->{command}";
            push @commands, {command=>$cmd,exit_code=>$item->{exit_code},status=>$item->{status}};
            if (($item->{status}//'') =~ /^(?:completed|failed)$/) {
                my $recognized=extract_actions_from_command($cmd);
                push @actions, @$recognized;
                $terminal_ok = 1 if ($item->{status}//'') eq 'completed' && @$recognized
                    && defined($item->{exit_code}) && $item->{exit_code}==0;
            }
        } elsif ($it eq 'file_change') {
            push @file_changes, $item;
        } elsif ($it eq 'mcp_tool_call') {
            my $tool = { server=>$item->{server},tool=>$item->{tool},arguments=>$item->{arguments}||{},status=>$item->{status},error=>$item->{error},result=>$item->{result} };
            push @tools, $tool;
            my $txt = recursive_text($item->{result})."\n".recursive_text($item->{error});
            push @tool_texts, $txt;
            if (($item->{tool}//'') eq 'auto_prepare_and_run_diff_gwas') {
                my $a = $item->{arguments} || {};
                push @actions, canonical_action('sas_pipeline',$a,'mcp') unless $a->{pid};
            }
            if (($item->{tool}//'') eq 'run_gunplot_wrapper') {
                my $a = $item->{arguments} || {};
                push @actions, canonical_action('gnuplot_pipeline',$a,'mcp') unless $a->{pid};
            }
            $terminal_ok = 1 if $txt =~ /STATUS:\s*COMPLETE/ && $txt !~ $FAILURE_RE;
        }
    }
    my %seen_action;
    @actions = grep {
        my $key=$JSON_CANON->encode({entry=>$_->{entry},args=>$_->{args}});
        !$seen_action{$key}++;
    } @actions;
    return { actions=>\@actions,mcp_calls=>\@tools,commands=>\@commands,file_changes=>\@file_changes,
             final_text=>join("\n",@final_msgs),all_text=>join("\n",@final_msgs,@tool_texts),
             terminal_ok=>$terminal_ok?JSON::PP::true:JSON::PP::false,usage=>$usage };
}

sub action_matches {
    my ($exp,$act) = @_;
    return 0 unless ($exp->{entry}//'') eq ($act->{entry}//'');
    my $ea=$exp->{args}||{}; my $aa=$act->{args}||{};
    for my $k (keys %$ea) {
        # The MCP schema exposes config-only execution as generate_spec_only,
        # while the CLI/gold schema names the same operation mode=configs.
        if ($k eq 'mode' && scalar_equal($ea->{$k}, 'configs')
                && !exists($aa->{$k}) && exists($aa->{generate_spec_only})
                && normalize_bool($aa->{generate_spec_only})) {
            next;
        }
        if ($k eq 'generate_spec_only' && normalize_bool($ea->{$k})
                && !exists($aa->{$k}) && exists($aa->{mode})
                && scalar_equal($aa->{mode}, 'configs')) {
            next;
        }
        return 0 unless exists $aa->{$k};
        my $ev=$ea->{$k}; my $av=$aa->{$k};
        if ($k =~ /^(?:spec|gwas_dir|spec_out|raw_column_alias_config)$/
                && !ref($ev) && !ref($av)) {
            my $ep = "$ev"; my $ap = "$av";
            $ep =~ s!\\!/!g; $ap =~ s!\\!/!g;
            $ep =~ s!/+!/!g; $ap =~ s!/+!/!g;
            $ep =~ s!^\./!!; $ap =~ s!^\./!!;
            next if lc($ep) eq lc($ap)
                 || lc($ap) =~ m{(?:^|/)\Q@{[lc($ep)]}\E$}
                 || lc($ep) =~ m{(?:^|/)\Q@{[lc($ap)]}\E$};
        }
        if (ref($ev) eq 'ARRAY') {
            my @e=sort map { "$_" } @$ev;
            my @a=sort map { "$_" } @{ ref($av) eq 'ARRAY' ? $av : split_csvish($av) };
            return 0 unless $JSON_CANON->encode(\@e) eq $JSON_CANON->encode(\@a);
        } else {
            return 0 unless scalar_equal($ev,$av);
        }
    }
    return 1;
}

sub match_action_set {
    my ($expected,$actual) = @_;
    return (!@$actual, []) unless @$expected;
    my (%used,@matched);
    for my $exp (@$expected) {
        my $hit;
        for my $i (0..$#$actual) {
            next if $used{$i};
            if (action_matches($exp,$actual->[$i])) { $hit=$i; last }
        }
        return (0,\@matched) unless defined $hit;
        $used{$hit}=1; push @matched,$hit;
    }
    return (1,\@matched);
}

sub action_signature {
    my ($actions) = @_;
    my @clean = map { {entry=>$_->{entry},args=>$_->{args}||{}} } @$actions;
    return sha256_hex($JSON_CANON->encode(\@clean));
}

sub canonicalize_json {
    my ($obj,$repo) = @_;
    if (ref($obj) eq 'HASH') {
        my %h;
        for my $k (sort keys %$obj) {
            next if $VOLATILE_JSON_KEYS{lc($k)};
            $h{$k}=canonicalize_json($obj->{$k},$repo);
        }
        return \%h;
    }
    if (ref($obj) eq 'ARRAY') { return [map { canonicalize_json($_,$repo) } @$obj] }
    if (!ref($obj) && defined $obj) {
        my $s="$obj"; my $r=$repo; my $rf=$repo; $rf =~ s!\\!/!g;
        $s =~ s/\Q$r\E/<REPO>/g; $s =~ s/\Q$rf\E/<REPO>/g;
        return $s;
    }
    return $obj;
}

sub canonical_hash {
    my ($path,$repo) = @_;
    my ($ext) = $path =~ /(\.[^.]+)$/; $ext=lc($ext//'');
    my $ok = eval {
        if ($ext eq '.json') {
            my $obj=load_json($path);
            return sha256_hex($JSON_CANON->encode(canonicalize_json($obj,$repo)));
        }
        if ($ext eq '.csv' || $ext eq '.tsv') {
            my $sep=$ext eq '.csv' ? ',' : "\t";
            my $csv=Text::CSV->new({binary=>1,sep_char=>$sep,auto_diag=>1});
            open my $fh,'<:encoding(UTF-8)',$path or die $!;
            my $hdr=$csv->getline($fh); return sha256_file($path) unless $hdr;
            $csv->column_names(@$hdr);
            my @rows;
            while (my $r=$csv->getline_hr($fh)) {
                my %c;
                for my $k (@$hdr) {
                    my $lk=lc($k);
                    next if $VOLATILE_JSON_KEYS{$lk} || $lk =~ /timestamp|elapsed|runtime|output_path/;
                    my $v=defined($r->{$k})?"$r->{$k}":'';
                    $v =~ s/\Q$repo\E/<REPO>/g;
                    $c{$k}=$v;
                }
                push @rows,\%c;
            }
            close $fh;
            @rows=sort { $JSON_CANON->encode($a) cmp $JSON_CANON->encode($b) } @rows;
            return sha256_hex($JSON_CANON->encode(\@rows));
        }
        if ($ext =~ /^\.(?:txt|log|md|html)$/) {
            my $txt=slurp($path); my $rf=$repo; $rf =~ s!\\!/!g;
            $txt =~ s/\Q$repo\E/<REPO>/g; $txt =~ s/\Q$rf\E/<REPO>/g;
            return sha256_hex($txt);
        }
        return sha256_file($path);
    };
    return $ok if defined $ok && !$@;
    return sha256_file($path);
}

sub artifact_snapshot {
    my ($repo,$patterns) = @_;
    my %snap;
    for my $pat (@$patterns) {
        $pat=resolve_placeholders($pat,{REPO=>$repo});
        my $abs=File::Spec->file_name_is_absolute($pat) ? $pat : File::Spec->catfile($repo,$pat);
        for my $f (bsd_glob($abs, GLOB_NOSORT)) {
            next unless -f $f;
            my $af=abs_path($f) || $f;
            my $rel=File::Spec->abs2rel($af,$repo); $rel =~ s!\\!/!g;
            $snap{$rel}={path=>$af,size=>(-s $af),sha256=>sha256_file($af),canonical_sha256=>canonical_hash($af,$repo)};
        }
    }
    return \%snap;
}

sub compare_artifacts {
    my ($ref,$obs) = @_;
    return (undef,0,0,[]) unless $ref && %$ref;
    my ($compared,$matched)=(0,0); my @mis;
    for my $k (keys %$ref) {
        if (!exists $obs->{$k}) { push @mis,"missing:$k"; next }
        $compared++;
        if (($ref->{$k}{canonical_sha256}//'') eq ($obs->{$k}{canonical_sha256}//'')) { $matched++ }
        else { push @mis,"mismatch:$k" }
    }
    my $ok = ($matched == scalar(keys %$ref) && !@mis) ? JSON::PP::true : JSON::PP::false;
    return ($ok,$matched,scalar(keys %$ref),\@mis);
}

sub prompt_preamble {
    my ($run_id,$preferred) = @_;
    my %line=(
        mcp=>'For this task, use the configured MultiGWAS-Explorer MCP tool auto_prepare_and_run_diff_gwas as the execution interface. Do not bypass it with direct shell execution unless the MCP tool itself cannot express the request.',
        native=>'For this task, use the repository-provided scripts/skill directly when the requested function is not exposed by the main MCP automation tool. Do not reimplement the scientific logic.',
        mixed=>'For this task, use the MCP automation path for the SAS workflow and the repository-provided local script for the gnuplot workflow.',
        clarify=>'This task may contain incompatible instructions. Do not execute anything until the contradiction is resolved.'
    );
    my $iface=$line{$preferred}//"Use the repository's versioned tools rather than reimplementing the analysis.";
    return "BENCHMARK RUN $run_id\n\nYou are being evaluated on faithful orchestration of the existing MultiGWAS-Explorer repository.\n".
           "- Do not edit source code, benchmark files, or manuscript files.\n- Do not install or upgrade software.\n".
           "- Do not invent or silently substitute a missing spec, SNP, GWAS track, path, or column mapping.\n".
           "- Use only the existing repository pipeline and its configured MCP server/skill.\n".
           "- If an MCP call starts a background job, follow the server's polling guidance until a terminal status is obtained.\n".
           "- Treat intentionally invalid inputs as safe failures, not as requests to repair or replace them.\n".
           "- At the end, state only the observed status and artifact/log paths; do not claim execution that did not occur.\n$iface\n\nUSER BENCHMARK PROMPT:\n";
}

sub codex_cli_path {
    my ($path)=@_;
    return $path unless $^O eq 'cygwin';
    if (open my $fh, '-|', 'cygpath', '-w', $path) {
        my $converted=<$fh>;
        close $fh;
        if (defined $converted) {
            $converted =~ s/[\r\n]+$//;
            return $converted if length $converted;
        }
    }
    return $path;
}

sub build_codex_cmd {
    my ($command,$repo,$prompt,$model,$sandbox,$approval_mode,$writable_dirs,$reasoning_effort)=@_;
    $sandbox ||= 'workspace-write';
    $approval_mode ||= 'sandbox';
    my $codex_repo=codex_cli_path($repo);
    my @cmd=($command,'exec','--json','--ephemeral','--skip-git-repo-check');
    if ($approval_mode eq 'approve-for-me') {
        push @cmd,'--approve-for-me';
    } elsif ($approval_mode eq 'sandbox') {
        push @cmd,('--sandbox',$sandbox);
    } else {
        die "Unknown Codex approval_mode '$approval_mode' (expected sandbox or approve-for-me)\n";
    }
    push @cmd,('-C',$codex_repo);
    for my $dir (@{$writable_dirs||[]}) {
        next unless defined($dir) && length($dir);
        push @cmd,('--add-dir',codex_cli_path($dir));
    }
    push @cmd,('--model',$model) if defined($model) && length($model);
    push @cmd,('-c',qq{model_reasoning_effort="$reasoning_effort"})
        if defined($reasoning_effort) && $reasoning_effort =~ /^(?:low|medium|high|xhigh)$/;
    push @cmd,$prompt;
    return \@cmd;
}

sub codex_unavailable_message {
    my ($text)=@_;
    return '' unless defined($text) && $text =~ $CODEX_UNAVAILABLE_RE;
    my @lines=grep { /$CODEX_UNAVAILABLE_RE/ } split /\r?\n/, $text;
    my $message=@lines ? $lines[0] : 'Codex is unavailable.';
    $message =~ s/^\s+|\s+$//g;
    $message=substr($message,0,1000) if length($message)>1000;
    return $message;
}

sub check_codex_access {
    my (%a)=@_;
    my $dir=File::Spec->catdir($a{out},'codex_access_check');
    make_path($dir);
    my $stdout=File::Spec->catfile($dir,'stdout.jsonl');
    my $stderr=File::Spec->catfile($dir,'stderr.txt');
    my $command=$a{agent_cfg}{command}//'codex';
    die "Required command not found for codex: $command\n" unless which_cmd($command);
    my $cmd=build_codex_cmd(
        $command,$a{repo},
        'Availability check only. Do not use tools. Reply exactly CODEX_ACCESS_OK.',
        $a{agent_cfg}{model}//'',
        $a{agent_cfg}{sandbox}//'workspace-write',
        $a{agent_cfg}{approval_mode}//'sandbox',
        $a{agent_cfg}{writable_dirs}||[],
        $a{agent_cfg}{reasoning_effort}//''
    );
    my ($code,$elapsed,$timed_out)=run_capture(
        cmd=>$cmd,cwd=>$a{repo},timeout=>$a{timeout}//120,
        stdout=>$stdout,stderr=>$stderr
    );
    my $text=(-e $stdout?slurp($stdout):'')."\n".(-e $stderr?slurp($stderr):'');
    my $unavailable=codex_unavailable_message($text);
    my $accessible=($code==0 && !$timed_out && !$unavailable && $text=~/CODEX_ACCESS_OK/)?1:0;
    my $result={
        checked_at=>utcnow(),accessible=>$accessible?JSON::PP::true:JSON::PP::false,
        exit_code=>$code,timed_out=>$timed_out,elapsed_seconds=>$elapsed,
        unavailable_message=>$unavailable,stdout_path=>$stdout,stderr_path=>$stderr
    };
    write_json(File::Spec->catfile($dir,'result.json'),$result);
    unless ($accessible) {
        my $reason=$unavailable || ($timed_out ? 'Codex availability check timed out.' : "Codex availability check failed (exit code $code or missing CODEX_ACCESS_OK response).");
        die "WARNING: Codex is not accessible; benchmark execution has been stopped before any scientific test.\nReason: $reason\nDetails: ".File::Spec->catfile($dir,'result.json')."\n";
    }
    print "[access] Codex availability check passed.\n";
    return $result;
}

sub maybe_configure_mcp {
    my ($cfg,$repo,$out)=@_;
    my $mcp=$cfg->{mcp}||{};
    return {attempted=>JSON::PP::false,note=>'configure_clients=false'} unless $mcp->{configure_clients};
    my $name=$mcp->{name}//'perl-bio'; my $url=$mcp->{url}//'http://127.0.0.1:8080/mcp';
    my %cmds=(codex=>[$cfg->{agents}{codex}{command}//'codex','mcp','add',$name,'--url',$url]);
    my %res;
    for my $agent (qw(codex)) {
        my $cmd=$cmds{$agent};
        if (!which_cmd($cmd->[0])) { $res{$agent}={exit_code=>127,output=>'command not found'}; next }
        my ($code,$txt)=capture_text(cmd=>$cmd,cwd=>$repo,timeout=>60);
        $txt=substr($txt,-5000) if length($txt)>5000;
        $res{$agent}={exit_code=>$code,output=>$txt};
    }
    write_json(File::Spec->catfile($out,'mcp_configuration.json'),\%res);
    return {attempted=>JSON::PP::true,%res};
}

sub start_mcp_if_needed {
    my ($cfg,$repo,$out)=@_;
    my $mcp=$cfg->{mcp}||{};
    return undef if exists($mcp->{start_server}) && !$mcp->{start_server};
    my $url=$mcp->{url}//'http://127.0.0.1:8080/mcp';
    my ($host,$port)=parse_mcp_url($url);
    return undef if wait_for_port($host,$port,1);
    my $server_dir=$mcp->{server_dir}//$repo;
    $server_dir=File::Spec->rel2abs($server_dir,$repo)
        unless File::Spec->file_name_is_absolute($server_dir);
    die "MCP server directory does not exist: $server_dir\n" unless -d $server_dir;
    my $server_script=File::Spec->catfile($server_dir,'server.pl');
    die "MCP server script does not exist: $server_script\n" unless -f $server_script;
    my $log=File::Spec->catfile($out,'mcp_server.log');
    my $pid=fork(); die "fork failed starting MCP: $!\n" unless defined $pid;
    if ($pid==0) {
        chdir $server_dir or exit 127;
        open STDOUT,'>>:raw',$log or exit 127; open STDERR,'>&STDOUT' or exit 127;
        exec 'perl','server.pl','daemon','-m','production','-l',"http://$host:$port";
        exit 127;
    }
    my $startup_timeout=$mcp->{startup_timeout_seconds}//60;
    if (!wait_for_port($host,$port,$startup_timeout)) { kill 'TERM',$pid; die "MCP server did not open $host:$port within ${startup_timeout}s; see $log\n" }
    return $pid;
}

# Deterministic, language-independent pseudo-random ordering based on SHA-256.
sub stable_order {
    my ($rows,$seed,$salt)=@_;
    return [ sort { sha256_hex("$seed|$salt|$a->{prompt_id}") cmp sha256_hex("$seed|$salt|$b->{prompt_id}") } @$rows ];
}

sub select_prompts {
    my (%a)=@_;
    my @rows=@{$a{rows}};
    if ($a{prompt_ids}) { my %s=map {$_=>1}@{$a{prompt_ids}}; @rows=grep {$s{$_->{prompt_id}}}@rows }
    if ($a{families}) { my %s=map {$_=>1}@{$a{families}}; @rows=grep {$s{$_->{id}}}@rows }
    return \@rows if $a{mode} eq 'all';
    my %fam; push @{$fam{$_->{id}}},$_ for @rows;
    my @sel;
    if ($a{mode} eq 'sample') {
        for my $fid (sort keys %fam) {
            my $v=stable_order($fam{$fid},$a{seed},$fid);
            my $n=$a{per_family}; $n=@$v if $n>@$v;
            push @sel,@$v[0..$n-1] if $n>0;
        }
        return stable_order(\@sel,$a{seed},'sample-final');
    }
    if ($a{mode} eq 'fold') {
        die "fold mode requires --fold 1..4\n" unless $a{fold} && $a{fold}>=1 && $a{fold}<=4;
        for my $fid (sort keys %fam) {
            my $v=stable_order($fam{$fid},$a{seed},$fid);
            die "$fid has fewer than four paraphrases\n" if @$v<4;
            push @sel,$v->[$a{fold}-1];
        }
        return stable_order(\@sel,$a{seed},'fold-final-'.$a{fold});
    }
    die "Unknown selection mode: $a{mode}\n";
}

sub gold_to_cli {
    my ($action)=@_;
    my $entry=$action->{entry}//''; my $a=$action->{args}||{};
    return ['bash','install/check_pipeline_install.sh'] if $entry eq 'install_check';
    return ['perl','DiffGWASDeps/test_top_hit_maf_filter.pl','--no-real','--keep-workdir'] if $entry eq 'maf_test';
    return ['perl','./DiffGWASDeps/run_sas_codes_or_script_in_ODA.pl','--check-sas-oda-login-only'] if $entry eq 'sas_login_check';
    return undef unless $entry eq 'sas_pipeline' || $entry eq 'gnuplot_pipeline';
    my $script=$entry eq 'sas_pipeline'?'./auto_prepare_and_run_diff_gwas.pl':'./auto_prepare_and_run_diff_gwas_with_gunplot.pl';
    my @cmd=('perl',$script);
    my %map=(spec=>'--spec',gwas_dir=>'--gwas-dir',spec_out=>'--spec-out',raw_column_alias_config=>'--raw-column-alias-config',
             mode=>'--mode',target_snps=>'--target-snps',target_snp_genes=>'--target-snp-genes',display_gwas=>'--display-gwas',
             local_gtf_window_bp=>'--local-gtf-window-bp');
    for my $k (qw(spec gwas_dir spec_out raw_column_alias_config mode target_snps target_snp_genes display_gwas local_gtf_window_bp)) {
        push @cmd,($map{$k},"$a->{$k}") if exists $a->{$k};
    }
    push @cmd,'--generate-spec-only' if $a->{generate_spec_only};
    push @cmd,'--preview-spec' if $a->{preview_spec};
    push @cmd,'--force' if $a->{force};
    push @cmd,'--get-common-associations' if $a->{common};
    for my $s (@{$a->{steps}||[]}) { push @cmd,('--step',$s) }
    my $plots=$a->{plots}||[];
    if ($entry eq 'gnuplot_pipeline' && @$plots) { push @cmd,('--plots',join(',',@$plots)) }
    elsif ($entry eq 'sas_pipeline') { for my $p (@$plots) { push @cmd,('--step','plot_'.$p) } }
    return \@cmd;
}

sub run_cli_references {
    my (%a)=@_;
    my %families=map { $_->{id}=>1 } grep { ($_->{gold_status}//'') ne 'clarify' } @{$a{selected}};
    my %refs;
    for my $fid (sort keys %families) {
        my $spec=$a{gold}{$fid}||{}; my $actions=$spec->{actions}||[];
        my $fdir=File::Spec->catdir($a{out},'cli_reference',$fid); make_path($fdir);
        my (@results,$all_ok); $all_ok=1; my $total_elapsed=0;
        my $j=0;
        for my $action (@$actions) {
            $j++; my $cmd=gold_to_cli($action);
            if (!$cmd) { push @results,{skipped=>JSON::PP::true,reason=>'no executable gold command',action=>$action}; next }
            my $stdout=File::Spec->catfile($fdir,sprintf('action_%02d.stdout.txt',$j));
            my $stderr=File::Spec->catfile($fdir,sprintf('action_%02d.stderr.txt',$j));
            my ($code,$elapsed,$to)=run_capture(cmd=>$cmd,cwd=>$a{repo},timeout=>$a{timeout},stdout=>$stdout,stderr=>$stderr);
            $total_elapsed+=$elapsed;
            my ($row)=grep { $_->{id} eq $fid } @{$a{selected}};
            my $safe=(($row->{gold_status}//'success') eq 'safe_failure');
            my $txt=(-e $stdout?slurp($stdout):'').(-e $stderr?slurp($stderr):'');
            my $ok=$safe ? ($code!=0 || $txt=~$FAILURE_RE) : ($code==0);
            $all_ok &&= $ok;
            push @results,{command=>$cmd,exit_code=>$code,elapsed_seconds=>$elapsed,timed_out=>$to,reference_ok=>$ok?JSON::PP::true:JSON::PP::false};
        }
        my $patterns=$a{artifact_map}{$fid}||[];
        my $snap=@$patterns ? artifact_snapshot($a{repo},$patterns) : {};
        write_json(File::Spec->catfile($fdir,'artifacts.json'),$snap);
        $refs{$fid}={reference_ok=>$all_ok?JSON::PP::true:JSON::PP::false,elapsed_seconds=>$total_elapsed,actions=>\@results,artifacts=>$snap};
        write_json(File::Spec->catfile($fdir,'meta.json'),$refs{$fid});
    }
    return \%refs;
}

sub load_existing_cli_refs {
    my ($out,$selected)=@_; my %refs; my %f=map {$_->{id}=>1} @$selected;
    for my $fid (sort keys %f) {
        my $p=File::Spec->catfile($out,'cli_reference',$fid,'meta.json');
        $refs{$fid}=load_json($p) if -e $p;
    }
    return \%refs;
}

sub collect_existing_results {
    my ($out)=@_; my @res;
    my $root=File::Spec->catdir($out,'runs'); return [] unless -d $root;
    find(sub {
        return unless $_ eq 'meta.json' && $File::Find::name =~ m![\\/]rep_\d+[\\/]meta\.json$!;
        my $r=eval { load_json($File::Find::name) };
        push @res,$r if ref($r) eq 'HASH' && $r->{finished};
    },$root);
    @res=sort { ($a->{agent}//'') cmp ($b->{agent}//'') || ($a->{prompt_id}//'') cmp ($b->{prompt_id}//'') || ($a->{replicate}//0)<=>($b->{replicate}//0) } @res;
    return \@res;
}

sub run_one_agent {
    my (%a)=@_;
    my $pid=$a{prompt_row}{prompt_id};
    my $rdir=File::Spec->catdir($a{out},'runs',$a{agent},$pid,sprintf('rep_%02d',$a{rep})); make_path($rdir);
    my $meta=File::Spec->catfile($rdir,'meta.json');
    if (-e $meta) { my $old=eval{load_json($meta)}; return $old if ref($old) eq 'HASH' && $old->{finished} }
    my $preferred=$a{gold_spec}{preferred_interface}//'native';
    my $prompt=prompt_preamble("$a{agent}-$pid-r$a{rep}",$preferred).resolve_placeholders($a{prompt_row}{prompt_text},$a{placeholders});
    spew(File::Spec->catfile($rdir,'prompt.txt'),$prompt);
    my $stdout=File::Spec->catfile($rdir,'stdout.jsonl'); my $stderr=File::Spec->catfile($rdir,'stderr.txt');
    my $command=$a{agent_cfg}{command}//$a{agent}; my $model=$a{agent_cfg}{model}//'';
    die "Unknown agent $a{agent}\n" unless $a{agent} eq 'codex';
    my $cmd=build_codex_cmd($command,$a{repo},$prompt,$model,
                            $a{agent_cfg}{sandbox}//'workspace-write',
                            $a{agent_cfg}{approval_mode}//'sandbox',
                            $a{agent_cfg}{writable_dirs}||[],
                            $a{agent_cfg}{reasoning_effort}//'');
    my $tracked_before=git_diff_fingerprint($a{repo}); my $start=utcnow();
    my ($code,$elapsed,$timed_out)=run_capture(cmd=>$cmd,cwd=>$a{repo},timeout=>$a{timeout},stdout=>$stdout,stderr=>$stderr);
    my $end=utcnow(); my $tracked_after=git_diff_fingerprint($a{repo});
    my $events=parse_jsonl_file($stdout); my $parsed=codex_parse($events);
    my ($task_match)=match_action_set($a{gold_spec}{actions}||[],$parsed->{actions});
    my $expected=$a{prompt_row}{gold_status}//'success';
    my $all_text=$parsed->{all_text}."\n".(-e $stderr?slurp($stderr):'');
    my $agent_unavailable=codex_unavailable_message($all_text);
    my $failure=($all_text=~$FAILURE_RE)?1:0;
    for my $c (@{$parsed->{commands}}) { $failure=1 if defined($c->{exit_code}) && $c->{exit_code}!=0 }
    my $clarify=($parsed->{final_text}=~$CLARIFY_RE)?1:0;
    my $tracked_changed=($tracked_before->{sha256}//'') ne ($tracked_after->{sha256}//'');
    my $source_edit=(@{$parsed->{file_changes}} || $tracked_changed)?1:0;
    my ($scored,$safe_correct)=(0,undef);
    if ($expected eq 'success') { $scored=($code==0 && !$timed_out && $task_match && $parsed->{terminal_ok} && !$source_edit)?1:0 }
    elsif ($expected eq 'safe_failure') { $scored=($code==0 && !$timed_out && $task_match && $failure && !$source_edit)?1:0; $safe_correct=$scored }
    elsif ($expected eq 'clarify') { $scored=($code==0 && !$timed_out && !@{$parsed->{actions}} && $clarify && !$source_edit)?1:0 }
    my $obs=@{$a{artifact_patterns}}?artifact_snapshot($a{repo},$a{artifact_patterns}):{};
    my $ref=($a{cli_ref}||{})->{artifacts}||{};
    my ($art_parity,$art_match,$art_total,$art_mis)=compare_artifacts($ref,$obs);
    my $result={
        finished=>JSON::PP::true,timestamp_start=>$start,timestamp_end=>$end,agent=>$a{agent},model=>length($model)?$model:'default',
        prompt_id=>$pid,family_id=>$a{prompt_row}{id},family_name=>$a{prompt_row}{name}//'',variant=>$a{prompt_row}{variant},
        category=>$a{prompt_row}{category}//'',difficulty=>$a{prompt_row}{difficulty}//'',dataset=>$a{prompt_row}{dataset}//'',backend=>$a{prompt_row}{backend}//'',
        preferred_interface=>$preferred,expected_status=>$expected,replicate=>$a{rep},exit_code=>$code,timed_out=>$timed_out,
        elapsed_seconds=>$elapsed,process_ok=>($code==0&&!$timed_out)?JSON::PP::true:JSON::PP::false,task_match=>$task_match?JSON::PP::true:JSON::PP::false,
        terminal_execution_seen=>$parsed->{terminal_ok},failure_detected=>$failure?JSON::PP::true:JSON::PP::false,
        clarification_detected=>$clarify?JSON::PP::true:JSON::PP::false,source_edit_detected=>$source_edit?JSON::PP::true:JSON::PP::false,
        scored_success=>$scored?JSON::PP::true:JSON::PP::false,safe_failure_correct=>defined($safe_correct)?($safe_correct?JSON::PP::true:JSON::PP::false):undef,
        action_signature=>action_signature($parsed->{actions}),actual_actions=>$parsed->{actions},mcp_call_count=>scalar(@{$parsed->{mcp_calls}}),
        shell_command_count=>scalar(@{$parsed->{commands}}),file_change_count=>scalar(@{$parsed->{file_changes}}),
        tracked_diff_changed=>$tracked_changed?JSON::PP::true:JSON::PP::false,tracked_changed_files_before=>$tracked_before->{tracked_changed_files},
        tracked_changed_files_after=>$tracked_after->{tracked_changed_files},final_text=>$parsed->{final_text},usage=>$parsed->{usage},
        artifact_parity=>$art_parity,artifact_matched=>$art_match,artifact_reference_count=>$art_total,artifact_mismatches=>$art_mis,
        artifact_snapshot=>$obs,stdout_path=>$stdout,stderr_path=>$stderr
    };
    $result->{agent_unavailable}=$agent_unavailable?JSON::PP::true:JSON::PP::false;
    $result->{agent_unavailable_message}=$agent_unavailable;
    write_json($meta,$result);
    write_json(File::Spec->catfile($rdir,'tool_calls.json'),{actions=>$parsed->{actions},mcp_calls=>$parsed->{mcp_calls},commands=>$parsed->{commands},file_changes=>$parsed->{file_changes}});
    write_json(File::Spec->catfile($rdir,'artifacts.json'),$obs);
    return $result;
}

sub flatten_csv_value {
    my ($v)=@_; return '' unless defined $v;
    return $JSON_CANON->encode($v) if ref($v);
    return "$v";
}

sub write_csv_rows {
    my ($path,$rows,$cols,$sep)=@_; $sep//=',';
    my $csv=Text::CSV->new({binary=>1,eol=>"\n",sep_char=>$sep});
    open my $fh,'>:encoding(UTF-8)',$path or die "Cannot write $path: $!\n";
    $csv->print($fh,$cols);
    for my $r (@$rows) { $csv->print($fh,[map { flatten_csv_value($r->{$_}) } @$cols]) }
    close $fh;
}

sub write_results {
    my ($out,$results)=@_;
    my $jsonl=File::Spec->catfile($out,'run_results.jsonl');
    open my $fh,'>:encoding(UTF-8)',$jsonl or die $!;
    print {$fh} $JSON->encode($_)."\n" for @$results; close $fh;
    my @cols=qw(agent model prompt_id family_id family_name variant category difficulty dataset backend preferred_interface expected_status replicate exit_code timed_out elapsed_seconds process_ok task_match terminal_execution_seen failure_detected clarification_detected source_edit_detected scored_success safe_failure_correct action_signature mcp_call_count shell_command_count file_change_count tracked_diff_changed tracked_changed_files_before tracked_changed_files_after artifact_parity artifact_matched artifact_reference_count artifact_mismatches timestamp_start timestamp_end stdout_path stderr_path);
    write_csv_rows(File::Spec->catfile($out,'run_results.csv'),$results,\@cols,',');
}

sub agent_summary {
    my ($results,$agent)=@_; my @rs=grep {$_->{agent} eq $agent} @$results;
    my @pos=grep {$_->{expected_status} eq 'success'} @rs; my @neg=grep {$_->{expected_status} eq 'safe_failure'} @rs; my @clar=grep {$_->{expected_status} eq 'clarify'} @rs;
    my @lat=map {0+$_->{elapsed_seconds}} grep {!$_->{timed_out}} @rs;
    my $task=scalar grep {$_->{task_match}} @rs; my $succ=scalar grep {$_->{scored_success}} @pos; my $safe=scalar grep {$_->{safe_failure_correct}} @neg; my $clarok=scalar grep {$_->{scored_success}} @clar;
    # Artifact hashes are scientifically interpretable only after the agent
    # selected the preregistered action.  Otherwise an unchanged reference
    # artifact can create a false parity result when the agent ran a different
    # specification and wrote different output paths.
    my @art=grep {$_->{task_match} && defined $_->{artifact_parity}} @rs;
    my $artok=scalar grep {$_->{artifact_parity}} @art;
    my (%byp,%byf); push @{$byp{$_->{prompt_id}}},$_ for @rs; push @{$byf{$_->{family_id}}},$_ for @rs;
    my $repeat=0;
    my $repeat_assessed=0;
    for my $rows (values %byp) {
        next unless @$rows >= 2;
        $repeat_assessed++;
        my %s=map {(($_->{action_signature}//'').'|'.($_->{scored_success}?1:0))=>1} @$rows;
        $repeat++ if keys(%s)==1;
    }
    my ($pa,$pc)=(0,0);
    for my $rows (values %byf) { my %p=map {$_->{prompt_id}=>1} @$rows; if (keys(%p)>=4) { $pa++; my $ok=1; $ok&&=($_->{task_match}&&$_->{scored_success}) for @$rows; $pc++ if $ok } }
    my ($q1,$q3)=iqr(\@lat);
    return {
        agent=>$agent,runs=>scalar(@rs),overall_success=>scalar(grep {$_->{scored_success}} @rs),positive_runs=>scalar(@pos),positive_success=>$succ,completion_rate=>@pos?$succ/@pos:undef,failure_rate=>@pos?1-$succ/@pos:undef,
        task_match=>$task,task_match_rate=>@rs?$task/@rs:undef,safe_failure_runs=>scalar(@neg),safe_failure_correct=>$safe,safe_failure_rate=>@neg?$safe/@neg:undef,
        clarification_runs=>scalar(@clar),clarification_correct=>$clarok,prompts_repeated=>$repeat_assessed,repeat_consistent_prompts=>$repeat,repeat_consistency_rate=>$repeat_assessed?$repeat/$repeat_assessed:undef,
        families=>scalar(keys%byf),paraphrase_assessed_families=>$pa,paraphrase_consistent_families=>$pc,paraphrase_consistency_rate=>$pa?$pc/$pa:undef,
        artifact_assessed_runs=>scalar(@art),artifact_parity_runs=>$artok,artifact_parity_rate=>@art?$artok/@art:undef,
        median_latency=>median(\@lat),q1_latency=>$q1,q3_latency=>$q3,p95_latency=>percentile(\@lat,.95),
        source_edit_runs=>scalar(grep {$_->{source_edit_detected}} @rs),timeouts=>scalar(grep {$_->{timed_out}} @rs),
        mcp_calls=>do{my $z=0;$z+=0+($_->{mcp_call_count}//0) for @rs;$z},shell_commands=>do{my $z=0;$z+=0+($_->{shell_command_count}//0) for @rs;$z}
    };
}

sub summarize {
    my ($out,$results,$pf)=@_;
    my %agents=map {$_->{agent}=>1} @$results; my @sums=map {agent_summary($results,$_)} sort keys %agents;
    my @agent_cols=qw(agent runs overall_success positive_runs positive_success completion_rate failure_rate task_match task_match_rate safe_failure_runs safe_failure_correct safe_failure_rate clarification_runs clarification_correct prompts_repeated repeat_consistent_prompts repeat_consistency_rate families paraphrase_assessed_families paraphrase_consistent_families paraphrase_consistency_rate artifact_assessed_runs artifact_parity_runs artifact_parity_rate median_latency q1_latency q3_latency p95_latency source_edit_runs timeouts mcp_calls shell_commands);
    write_csv_rows(File::Spec->catfile($out,'summary_by_agent.csv'),\@sums,\@agent_cols,',') if @sums;
    my (%by,@fam_rows);
    push @{$by{$_->{family_id}.'\0'.$_->{agent}}},$_ for @$results;
    for my $key (sort keys %by) {
        my $rs=$by{$key}; my @lat=map {0+$_->{elapsed_seconds}} grep {!$_->{timed_out}} @$rs;
        push @fam_rows,{family_id=>$rs->[0]{family_id},family_name=>$rs->[0]{family_name},agent=>$rs->[0]{agent},n_runs=>scalar(@$rs),task_match=>scalar(grep{$_->{task_match}}@$rs),scored_success=>scalar(grep{$_->{scored_success}}@$rs),failures=>scalar(grep{!$_->{scored_success}}@$rs),median_latency=>median(\@lat)//'',expected_status=>$rs->[0]{expected_status},preferred_interface=>$rs->[0]{preferred_interface}};
    }
    write_csv_rows(File::Spec->catfile($out,'summary_by_family.csv'),\@fam_rows,[qw(family_id family_name agent n_runs task_match scored_success failures median_latency expected_status preferred_interface)],',') if @fam_rows;

    my (%bi,@interface_rows);
    push @{$bi{$_->{agent}.'\0'.$_->{preferred_interface}}},$_ for @$results;
    for my $key (sort keys %bi) {
        my $rs=$bi{$key}; my @pos=grep{$_->{expected_status} eq 'success'}@$rs; my @neg=grep{$_->{expected_status} eq 'safe_failure'}@$rs; my @lat=map{0+$_->{elapsed_seconds}}grep{!$_->{timed_out}}@$rs; my @art=grep{$_->{task_match} && defined$_->{artifact_parity}}@$rs; my($q1,$q3)=iqr(\@lat);
        push @interface_rows,{agent=>$rs->[0]{agent},interface=>$rs->[0]{preferred_interface},runs=>scalar(@$rs),positive_runs=>scalar(@pos),positive_success=>scalar(grep{$_->{scored_success}}@pos),task_match=>scalar(grep{$_->{task_match}}@$rs),safe_failure_runs=>scalar(@neg),safe_failure_correct=>scalar(grep{$_->{safe_failure_correct}}@neg),artifact_assessed_runs=>scalar(@art),artifact_parity_runs=>scalar(grep{$_->{artifact_parity}}@art),median_latency=>median(\@lat),q1_latency=>$q1,q3_latency=>$q3,p95_latency=>percentile(\@lat,.95)};
    }
    write_csv_rows(File::Spec->catfile($out,'summary_by_interface.csv'),\@interface_rows,[qw(agent interface runs positive_runs positive_success task_match safe_failure_runs safe_failure_correct artifact_assessed_runs artifact_parity_runs median_latency q1_latency q3_latency p95_latency)],',') if @interface_rows;

    my @mcp=grep{$_->{interface} eq 'mcp'}@interface_rows;
    my @mcp_cols=('Agent','MCP runs','Positive tasks completed','Completion rate (95% CI)','Prespecified invocation concordance (95% CI)','Safe-failure detection (95% CI)','Canonical artifact parity','Median end-to-end latency (IQR), s','95th percentile latency, s');
    my @mcp_out;
    for my $s (@mcp) { push @mcp_out,{ 'Agent'=>$s->{agent},'MCP runs'=>$s->{runs},'Positive tasks completed'=>"$s->{positive_success}/$s->{positive_runs}",'Completion rate (95% CI)'=>fmt_pct($s->{positive_success},$s->{positive_runs}),'Prespecified invocation concordance (95% CI)'=>fmt_pct($s->{task_match},$s->{runs}),'Safe-failure detection (95% CI)'=>fmt_pct($s->{safe_failure_correct},$s->{safe_failure_runs}),'Canonical artifact parity'=>$s->{artifact_assessed_runs}?"$s->{artifact_parity_runs}/$s->{artifact_assessed_runs}":'not assessed','Median end-to-end latency (IQR), s'=>fmt_num($s->{median_latency}).' ('.fmt_num($s->{q1_latency}).'-'.fmt_num($s->{q3_latency}).')','95th percentile latency, s'=>fmt_num($s->{p95_latency})} }
    write_csv_rows(File::Spec->catfile($out,'supplementary_table_mcp_core.tsv'),\@mcp_out,\@mcp_cols,"\t");

    my @perf_cols=('Agent','Runs','All criteria met','Overall success rate (95% CI)','Positive tasks completed','Completion rate (95% CI)','Task/invocation concordance (95% CI)','Repeated-run consistency','Paraphrase consistency','Safe-failure detection (95% CI)','Clarification handling (95% CI)','Canonical artifact parity','Median latency (IQR), s','95th percentile latency, s');
    my @perf;
    for my $s (@sums) { push @perf,{'Agent'=>$s->{agent},'Runs'=>$s->{runs},'All criteria met'=>"$s->{overall_success}/$s->{runs}",'Overall success rate (95% CI)'=>fmt_pct($s->{overall_success},$s->{runs}),'Positive tasks completed'=>"$s->{positive_success}/$s->{positive_runs}",'Completion rate (95% CI)'=>fmt_pct($s->{positive_success},$s->{positive_runs}),'Task/invocation concordance (95% CI)'=>fmt_pct($s->{task_match},$s->{runs}),'Repeated-run consistency'=>$s->{prompts_repeated}?"$s->{repeat_consistent_prompts}/$s->{prompts_repeated}":'not assessed','Paraphrase consistency'=>$s->{paraphrase_assessed_families}?"$s->{paraphrase_consistent_families}/$s->{paraphrase_assessed_families}":'not assessed','Safe-failure detection (95% CI)'=>fmt_pct($s->{safe_failure_correct},$s->{safe_failure_runs}),'Clarification handling (95% CI)'=>fmt_pct($s->{clarification_correct},$s->{clarification_runs}),'Canonical artifact parity'=>$s->{artifact_assessed_runs}?"$s->{artifact_parity_runs}/$s->{artifact_assessed_runs}":'not assessed','Median latency (IQR), s'=>fmt_num($s->{median_latency}).' ('.fmt_num($s->{q1_latency}).'-'.fmt_num($s->{q3_latency}).')','95th percentile latency, s'=>fmt_num($s->{p95_latency})} }
    write_csv_rows(File::Spec->catfile($out,'supplementary_table_agent_performance.tsv'),\@perf,\@perf_cols,"\t");

    my $commit=$pf->{git_commit}//'unknown';
    my $working_tree_modified=length($pf->{git_status}//'') ? 1 : 0;
    my $revision_note=$working_tree_modified
      ? "repository commit `$commit` plus recorded working-tree modifications (tracked-diff SHA-256 `$pf->{git_diff_sha256}`)"
      : "frozen repository commit `$commit`";
    my %p=map{$_->{prompt_id}=>1}@$results; my %f=map{$_->{family_id}=>1}@$results; my %r=map{$_->{replicate}=>1}@$results;
    my $prompt_n=scalar(keys%p); my $family_n=scalar(keys%f); my $rep_n=scalar(keys%r);
    my $safe_n=scalar grep {($_->{expected_status}//'') eq 'safe_failure'} @$results;
    my $clarify_n=scalar grep {($_->{expected_status}//'') eq 'clarify'} @$results;
    my $family_label=$family_n==1?'task family':'task families';
    my @lines=('# Reviewer-response benchmark summary','',"Repository commit: `$commit`  ",($working_tree_modified ? "Tracked working-tree diff SHA-256: `$pf->{git_diff_sha256}`  " : ()),"Evaluated prompts: $prompt_n across $family_n $family_label  ",'Replicates represented: '.join(', ',sort{$a<=>$b}keys%r),'','## Ready-to-adapt response to the reviewer','',"We conducted a focused AI-interface evaluation using $revision_note. We evaluated $prompt_n natural-language prompts representing $family_n $family_label with fresh Codex sessions. The command-line workflow was treated as the scientific reference. We recorded whether Codex selected the prespecified repository action, whether execution reached a terminal outcome, source-edit violations, and end-to-end wall time. For families configured with explicit scientific artifact paths, task-concordant outputs were additionally canonicalized and compared with the corresponding CLI reference hashes.",'');
    for my $s (@sums) {
        my $line='**'.ucfirst($s->{agent}).'.** Overall benchmark success was '.fmt_pct($s->{overall_success},$s->{runs}).'; positive-task completion was '.fmt_pct($s->{positive_success},$s->{positive_runs}).'; prespecified action/invocation concordance was '.fmt_pct($s->{task_match},$s->{runs}).'; ';
        $line.='repeated-run consistency was '.fmt_pct($s->{repeat_consistent_prompts},$s->{prompts_repeated}).'; ' if $s->{prompts_repeated};
        $line.='paraphrase consistency was '.fmt_pct($s->{paraphrase_consistent_families},$s->{paraphrase_assessed_families}).'; ' if $s->{paraphrase_assessed_families};
        $line.='safe-failure detection was '.fmt_pct($s->{safe_failure_correct},$s->{safe_failure_runs}).'; ' if $s->{safe_failure_runs};
        $line.='clarification handling was '.fmt_pct($s->{clarification_correct},$s->{clarification_runs}).'; ' if $s->{clarification_runs};
        $line.='median end-to-end latency was '.fmt_num($s->{median_latency}).' s (IQR '.fmt_num($s->{q1_latency}).'-'.fmt_num($s->{q3_latency}).' s; p95 '.fmt_num($s->{p95_latency}).' s). ';
        $line.=$s->{artifact_assessed_runs}?"Canonical artifact parity was observed in $s->{artifact_parity_runs}/$s->{artifact_assessed_runs} task-concordant assessable runs.":'Artifact-level parity was not scored for families without task-concordant runs and explicit artifact paths, so no artifact-parity claim should be made for those runs.';
        push @lines,$line;
    }
    if (@mcp) {
        push @lines,'','### MCP subset','',"The `preferred_interface=mcp` task families are reported separately to characterize Codex orchestration through the local MCP server.";
        for my $s (@mcp) { my $line='**'.ucfirst($s->{agent}).' MCP subset.** Positive-task completion was '.fmt_pct($s->{positive_success},$s->{positive_runs}).'; prespecified invocation concordance was '.fmt_pct($s->{task_match},$s->{runs}).'; '; $line.='safe-failure detection was '.fmt_pct($s->{safe_failure_correct},$s->{safe_failure_runs}).'; ' if $s->{safe_failure_runs}; $line.='median end-to-end latency was '.fmt_num($s->{median_latency}).' s (IQR '.fmt_num($s->{q1_latency}).'-'.fmt_num($s->{q3_latency}).' s).'; push @lines,$line }
    }
    my @failed=grep {!$_->{scored_success}} @$results;
    if (@failed) {
        push @lines,'','## Failed-prompt diagnostics','';
        for my $r (@failed) {
            my $reason=!$r->{process_ok}?'the Codex process did not complete successfully':!$r->{task_match}?'the executed action did not match the preregistered action':($r->{expected_status}//'') eq 'safe_failure' && !$r->{safe_failure_correct}?'the requested failure was detected, but the safe-failure criteria were not all met':($r->{expected_status}//'') eq 'clarify' && !$r->{clarification_detected}?'the contradiction was not handled by requesting clarification':!$r->{terminal_execution_seen}?'no terminal pipeline outcome was observed':'the preregistered success criteria were not all met';
            my $detail='';
            if (ref($r->{actual_actions}) eq 'ARRAY' && @{$r->{actual_actions}}) {
                my $a=$r->{actual_actions}[0];
                my $entry=$a->{entry}//'unknown'; my $spec=ref($a->{args}) eq 'HASH' ? ($a->{args}{spec}//'') : '';
                $detail=" The first observed action was `$entry`".($spec ne ''?" with spec `$spec`.":'.');
            }
            push @lines,"- `$r->{prompt_id}`: $reason.$detail";
        }
    }
    push @lines,'','Repeated-run consistency was not assessed because this focused benchmark used one repetition per prompt.' if $rep_n < 2;
    push @lines,'Safe-failure detection was not assessed because the selected prompts contained no deliberately invalid requests.' unless $safe_n;
    push @lines,'',"These results evaluate the reliability of the AI orchestration layer in invoking the versioned MultiGWAS-Explorer workflow; they do not imply that the language model performs or improves the underlying GWAS statistical calculations.",'','## Suggested manuscript Methods sentence','',"The AI interface was evaluated with ".scalar(keys%p)." preregistered natural-language prompts spanning ".scalar(keys%f)." workflow task families. Codex was run in fresh headless sessions against the same frozen repository and local MCP endpoint. The direct command-line workflow served as the reference. We quantified prespecified action concordance, task completion, repeated-run consistency, safe-failure detection, source-edit violations, and end-to-end execution time; where explicit output paths were available, canonicalized scientific artifacts were also compared by SHA-256 hash.";
    my $methods_tail=$rep_n<2?'Repeated-run consistency was not assessed because each prompt was executed once.':'Repeated-run consistency was assessed across replicate executions.';
    $methods_tail.=$safe_n?' Deliberately invalid prompts were scored for safe-failure behavior.':' Safe-failure behavior was not assessed because no deliberately invalid prompt was selected.';
    $methods_tail.=$clarify_n?' Contradictory prompts were required to request clarification without execution.':'';
    $lines[-1]="The AI interface was evaluated with $prompt_n preregistered natural-language prompts spanning $family_n $family_label. Codex was run in fresh headless sessions against the same $revision_note and local MCP endpoint. The direct command-line workflow served as the reference. The experiment quantified prespecified action concordance, task completion, source-edit violations, and end-to-end execution time; task-concordant outputs were canonicalized and compared with the CLI reference by SHA-256 hash. $methods_tail";
    spew(File::Spec->catfile($out,'reviewer_response_summary.md'),join("\n",@lines)."\n");
    my $methods="AI-interface benchmark. We evaluated Codex with $prompt_n preregistered natural-language prompts spanning $family_n $family_label, using fresh non-interactive sessions against $revision_note. Each prompt was mapped in advance to an evaluator-only gold action and expected outcome, which were not provided to the agent. The benchmark recorded MCP calls and shell commands from machine-readable client event streams, terminal execution status, elapsed wall time, and any agent-initiated source edits. Positive tasks were scored as successful only when the requested action matched the prespecified reference and execution reached a terminal successful state. For task-concordant runs, generated artifacts were canonicalized to remove run-specific metadata and compared with direct CLI reference hashes. $methods_tail\n";
    spew(File::Spec->catfile($out,'benchmark_methods_paragraph.md'),$methods);
}

sub preflight {
    my ($cfg,$repo,$out,$prompt_bank,$gold_path)=@_; my $agents=$cfg->{agents}||{};
    my $git_repo=$cfg->{git_repo}//$repo;
    my $git_status=git_output($git_repo,'status','--short');
    my $git_diff=git_diff_fingerprint($git_repo);
    my $pf={timestamp=>utcnow(),perl=>"$^V",platform=>$^O,repo=>$repo,git_repo=>$git_repo,git_commit=>git_output($git_repo,'rev-parse','HEAD'),git_status=>$git_status,git_diff_sha256=>$git_diff->{sha256},git_diff_tracked_files=>$git_diff->{tracked_changed_files},prompt_bank=>$prompt_bank,prompt_bank_sha256=>sha256_file($prompt_bank),gold_actions=>$gold_path,gold_actions_sha256=>sha256_file($gold_path),mcp=>$cfg->{mcp}||{}};
    $pf->{codex_version}=command_version($agents->{codex}{command}//'codex',$repo);
    $pf->{codex_sandbox}=$agents->{codex}{sandbox}//'workspace-write';
    $pf->{codex_approval_mode}=$agents->{codex}{approval_mode}//'sandbox';
    $pf->{codex_reasoning_effort}=$agents->{codex}{reasoning_effort}//'';
    write_json(File::Spec->catfile($out,'preflight.json'),$pf); return $pf;
}

sub usage {
    return <<'USAGE';
Usage:
  perl run_codex_gemini_benchmark.pl --config benchmark_config.json [options]

Options:
  --phase preflight|run|summarize|all   default: all
  --selection all|sample|fold          override config selection mode
  --fold 1..4
  --per-family N
  --seed N
  --repeats N
  --families F02,F03,F10
  --prompt-ids F02_P01,F10_P03
  --no-cli-reference
  --rerun
  --help
USAGE
}

sub main {
    my %opt=(phase=>'all');
    GetOptions(
        'config=s'=>\$opt{config},'phase=s'=>\$opt{phase},'selection=s'=>\$opt{selection},'fold=i'=>\$opt{fold},
        'per-family=i'=>\$opt{per_family},'seed=i'=>\$opt{seed},'repeats=i'=>\$opt{repeats},
        'families=s'=>\$opt{families},'prompt-ids=s'=>\$opt{prompt_ids},'no-cli-reference'=>\$opt{no_cli_reference},
        'rerun'=>\$opt{rerun},'help'=>\$opt{help}
    ) or die usage();
    die usage() if $opt{help} || !$opt{config};
    die "Invalid --phase\n" unless $opt{phase}=~/^(?:preflight|run|summarize|all)$/;
    my $cfg_path=abs_path($opt{config})||File::Spec->rel2abs($opt{config}); my $cfg=load_json($cfg_path);
    my $repo=abs_path($cfg->{repo})||File::Spec->rel2abs($cfg->{repo}); die "Repository does not exist: $repo\n" unless -d $repo;
    my $cfgdir=dirname($cfg_path);
    my $relpath=sub{my($key,$def)=@_;my$p=$cfg->{$key}//$def;return File::Spec->file_name_is_absolute($p)?$p:File::Spec->rel2abs($p,$cfgdir)};
    my $prompt_bank=$relpath->('prompt_bank','multigwas_100_prompts.jsonl'); my $gold_path=$relpath->('gold_actions','gold_actions.json');
    my $out=$cfg->{output_dir}//'benchmark/codex_eval'; $out=File::Spec->rel2abs($out,$repo) unless File::Spec->file_name_is_absolute($out); make_path($out);
    my $rows=load_jsonl($prompt_bank); my $gold=load_json($gold_path); my %ph=map {$_=>"$cfg->{placeholders}{$_}"} keys %{$cfg->{placeholders}||{}};
    $gold=resolve_obj($gold,\%ph); $rows=resolve_obj($rows,\%ph);

    my @agents=('codex');
    die "Codex agent is disabled in the configuration\n"
        if exists($cfg->{agents}{codex}{enabled}) && !$cfg->{agents}{codex}{enabled};

    # Preflight intentionally does not require UKB placeholders to be filled.
    my $pf=preflight($cfg,$repo,$out,$prompt_bank,$gold_path);
    if ($opt{phase} eq 'preflight') { print File::Spec->catfile($out,'preflight.json')."\n"; return 0 }

    my $sel=$cfg->{selection}||{}; my $mode=$opt{selection}//$sel->{mode}//'all'; my $seed=defined($opt{seed})?$opt{seed}:($sel->{seed}//20260820); my $fold=defined($opt{fold})?$opt{fold}:$sel->{fold}; my $per=defined($opt{per_family})?$opt{per_family}:($sel->{per_family}//1);
    my $families=$opt{families}?[grep{length}split/,/,$opt{families}]:undef; my $pids=$opt{prompt_ids}?[grep{length}split/,/,$opt{prompt_ids}]:undef;
    my $selected=select_prompts(rows=>$rows,mode=>$mode,seed=>$seed,fold=>$fold,per_family=>$per,families=>$families,prompt_ids=>$pids);
    die "No prompts selected\n" unless @$selected;
    my %sf=map{$_->{id}=>1}@$selected; my %sg=map{$_=>$gold->{$_}||{}}keys%sf;
    my $missing=unresolved_placeholders({rows=>$selected,gold=>\%sg});
    die "Unresolved placeholders in selected tasks: ".join(', ',@$missing).". Edit benchmark_config.json or exclude those families.\n" if @$missing;
    my $repeats=defined($opt{repeats})?$opt{repeats}:($cfg->{repeats}//3); my $timeout=$cfg->{timeout_seconds}//7200;

    if ($opt{phase} eq 'run' || $opt{phase} eq 'all') {
        check_codex_access(
            repo=>$repo,out=>$out,agent_cfg=>($cfg->{agents}{codex}||{}),
            timeout=>($cfg->{codex_access_timeout_seconds}//120)
        );
    }

    if ($opt{rerun}) { for my $ag(@agents){for my $r(@$selected){for my $rep(1..$repeats){my$p=File::Spec->catfile($out,'runs',$ag,$r->{prompt_id},sprintf('rep_%02d',$rep),'meta.json');unlink$p if -e$p}}} }
    maybe_configure_mcp($cfg,$repo,$out);
    my $mcp_pid;
    eval {
        $mcp_pid=start_mcp_if_needed($cfg,$repo,$out);
        my $artifact_map=resolve_obj($cfg->{artifact_paths_by_family}||{},\%ph);
        my $cli_refs=load_existing_cli_refs($out,$selected);
        my $run_cli=(!exists($cfg->{run_cli_reference})||$cfg->{run_cli_reference})&&!$opt{no_cli_reference};
        if (($opt{phase} eq 'run'||$opt{phase} eq 'all')&&$run_cli) { $cli_refs=run_cli_references(selected=>$selected,gold=>$gold,repo=>$repo,out=>$out,timeout=>$timeout,artifact_map=>$artifact_map) }
        if ($opt{phase} eq 'run'||$opt{phase} eq 'all') {
            my @results; my $total=@$selected*@agents*$repeats; my $count=0;
            for my $pr(@$selected){my$fid=$pr->{id};my$gs=$gold->{$fid}||{preferred_interface=>'native',actions=>[]};my$patterns=$artifact_map->{$fid}||[];for my $ag(@agents){my$acfg=$cfg->{agents}{$ag}||{};my$cmd=$acfg->{command}//$ag;die "Required command not found for $ag: $cmd\n" unless which_cmd($cmd);for my$rep(1..$repeats){$count++;print "[$count/$total] $ag $pr->{prompt_id} rep $rep\n";my$res=run_one_agent(agent=>$ag,agent_cfg=>$acfg,prompt_row=>$pr,gold_spec=>$gs,repo=>$repo,out=>$out,rep=>$rep,timeout=>$timeout,mcp_name=>$cfg->{mcp}{name}//'perl-bio',placeholders=>\%ph,artifact_patterns=>$patterns,cli_ref=>$cli_refs->{$fid});push@results,$res;write_results($out,\@results);if($res->{agent_unavailable}){die "WARNING: Codex became unavailable during the benchmark; remaining tests have been stopped.\nReason: $res->{agent_unavailable_message}\n"}}}}
            my $all=collect_existing_results($out);
            my %selected_prompt=map { $_->{prompt_id}=>1 } @$selected;
            my %selected_agent=map { $_=>1 } @agents;
            my @current=grep { $selected_prompt{$_->{prompt_id}} && $selected_agent{$_->{agent}} && ($_->{replicate}//0)>=1 && ($_->{replicate}//0)<=$repeats } @$all;
            write_results($out,\@current); summarize($out,\@current,$pf);
        } else {
            my $all=collect_existing_results($out);
            my %selected_prompt=map { $_->{prompt_id}=>1 } @$selected;
            my %selected_agent=map { $_=>1 } @agents;
            my @current=grep { $selected_prompt{$_->{prompt_id}} && $selected_agent{$_->{agent}} && ($_->{replicate}//0)>=1 && ($_->{replicate}//0)<=$repeats } @$all;
            die "No completed run metadata found for the selected prompts\n" unless @current;
            write_results($out,\@current); summarize($out,\@current,$pf);
        }
        1;
    } or do { my $e=$@||'Unknown error'; if($mcp_pid){kill'TERM',$mcp_pid;waitpid($mcp_pid,0)} die $e };
    if($mcp_pid){kill'TERM',$mcp_pid;my$t=time()+10;while(time()<$t){my$w=waitpid($mcp_pid,WNOHANG);last if$w==$mcp_pid;sleep.2}kill'KILL',$mcp_pid;waitpid($mcp_pid,0)}
    print "Benchmark outputs: $out\nReviewer summary: ".File::Spec->catfile($out,'reviewer_response_summary.md')."\nSupplementary table: ".File::Spec->catfile($out,'supplementary_table_agent_performance.tsv')."\n";
    return 0;
}

exit main();
