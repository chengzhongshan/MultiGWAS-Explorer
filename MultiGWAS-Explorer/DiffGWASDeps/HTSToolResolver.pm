package HTSToolResolver;

use strict;
use warnings;
use Exporter 'import';
use Cwd qw(abs_path);
use File::Basename qw(dirname);
use File::Spec;

our @EXPORT_OK = qw(resolve_hts_tool external_tool_path hts_tool_is_native);

sub resolve_hts_tool {
    my ($name, %args) = @_;
    die "HTS tool name is required\n" unless defined($name) && length($name);
    my $explicit = $args{explicit} // '';
    my $start_dir = abs_path($args{start_dir} || '.') || ($args{start_dir} || '.');
    my @names = $^O =~ /^(?:cygwin|MSWin32)$/i ? ($name, "$name.exe") : ($name);
    my @candidates;

    if (length $explicit) {
        if (-d $explicit) {
            push @candidates, map { File::Spec->catfile($explicit, $_) } @names;
        }
        else {
            push @candidates, $explicit;
        }
    }

    my $env_name = uc($name) . '_BIN';
    if (!$explicit && defined($ENV{$env_name}) && length($ENV{$env_name})) {
        my $value = $ENV{$env_name};
        push @candidates, -d $value
          ? map { File::Spec->catfile($value, $_) } @names
          : $value;
    }

    my $ancestor = $start_dir;
    for (0 .. 6) {
        my $bin = File::Spec->catdir($ancestor, 'local', 'bin');
        push @candidates, map { File::Spec->catfile($bin, $_) } @names;
        my $parent = abs_path(File::Spec->catdir($ancestor, File::Spec->updir()));
        last unless defined($parent) && length($parent) && $parent ne $ancestor;
        $ancestor = $parent;
    }

    if ($^O =~ /cygwin/i) {
        push @candidates, map { File::Spec->catfile('/usr/bin', $_) } @names;
        push @candidates, map { File::Spec->catfile('/usr/local/bin', $_) } @names;
    }
    for my $dir (File::Spec->path()) {
        push @candidates, map { File::Spec->catfile($dir, $_) } @names;
    }

    my %seen;
    for my $candidate (@candidates) {
        next unless defined($candidate) && length($candidate);
        my $key = lc($candidate);
        next if $seen{$key}++;
        next unless -f $candidate && -x $candidate;
        next unless hts_tool_is_native($candidate);
        return abs_path($candidate) || $candidate;
    }

    return;
}

sub hts_tool_is_native {
    my ($tool) = @_;
    return 0 unless defined($tool) && -f $tool && -x $tool;
    return 1 unless $^O =~ /cygwin/i;
    return 1 if _truthy($ENV{PIPELINE_ALLOW_WINDOWS_HTSLIB});

    my $pid = open my $fh, '-|';
    return 0 unless defined $pid;
    if ($pid == 0) {
        open STDERR, '>', File::Spec->devnull();
        exec 'cygcheck', $tool;
        exit 127;
    }
    my $native = 0;
    while (my $line = <$fh>) {
        $native = 1 if $line =~ /cygwin1\.dll/i;
    }
    close $fh;
    return $native;
}

sub external_tool_path {
    my ($tool, $path) = @_;
    return $path unless defined($path) && length($path) && $^O =~ /cygwin/i;
    return $path if hts_tool_is_native($tool);
    my $converted = '';
    if (open my $fh, '-|', 'cygpath', '-m', $path) {
        $converted = <$fh> // '';
        chomp $converted;
        close $fh;
    }
    return length($converted) ? $converted : $path;
}

sub _truthy {
    my ($value) = @_;
    return defined($value) && $value =~ /^(?:1|true|yes|y|on)$/i;
}

1;
