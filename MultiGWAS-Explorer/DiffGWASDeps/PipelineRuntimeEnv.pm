package PipelineRuntimeEnv;

use strict;
use warnings;
use Config ();
use Cwd ();
use File::Basename ();
use File::Spec ();

sub bootstrap_local_perl {
    my (%args) = @_;
    my $script_dir = Cwd::abs_path($args{script_dir} || '.') || ($args{script_dir} || '.');
    my @required = @{ $args{required_modules} || [] };
    my $platform_tag = lc($^O || '');
    $platform_tag =~ s/[^a-z0-9]+/_/g;

    my @roots;
    if (defined($ENV{PIPELINE_PERL_LOCAL_DIR}) && length($ENV{PIPELINE_PERL_LOCAL_DIR})) {
        push @roots, [ $ENV{PIPELINE_PERL_LOCAL_DIR}, 'PIPELINE_PERL_LOCAL_DIR' ];
    }

    my $ancestor = $script_dir;
    for my $depth (0 .. 6) {
        push @roots,
          [ File::Spec->catdir($ancestor, 'local', "perl5-$platform_tag"), $depth ? "ancestor level $depth" : 'current checkout' ],
          [ File::Spec->catdir($ancestor, 'local', 'perl5'), $depth ? "ancestor level $depth" : 'current checkout' ];
        my $parent = Cwd::abs_path(File::Spec->catdir($ancestor, File::Spec->updir()));
        last unless defined($parent) && length($parent) && $parent ne $ancestor;
        $ancestor = $parent;
    }

    my %seen;
    my @rejected;
    my @rejected_libs;
    for my $entry (@roots) {
        my ($root, $source) = @$entry;
        next unless defined($root) && length($root);
        my $canonical = Cwd::abs_path($root) || $root;
        next if $seen{$canonical}++;
        my $base = File::Spec->catdir($canonical, 'lib', 'perl5');
        next unless -d $base;
        my @arch = grep { _arch_dir_matches($_) } glob(File::Spec->catdir($base, '*'));
        my @libs = ($base, @arch);
        unless (_perl_tree_is_compatible(\@libs, \@required)) {
            push @rejected, $canonical;
            push @rejected_libs, @libs;
            next;
        }

        _remove_rejected_libraries(\@rejected_libs) if @rejected_libs;
        require lib;
        lib->import(@libs);
        _prepend_env_list('PERL5LIB', @libs);
        $ENV{PIPELINE_PERL_LOCAL_DIR} = $canonical;
        my $bin = File::Spec->catdir($canonical, 'bin');
        _prepend_env_list('PATH', $bin) if -d $bin;
        if ($source ne 'current checkout' && $source ne 'PIPELINE_PERL_LOCAL_DIR') {
            print STDERR "[runtime] Using compatible local Perl environment from $canonical ($source).\n";
        }
        return $canonical;
    }

    if (@rejected) {
        _remove_rejected_libraries(\@rejected_libs);
        if (defined($ENV{PIPELINE_PERL_LOCAL_DIR}) && length($ENV{PIPELINE_PERL_LOCAL_DIR})) {
            my $configured = Cwd::abs_path($ENV{PIPELINE_PERL_LOCAL_DIR})
              || $ENV{PIPELINE_PERL_LOCAL_DIR};
            delete $ENV{PIPELINE_PERL_LOCAL_DIR}
              if grep { $_ eq $configured } @rejected;
        }
        print STDERR "[runtime] Skipped incompatible local Perl environment(s): "
          . join(', ', @rejected) . "\n";
    }
    return '';
}

sub _remove_rejected_libraries {
    my ($rejected) = @_;
    my %bad = map {
        my $path = Cwd::abs_path($_) || $_;
        $path =~ s{[\\/]+$}{};
        (lc($path) => 1)
    } @$rejected;
    return unless %bad;

    @INC = grep { !_is_rejected_path($_, \%bad) } @INC;
    my $separator = $^O eq 'MSWin32' ? ';' : ':';
    my @perl5lib = grep { length && !_is_rejected_path($_, \%bad) }
      split /\Q$separator\E/, ($ENV{PERL5LIB} // '');
    if (@perl5lib) {
        $ENV{PERL5LIB} = join($separator, @perl5lib);
    } else {
        delete $ENV{PERL5LIB};
    }
}

sub _is_rejected_path {
    my ($path, $bad) = @_;
    return 0 if ref($path);
    my $canonical = Cwd::abs_path($path) || $path;
    $canonical =~ s{[\\/]+$}{};
    my $lower = lc($canonical);
    return 1 if $bad->{$lower};
    for my $root (keys %$bad) {
        return 1 if index($lower, "$root/") == 0 || index($lower, "$root\\") == 0;
    }
    return 0;
}

sub _arch_dir_matches {
    my ($dir) = @_;
    return 0 unless -d $dir;
    my $name = File::Basename::basename($dir);
    my $arch = lc($Config::Config{archname} || '');
    my $os = lc($^O || '');
    my $lower = lc($name);
    return 1 if $arch && ($lower eq $arch || index($arch, $lower) >= 0 || index($lower, $arch) >= 0);
    return 1 if $os eq 'cygwin' && $lower =~ /cygwin/;
    return 1 if $os =~ /linux/ && $lower =~ /(?:linux|gnu)/;
    return 1 if $os =~ /darwin/ && $lower =~ /darwin/;
    return 1 if $os =~ /mswin32/ && $name =~ /MSWin32/i;
    return 0;
}

sub _perl_tree_is_compatible {
    my ($libs, $required) = @_;
    my @probe_modules = @$required;
    @probe_modules = ('Compress::Raw::Zlib') unless @probe_modules;
    my @cmd = (
        $^X,
        (map { ('-I', $_) } @$libs),
        (map { "-M$_" } @probe_modules),
        '-e', '1',
    );
    open my $saved_stdout, '>&', \*STDOUT or return 0;
    open my $saved_stderr, '>&', \*STDERR or return 0;
    open STDOUT, '>', File::Spec->devnull() or return 0;
    open STDERR, '>', File::Spec->devnull() or return 0;
    my $status;
    {
        local %ENV = %ENV;
        delete @ENV{qw(PERL5LIB PERL_LOCAL_LIB_ROOT PERL_MB_OPT PERL_MM_OPT)};
        $status = system { $^X } @cmd;
    }
    open STDOUT, '>&', $saved_stdout;
    open STDERR, '>&', $saved_stderr;
    return $status == 0 ? 1 : 0;
}

sub _prepend_env_list {
    my ($name, @values) = @_;
    my $separator = $^O eq 'MSWin32' ? ';' : ':';
    my @current = grep { length } split /\Q$separator\E/, ($ENV{$name} // '');
    my %seen = map { $_ => 1 } @current;
    my @prefix = grep { defined($_) && length($_) && !$seen{$_}++ } @values;
    $ENV{$name} = join($separator, @prefix, @current);
}

1;
