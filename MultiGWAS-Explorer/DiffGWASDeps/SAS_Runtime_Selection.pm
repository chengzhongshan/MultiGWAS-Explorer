package SAS_Runtime_Selection;
use strict;
use warnings;
use File::Spec;

# Return only public routing settings. Credentials remain in SASPy configuration.
sub resolve {
    my (%args) = @_;
    my $cli = $args{cli} || {};
    my $spec = $args{spec} || {};
    my $env = $args{env} || {};
    my %selected;
    for my $field ([cfgname => 'saspy_cfgname', 'SASPY_CFGNAME'],
                   [cfgfile => 'saspy_cfgfile', 'SASPY_CFGFILE'],
                   [home => 'dir4local_sas_home', 'SAS_LOCAL_HOME_DIR']) {
        my ($key, $spec_key, $env_key) = @$field;
        $selected{$key} = defined($cli->{$key}) ? $cli->{$key}
            : defined($spec->{$spec_key}) ? $spec->{$spec_key}
            : ($env->{$env_key} // ($key eq 'cfgname' ? $env->{SASPY_CONFIG_NAME} : '') // '');
    }
    $selected{cfgname} ||= 'oda';
    my $local = $selected{cfgname} =~ /^(?:winlocal|winlocal_com|linuxlocal|local|stdio)$/i;
    if (lc($selected{cfgname}) eq 'oda') {
        die "A local SAS home cannot be used with --saspy-cfgname oda\n"
            if length($cli->{home} // $spec->{dir4local_sas_home} // '');
        # An explicit return to ODA must not carry a previous local staging root.
        $selected{home} = '';
    }
    elsif ($local && !length($selected{home})) {
        my $tag = $args{artifact_stem} || 'diff_gwas';
        $tag =~ s/[^A-Za-z0-9_.-]/_/g;
        $selected{home} = File::Spec->rel2abs(
            File::Spec->catdir($args{output_dir} || '.', "SAS_Local_Home_$tag"));
        $selected{home} =~ s{\\}{/}g;
    }
    return {
        SASPY_CFGNAME => $selected{cfgname},
        SASPY_CFGFILE => $selected{cfgfile},
        SAS_LOCAL_HOME_DIR => $selected{home},
    };
}
1;
