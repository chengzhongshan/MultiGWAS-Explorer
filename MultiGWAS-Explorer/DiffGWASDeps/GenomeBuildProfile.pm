package GenomeBuildProfile;
use strict;
use warnings;

use Exporter 'import';

our @EXPORT_OK = qw(
  canonicalize_reference_build
  genome_build_profiles
  build_profile_for
  detect_reference_build_profile
  detect_reference_build_token
);

sub genome_build_profiles {
    return {
        hg19 => {
            build      => 'hg19',
            label      => 'hg19 / GRCh37',
            gtf_url    => 'https://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_49/GRCh37_mapping/gencode.v49lift37.annotation.gtf.gz',
            local_dsd  => 'gtf_hg19',
            shared_dsd => 'FM.GTF_HG19',
            aliases    => [qw(hg19 grch37 b37 build37 lift37 37lift37 release_49lift37)],
        },
        hg38 => {
            build      => 'hg38',
            label      => 'hg38 / GRCh38',
            gtf_url    => 'https://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_49/gencode.v49.annotation.gtf.gz',
            local_dsd  => 'gtf_hg38',
            shared_dsd => 'FM.GTF_HG38',
            aliases    => [qw(hg38 grch38 b38 build38)],
        },
        t2t => {
            build      => 't2t',
            label      => 'T2T / hs1 / CHM13v2.0',
            gtf_url    => 'https://hgdownload.soe.ucsc.edu/goldenPath/hs1/bigZips/genes/hs1.ncbiRefSeq.gtf.gz',
            local_dsd  => 'gtf_t2t',
            shared_dsd => 'FM.GTF_T2T',
            aliases    => [qw(t2t hs1 chm13 chm13v2 chm13v2.0 t2t-chm13)],
        },
    };
}

sub canonicalize_reference_build {
    my ($value) = @_;
    return '' unless defined $value;
    my $norm = lc $value;
    $norm =~ s/^\s+|\s+$//g;
    return '' unless length $norm;

    my $profiles = genome_build_profiles();
    for my $build (sort keys %{$profiles}) {
        return $build if $norm eq $build;
        for my $alias (@{ $profiles->{$build}{aliases} || [] }) {
            return $build if $norm eq lc($alias);
        }
    }
    return '';
}

sub build_profile_for {
    my (%args) = @_;
    my $build = canonicalize_reference_build($args{build});
    return undef unless length $build;
    my $profiles = genome_build_profiles();
    my %copy = %{ $profiles->{$build} || {} };
    $copy{build} = $build;
    $copy{source} = $args{source} if defined $args{source};
    $copy{evidence} = $args{evidence} if defined $args{evidence};
    return \%copy;
}

sub detect_reference_build_profile {
    my (%args) = @_;
    my $default_build = canonicalize_reference_build($args{default_build}) || 'hg38';
    if (my $explicit = canonicalize_reference_build($args{explicit_build})) {
        return build_profile_for(
            build    => $explicit,
            source   => 'explicit_override',
            evidence => $args{explicit_build},
        );
    }

    for my $header (@{ $args{header_lines} || [] }) {
        my ($build, $evidence) = detect_reference_build_token($header);
        next unless $build;
        return build_profile_for(
            build    => $build,
            source   => 'header_token',
            evidence => $evidence,
        );
    }

    for my $path (@{ $args{file_paths} || [] }) {
        my ($build, $evidence) = detect_reference_build_token($path);
        next unless $build;
        return build_profile_for(
            build    => $build,
            source   => 'path_token',
            evidence => $evidence,
        );
    }

    return build_profile_for(
        build    => $default_build,
        source   => 'fallback_default',
        evidence => $default_build,
    );
}

sub detect_reference_build_token {
    my ($text) = @_;
    return ('', '') unless defined $text;
    my $norm = lc $text;
    return ('t2t', $1) if $norm =~ /(?:^|[^a-z0-9])(t2t(?:-?chm13)?|hs1|chm13(?:v2(?:\.0)?)?)(?:[^a-z0-9]|$)/;
    return ('hg38', $1) if $norm =~ /(?:^|[^a-z0-9])(hg38|grch38|b38|build38|pos_hg38|bp_hg38)(?:[^a-z0-9]|$)/;
    return ('hg19', $1) if $norm =~ /(?:^|[^a-z0-9])(hg19|grch37|b37|build37|lift37|pos_hg19|bp_hg19)(?:[^a-z0-9]|$)/;
    return ('', '');
}

1;
