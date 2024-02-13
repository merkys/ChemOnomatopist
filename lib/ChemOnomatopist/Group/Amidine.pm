package ChemOnomatopist::Group::Amidine;

# ABSTRACT: Amidine group
# VERSION

use strict;
use warnings;

use parent ChemOnomatopist::Group::, ChemOnomatopist::Chain::;

use Scalar::Util qw( blessed );

sub new
{
    my( $class, $graph, @vertices ) = @_;
    return bless { graph => $graph, vertices => \@vertices }, $class;
}

sub nonstandard_valence_positions()
{
    my( $self ) = @_;

    return @{$self->{nonstandard_valence_positions}} if $self->{nonstandard_valence_positions};

    my @vertices = $self->vertices;
    my @nonstandard_valence_positions;
    for (1..$#vertices) { # Nonstandard valence of the central atom is not important, hence skipped
        next if blessed $vertices[$_];
        next if ChemOnomatopist::is_element( $vertices[$_], 'C' );
        next unless exists $vertices[$_]->{valence};
        push @nonstandard_valence_positions, $_;
    }

    $self->{nonstandard_valence_positions} = \@nonstandard_valence_positions;
    return @nonstandard_valence_positions;
}

sub needs_heteroatom_locants() { '' }
sub needs_heteroatom_names() { '' }
sub needs_substituent_locants { '' }

sub prefix() { 'carbamimidoyl' }
sub suffix()
{
    my( $self ) = @_;

    # FIXME: Not fully implemented
    my( $central_atom, @others ) = $self->vertices;
    if( $central_atom->{symbol} eq 'S' ) {
        return 'sulfinimidamide';
    }

    return 'imidamide';
}

1;
