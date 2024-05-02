package ChemOnomatopist::Group::Hydrazine;

# ABSTRACT: Hydrazine group
# VERSION

use strict;
use warnings;

use parent ChemOnomatopist::Group::, ChemOnomatopist::Chain::;

use Chemistry::OpenSMILES qw( is_double_bond );
use List::Util qw( any first );

sub new
{
    my( $class, $graph, @vertices ) = @_;
    return bless { graph => $graph, vertices => \@vertices }, $class;
}

sub candidates()
{
    my( $self ) = @_;

    my @chains = ( $self, ChemOnomatopist::Group::Hydrazine->new( $self->graph, reverse $self->vertices ) );
    $chains[1]->{candidate_for} = $self;

    return @chains;
}

sub needs_heteroatom_locants() { '' }
sub needs_heteroatom_names() { '' }
sub needs_substituent_locants()
{
    my( $self ) = @_;

    my $graph = $self->graph;
    my $parent = $self->{parent};
    if( $parent ) {
        my $vertex = first { $graph->has_edge( $parent, $_ ) } $self->vertices;
        return '' if is_double_bond( $self->graph, $vertex, $parent );
    }

    return $self->number_of_branches > 1 && $self->number_of_branches < $self->max_valence;
}

sub prefix() { 'hydrazinyl' }
sub suffix() { 'hydrazine' }

1;
