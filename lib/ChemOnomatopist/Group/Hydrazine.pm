package ChemOnomatopist::Group::Hydrazine;

use strict;
use warnings;

# ABSTRACT: Hydrazine group
# VERSION

use parent ChemOnomatopist::Group::, ChemOnomatopist::Chain::;

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

sub is_part_of_chain() { return 1 }

sub prefix(;$)
{
    my( $self, $parent ) = @_;
    return 'hydrazine';
}

sub suffix()
{
    my( $self ) = @_;
    return 'hydrazine';
}

1;
