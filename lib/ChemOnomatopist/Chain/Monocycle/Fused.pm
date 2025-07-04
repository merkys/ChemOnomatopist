package ChemOnomatopist::Chain::Monocycle::Fused;

use strict;
use warnings;

# ABSTRACT: Monocyclic group which is fused to some other cycle
# VERSION

use parent ChemOnomatopist::Chain::Monocycle::;

use ChemOnomatopist;
use ChemOnomatopist::Util::SMILES qw( cycle_SMILES );

sub new
{
    my( $class, $graph, $system, @vertices ) = @_;
    return bless { graph => $graph, system => $system, vertices => \@vertices }, $class;
}

sub system()
{
    my( $self ) = @_;
    return $self->{system};
}

sub backbone_SMILES()
{
    my( $self ) = @_;
    return cycle_SMILES( $self->graph, $self->vertices );
}

# Returns a copy of the monocycle flipped around the bridge
sub flipped()
{
    my( $self ) = @_;
    my @vertices = $self->vertices;
    my @bridge = splice @vertices, -2;
    my @flipped = ( reverse( @vertices ), reverse( @bridge ) );
    return ChemOnomatopist::Chain::Monocycle::Fused->new( $self->graph,
                                                          $self->system,
                                                          @flipped );
}

# Find the best orientation between self and flipped self.
# Return 1 if flipped.
sub orient()
{
    my( $self ) = @_;
    my( $chain ) = ChemOnomatopist::filter_chains( $self, $self->flipped );
    return '' if $chain == $self;

    $self->vertices( $chain->vertices );
    return 1;
}

1;
