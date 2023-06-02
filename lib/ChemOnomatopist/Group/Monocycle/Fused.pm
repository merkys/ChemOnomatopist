package ChemOnomatopist::Group::Monocycle::Fused;

use strict;
use warnings;

# ABSTRACT: Monocyclic group which is fused to some other cycle
# VERSION

use parent ChemOnomatopist::Group::Monocycle::;

use ChemOnomatopist;
use List::Util qw( all );

sub new
{
    my( $class, $graph, $system, @vertices ) = @_;

    return bless { graph => $graph, system => \$system, vertices => \@vertices }, $class;
}

sub is_aromatic()
{
    my( $self ) = @_;
    return 1 if $self->SUPER::is_aromatic;
    return $self->{system}->is_aromatic;
}

# Returns a copy of the monocycle flipped around the bridge
sub flipped()
{
    my( $self ) = @_;
    my @vertices = $self->vertices;
    my @flipped = ( reverse( @vertices[0..-2] ), reverse( @vertices[-2..-1] ) );
    return ChemOnomatopist::Group::Monocycle::Fused->new( $self->graph,
                                                          $self->{system},
                                                          @flipped );
}

# Find the best orientation between self and flipped self.
# Return 1 if flipped.
sub orient()
{
    my( $self ) = @_;
    my( $chain ) = ChemOnomatopist::filter_chains( $self, $self->flipped );
    return '' if $chain == $self;

    $self->{vertices} = $chain->{vertices};
    return 1;
}

1;
