package ChemOnomatopist::Chain::Circular::Numbering;

use strict;
use warnings;

use ChemOnomatopist::Chain::Circular; # FIXME: Not sure why it is needed

use parent ChemOnomatopist::Chain::Circular::;

sub new
{
    my( $class, $chain, $offset, $direction ) = @_;
    return bless { chain => $chain, offset => $offset, direction => $direction }, $class;
}

sub vertices()
{
    my( $self ) = @_;
    my @vertices = $self->vertices;
    @vertices = shift @vertices, reverse @vertices if $self->{direction} == -1;
    for (1..$self->{offset}) {
        push @vertices, shift @vertices;
    }
    return @vertices;
}

# Properties

# TODO: Rewrite

1;
