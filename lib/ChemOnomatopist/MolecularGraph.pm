package ChemOnomatopist::MolecularGraph;

use strict;
use warnings;

# ABSTRACT: Graph extension for molecular graphs
# VERSION

use parent Graph::Undirected::;

sub new
{
    my( $class, $what ) = @_;
    return bless Graph::Undirected->new( refvertexed => 1 ), $class unless $what;

    die "molecular graphs can only be made from graphs\n" unless $what->isa( Graph::Undirected:: );

    return bless $what, $class;
}

1;
