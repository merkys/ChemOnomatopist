package ChemOnomatopist::Group::Monospiro;

use strict;
use warnings;

# ABSTRACT: Monospiro compound
# VERSION

use parent ChemOnomatopist::Group::;

use ChemOnomatopist::Chain::Circular;

sub new
{
    my( $class, $graph, @vertices ) = @_;

    my $subgraph = $graph->subgraph( \@vertices );
    my( $spiro_atom ) = grep { $subgraph->degree( $_ ) == 4 } @vertices;
    $subgraph->delete_vertex( $spiro_atom );
    my @components = sort { scalar @$a <=> scalar @$b } $subgraph->connected_components;

    return bless { graph => $graph, vertices => \@vertices, spiro_atom => $spiro_atom, components => \@components }, $class;
}

1;
