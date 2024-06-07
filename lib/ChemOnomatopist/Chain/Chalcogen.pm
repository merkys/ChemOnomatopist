package ChemOnomatopist::Chain::Chalcogen;

# ABSTRACT: Chalcogen chain
# VERSION

use strict;
use warnings;

use ChemOnomatopist::Name;
use Graph::Traversal::DFS;
use List::Util qw( first );

use parent ChemOnomatopist::Chain::;

sub new
{
    my( $class, $graph, $parent, @vertices ) = @_;
    my $subgraph = $graph->subgraph( @vertices );
    @vertices = Graph::Traversal::DFS->new( $subgraph,
                                            start => first { $subgraph->degree( $_ ) == 1 } @vertices )->dfs;
    return bless { graph => $graph, vertices => \@vertices }, $class;
}

sub needs_suffix_locant() { '' }

sub suffix() { ChemOnomatopist::Name->new( 'idane' ) }

1;
