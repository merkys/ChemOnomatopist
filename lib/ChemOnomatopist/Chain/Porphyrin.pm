package ChemOnomatopist::Chain::Porphyrin;

# ABSTRACT: Porphyrin compound
# VERSION

use strict;
use warnings;

use parent ChemOnomatopist::Chain::Circular::;

use ChemOnomatopist::Name;
use ChemOnomatopist::Util;
use ChemOnomatopist::Util::Graph qw( graph_without_edge_attributes );
use Graph::Nauty qw( are_isomorphic );
use Graph::Traversal::DFS;
use Graph::Undirected;
use List::Util qw( all first );

sub new
{
    my( $class, $graph, @vertices ) = @_;

    my $subgraph = $graph->subgraph( \@vertices );

    my $first = first { ChemOnomatopist::is_element( $_, 'C' ) &&
                      $subgraph->degree( $_ ) == 3 } @vertices;
    my $last = first { ChemOnomatopist::is_element( $_, 'C' ) &&
                       all { $subgraph->degree( $_ ) == 3 }
                           $subgraph->neighbours( $_ ) }
                     $subgraph->neighbours( $first );

    my @N = grep { ChemOnomatopist::is_element( $_, 'N' ) } @vertices;
    $subgraph->delete_edge( $first, $last );
    $subgraph->delete_vertices( @N );
    my @order = reverse Graph::Traversal::DFS->new( $subgraph, start => $first )->dfs;
    @vertices = ( @order, @N );

    # FIXME: Order N atoms
    # FIXME: Need to do better to achieve minimal locants

    return bless { graph => $graph, vertices => \@vertices }, $class;
}

sub needs_heteroatom_locants() { '' }
sub needs_heteroatom_names() { '' }

sub locants(@) { map { $_ + 1 } @_[1..$#_] }

sub has_form($$)
{
    my( $class, $graph ) = @_;
    my @vertices = $graph->vertices;

    return '' unless @vertices == 24;
    return '' unless (grep { ChemOnomatopist::is_element( $_, 'C' ) } @vertices) == 20;
    return '' unless (grep { ChemOnomatopist::is_element( $_, 'N' ) } @vertices) ==  4;

    return are_isomorphic( graph_without_edge_attributes( $graph ),
                           $class->ideal_graph,
                           sub { ChemOnomatopist::Util::element( $_[0] ) } );
}

sub ideal_graph($)
{
    my( $class ) = @_;
    my $graph = Graph::Undirected->new( refvertexed => 1 );
    my @vertices = map { { symbol => 'C' } } 1..20;
    $graph->add_cycle( @vertices );
    for (0..3) {
        $graph->add_path( $vertices[$_ * 5], { symbol => 'N' }, $vertices[$_ * 5 + 3] );
    }
    return $graph;
}

sub number_of_rings() { 5 }

sub prefix() { ChemOnomatopist::Name->new( 'porphyrin' ) }
sub suffix() { ChemOnomatopist::Name->new( 'porphyrin' ) }

1;
