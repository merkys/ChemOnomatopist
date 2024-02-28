package ChemOnomatopist::Chain::VonBaeyer;

# ABSTRACT: Von Baeyer hydrocarbon
# VERSION

use strict;
use warnings;

use Graph::Traversal::DFS;
use List::Util qw( first );

use parent ChemOnomatopist::Chain::Circular::;

sub new
{
    my( $class, $graph, @vertices ) = @_;

    my $subgraph = $graph->subgraph( \@vertices );
    my @d3 = grep { $subgraph->degree( $_ ) == 3 } @vertices;

    $subgraph->delete_vertices( @d3 );
    my @components = sort { @$b <=> @$a } $graph->connected_components;

    $subgraph = $graph->subgraph( \@vertices );
    my $first_of_bridge = first { $subgraph->has_edge( $d3[0], $_ ) }
                                 @{$components[-1]};
    my $last_of_bridge  = first { $subgraph->has_edge( $d3[1], $_ ) }
                                 @{$components[-1]};
    # Disconnect the main bridge
    $subgraph->delete_edge( $d3[0], $first_of_bridge );
    $subgraph->delete_edge( $d3[1], $last_of_bridge );

    # Find the last atom of the main ring
    my $last_of_the_main = first { $subgraph->has_edge( $d3[0], $_ ) }
                                 @{$components[1]};
    $subgraph->delete_edge( $d3[0], $last_of_the_main ); # Disconnect
    # Connect in order to get the correct numbering
    $subgraph->add_edge( $last_of_the_main, $first_of_bridge );

    @vertices = Graph::Traversal::DFS->new( $subgraph, start => $d3[0] )->dfs;

    return bless { graph => $graph,
                   vertices => \@vertices,
                   sizes => [ map { scalar @$_ } @components ] }, $class;
}

sub has_form($$)
{
    my( $class, $graph ) = @_;

    my @d2 = grep { $graph->degree( $_ ) == 2 } $graph->vertices;
    my @d3 = grep { $graph->degree( $_ ) == 3 } $graph->vertices;

    return '' unless @d3 == 2;
    return '' unless @d2 + @d3 == scalar $graph->vertices;

    $graph = $graph->copy->delete_vertices( @d3 );
    return scalar( $graph->connected_components ) == 3;
}

sub suffix()
{
    my( $self ) = @_;
    return 'bicyclo' . $self->length;
}

1;
