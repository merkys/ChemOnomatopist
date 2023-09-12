package ChemOnomatopist::Group::Polyaphene;

use strict;
use warnings;

# ABSTRACT: Polyaphene
# VERSION

use parent ChemOnomatopist::Group::, ChemOnomatopist::Chain::Circular::;

use ChemOnomatopist::Util::Graph qw( merge_graphs );
use Graph::Undirected;

sub ideal_graph($$)
{
    my( $class, $N ) = @_;
    my @sizes;
    if( ($N - 6) / 4 % 2 ) {
        # Two unequal branches
        @sizes = map { ( ( ($N - 6) / 4 - 1 ) / 2 + $_ ) * 4 + 2 } (0, 1);
    } else {
        # Two equal-sized branches
        @sizes = ( ($N - 2)/2, ($N - 2)/2 );
    }
    
    my @graphs = map { ChemOnomatopist::Group::Polyacene->ideal_graph( $_ ) } @sizes;
    my $graph = merge_graphs( @graphs );

    # Locate the terminal edges
    my @termini;
    for my $g (@graphs) {
        $g->delete_vertices( map  { $g->neighbours( $_ ) }
                             grep { $g->degree( $_ ) == 3 } $g->vertices );
        $g->delete_vertices( grep { $g->degree( $_ ) == 3 } $g->vertices );
        my( $terminus ) = $g->edges;
        push @termini, $terminus;
    }

    # Create the common ring (BBv2 P-25.1.2.2)
    $graph->add_edge( map { $termini[$_]->[0] } (0, 1) );
    $graph->add_path( $termini[0]->[1],
                      { symbol => 'C' },
                      { symbol => 'C' },
                      $termini[1]->[1] );

    return $graph;
}

1;
