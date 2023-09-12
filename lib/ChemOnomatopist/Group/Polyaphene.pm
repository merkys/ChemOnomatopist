package ChemOnomatopist::Group::Polyaphene;

use strict;
use warnings;

# ABSTRACT: Polyaphene
# VERSION

use parent ChemOnomatopist::Group::, ChemOnomatopist::Chain::Circular::;

use ChemOnomatopist::Util::Graph qw(
    merge_graphs
    subgraph
);

sub new
{
    my( $class, $graph, @cycles ) = @_;

    my $subgraph = subgraph( $graph, map { $_->vertices } @cycles );
    my $common_ring; # The common ring has 3 edges with degree 3 at their ends
    # Terminal rings have 4 vertices of degree 2
    my @termini;
    for my $cycle (@cycles) {
        my $nchords = scalar grep { $subgraph->degree( $_->[0] ) == 3 &&
                                    $subgraph->degree( $_->[1] ) == 3 }
                                  subgraph( $subgraph, $cycle->vertices )->edges;
        $common_ring = $cycle if $nchords == 3;

        if( scalar( grep { $subgraph->degree( $_ ) == 2 } $cycle->vertices ) == 4 ) {
            push @termini, $cycle;
        }
    }

    # Finding out the smaller of the branches
    $subgraph->delete_cycle( $common_ring->vertices );
    my( $smaller ) = sort { @$a <=> @$b } $subgraph->connected_components;

    $subgraph = subgraph( $graph, map { $_->vertices } @cycles ); # Restore
    $subgraph->delete_vertices( @$smaller );
    # Junction atom is the one from the common ring which has the shortest distance to the atom number 1.
    my( $junction ) = map  { $subgraph->degree( $_->[0] ) == 1 ? $_->[0] : $_->[1] }
                      grep { ($subgraph->degree( $_->[0] ) == 1 &&
                              $subgraph->degree( $_->[1] ) == 3) ||
                             ($subgraph->degree( $_->[0] ) == 3 &&
                              $subgraph->degree( $_->[1] ) == 1) }
                           $subgraph->edges;

    $subgraph = subgraph( $graph, map { $_->vertices } @smaller );
    my @distances = map { scalar $subgraph->SP_Dijkstra( $junction, $_ ) }
}

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

sub suffix()
{
    my( $self ) = @_;
    return ChemOnomatopist::IUPAC_numerical_multiplier( ($self->length - 2) / 4 ) . 'aphene';
}

1;
