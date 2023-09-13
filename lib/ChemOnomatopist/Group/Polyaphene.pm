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
use Graph::Undirected;
use Set::Object qw( set );

sub new
{
    my( $class, $graph, @cycles ) = @_;

    my $subgraph = subgraph( $graph, map { $_->vertices } @cycles );

    # Constructing the connectivity graph for cycles
    my $connectivity_graph = Graph::Undirected->new( refvertexed => 1 );
    for my $vertex ($subgraph->vertices) {
        next unless $subgraph->degree( $vertex ) == 3;
        $connectivity_graph->add_edge( grep { set( $_->vertices )->has( $vertex ) } @cycles );
    }

    # Detecting the common ring which has three edges having degree 3 at their ends
    my $common_ring;
    for my $cycle (@cycles) {
        next unless scalar( grep { $subgraph->degree( $_->[0] ) == 3 &&
                                   $subgraph->degree( $_->[1] ) == 3 }
                                 subgraph( $graph, $cycle->vertices )->edges ) == 3;
        $common_ring = $cycle;
    }

    # Finding the correct order of cycles, flipping if needed
    my $common_ring_pos = @cycles % 2 ? (@cycles - 1) / 2 : @cycles / 2 - 1;
    my( $start ) = grep { $connectivity_graph->degree( $_ ) == 1 } @cycles;
    my @cycles_in_order = Graph::Traversal::DFS->new( $connectivity_graph,
                                                      start => $start )->dfs;
    if( !@cycles % 2 && $cycles_in_order[$common_ring_pos] != $common_ring ) {
        @cycles_in_order = reverse @cycles_in_order;
    }

    # Finding the atom in the common ring which will get the lowest number
    $subgraph = subgraph( $graph, $common_ring->vertices );
    $subgraph->delete_edges( map { @$_ }
                             map { subgraph( $graph, $_->vertices )->edges }
                                 ( $cycles_in_order[$common_ring_pos-1],
                                   $cycles_in_order[$common_ring_pos+1] ) );
    my( $short_edge ) = grep { $subgraph->degree( $_->[0] ) == 1 &&
                               $subgraph->degree( $_->[1] ) == 1 }
                               $subgraph->edges;
    my( $junction ) = (set( $cycles_in_order[$common_ring_pos-1]->vertices ) *
                       set( @$short_edge ))->members;

    # Finding the candidates of the starting atom
    $subgraph = subgraph( $graph, $cycles_in_order[0]->vertices );
    $subgraph->delete_vertices(   $cycles_in_order[1]->vertices );
    my @candidates = grep { $subgraph->degree( $_ ) == 1 } $subgraph->vertices;

    # Finding the first and the last atom in the enumeration order
    $subgraph = subgraph( $graph, map { $_->vertices } @cycles );
    my $shortest_paths = $subgraph->single_source_shortest_paths( $junction );
    my $min_length;
    my $first;
    for my $vertex (@candidates) {
        my $length = 0;
        my $v = $vertex;
        while(   $shortest_paths->has_vertex_attribute( $v, 'p' ) ) {
            $v = $shortest_paths->get_vertex_attribute( $v, 'p' );
            $length++;
        }
        if( !defined $min_length || $min_length > $length ) {
            $min_length = $length;
            $first = $vertex;
        }
    }

    my( $last ) = grep { $subgraph->degree( $_ ) == 3 }
                       $subgraph->neighbours( $first );

    # Deleting chords and connection between first and last atoms
    $subgraph->delete_edges( map { (set( $cycles_in_order[$_] ) * set( $cycles_in_order[$_+1] ))->members }
                                 0..$#cycles-1 );
    $subgraph->delete_edge( $first, $last );

    my @vertices = Graph::Traversal::DFS->new( $subgraph, start => $last )->dfs;
    return bless { graph => $graph, vertices => \@vertices }, $class;
}

sub candidates()
{
    my( $self ) = @_;
    return ( $self ); # FIXME: For now
}

sub locants(@) # FIXME: Complete
{
    my $self = shift;
    my $N = ($self->length - 6) / 4;
    my @rings = $N % 2 ? ( ($N-1)/2, ($N-1)/2+1 ) : ( $N/2, $N/2 );
    my @locants = ( 1..3,
                    (map { 3 + $_, (3 + $_) . 'a' } 1..$rings[0]),
                    4+$rings[0],
                    (map { 3+$rings[0]+1 + $_, (3+$rings[0]+1 + $_) . 'a' } 1..$rings[1]),
                    (4 + $rings[0] + $rings[1]) . 'b',
                    (map { 4 + $rings[0] + $rings[1] + $_, 4 + $rings[0] + $rings[1] + $_ . 'a' } 1..$rings[0]) );
    # print join ',', @locants;
    return @_;
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

sub prefix()
{
    my( $self ) = @_;
    return ChemOnomatopist::IUPAC_numerical_multiplier( ($self->length - 2) / 4 ) . 'aphene';
}

sub suffix() { return $_[0]->prefix }

1;