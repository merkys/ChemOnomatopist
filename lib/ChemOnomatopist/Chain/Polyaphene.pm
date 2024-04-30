package ChemOnomatopist::Chain::Polyaphene;

use strict;
use warnings;

# ABSTRACT: Polyaphenes, including phenanthrene
# VERSION

use parent ChemOnomatopist::Chain::Circular::;

use ChemOnomatopist::Util::Graph qw(
    graph_without_edge_attributes
    merge_graphs
    subgraph
);
use Graph::Nauty qw( are_isomorphic );
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
    $subgraph->delete_edges( map { (set( $cycles_in_order[$_  ]->vertices ) *
                                    set( $cycles_in_order[$_+1]->vertices ))->members }
                                 0..$#cycles-1 );
    $subgraph->delete_edge( $first, $last );

    my @vertices = Graph::Traversal::DFS->new( $subgraph, start => $last )->dfs;
    return bless { graph => $graph, vertices => \@vertices }, $class;
}

sub candidates()
{
    my( $self ) = @_;
    my @candidates = ( $self );

    if( (($self->length - 6) / 4) % 2 == 0 ) {
        my @vertices = reverse $self->vertices;
        my $N = ($self->length - 6) / 4 / 2; # Number of rings in one "line"
        for (1..($N-1)*4+2) {
            push @vertices, shift @vertices;
        }
        push @candidates,
             bless { graph => $self->graph,
                     vertices => \@vertices,
                     candidate_for => $self };
    }

    return @candidates;
}

sub needs_substituent_locants() { 1 }

sub has_form($$)
{
    my( $class, $graph ) = @_;
    my %degrees = map { $graph->degree( $_ ) => 1 } $graph->vertices;
    return '' unless join( ',', sort keys %degrees ) eq '2,3';

    my $N = scalar $graph->vertices;
    return '' if  $N < 6 + 3 * 4;
    return '' if ($N - 6) % 4;

    return are_isomorphic( graph_without_edge_attributes( $graph ),
                           $class->ideal_graph( $N ),
                           sub { 'C' } );
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

    my @graphs = map { ChemOnomatopist::Chain::Polyacene->ideal_graph( $_ ) } @sizes;
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

sub number_of_rings() { ($_[0]->length - 2) / 4 }

sub prefix()
{
    my( $self ) = @_;
    return ChemOnomatopist::IUPAC_numerical_multiplier( $self->number_of_rings ) . 'aphene';
}

sub suffix() { $_[0]->prefix }

1;
