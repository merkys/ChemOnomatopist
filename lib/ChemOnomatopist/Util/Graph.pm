package ChemOnomatopist::Util::Graph;

use strict;
use warnings;

# ABSTRACT: Generic graph utilities
# VERSION

use ChemOnomatopist::Util qw( copy );
use Exporter;
use Graph::Traversal::BFS;
use List::Util qw( any sum0 );
use Set::Object qw( set );

use parent Exporter::;

our @EXPORT_OK = qw(
    BFS_calculate_chain_length
    BFS_is_chain_branched
    cyclic_components
    graph_center
    graph_cycle_core
    graph_cycles
    graph_has_cycle
    graph_longest_paths
    graph_longest_paths_from_vertex
    graph_path_between_vertices
    graph_replace
    merge_graphs
    subgraph
    tree_branch_positions
    tree_number_of_branches
);

# Calculates length of given graph (vertices count)
sub BFS_calculate_chain_length
{
    my( $graph, $start ) = @_;
    my $bfs = Graph::Traversal::BFS->new( $graph, start => $start );
    return scalar $bfs->bfs;
}

# Returns 1 if there is any braches in the given graph and if there is none
sub BFS_is_chain_branched
{
    my( $graph, $start ) = @_;

    # FIXME: Not entirely sure why visited vertices are removed (A.M.)
    my $graph_copy = $graph->copy;

    my $branched = 0;
    my $bfs = Graph::Traversal::BFS->new(
        $graph,
        pre => sub {
            my @neighbours = $graph_copy->neighbours( $_[0] );
            $branched = 1 if scalar @neighbours > 1;
            $graph_copy->delete_vertex( $_[0] );
        },
        start => $start );
    $bfs->bfs;

    return $branched;
}

sub cyclic_components
{
    my( $graph ) = @_;

    $graph = copy $graph;

    # Due to the issue in Graph, bridges() returns strings instead of real objects.
    # Graph issue: https://github.com/graphviz-perl/Graph/issues/29
    my %vertices_by_name = map { $_ => $_ } $graph->vertices;
    $graph->delete_edges( map { map { $vertices_by_name{$_} } @$_ } $graph->bridges );
    $graph->delete_vertices( grep { !$graph->degree( $_ ) } $graph->vertices );

    return () unless $graph->vertices; # No vertices = no cycles

    return map { $graph->subgraph( $_ ) } $graph->connected_components;
}

# Find how many side attachments are at every position of the given path.
sub graph_attachment_positions
{
    my( $graph, @vertices ) = @_;

    $graph = $graph->copy;
    $graph->delete_path( @vertices );

    return map { $graph->degree( $_ ) } @vertices;
}

# Finds center (or two centers) of a tree graph.
# Returns one or two vertices constituting the center.
sub graph_center
{
    my( $graph ) = @_;

    $graph = $graph->copy;
    my $nvertices = scalar $graph->vertices;
    while( $graph->vertices > 2 ) {
        $graph->delete_vertices( grep { $graph->degree( $_ ) == 1 }
                                      $graph->vertices );
        my $nvertices_now = scalar $graph->vertices;
        if( $nvertices_now == $nvertices ) {
            # Safeguard for cycles and/or isolated vertices
            die 'cannot find center of cyclic or isolated graphs';
        }
        $nvertices = $nvertices_now;
    }
    return $graph->vertices;
}

# Iteratively removes leaves of a graph until cycle core remains.
sub graph_cycle_core
{
    my( $graph ) = @_;

    $graph = $graph->copy;
    while( my @leaves = grep { $graph->degree( $_ ) == 1 } $graph->vertices ) {
        $graph->delete_vertices( @leaves );
    }
    return $graph;
}

# Decomposes bridgeless graph into cycles
# FIXME: Experimental
sub graph_cycles
{
    my( $graph ) = @_;

    $graph = copy $graph;

    my @cycles;
    while( $graph->vertices ) {
        my $triple_connected = set( grep { $graph->degree( $_ ) == 3 } $graph->vertices );
        if( !@$triple_connected ) {
            push @cycles, $graph;
            last;
        }

        my @chords = grep { $triple_connected->has( $_->[0] ) &&
                            $triple_connected->has( $_->[1] ) } $graph->edges;
        my $wo_chords = copy( $graph )->delete_vertices( @$triple_connected );
        for my $component ($wo_chords->connected_components) {
            next if @$component == 1;
            my @ends = grep { $wo_chords->degree( $_ ) == 1 } @$component; # Should be two
            my( $chord ) = grep { ($graph->has_edge( $_->[0], $ends[0] ) && $graph->has_edge( $_->[1], $ends[1] )) ||
                                  ($graph->has_edge( $_->[0], $ends[1] ) && $graph->has_edge( $_->[1], $ends[0] )) } @chords;
            next unless $chord;
            # Found a chord which completes a cycle
            my $vertices = set( @$component, @$chord );
            my $cycle = copy( $graph )->delete_vertices( grep { !$vertices->has( $_ ) } $graph->vertices );
            push @cycles, $cycle;
            $graph->delete_vertices( @$component );
        }
    }
    return @cycles;
}

sub graph_has_cycle
{
    my( $graph ) = @_;
    return $graph->edges != $graph->vertices - 1;
}

# Finds longest paths in a tree graph.
# The subroutine finds all longest paths originating at graph center(s) and produces all their combinations.
# No two paths containing the same vertices are returned.
sub graph_longest_paths
{
    my( $graph ) = @_;

    my @centers = graph_center( $graph );
    my @longest_paths;
    if( @centers == 1 ) {
        # Single-centered graph
        # Removing the center from longest path parts, to be added later
        my @longest_path_parts = map { [ @{$_}[1..$#$_] ] }
                                     graph_longest_paths_from_vertex( $graph, $centers[0] );
        for my $i (0..$#longest_path_parts) {
            for my $j (0..$#longest_path_parts) {
                next if $i >= $j;
                # Ensure that two paths do not start at the same vertex.
                next if $longest_path_parts[$i]->[0] eq $longest_path_parts[$j]->[0];
                push @longest_paths, [ reverse( @{$longest_path_parts[$i]} ),
                                       $centers[0],
                                       @{$longest_path_parts[$j]} ];
            }
        }
    } else {
        # Double-centered graph
        $graph = $graph->copy;
        $graph->delete_edge( @centers );
        my @longest_path_parts1 = graph_longest_paths_from_vertex( $graph, $centers[0] );
        my @longest_path_parts2 = graph_longest_paths_from_vertex( $graph, $centers[1] );
        for my $i (0..$#longest_path_parts1) {
            for my $j (0..$#longest_path_parts2) {
                push @longest_paths, [ reverse( @{$longest_path_parts1[$i]} ),
                                                @{$longest_path_parts2[$j]} ];
            }
        }
    }

    return @longest_paths;
}

# Finds all the longest paths from given vertex
sub graph_longest_paths_from_vertex
{
    my( $graph, $vertex ) = @_;

    my %from   = ( $vertex => undef );
    my %length = ( $vertex => 0 );
    my $bfs = Graph::Traversal::BFS->new(
        $graph,
        tree_edge =>
            sub {
                my( $u, $v ) = @_;
                ( $u, $v ) = ( $v, $u ) if exists $from{$v};
                $from{$v} = $u;
                $length{$v} = $length{$u} + 1;
            },
        start => $vertex,
    );
    $bfs->bfs;

    my @furthest_leaves;
    my $furthest_distance = 0;
    for my $vertex ( $graph->vertices ) {
        next unless exists $length{$vertex}; # May happen in disconnected graphs
        if(      $length{$vertex} < $furthest_distance ) {
            next;
        } elsif( $length{$vertex} == $furthest_distance ) {
            push @furthest_leaves, $vertex;
        } else {
            @furthest_leaves = ( $vertex );
            $furthest_distance = $length{$vertex};
        }
    }

    # Backtrack starting from the furthest leaves to collect the longest
    # paths. In the returned result path, starting vertex is the first.
    my @longest_paths;
    for my $vertex ( @furthest_leaves ) {
        my @path;
        while( defined $vertex ) {
            push @path, $vertex;
            $vertex = $from{$vertex};
        }
        push @longest_paths, [ reverse @path ];
    }

    return @longest_paths;
}

# Replace one or more old vertices with a single new one
sub graph_replace
{
    my( $graph, $new, @old ) = @_;

    $graph->add_vertex( $new );

    my $old = set( @old );
    for my $edge (grep { ($old->has( $_->[0] ) && !$old->has( $_->[1] )) ||
                         ($old->has( $_->[1] ) && !$old->has( $_->[0] )) }
                       $graph->edges) {
        my( $vertex, $neighbour ) = $old->has( $edge->[0] ) ? @$edge : reverse @$edge;
        next if $graph->has_edge( $new, $neighbour );
        $graph->add_edge( $new, $neighbour );
        next unless $graph->has_edge_attributes( @$edge );
        $graph->set_edge_attributes( $new, $neighbour, $graph->get_edge_attributes( @$edge ) );
    }
    $graph->delete_vertices( @old );

    return $graph;
}

sub merge_graphs
{
    my( $A, $B, $A_vertex, $B_vertex ) = @_;

    my $merged = copy $A;
    for my $edge ($B->edges) {
        $merged->add_edge( @$edge );
        next unless $B->has_edge_attributes( @$edge );
        $merged->set_edge_attributes( @$edge, $B->get_edge_attributes( @$edge ) );
    }
    $merged->add_edge( $A_vertex, $B_vertex ) if $A_vertex && $B_vertex;

    return $merged;
}

sub subgraph
{
    my( $graph, @vertices ) = @_;
    my $vertices = set( @vertices );
    $graph = copy $graph;
    $graph->delete_vertices( grep { !$vertices->has( $_ ) } $graph->vertices );
    return $graph;
}

# Given a tree and a path, finds the number of branches branching off the given path.
# It is equal to the sum of all degrees minus the edges between vertices in a path.
sub tree_number_of_branches
{
    my( $tree, @vertices ) = @_;
    return sum0( map { $tree->degree( $_ ) } @vertices ) - 2 * (scalar @vertices - 1);
}

# Returns a list of 0-based indices of branch positions.
sub tree_branch_positions
{
    my( $tree, @vertices ) = @_;
    return map  { ( $_ ) x ( $tree->degree( $vertices[$_] ) - 2 ) }
           grep { $tree->degree( $vertices[$_] ) > 2 }
                0..$#vertices;
}

# Find a path between two vertices in an acyclic graph.
sub graph_path_between_vertices
{
    my( $graph, $A, $B ) = @_;

    if( graph_has_cycle( $graph ) ) {
        die "cannot call graph_path_between_vertices() on graph with cycles\n";
    }

    $graph = $graph->copy;
    while( my @leaves = grep { $graph->degree( $_ ) == 1 && $_ != $A && $_ != $B } $graph->vertices ) {
        $graph->delete_vertices( @leaves );
    }

    my @path;
    my $vertex = $A;
    while( $vertex ) {
        push @path, $vertex;
        my( $vertex_now ) = $graph->neighbours( $vertex );
        $graph->delete_vertex( $vertex );
        $vertex = $vertex_now;
    }

    return @path;
}

1;
