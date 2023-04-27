package ChemOnomatopist::Util::Graph;

use strict;
use warnings;

# ABSTRACT: Generic graph utilities
# VERSION

use Exporter;
use Graph::Traversal::BFS;
use List::Util qw( sum0 );

use parent Exporter::;

our @EXPORT_OK = qw(
    BFS_calculate_chain_length
    BFS_is_chain_branched
    graph_center
    graph_has_cycle
    graph_longest_paths
    graph_longest_paths_from_vertex
    graph_path_between_vertices
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
