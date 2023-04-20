package ChemOnomatopist::ChainHalf;

use strict;
use warnings;

use ChemOnomatopist;
use ChemOnomatopist::Util::Graph qw(
    tree_branch_positions
);
use Graph::Traversal::DFS;
use List::Util qw( sum0 );
use Memoize;

# ABSTRACT: Half of a longest chain
# VERSION

# memoize( 'branch_positions' );

sub new
{
    my( $class, $graph, $other_center, @vertices ) = @_;
    my $self = { vertices => \@vertices, graph => $graph, other_center => $other_center };
    return bless $self, $class;
}

# Accessors

# Groups are used to check which halves of a chain can be combined together.
# If a graph contains single center, all halves will share the center.
sub group()
{
    my( $self ) = @_;
    return $self->{vertices}[1 - defined $self->{other_center}];
}

sub vertices()
{
    my( $self ) = @_;
    return @{$self->{vertices}};
}

# Properties

sub branch_positions()
{
    my( $self ) = @_;

    my $graph = $self->_disconnected_chain_graph;
    my @vertices = $self->vertices;

    return map  { ( $_ ) x $graph->degree( $vertices[$_] ) }
           grep { $graph->degree( $vertices[$_] ) }
                0..$#vertices;
}

sub length()
{
    my( $self ) = @_;
    return scalar $self->vertices;
}

sub locant_names()
{
    my( $self ) = @_;

    my $graph = $self->_disconnected_chain_graph;

    my @locants;
    for my $vertex ($self->vertices) {
        my @current_locants;
        for my $neighbour ($graph->neighbours( $vertex )) {
            $graph->delete_edge( $vertex, $neighbour );
            push @current_locants, ChemOnomatopist::get_sidechain_name( $graph, $neighbour );
        }
        push @locants, \@current_locants;
    }

    return @locants;
}

# Sum of locant positions (1-based) as calculated as if this chain half is terminal.
sub locant_positions_forward()
{
    my( $self ) = @_;
    my @locants = $self->branch_positions;
    return ($self->length + (defined $self->{other_center})) * @locants + sum0 @locants;
}

# Sum of locant positions (1-based) as calculated as if this chain half is the first part.
sub locant_positions_backward()
{
    my( $self ) = @_;
    my @locants = $self->branch_positions;
    return $self->length * @locants - sum0 @locants;
}

sub number_of_branches_in_sidechains()
{
    my( $self ) = @_;

    my $graph = $self->_disconnected_chain_graph;
    my @vertex_neighbours = map { $graph->neighbours( $_ ) } $self->vertices;
    $graph->delete_vertices( $self->vertices );

    return sum0 map { $_ > 2 ? $_ - 2 : 0 }
                map { $graph->degree( $_ ) }
                map { Graph::Traversal::DFS->new( $graph, start => $_ )->dfs }
                    @vertex_neighbours;
}

sub number_of_carbons()
{
    my( $self ) = @_;

    my $graph = $self->_disconnected_chain_graph;

    my $C = grep { ChemOnomatopist::is_element( $_, 'C' ) }
            map  { Graph::Traversal::DFS->new( $graph, start => $_ )->dfs }
                 $self->vertices;

    # Since main chain carbons are included in the count, they have to be subtracted.
    return $C - $self->length;
}

sub number_of_branches()
{
    my( $self ) = @_;
    return scalar $self->branch_positions;
}

sub _disconnected_chain_graph()
{
    my( $self ) = @_;

    my $graph = $self->{graph}->copy;
    my @vertices = $self->vertices;

    if( $self->{other_center} ) {
        # Cut the edge to the other center
        $graph->delete_edge( $vertices[0], $self->{other_center} );
    } else {
        # Cut the edges to the other candidates
        for ($graph->neighbours( $vertices[0] )) {
            $graph->delete_edge( $vertices[0], $_ );
        }
    }
    $graph->delete_path( @vertices );

    return $graph;
}

1;
