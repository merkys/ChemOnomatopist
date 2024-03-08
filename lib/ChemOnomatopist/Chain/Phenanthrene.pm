package ChemOnomatopist::Chain::Phenanthrene;

# ABSTRACT: Phenanthrene or its derivative
# VERSION

use strict;
use warnings;

use ChemOnomatopist::Chain::Polyaphene;
use ChemOnomatopist::Name;
use ChemOnomatopist::Util::Graph qw( merge_graphs );
use Graph::Undirected;
use List::Util qw( first all any );

use parent ChemOnomatopist::Chain::Polyaphene::;

sub new
{
    my( $class, $graph, @cycles ) = @_;

    my $subgraph = $graph->subgraph( map { $_->vertices } @cycles );

    # Deleting all edges having degree 3 vertices at both ends
    $subgraph->delete_edges( map  { @$_ }
                             grep { $subgraph->degree( $_->[0] ) == 3 &&
                                    $subgraph->degree( $_->[1] ) == 3 }
                                  $subgraph->edges );

    # Find an order
    my( $start ) = grep { $subgraph->degree( $_ ) == 1 } $subgraph->vertices;
    my @vertices = Graph::Traversal::DFS->new( $subgraph, start => $start )->dfs;

    # Adjust the order
    if( any { ChemOnomatopist::element( $_ ) eq 'N' } @vertices ) {
        # Find the order so as N is closest to the begining of the chain
        # CHECKME: This might not be correct due to offset
        my $first = first { ChemOnomatopist::element( $vertices[$_] ) eq 'N' }
                          0..$#vertices;
        my $last  = first { ChemOnomatopist::element( $vertices[-1-$_] ) eq 'N' }
                          0..$#vertices;
        @vertices = reverse @vertices if $last < $first;
        push @vertices, shift @vertices;

        # Phenanthridine has a strict order
        if( (grep { ChemOnomatopist::element( $_ ) eq 'N' } @vertices) == 1 &&
            ChemOnomatopist::element( $vertices[5] ) ne 'N' ) {
            die "cannot handle complicated cyclic compounds\n";
        }
    } else {
        for (1..5) {
            push @vertices, shift @vertices;
        }
        @vertices = reverse @vertices;
    }

    return bless { graph => $graph, vertices => \@vertices }, $class;
}

sub candidates()
{
    my( $self ) = @_;
    my @candidates = ( $self );

    if( $self->is_hydrocarbon ) {
        my @vertices = $self->vertices;
        push @candidates,
             bless { graph => $self->graph,
                     vertices => [ reverse map { $vertices[$_] } ( 10..13, 0..9 ) ],
                     candidate_for => $self };
    }

    return @candidates;
}

sub ideal_graph($)
{
    my( $class ) = @_;

    my @graphs;
    for (0..1) {
        my $graph = Graph::Undirected->new( refvertexed => 1 );
        $graph->add_cycle( map { { symbol => 'C', number => $_-1 } } 1..6 );
        push @graphs, $graph;
    }
    my $graph = merge_graphs( @graphs );

    # Pick an edge from each graph
    my( $A ) = $graphs[0]->edges;
    my( $B ) = $graphs[1]->edges;

    # Join a pair of atoms with an edge
    $graph->add_edge( $A->[0], $B->[0] );

    # Add a longer arc between other two atoms
    $graph->add_path( $A->[1], { symbol => 'C' }, { symbol => 'C' }, $B->[1] );

    return $graph;
}

sub needs_heteroatom_locants()
{
    my( $self ) = @_;
    return $self->number_of_heteroatoms == 2;
}

sub needs_heteroatom_names() { '' }

sub prefix()
{
    my( $self ) = @_;

    if( all { $_ eq 'N' } $self->heteroatoms ) {
        return ChemOnomatopist::Name->new( 'phenanthridine' ) if $self->number_of_heteroatoms == 1;
        return ChemOnomatopist::Name->new( 'phenanthroline' ) if $self->number_of_heteroatoms == 2;
    }

    if( $self->number_of_heteroatoms == 1 ) {
        return ChemOnomatopist::Name->new( 'arsanthridine' )    if all { $_ eq 'As' } $self->heteroatoms;
        return ChemOnomatopist::Name->new( 'phosphanthridine' ) if all { $_ eq 'P'  } $self->heteroatoms;
    }

    return ChemOnomatopist::Name->new( 'phenanthrene' );
}

sub suffix() { $_[0]->prefix }

1;
