package ChemOnomatopist::Chain::Polyhelicene;

# ABSTRACT: Polyhelicenes
# VERSION

use strict;
use warnings;

use parent ChemOnomatopist::Chain::Circular::;

use ChemOnomatopist::Util::Graph qw(
    graph_without_edge_attributes
    subgraph
);
use Graph::Nauty qw( are_isomorphic );
use Graph::Undirected;
use List::Util qw( first );
use Set::Object qw( set );

sub new
{
    my( $class, $graph, @cycles ) = @_;

    my @vertices = map { $_->vertices } @cycles;

    # Take a subgraph of only triple-connected vertices
    my $g3 = subgraph( $graph, @vertices );
    # Leave only triple-connected vertices
    $g3->delete_vertices( grep { $g3->degree( $_ ) < 3 } $g3->vertices );
    # There will be only two two-connected vertices now
    my $end = first { $g3->degree( $_ ) == 2 } $g3->vertices;

    # Take a subgraph of all participating vertices
    my $subgraph = subgraph( $graph, @vertices );
    # Find the first vertex
    my $start = first { $subgraph->degree( $_ ) == 2 } $subgraph->neighbours( $end );

    # Delete all bridges from the subgraph and the last edge as well
    $subgraph->delete_edges( $start, $end,
                             map  { @$_ }
                             grep { $subgraph->degree( $_->[0] ) == 1 ||
                                    $subgraph->degree( $_->[1] ) == 1 }
                                    $subgraph->edges );

    @vertices = reverse( Graph::Traversal::DFS->new( $subgraph, start => $start )->dfs );

    return bless { graph => $graph, vertices => \@vertices }, $class;
}

sub candidates()
{
    my( $self ) = @_;
    my @candidates = ( $self );

    my @vertices = reverse $self->vertices;
    for (1..5) {
        push @vertices, shift @vertices;
    }
    push @candidates, bless { graph => $self->graph,
                              vertices => \@vertices,
                              candidate_for => $self };

    return @candidates;
}

sub needs_substituent_locants() { 1 }

sub has_form($$)
{
    my( $class, $graph ) = @_;
    my %degrees = map { $graph->degree( $_ ) => 1 } $graph->vertices;
    return '' unless join( ',', sort keys %degrees ) eq '2,3';

    my $N = scalar $graph->vertices;
    return '' if ($N-2) % 4;
    return '' if ($N-2) / 4 < 6;

    return are_isomorphic( graph_without_edge_attributes( $graph ),
                           $class->ideal_graph( $N ),
                           sub { 'C' } );
}

sub ideal_graph($$)
{
    my( $class, $N ) = @_;
    die "cannot construct helicene with $N vertices\n" if ($N-2) % 4 || ($N-2) / 4 < 6;

    my $graph = Graph::Undirected->new( refvertexed => 1 );
    $graph->add_cycle( map { { symbol => 'C' } } 1..6 );
    my( $edge ) = $graph->edges;
    for (2..($N-2) / 4) {
        my @vertices = map { { symbol => 'C' } } 1..4;
        $graph->add_path( $edge->[0], @vertices, $edge->[1] );
        $edge = [ @vertices[0..1] ];
    }
    return $graph;
}

sub number_of_rings() { ($_[0]->length - 2) / 4 }

sub suffix
{
    my( $self ) = @_;
    return ChemOnomatopist::IUPAC_numerical_multiplier( $self->number_of_rings, 1 ) . 'helicene';
}

1;
