package ChemOnomatopist::Chain::Polyhelicene;

use strict;
use warnings;

# ABSTRACT: Polyhelicenes
# VERSION

use parent ChemOnomatopist::Chain::Circular::;

use ChemOnomatopist::Util::Graph qw(
    graph_without_edge_attributes
    subgraph
);
use Graph::Nauty qw( are_isomorphic );
use Graph::Undirected;
use Set::Object qw( set );

sub new
{
    my( $class, $graph, @cycles ) = @_;

    my @vertices = map { $_->vertices } @cycles;
    my $subgraph = subgraph( $graph, @vertices );
    # Identify the triple-connected vertices
    my $d3 = set( grep { $subgraph->degree( $_ ) == 3 } $subgraph->vertices );
    # Remove them and intermediate atoms
    $subgraph->delete_vertices( @$d3,
                                grep { set( $subgraph->neighbours( $_ ) ) <= $d3 }
                                     $subgraph->vertices );
    # Find the first vertex
    my( $start ) = grep { $subgraph->degree( $_ ) == 1 }
                        $subgraph->vertices;
    $subgraph = subgraph( $graph, @vertices ); # Restore the subgraph
    # Delete the edge which closes the all-encompassing cycle
    $subgraph->delete_edge( $start, grep { $subgraph->degree( $_ ) == 3 }
                                         $subgraph->neighbours( $start ) );
    $subgraph->delete_edges( map  { @$_ }
                             grep { $subgraph->degree( $_->[0] ) == 3 &&
                                    $subgraph->degree( $_->[1] ) == 3 }
                                  $subgraph->edges );
    @vertices = reverse( Graph::Traversal::DFS->new( $subgraph, start => $start )->dfs );

    return bless { graph => $graph, vertices => \@vertices }, $class;
}

sub candidates()
{
    my( $self ) = @_;
    my @candidates = ( $self );

    # TODO: Add the mirror copy

    return @candidates;
}

sub needs_substituent_locants() { return 1 }

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
                           sub { return 'C' } );
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

sub suffix
{
    my( $self ) = @_;
    return ChemOnomatopist::IUPAC_numerical_multiplier( ($self->length - 2) / 4 ) . 'helicene';
}

1;
