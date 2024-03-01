package ChemOnomatopist::Chain::Polyacene;

# ABSTRACT: Polyacenes, including anthracene
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
    my @candidates = ( $self,
                       $self->flipped_vertically,
                       $self->flipped_horizontally,
                       $self->flipped_vertically->flipped_horizontally );
    for (1..$#candidates) {
        $candidates[$_]->{candidate_for} = $self;
    }
    return @candidates;
}

sub flipped_horizontally()
{
    my( $self ) = @_;

    my $subgraph = subgraph( $self->graph, $self->vertices );
    my @chords = grep { $subgraph->degree( $_->[0] ) == 3 &&
                        $subgraph->degree( $_->[1] ) == 3 } $subgraph->edges;
    $subgraph->delete_edges( $self->{vertices}[$self->length / 2 + 3],
                             $self->{vertices}[$self->length / 2 + 4],
                             map { @$_ } @chords );
    return bless { graph => $self->graph,
                   vertices => [ reverse Graph::Traversal::DFS->new( $subgraph, start => $self->{vertices}[3] )->dfs ] };
}

sub flipped_vertically()
{
    my( $self ) = @_;

    my $subgraph = subgraph( $self->graph, $self->vertices );
    my @chords = grep { $subgraph->degree( $_->[0] ) == 3 &&
                        $subgraph->degree( $_->[1] ) == 3 } $subgraph->edges;
    $subgraph->delete_edges( $self->{vertices}[3],
                             $self->{vertices}[4],
                             map { @$_ } @chords );
    return bless { graph => $self->graph,
                   vertices => [ reverse Graph::Traversal::DFS->new( $subgraph, start => $self->{vertices}[3] )->dfs ] };
}

sub needs_substituent_locants() { 1 }

sub has_form($$)
{
    my( $class, $graph ) = @_;
    my %degrees = map { $graph->degree( $_ ) => 1 } $graph->vertices;
    return '' unless join( ',', sort keys %degrees ) eq '2,3';

    my $N = scalar $graph->vertices;
    return '' if $N < 14;
    return '' if ($N - 14) % 4;

    return are_isomorphic( graph_without_edge_attributes( $graph ),
                           $class->ideal_graph( $N ),
                           sub { return 'C' } );
}

sub ideal_graph($$)
{
    my( $class, $N ) = @_;
    die "cannot construct polyacene with $N vertices\n" if $N < 10 || ($N-2) % 4;

    my $graph = Graph::Undirected->new( refvertexed => 1 );
    my @vertices = map { { symbol => 'C', number => $_-1 } } 1..$N;
    $graph->add_cycle( @vertices );
    for (0..($N-6) / 4 -1) {
        $graph->add_edge( map { $vertices[$_] } ( 4 + 2*$_, $N - 1 - 2*$_ ) );
    }
    return $graph;
}

sub suffix
{
    my( $self ) = @_;
    return 'anthracene' if $self->length == 14;
    return ChemOnomatopist::IUPAC_numerical_multiplier( ($self->length - 2) / 4 ) . 'acene';
}

1;
