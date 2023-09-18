package ChemOnomatopist::Group::Polyacene;

use strict;
use warnings;

# ABSTRACT: Polyacene
# VERSION

use parent ChemOnomatopist::Group::, ChemOnomatopist::Chain::Circular::;

use ChemOnomatopist::Util::Graph qw( subgraph );
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

    my $subgraph = subgraph( $self->graph, $self->vertices );
    my @chords = grep { $subgraph->degree( $_->[0] ) == 3 &&
                        $subgraph->degree( $_->[1] ) == 3 } $subgraph->edges;
    $subgraph->delete_edges( $self->{vertices}[3],
                             $self->{vertices}[4],
                             map { @$_ } @chords );

    push @candidates,
         bless { graph => $self->graph,
                 vertices => [ reverse Graph::Traversal::DFS->new( $subgraph, start => $self->{vertices}[3] )->dfs ],
                 candidate_for => $self };

    return @candidates;
}

sub locants(@)
{
    my $self = shift;
    my $N = $self->length;
    my @locant_map = ( 1..3,
                       (map { 4 + $_, (4 + $_) . 'a' } 0..($N-6) / 4 - 1),
                       ($N / 2)..($N / 2 + 2),
                       (map { $N / 2 + 3 + $_, ($N / 2 + 3 + $_) . 'a' } 0..($N-6) / 4 - 1) );
    return map { $locant_map[$_] } @_;
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
    return ChemOnomatopist::IUPAC_numerical_multiplier( ($self->length - 2) / 4 ) . 'acene';
}

1;
