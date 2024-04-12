package ChemOnomatopist::Chain::Aceylene;

# ABSTRACT: Ace...ylene chain, as per BBv3 P-25.1.2.7
# VERSION

use strict;
use warnings;

use parent ChemOnomatopist::Chain::Circular::;

use ChemOnomatopist::Chain::Phenanthrene;
use ChemOnomatopist::Chain::Polyacene;
use ChemOnomatopist::Util::Graph qw(
    graph_without_edge_attributes
    subgraph
);
use Graph::Nauty qw( are_isomorphic );
use Graph::Traversal::DFS;
use Graph::Undirected;
use List::Util qw( first );

sub new
{
    my( $class, $graph, @cycles ) = @_;

    my $pentane = first { $_->length == 5 } @cycles;
    my $subgraph = $graph->subgraph( map { $_->vertices } @cycles );
    my %cycles_per_atom;
    for (@cycles) {
        for ($_->vertices) {
            $cycles_per_atom{$_}++;
        }
    }
    my $center = first { $cycles_per_atom{$_} == 3 } $subgraph->vertices;

    my @vertices;
    if( @cycles == 3 ) {
        my $first = first { $cycles_per_atom{$_} == 1 } $pentane->vertices;
        my $last  = first { $cycles_per_atom{$_} == 2 } $subgraph->neighbours( $first );
        $subgraph->delete_edge( $first, $last );
        $subgraph->delete_vertex( $center );
        @vertices = reverse Graph::Traversal::DFS->new( $subgraph, start => $first )->dfs;
        @vertices = ( @vertices[0..1], $center, @vertices[2..$#vertices] );
    } else {
        my $d3_subgraph = $subgraph->subgraph( grep { $subgraph->degree( $_ ) == 3 } $subgraph->vertices );
        my $last = first { $d3_subgraph->degree( $_ ) == 2 }
                         $subgraph->neighbours( $center );
        my $first = first { $cycles_per_atom{$_} == 1 } $subgraph->neighbours( $last );
        $subgraph->delete_edge( $first, $last );
        $subgraph->delete_vertex( $center );
        $subgraph->delete_edges( map  { @$_ }
                                 grep { $subgraph->degree( $_->[0] ) == 3 &&
                                        $subgraph->degree( $_->[1] ) == 3 }
                                      $subgraph->edges );
        @vertices = reverse Graph::Traversal::DFS->new( $subgraph, start => $first )->dfs;

        # Restore the original subgraph
        $subgraph = $graph->subgraph( map { $_->vertices } @cycles );
        my $first_d3 = first { $subgraph->degree( $vertices[$_] ) == 3 }
                             0..$#vertices;
        @vertices = ( @vertices[0..$first_d3-1], $center, @vertices[$first_d3..$#vertices] );
    }

    return bless { graph => $graph, vertices => \@vertices }, $class;
}

sub is_acenaphthylene()    { $_[0]->length == 12 }
sub is_aceanthrylene()     { $_[0]->length == 16 && $_[0]->graph->degree( $_[0]->{vertices}[2] ) == 3 } # CHECKME: This should probably use subgraph
sub is_acephenanthrylene() { $_[0]->length == 16 && $_[0]->graph->degree( $_[0]->{vertices}[3] ) == 3 } # CHECKME: This should probably use subgraph

sub has_form($$)
{
    my( $class, $graph ) = @_;

    return '' unless $graph->vertices == 12 || $graph->vertices == 16;

    return 1 if are_isomorphic( graph_without_edge_attributes( $graph ),
                                $class->ideal_graph_acenaphthylene,
                                sub { ChemOnomatopist::element( $_[0] ) } );
    return 1 if are_isomorphic( graph_without_edge_attributes( $graph ),
                                $class->ideal_graph_aceanthrylene,
                                sub { ChemOnomatopist::element( $_[0] ) } );
    return 1 if are_isomorphic( graph_without_edge_attributes( $graph ),
                                $class->ideal_graph_acephenanthrylene,
                                sub { ChemOnomatopist::element( $_[0] ) } );
    return '';
}

sub ideal_graph_acenaphthylene()
{
    my( $class ) = @_;

    my $graph = Graph::Undirected->new( refvertexed => 1 );
    my @vertices = map { { symbol => 'C' } } 1..11;
    $graph->add_cycle( @vertices );
    my $center = { symbol => 'C' };
    $graph->add_path( $vertices[2], $center, $vertices[6] );
    $graph->add_path( $vertices[2], $center, $vertices[10] );

    return $graph;
}

sub ideal_graph_aceanthrylene()
{
    my( $class ) = @_;

    my $graph = ChemOnomatopist::Chain::Polyacene->ideal_graph( 14 );
    my $d3 = first { $graph->degree( $_ ) == 3 } $graph->vertices;
    my @d2 = grep  { $graph->degree( $_ ) == 2 } $graph->neighbours( $d3 );
    $graph->add_path( $d2[0], { symbol => 'C' }, { symbol => 'C' }, $d2[1] );

    return $graph;
}

sub ideal_graph_acephenanthrylene()
{
    my( $class ) = @_;

    my $graph = ChemOnomatopist::Chain::Phenanthrene->ideal_graph;
    my $d3_subgraph = subgraph( $graph, grep { $graph->degree( $_ ) == 3 } $graph->vertices );
    my $d3 = first { $d3_subgraph->degree( $_ ) == 1 } $d3_subgraph->vertices;
    my @d2 = grep { $graph->degree( $_ ) == 2 } $graph->neighbours( $d3 );
    $graph->add_path( $d2[0], { symbol => 'C' }, { symbol => 'C' }, $d2[1] );

    return $graph;
}

# TODO: Override locants() to give the internal vertex a number

sub prefix() { &suffix }
sub suffix()
{
    my( $self ) = @_;
    return 'acenaphthylene'    if $self->is_acenaphthylene;
    return 'aceanthrylene'     if $self->is_aceanthrylene;
    return 'acephenanthrylene' if $self->is_acephenanthrylene;
    die "unknown ace...ylene\n";
}

1;
