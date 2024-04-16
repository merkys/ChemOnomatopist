package ChemOnomatopist::Chain::Bicycle::Purine;

# ABSTRACT: Purine chain
# VERSION

use strict;
use warnings;

use Graph::Traversal::DFS;
use List::Util qw( first );

use parent ChemOnomatopist::Chain::Bicycle::;

sub new
{
    my( $class, $graph, @cycles ) = @_;

    my $pyrimidine = first { $_->length == 6 } @cycles;
    if( ChemOnomatopist::element( $pyrimidine->{vertices}[0] ) eq 'N' ) {
        $pyrimidine = $pyrimidine->flipped;
    }

    my $subgraph = $graph->subgraph( map { $_->vertices } @cycles );
    $subgraph->delete_edge( @{$pyrimidine->{vertices}}[0..1] );
    $subgraph->delete_edge( map  { @$_ }
                            grep { $subgraph->degree( $_->[0] ) == 3 &&
                                   $subgraph->degree( $_->[1] ) == 3 }
                                 $subgraph->edges );
    my @vertices = Graph::Traversal::DFS->new( $subgraph, start => $pyrimidine->{vertices}[0] )->dfs;
    return bless { graph => $graph,
                   cycles => \@cycles,
                   vertices => \@vertices }, $class;
}

sub is_purine() { 1 }

sub locants(@)
{
    my $self = shift;
    return map { $_ + 1 } @_;
}

sub suffix() { ChemOnomatopist::Name->new( 'purine' ) }

1;
