package ChemOnomatopist::Group::Polyacene;

use strict;
use warnings;

# ABSTRACT: Polyacene
# VERSION

use parent ChemOnomatopist::Group::, ChemOnomatopist::Chain::Circular::;

use Graph::Undirected;

sub ideal_graph($$)
{
    my( $class, $N ) = @_;
    die "cannot construct polyacene with $N vertices\n" if $N < 18 || ($N-2) % 4;

    my $graph = Graph::Undirected->new( refvertexed => 1 );
    my @vertices = map { { symbol => 'C', number => $_-1 } } 1..$N;
    $graph->add_cycle( @vertices );
    for (0..($N-6) / 4 -1) {
        $graph->add_edge( map { $vertices[$_] } ( 4 + 2*$_, $N - 1 - 2*$_ ) );
    }
    return $graph;
}

1;
