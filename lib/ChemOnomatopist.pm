package ChemOnomatopist;

use strict;
use warnings;

use Graph::Traversal::BFS;

sub get_name
{
    my( $graph ) = @_;

    my $last;
    my $bfs = Graph::Traversal::BFS->new( $graph,
                                          # next_root  => undef,
                                        );
    my @order = $bfs->bfs;
    $last = pop @order; # should return the last-visited node
    get_chain( $graph, $last );
    print "e\n";
}

sub get_chain
{
    my( $graph, $start ) = @_;

    my $bfs = Graph::Traversal::BFS->new( $graph,
                                          start => $start,
                                        );
    my @order = $bfs->bfs;
    my %order;
    for my $i (0..$#order) {
        $order{$order[$i]} = $i;
    }

    # Identify the main chain by backtracking
    my $end = pop @order;
    my @chain;
    while( $end ) {
        push @chain, $end;

        my $min = $end;
        for my $neighbour ($graph->neighbours( $end )) {
            next if $order{$neighbour} >= $order{$min};
            $min = $neighbour;
        }
        last if $min eq $end;
        $end = $min;
    }

    @chain = reverse @chain;
    $graph->delete_path( @chain );

    for my $i (0..$#chain) {
        my $atom = $chain[$i];
        for my $neighbour ($graph->neighbours( $atom )) {
            $graph->delete_edge( $atom, $neighbour );
            printf '%d-(', $i+1;
            get_chain( $graph, $neighbour );
            print ')';
        }
    }

    printf '%den-', scalar @chain;
}

1;
