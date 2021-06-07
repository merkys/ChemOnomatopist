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
}

sub get_chain
{
    my( $graph, $start ) = @_;

    my %order;
    my $bfs = Graph::Traversal::BFS->new( $graph,
                                          start => $start,
                                          pre  => sub { $order{$_[0]} = scalar %order },
                                        );
    my @order = $bfs->bfs;

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

    print STDERR sprintf "%d-mer\n", scalar @chain;
}

1;
