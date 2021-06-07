package ChemOnomatopist;

use strict;
use warnings;

use Graph::Traversal::BFS;

our @prefixes = qw(
    ?
    meth
    eth
    prop
    but
    bent
    hex
    hept
    oct
    non
    dec
    undec
    dodec
    tridec
    tetradec
    pentadec
    hexadec
    heptadec
    octadec
    nonadec
    icos
);

sub get_name
{
    my( $graph ) = @_;

    my $last;
    my $bfs = Graph::Traversal::BFS->new( $graph,
                                          # next_root => undef,
                                        );
    my @order = $bfs->bfs;
    get_chain( $graph, pop @order, { choose_direction => 1 } );
    print "ane\n";
}

sub get_chain
{
    my( $graph, $start, $options ) = @_;

    $options = {} unless $options;

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

    # Establishing a stable order similarly as suggested by the IUPAC
    # rules
    if( $options->{choose_direction} ) {
        for my $i (0..int(@chain/2)-1) {
            if( $graph->degree( $chain[$i] ) !=
                $graph->degree( $chain[$#chain-$i] ) ) {
                if( $graph->degree( $chain[$i] ) <
                    $graph->degree( $chain[$#chain-$i] ) ) {
                    @chain = reverse @chain;
                }
                last;
            }
        }
    }

    for my $i (0..$#chain) {
        my $atom = $chain[$i];
        for my $neighbour ($graph->neighbours( $atom )) {
            $graph->delete_edge( $atom, $neighbour );
            printf '%d-(', $i+1;
            get_chain( $graph, $neighbour );
            print 'yl)-';
        }
    }

    print $prefixes[scalar @chain];
}

1;
