package ChemOnomatopist;

use strict;
use warnings;

use ChemOnomatopist::Group;
use ChemOnomatopist::Group::Carbonyl;
use Graph::Traversal::BFS;
use Scalar::Util qw(blessed);

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
    henicos
    docos
    tricos
    tetracos
    pentacos
    hexacos
    heptacos
    octacos
    nonacos
    triacont
    hentriacont
    dotriacont
    tritriacont
    tetratriacont
    pentatriacont
    hexatriacont
    heptatriacont
    octatriacont
    nonatriacont
    tetracont
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

sub find_groups
{
    my( $graph ) = @_;

    for my $atom ($graph->vertices) {
        my @neighbours = $graph->neighbours( $atom );

        # Detecting carbonyl
        if( is_element( $atom, 'O' ) &&
            scalar @neighbours == 1 &&
            is_element( $neighbours[0], 'C' ) ) {
            my $carbonyl = ChemOnomatopist::Group::Carbonyl->new( $neighbours[0] );
            for ($graph->neighbours( $neighbours[0] )) {
                $graph->add_edge( $_, $carbonyl );
                $graph->delete_edge( $_, $neighbours[0] );
            }
            $graph->delete_vertex( $atom );

            # Carbonyl derivatives should be detected here
        # Detecting hydroxy
        } elsif( is_element( $atom, 'O' ) &&
                 scalar @neighbours == 2 &&
                 grep { is_element( $_, 'H' ) } @neighbours == 1 ) {
            my $hydroxy  = ChemOnomatopist::Group::Hydroxy->new( $atom );
            my $hydrogen = grep { is_element( $_, 'H' ) } @neighbours;
            for (@neighbours) {
                $graph->add_edge( $_, $hydroxy );
                $graph->delete_edge( $_, $atom );
            }
            $graph->delete_vertex( $hydrogen );
        }
    }

    # Second pass is needed to build on top of these trivial groups
    for my $atom ($graph->vertices) {
        my @neighbours = $graph->neighbours( $atom );

        # Detecging carboxyl
        if( blessed $atom &&
            $atom->isa( ChemOnomatopist::Group::Hydroxy:: ) &&
            scalar @neighbours == 1 &&
            blessed $neighbours[0] &&
            $neighbours[0]->isa( ChemOnomatopist::Group::Carbonyl:: ) ) {
            my( $hydroxy, $carbonyl ) = ( $atom, @neighbours );
            my $carboxyl = ChemOnomatopist::Group::Carbonyl->new( $carbonyl );
            for ($graph->neighbours( $carbonyl )) {
                $graph->add_edge( $_, $carboxyl );
                $graph->delete_edge( $_, $carbonyl );
            }
            $graph->delete_vertex( $carboxyl );
            $graph->delete_vertex( $hydroxyl );
        }
    }
}

sub is_element
{
    my( $atom, $element ) = @_;

    return unless ref $atom;

    if( blessed $atom ) {
        if( $atom->isa( ChemOnomatopist::Group:: ) ) {
            if(      $element eq 'C' ) {
                return $atom->is_carbon;
            } elsif( $element eq 'O' ) {
                return $atom->is_oxygen;
            }
            warn "cannot say whether $atom is $element\n";
        }
        return '';
    }

    return ref $atom eq 'HASH' &&
           exists $atom->{symbol} &&
           $atom->{symbol} eq $element;
}

1;
