package ChemOnomatopist;

use strict;
use warnings;

# ABSTRACT: Give molecule a name
# VERSION

use ChemOnomatopist::Group;
use ChemOnomatopist::Group::Carbonyl;
use ChemOnomatopist::Group::Carboxyl;
use ChemOnomatopist::Group::Hydroxy;
use Graph::Traversal::BFS;
use Graph::Undirected;
use Scalar::Util qw(blessed);

our @prefixes = qw(
    ?
    meth
    eth
    prop
    but
    pent
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

my @numbers = ( '?', '', 'di', 'tri', 'tetra',
                'penta', 'sexta', 'hepta', 'octa' );

sub get_name
{
    my( $what ) = @_;

    my( $graph );
    if( blessed $what && $what->isa( Graph::Undirected:: ) ) {
        $graph = $what;
    } else {
        # Assume SMILES string
        require Chemistry::OpenSMILES::Parser;
        my $parser = Chemistry::OpenSMILES::Parser->new;
        ( $graph ) = $parser->parse( $what, { raw => 1 } ); # Taking only the first graph
    }
    die "nothing supplied for get_name()\n" unless $graph;

    if( scalar $graph->edges != $graph->vertices - 1 ) {
        die "cannot handle graphs with cycles for now\n";
    }

    my $last;
    my $bfs = Graph::Traversal::BFS->new( $graph,
                                          # next_root => undef,
                                        );
    my @order = $bfs->bfs;
    return get_chain( $graph->copy,
                      pop @order,
                      { choose_direction => 1 } ) . 'ane';
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

    my %attachments;
    for my $i (0..$#chain) {
        my $atom = $chain[$i];
        for my $neighbour ($graph->neighbours( $atom )) {
            $graph->delete_edge( $atom, $neighbour );
            push @{$attachments{get_chain( $graph, $neighbour ) . 'yl'}},
                 $i;
        }
    }

    my $name = '';
    for my $attachment_name (sort { $numbers[scalar @{$attachments{$a}}] . $a cmp
                                    $numbers[scalar @{$attachments{$b}}] . $b }
                             keys %attachments) {
        $name = $name ? $name . '-' : $name;
        $name .= join( ',', map { $_ + 1 } @{$attachments{$attachment_name}} ) . '-' .
                 $numbers[scalar @{$attachments{$attachment_name}}] .
                 $attachment_name;
    }

    return $name . $prefixes[scalar @chain];
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
            my $carboxyl = ChemOnomatopist::Group::Carboxyl->new( $carbonyl );
            for ($graph->neighbours( $carbonyl )) {
                $graph->add_edge( $_, $carboxyl );
                $graph->delete_edge( $_, $carbonyl );
            }
            $graph->delete_vertex( $carboxyl );
            $graph->delete_vertex( $hydroxy );
        }
    }
}

sub remove_pendant_vertices
{
    my( $graph ) = @_;

    while( my @pendants = grep { $graph->degree( $_ ) == 1 }
                               $graph->vertices ) {
        $graph->delete_vertices( @pendants );
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
        } elsif( $atom->isa( 'Chemistry::Atom' ) ) {
            return $atom->symbol eq $element;
        }
        return '';
    }

    return ref $atom eq 'HASH' &&
           exists $atom->{symbol} &&
           $atom->{symbol} eq $element;
}

1;
