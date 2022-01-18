package ChemOnomatopist;

use strict;
use warnings;

# ABSTRACT: Give molecule a name
# VERSION

use ChemOnomatopist::Group;
use ChemOnomatopist::Group::Carbonyl;
use ChemOnomatopist::Group::Carboxyl;
use ChemOnomatopist::Group::Hydroxy;
use Chemistry::OpenSMILES::Writer qw( write_SMILES );
use Graph::Nauty qw( canonical_order );
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

my @numbers = ( '?', '', 'di', 'tri', 'tetra', 'penta',
                'sexta', 'hepta', 'octa', 'nona', 'deca',
                'undeca', 'dodeca', 'trideca', 'tetradeca', 'pentadeca',
                'hexadeca', 'octadeca' );

my %lengths;

sub get_name
{
    my( $what ) = @_;

    # Detect the type of the input data
    my( $graph );
    if( blessed $what && $what->isa( Graph::Undirected:: ) ) {
        $graph = $what;
    } else {
        # Assume SMILES string
        require Chemistry::OpenSMILES::Parser;
        my $parser = Chemistry::OpenSMILES::Parser->new;
        ( $graph ) = $parser->parse( $what ); # Taking only the first graph
    }
    die "nothing supplied for get_name()\n" unless $graph;

    # Check if graph is a tree as trees are easy to process
    if( scalar $graph->edges != $graph->vertices - 1 ) {
        # If it is not a tree, than the graph has cycles, and we have to
        # do our best to recognise them. To make it easier, hydrogen atoms
        # are removed here for now.
        $graph->delete_vertices( grep { $_->{symbol} eq 'H' } $graph->vertices );
        my $smiles = canonical_SMILES( $graph );
        while( $smiles =~ s/\(([^\()]+)\)\)/$1)/ ) {}; # need to simplify SMILES
        if( $smiles =~ /^C1\((C+)1\)$/ ) {
            # Cycloalkane detected
            return 'cyclo' . $prefixes[scalar $graph->vertices] . 'ane';
        } elsif( $smiles =~ /^c:1\(:c((:c)+):1\)$/ &&
                 ( length( $1 ) / 2 ) =~ /^(4|6|8|10|12|14|16)$/ ) {
            # Annulene detected
            return 'cyclo' . $numbers[scalar $graph->vertices] .
                   $numbers[scalar $graph->vertices / 2] . 'ene';
        }
        # No other types of graphs with cycles can be processed for now
        die "cannot handle graphs with cycles for now\n";
    }

    create_structure($graph);

    # Traverse the graph using breadth-first traversal and pick one of
    # the furthest vertices as a starting point for naming
    my @order = BFS_order_carbons_only($graph);

    return get_chain( $graph->copy,
                      pop @order,
                      { choose_direction => 1 } ) . 'ane';
}

sub get_chain
{
    my( $graph, $start, $options ) = @_;

    $options = {} unless $options;

    # As per https://www.geeksforgeeks.org/longest-path-undirected-tree/,
    # two BFSes are needed to find the longest path in a tree

    my @order = BFS_order_carbons_only($graph, $start);

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
            next if !exists $order{$neighbour};
            next if $order{$neighbour} >= $order{$min};
            $min = $neighbour;
        }
        last if $min eq $end;
        $end = $min;
    }
    @chain = reverse @chain;
    $graph->delete_path( @chain );

    # Establishing a stable direction similarly as suggested by the IUPAC
    # rules: choosing a direction which has more attachments to its
    # beginning than to its end
    if( $options->{choose_direction} ) {
        for my $i (0..int(@chain/2)-1) {
            my $forward_chain_degree =
                    scalar $graph->neighbours( $chain[$i] ) -
                    grep { is_element( $_, 'H' ) }
                    $graph->neighbours( $chain[$i] );
            my $backward_chain_degree =
                    scalar $graph->neighbours( $chain[$#chain-$i] ) -
                    grep { is_element( $_, 'H' ) }
                    $graph->neighbours( $chain[$#chain-$i] );

            if( $forward_chain_degree != $backward_chain_degree ) {
                if( $forward_chain_degree < $backward_chain_degree ) {
                    @chain = reverse @chain;
                }
                last;
            }
        }
    }

    # Examine the attachments to the main chain: delete the edges
    # connecting them to the main chain, at the same time giving them
    # names according to their lengths via calls to get_chain()
    my %attachments;
    my $attachment_name;
    for my $i (0..$#chain) {
        my $atom = $chain[$i];
        for my $neighbour ($graph->neighbours( $atom )) {
            $graph->delete_edge( $atom, $neighbour );
            unless (is_element( $neighbour, 'H' )){
                $attachment_name = get_chain( $graph, $neighbour );
                my $prefix = ($attachment_name =~ /^\(/) ? 'yl)' : 'yl';
                push @{$attachments{$attachment_name . $prefix}}, $i;
             }
        }
    }

    # Collecting names of all the attachments
    my $name = '';
    for my $attachment_name (sort { $a cmp $b } keys %attachments) {
        $name = $name ? $name . '-' : $name;
        $name .= join( ',', map { $_ + 1 } @{$attachments{$attachment_name}} )
                 . '-' . $numbers[scalar @{$attachments{$attachment_name}}] .
                 $attachment_name;
    }
    my $bracket =
        ($options->{choose_direction} || not ($name =~ /^[0-9]/)) ? '' : '(';

    return $bracket . $name . $prefixes[scalar @chain];
}

# FIXME: not used in the main code yet
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
                 grep { is_element( $_, 'H' ) == 1 } @neighbours ){
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

        # Detecting carboxyl
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

# Original source:
# URL: svn+ssh://saulius-grazulis.lt/home/andrius/svn-repositories/cod-smi-generation/trunk/bin/isomorphism.pl
# Relative URL: ^/trunk/bin/isomorphism.pl
# Repository Root: svn+ssh://saulius-grazulis.lt/home/andrius/svn-repositories/cod-smi-generation
# Repository UUID: 389e3913-ab09-4cbd-b281-3fb5d633c480
# Revision: 570
sub canonical_SMILES
{
    my( $graph, $color_sub ) = @_;

    my @order = canonical_order( $graph, $color_sub );
    my %order;
    for (0..$#order) {
        $order{$order[$_]} = $_;
    }

    # Ignoring most likely harmless warning emitted by Set::Object
    local $SIG{__WARN__} = sub {
        return if $_[0] =~ /^Reference found where even-sized list expected at \S+\/Set\/Object\.pm line [0-9]+\.\n$/;
        print STDERR $_[0];
    };

    my $smiles = write_SMILES(
        $graph,
        sub {
            my @sorted = sort { $order{$a} <=> $order{$b} }
                              keys %{$_[0]};
            return $_[0]->{shift @sorted};
        } );

    # A.M.: I cannot find a counter-example, thus the following seems
    # reasonable to me. In a SMILES descriptor, one can substitute all
    # '/' with '\' and vice versa, and retain correct cis/trans settings.
    if( $smiles =~ /([\/\\])/ && $1 eq '\\' ) {
        $smiles =~ tr/\/\\/\\\//;
    }

    return $smiles;
}

# Check if an object or Perl hash is of certain chemical element
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

sub BFS_order_carbons_only
{
    my $is_any_visited;
    %lengths = ();

    my ( $graph, $start ) = @_;
    my $carbon_graph = $graph->copy;

    $carbon_graph->delete_vertices( grep {!is_element( $_, 'C') } $carbon_graph->vertices );

    if ($start) {
        my $bfs = Graph::Traversal::BFS->new( $carbon_graph,
            pre => sub { if(!$is_any_visited) {
                            $lengths{$_[0]->{number}} = 0; $is_any_visited = 1  } },
            tree_edge => sub { if ( !defined $lengths{$_[0]->{number}} ) {
                ( $_[0], $_[1] ) = ( $_[1], $_[0] );
            }
            $lengths{$_[1]->{number}} = $lengths{$_[0]->{number}} + 1}, start => $start);
        my @order = $bfs->bfs;
        return @order;
    }
    else{
        my $bfs = Graph::Traversal::BFS->new( $carbon_graph);
        my @order = $bfs->bfs;
        return @order;
    }
}

sub create_structure
{
    my %tree;

    my ( $graph ) = @_;
    my @order = BFS_order_carbons_only($graph);

    my $start = pop @order;

    my @second_order = BFS_order_carbons_only($graph, $start);

    my $end = pop @second_order;

    my @farthest = grep { $lengths{$_} eq $lengths{$end->{number}} } keys %lengths;
    push(@farthest, $start->{number});

    for (my $i = 0; $i < scalar(@farthest); $i++) {
        %tree = ();
        my @value_array = 0;
        $tree{$farthest[$i]} = [@value_array];

        my $carbon_graph = $graph->copy;
        $carbon_graph->delete_vertices( grep {!is_element( $_, 'C') } $carbon_graph->vertices );

        my @vertice = grep { $_->{number} eq $farthest[$i] } $carbon_graph->vertices;

        create_tree($carbon_graph, $vertice[0], \%tree);
    }
}

sub create_tree
{
    my ( $graph, $atom, $tree ) = @_;

    my @neighbours = $graph->neighbours( $atom );

    my @array = @ { $tree->{$atom->{number}}};
    my @new_array;

    foreach my $arr (@array){
        push(@new_array, $arr + 1);
    }

    if (scalar @neighbours == 1) {
        unless (exists $tree->{$neighbours[0]->{number}}){
            $tree->{$neighbours[0]->{number}} = [@new_array];
            $graph->delete_vertex($atom);
            create_tree($graph, $neighbours[0], $tree)
        }
    }
    elsif(scalar @neighbours > 1){
        push(@new_array, 0);
        $graph->delete_vertex($atom);
        foreach my $neighbour (@neighbours) {
            $tree->{$neighbour->{number}} = [@new_array];
            create_tree($graph, $neighbour, $tree);
        }
    }
}

1;
