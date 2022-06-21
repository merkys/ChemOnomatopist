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
use Clone qw( clone );
use Graph::Nauty qw( canonical_order );
use Graph::Traversal::BFS;
use Graph::Undirected;
use List::Util qw( any );
use Scalar::Util qw( blessed );

no warnings 'recursion';

my @numbers = ( '?', '', 'di', 'tri', 'tetra', 'penta',
                'hexa', 'hepta', 'octa', 'nona', 'deca',
                'undeca', 'dodeca', 'trideca', 'tetradeca',
                'pentadeca', 'hexadeca', 'heptadeca', 'octadeca', 'nonadeca',
                'icosa', 'henicosa', 'docosa', 'tricosa',
                'tetracosa', 'pentacosa', 'hexacosa',
                'heptacosa', 'octacosa', 'nonacosa', 'triaconta',
                'hentriaconta', 'dotriaconta', 'tritriaconta', 'tetratriaconta',
                'pentatriaconta', 'hexatriaconta', 'heptatriaconta',
                'octatriaconta', 'nonatriaconta', 'tetraconta' );

my @numberskis = ( '?', '', 'bis', 'tris', 'tetrakis', 'pentakis',
                'hexakis', 'heptakis', 'octakis', 'nonakis', 'decakis' );

my %preferrable_names = ( 
                '(1-methylethyl)' => 'propan-2-yl',
                '(1-ethyl-1-methylpropyl)' => '(3-methylpentan-3-yl)',
                '(1,1-dimethylethyl)' => 'tert-butyl',
                '(1,3,3-trimethylbutyl)' => '(4,4-dimethylpentan-2-yl)',
                '(1-methylpropyl)' => 'butan-2-yl',
                '(1-ethylbutyl)' => 'hexan-3-yl',
                '(1-butylheptyl)' => 'undecan-5-yl',
                '(1-methyloctadecyl)' => 'nonadecan-2-yl',
                '(1-methylbutyl)' => 'pentan-2-yl',
                '(1-propylbutyl)' => 'heptan-4-yl',
                '(1-ethylpropyl)' =>  'pentan-3-yl',
                '(1-methylheptyl)' => 'octan-2-yl',
                '(1-methylpentyl)' => 'hexan-2-yl',
                '(1,2-dimethylpropyl)' => '(3-methylbutan-2-yl)',
                '(1,2-dimethylbutyl)' => '(3-methylpentan-2-yl)',
                '(1,1-dimethyldecyl)' => '(2-methylundecan-2-yl)',
                '(1,1-dimethylpentyl)' => '(2-methylhexan-2-yl)',
                '(1,1-dimethylbutyl)' => '(2-methylpentan-2-yl)',
                '(1-propylpentyl)' => 'octan-4-yl',
                '(1-ethyl-2-methylpropyl)' => '(2-methylpentan-3-yl)',
                '(1-butylhexyl)' => 'decan-5-yl',
                '(1-butylpentyl)' => 'nonan-5-yl',
                '(1,1-dimethylpropyl)' => '(2-methylbutan-2-yl)',
                '(1-(1-methylethyl)-2-methylpropyl)' => '(2,4-dimethylpentan-3-yl)',
                '(1-ethylpentyl)' => 'heptan-3-yl' );

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
            return 'cyclo' . alkane_chain_name( scalar $graph->vertices ) . 'ane';
        } elsif( $smiles =~ /^c:1\(:c((:c)+):1\)$/ &&
                 ( length( $1 ) / 2 ) =~ /^(4|6|8|10|12|14|16)$/ ) {
            # Annulene detected
            return 'cyclo' . $numbers[scalar $graph->vertices] .
                   $numbers[scalar $graph->vertices / 2] . 'ene';
        }
        # No other types of graphs with cycles can be processed for now
        die "cannot handle graphs with cycles for now\n";
    }
    my( $custom_order, $order ) = select_main_chain( $graph->copy );

    if( $custom_order ) {
        return get_chain_2( $graph->copy,
                            $order,
                            { choose_direction => 1 } ) . 'ane';
    } else {
        # Traverse the graph using breadth-first traversal and pick one of
        # the furthest vertices as a starting point for naming
        my @order = BFS_order_carbons_only($graph);

        return get_chain( $graph->copy,
                          pop @order,
                          { choose_direction => 1 } ) . 'ane';
    }
}

sub get_chain
{
    my( $graph, $start, $options ) = @_;

    $options = {} unless $options;

    # As per https://www.geeksforgeeks.org/longest-path-undirected-tree/,
    # two BFSes are needed to find the longest path in a tree
    my @order = BFS_order_carbons_only( $graph, $start );

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
        my $number = IUPAC_numerical_multiplier( scalar @{$attachments{$attachment_name}} );
        $number =~ s/^mono$//;
        $number .= 'a' unless $number =~ /^(|\?|.*i)$/;
        $name = $name ? $name . '-' : $name;
        $name .= join( ',', map { $_ + 1 } @{$attachments{$attachment_name}} ) .
                 '-' . $number . $attachment_name;
    }
    my $bracket =
        ($options->{choose_direction} || $name !~ /^[0-9]/) ? '' : '(';

    return $bracket . $name . alkane_chain_name( scalar @chain );
}

sub get_chain_2
{
    my( $graph, $main_chain, $options ) = @_;

    my @vertices = $graph->vertices;
    my @chain;

    # Recreate main chain order by the array in $main_chain
    for my $curr_vertex (@$main_chain) {
        my( $vertex ) = grep { $_->{number} == $curr_vertex } @vertices;
        push @chain, $vertex;
    }

    $graph->delete_path( @chain );

    # Examine the attachments to the main chain: delete the edges
    # connecting them to the main chain, at the same time giving them
    # names according to their lengths via calls to get_chain()
    my %attachments;
    my $attachment_name;
    for my $i (0..$#chain) {
        my $atom = $chain[$i];
        for my $neighbour ($graph->neighbours( $atom )) {
            $graph->delete_edge( $atom, $neighbour );
            unless (is_element( $neighbour, 'H' )) {
                $attachment_name = get_chain( $graph, $neighbour );
                my $prefix = ($attachment_name =~ /^\(/) ? 'yl)' : 'yl';
                push @{$attachments{$attachment_name . $prefix}}, $i;
             }
        }
    }

    # Replacing systematic IUPAC attachment names with their preferred ones
    for my $att_name (keys %attachments) {
        next unless exists $preferrable_names{$att_name};
        $attachments{$preferrable_names{$att_name}} = $attachments{$att_name};
        delete $attachments{$att_name};
        next unless $preferrable_names{$att_name} =~ /^[a-z0-9\-]+$/;

        # If there is more than one of the same attachment
        # so complete name would start with the prefix, brackets
        # should be added to the name
        if( scalar @{$attachments{$preferrable_names{$att_name}}} > 1 ) {
            $attachments{'(' . $preferrable_names{$att_name} . ')'} =
                            $attachments{$preferrable_names{$att_name}};
            delete $attachments{$preferrable_names{$att_name}};
        }
    }

    # Collecting names of all the attachments
    my $name = '';
    for my $attachment_name (sort compare_only_aphabetical keys %attachments) {
        $name = $name ? $name . '-' : $name;
        my $number;
        if( $attachment_name =~ /^\([0-9]/ ) {
            $number = $numberskis[scalar @{$attachments{$attachment_name}}];
        } else {
            $number = IUPAC_numerical_multiplier( scalar @{$attachments{$attachment_name}} );
            $number =~ s/^mono$//;
            $number .= 'a' unless $number =~ /^(|\?|.*i)$/;
        }
        $name .= join( ',', map { $_ + 1 } @{$attachments{$attachment_name}} ) .
                 '-' . $number . $attachment_name;
    }

    my $bracket =
        ($options->{choose_direction} || $name !~ /^[0-9]/) ? '' : '(';
    return $bracket . $name . alkane_chain_name( scalar @chain );
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
                 any { is_element( $_, 'H' ) == 1 } @neighbours ) {
            my $hydroxy  = ChemOnomatopist::Group::Hydroxy->new( $atom );
            my $hydrogen = any { is_element( $_, 'H' ) } @neighbours;
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

    # In a SMILES descriptor, one can substitute all '/' with '\'
    # and vice versa, retaining correct cis/trans settings.
    # Similar rule is explained in O'Boyle, 2012, Rule H.
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

# Subroutine gets an graph, removes all vertices that do not have C as their element.
# Performs BFS on that chain. During BFS, distance from start is calculated to each vertice
sub BFS_order_carbons_only_return_lengths
{
    my $is_any_visited;
    my %lengths;

    my ( $graph, $start ) = @_;
    my $carbon_graph = $graph->copy;

    $carbon_graph->delete_vertices( grep { !is_element( $_, 'C') } $carbon_graph->vertices );
    my $bfs = Graph::Traversal::BFS->new(
        $carbon_graph,
        pre => sub { if( !$is_any_visited ) {
                         $lengths{$_[0]->{number}} = 0;
                         $is_any_visited = 1;
                     } },
        tree_edge => sub { if( !defined $lengths{$_[0]->{number}} ) {
                               ( $_[0], $_[1] ) = ( $_[1], $_[0] );
                           }
                           $lengths{$_[1]->{number}} = $lengths{$_[0]->{number}} + 1;
                         },
        start => $start );

    return \%lengths, $bfs->bfs;
}

# BFS is performed for the given graph after all vertices that are not carbons
# removed
sub BFS_order_carbons_only
{
    my( $graph, $start ) = @_;
    my $carbon_graph = $graph->copy;
    $carbon_graph->delete_vertices( grep { !is_element( $_, 'C' ) } $carbon_graph->vertices );

    my $bfs;
    if( $start ) {
        $bfs = Graph::Traversal::BFS->new( $carbon_graph, start => $start );
    } else {
        $bfs = Graph::Traversal::BFS->new( $carbon_graph );
    }
    return $bfs->bfs;
}

# Calculates length of given graph (vertices count)
sub BFS_calculate_chain_length
{
    my( $graph, $start ) = @_;
    my $bfs = Graph::Traversal::BFS->new( $graph, start => $start );
    return scalar $bfs->bfs;
}

# Returns 1 if there is any braches in the given graph and if there is none
sub BFS_is_chain_branched
{
    my( $graph, $start ) = @_;

    # FIXME: Not entirely sure why visited vertices are removed (A.M.)
    my $graph_copy = $graph->copy;

    my $branched = 0;
    my $bfs = Graph::Traversal::BFS->new(
        $graph,
        pre => sub {
            my @neighbours = $graph_copy->neighbours( $_[0] );
            $branched = 1 if scalar @neighbours > 1;
            $graph_copy->delete_vertex( $_[0] );
        },
        start => $start );
    $bfs->bfs;

    return $branched;
}

# Returns main (parental) chain to be used during the naming
sub select_main_chain
{
    my( $graph ) = @_;
    my @order = BFS_order_carbons_only( $graph );

    my $start = $order[-1];

    my( $lengths, @second_order ) =
        BFS_order_carbons_only_return_lengths( $graph, $start );

    my $end = $second_order[-1];

    # Finding all farthest vertices from the starting point
    my @farthest = grep { $lengths->{$_} eq $lengths->{$end->{number}} }
                        keys %$lengths;

    # Also adding the first vertice to the array since it is farthest from other
    # ones
    push @farthest, $start->{number};

    # Going through every vertice in "farthest" array and creating tree-like structures
    my @all_trees;
    for (my $i = 0; $i < scalar @farthest; $i++) {
        my %tree = ( $farthest[$i] => [ $farthest[$i], 0 ] );

        my $carbon_graph = $graph->copy;
        $carbon_graph->delete_vertices( grep { !is_element( $_, 'C' ) } $carbon_graph->vertices );

        my( $vertice ) = grep { $_->{number} eq $farthest[$i] } $carbon_graph->vertices;
        # Start creation of the tree from all the starting vertices
        push @all_trees, \%{create_tree( $carbon_graph, $vertice, \%tree )};
    }

    my $trees;

    # Extracts arrays of all longest chains with numbers of vertices (in order)
    # from each tree-like structure
    my @main_chains = prepare_paths( @all_trees );

    my $carbon_graph = $graph->copy;
    $carbon_graph->delete_vertices(
        grep { !is_element( $_, 'C' ) } $carbon_graph->vertices
    );
    # From all possible main chains in tree-like structures, subroutine returns
    # the ones that has the greatest number of side chains. Also returns only the
    # trees that still have some possible main chains after selection
    ( $trees, @main_chains ) = rule_greatest_number_of_side_chains(
                                    $carbon_graph,
                                    \@main_chains,
                                    @all_trees
                               );
    return ( 1, @{$main_chains[0]} ) if scalar @{$main_chains[0]} == 1;

    # If more than one chain is left, second rule is applied
    $carbon_graph = $graph->copy;
    $carbon_graph->delete_vertices(
        grep { !is_element( $_, 'C' ) } $carbon_graph->vertices
    );
    # From all main chains left, subroutine selects all that have the
    # lowest numbered locants. Also the trees that have possible main
    # chains returned
    ( $trees, @main_chains ) = rule_lowest_numbered_locants(
                                $carbon_graph, @main_chains, @$trees
                               );

    return ( 1, @{$main_chains[0]} ) if scalar @{$main_chains[0]} == 1;

    # If more than one chain is left, third rule is applied
    $carbon_graph = $graph->copy;
    $carbon_graph->delete_vertices(
        grep { !is_element( $_, 'C' ) } $carbon_graph->vertices
    );
    # From all main chains left, subroutine selects all that have the
    # most carbons in side chains. Also the trees that have possible main
    # chains returned
    ( $trees, @main_chains ) = rule_most_carbon_in_side_chains(
                                $carbon_graph, @main_chains, @$trees
                               );

    return ( 1, @{$main_chains[0]} ) if scalar @{$main_chains[0]} == 1;

    # If more than one chain is left, fourth rule is applied
    $carbon_graph = $graph->copy;
    $carbon_graph->delete_vertices(
        grep { !is_element( $_, 'C' ) } $carbon_graph->vertices
    );
    # From all main chains left, subroutine selects all that have
    # the least branched side chains. Also the trees that have
    # possible main chains returned
    ( $trees, @main_chains ) = rule_least_branched_side_chains(
                                $carbon_graph, @main_chains, @$trees
                               );

    return ( 1, @{$main_chains[0]} ) if scalar @{$main_chains[0]} == 1;

    # If more than one chain is left, program states that there are
    # few chains that are identical by all the rules and selects
    # one from all that are left
    $carbon_graph = $graph->copy;
    $carbon_graph->delete_vertices(
        grep { !is_element( $_, 'C' ) } $carbon_graph->vertices
    );
    # One main chain is picked from all that are left
    my $main_chain = rule_pick_chain_from_valid(
                        $carbon_graph,
                        @main_chains,
                        @$trees
                     );
    return ( 1, $main_chain );
}

# Creating tree like structure for all the longest paths in molecule
sub create_tree
{
    my( $graph, $atom, $tree ) = @_;

    my @neighbours = $graph->neighbours( $atom );

    my @array = @{$tree->{$atom->{number}}};
    my @new_array = ($atom->{number});

    # Since first number in the stack-of-array boxes represents the parent of
    # the vertice, array with box information is shifted
    shift @array;

    # Box for the vertice is created by increasing box information in parental
    # box by one
    push @new_array, map { $_ + 1 } @array;

    # If there is one neighbour, it means that vertice do not have any branching.
    # Analysis of next vertice (the neighbour) is started
    if( scalar @neighbours == 1 ) {
        unless (exists $tree->{$neighbours[0]->{number}}) {
            $tree->{$neighbours[0]->{number}} = [ @new_array ];
            $graph->delete_vertex( $atom );
            create_tree( $graph, $neighbours[0], $tree );
        }
    }
    # If there is more than one neighour for the current vertice, analysis of
    # each of them is started independently
    elsif( scalar @neighbours > 1 ) {
        push @new_array, 0;
        $graph->delete_vertex( $atom );
        foreach my $neighbour ( @neighbours ) {
            $tree->{$neighbour->{number}} = [ @new_array ];
            create_tree( $graph, $neighbour, $tree );
        }
    }
    # If there is no neighbour or other vertices, all graph has been analyzed.
    # Created structure is returned
    return $tree;
}

# Create arrays of vertice numbers of all possible longest chains in the tree
sub prepare_paths
{
    my( @trees ) = @_;

    my $trees_copy = clone \@trees;
    my @all_chains;
    foreach my $tree ( @$trees_copy ) {
        my %structure = %$tree;

        # Tree-like structure is sorted by parental vertex number
        my @sorted = sort { $structure{$a}->[1] <=> $structure{$b}->[1] }
                          keys %structure;

        my $last = $sorted[-1];
        my @chain_ending = grep { $structure{$_}->[1] == $structure{$last}->[1] }
                                keys %structure;

        foreach my $ending ( @chain_ending ) {
            my @vertex_array = ($ending);
            push @all_chains, save_main_chain_vertices_in_array(
                                $ending,
                                \@vertex_array,
                                $tree
                              );
        }
    }

    # Adds reverted chains if they are not present yet as the longest chains
    for my $chain ( @all_chains ) {
        next if array_exists( [reverse @$chain], @all_chains );
        push @all_chains, [reverse @$chain];
    }

    return @all_chains;
}

# Checks if array exists in array of arrays
sub array_exists
{
    my( $chain, @all_chains ) = @_;

    my $same;
    for my $curr_chain ( @all_chains ) {
        $same = 1;
        for my $index ( 0..scalar @{$curr_chain}-1 ) {
            if( $chain->[$index] != $curr_chain->[$index] ) {
                $same = 0;
                last;
            }
        }
        return 1 if $same;
    }
    return 0;
}

# Tries to find the chain which has the greatest number of side chains
sub rule_greatest_number_of_side_chains
{
    my( $graph, $chains, @trees ) = @_;

    my $trees_copy = clone \@trees;
    my $index = 0;
    my @number_of_side_chains;

    foreach my $tree ( @$trees_copy ) {
        my %structure = %{clone $tree};

        # Reference to parental chain is removed from the boxes
        foreach my $key ( keys %structure ) {
            shift @{$structure{$key}};
        }

        # Beginning of the structure is found. Then all chains that belongs to
        # the current tree are selected
        my @first = grep { $structure{$_}->[0] == 0 } keys %structure;
        my @chains_in_the_tree =
            grep { $_->[0] == $first[0] || $_->[-1] == $first[0] } @$chains;

        # Structure with index of the tree, beginning and ending of the chain,
        # number of side chains in the chain created for each chain
        for my $chain ( @chains_in_the_tree ) {
            push @number_of_side_chains,
                 [$index, $chain->[0], $chain->[-1],
                    find_number_of_side_chains(
                        $graph,
                        \@{$chain},
                        $tree
                    )
                 ];
        }
        $index++;
    }

    # All chains that have the biggest number of side chains are selected and returned
    my @sorted_numbers = sort { $a->[3] <=> $b->[3] } @number_of_side_chains;

    my $path_length = $sorted_numbers[-1][3];
    my @biggest_number_of_side_chains = grep {$_->[3] == $path_length} @number_of_side_chains;
    my %seen;
    my @uniq_biggest_number_of_side_chains = grep { !$seen{$_->[0]}++ } @biggest_number_of_side_chains;
    my @result = @trees[map {$_->[0]} @uniq_biggest_number_of_side_chains];

    my @eligible_chains;
    for my $chain (@$chains) {
        if( any {$_->[1] == $chain->[0] && $_->[2] == $chain->[-1]} @biggest_number_of_side_chains ) {
            push @eligible_chains, $chain;
        }
    }
    return \@result, \@eligible_chains;
}

# Tries to find the chain which has the lowest-numbered locants
sub rule_lowest_numbered_locants
{
    my ( $graph, $chains, @trees ) = @_;

    my $trees_copy = clone \@trees;
    my $index = 0;
    my @locant_placing;

    foreach my $tree ( @$trees_copy ) {
        my %structure = %{clone $tree};

        # Reference to parental chain is removed from the boxes
        foreach my $key ( keys %structure ) {
            shift @{$structure{$key}};
        }

        # Beginning of the structure is found. Then all chains that belongs to
        # the current tree are selected
        my @first = grep { $structure{$_}->[0] == 0 } keys %structure;
        my @chains_in_the_tree = grep { @{$_}[0] == $first[0] || @{$_}[-1] == $first[0] } @$chains;

        # Structure with index of the tree, beginning and ending of the chain,
        # places of the locants in the chain created for each tree
        for my $chain ( @chains_in_the_tree ) {
            push @locant_placing,
                 [$index, $chain->[0], $chain->[-1],
                    [find_locant_placing(
                        $graph,
                        $chain,
                        $tree
                    )]
                 ];
        }
        $index++;
    }

    # All chains that have the lowest numbers of locants are selected and returned
    my @sorted_paths = sort compare_locant_placings reverse @locant_placing;

    my $lowest_locants = $sorted_paths[0][3];
    my @lowest_locants_paths =
        grep { join( '', @{$_->[3]} ) eq join( '', @$lowest_locants ) }
            @locant_placing;
    my %seen;
    my @uniq_lowest_locants_paths = grep { !$seen{$_->[0]}++ } @lowest_locants_paths;
    my @result = @trees[map {$_->[0]} @uniq_lowest_locants_paths];

    my @eligible_chains;
    for my $chain (@$chains) {
        if( any { $_->[1] == $chain->[0] and $_->[2] == $chain->[-1] } @lowest_locants_paths) {
            push @eligible_chains, $chain;
        }
    }
    return \@result, \@eligible_chains;
}

# Tries to find chain that has the greatest number of carbon atoms in the smaller
# side chains
sub rule_most_carbon_in_side_chains
{
    my ( $graph, $chains, @trees ) = @_;

    my $trees_copy = clone \@trees;
    my @side_chain_lengths;
    my $index = 0;

    foreach my $tree (@$trees_copy) {
        my %structure = %{clone $tree};
        my @all_vertices = keys %structure;

        # Reference to parental chain is removed from the boxes
        foreach my $key (keys %structure) {
            shift @{$structure{$key}};
        }

        my @sorted = sort {
                            @{$structure{$a}} <=> @{$structure{$b}}
                           or
                            $structure{$a}->[0] cmp $structure{$b}->[0]
                          } keys %structure;
        my $last = $sorted[-1];

        # Beginning of the structure is found. Then all chains that belongs to
        # the current tree are selected
        my @first = grep{ $structure{$_}->[0] == 0 } keys %structure;
        my @structure_chains = grep { @{$_}[0] == $first[0] || @{$_}[-1] == $first[0] } @$chains;

        # Structure with index of the tree, beginning and ending of the chain,
        # lengths of side chains of the chain created for each tree
        for my $chain (@structure_chains) {
            my $graph_copy = $graph->copy;
            my @side_chain_length;
            push ( @side_chain_lengths,
                [$index, $chain->[0], $chain->[-1],
                    [find_lengths_of_side_chains(
                        $graph_copy,
                        $chain->[-1],
                        \@{$chain},
                        \@side_chain_length,
                        $tree,
                        scalar @{$chain}
                    )]
                ]
            )
        }
        $index++;
    }

    # All chains that have the highest number of carbons in side chains are selected and returned
    my @sorted_final = sort compare_side_chain_lengths @side_chain_lengths;

    my $last = $sorted_final[-1][3];
    my @greatest_no_of_side_chains_paths =
        grep { join( '', @{$_->[3]} ) eq join( '', @$last ) }
            @sorted_final;
    my @eligible_chains;
    for my $chain (@{$chains}){
        if( any { $_->[1] == $chain->[0] && $_->[2] == $chain->[-1] }
                @greatest_no_of_side_chains_paths ) {
            push @eligible_chains, $chain;
        }
    }
    my %seen;
    my @uniq_side_chain_paths =
                grep { !$seen{$_->[0]}++ } @greatest_no_of_side_chains_paths;
    my @result = @trees[map { $_->[0] } @uniq_side_chain_paths];

    return \@result, \@eligible_chains;
}

# Tries to find chain that have the least branched side chains
sub rule_least_branched_side_chains
{
    my ( $graph, $chains, @trees ) = @_;

    my $trees_copy = clone \@trees;
    my $index = 0;
    my @number_of_branched_side_chains;

    foreach my $tree (@$trees_copy) {
        my %structure = %{clone $tree};

        # Reference to parental chain is removed from the boxes
        foreach my $key (keys %structure) {
            shift @{$structure{$key}};
        }

        # Beginning of the structure is found. Then all chains that belongs to
        # the current tree are selected
        my @first = grep { $structure{$_}->[0] == 0 } keys %structure;
        my @chains_in_the_tree = grep { $_->[0] == $first[0] || $_->[-1] == $first[0] } @$chains;

        for my $chain (@chains_in_the_tree) {
            push @number_of_branched_side_chains,
                 [$index, $chain->[0], $chain->[-1],
                    find_number_of_branched_side_chains(
                        $graph,
                        \@{$chain},
                        $tree
                    )
                 ];
        }
        $index++;
    }

    # All chains that have the least amount of branches side chains are selected and returned
    my @sorted_paths = sort { $a->[3] <=> $b->[3] }
                            @number_of_branched_side_chains;
    my $path_length = $sorted_paths[0][3];
    my @longest_paths = grep { $_->[3] == $path_length } @number_of_branched_side_chains;
    my %seen;
    my @uniq_longest_paths = grep { !$seen{$_->[0]}++ } @longest_paths;
    my @result = @trees[map { $_->[0] } @uniq_longest_paths];
    my @eligible_chains;
    for my $chain (@$chains) {
        if( any { $_->[1] == $chain->[0] && $_->[2] == $chain->[-1] } @longest_paths ){
            push @eligible_chains, $chain;
        }
    }
    return \@result, \@eligible_chains;
}

# Subroutine sorts all valid chains that are left and returns the first one -
# the one that have carbons with lowest indexes if there is no differences
# regarding the attachment names. If there is, then selects the ones that have
# lowest attachment indexes with lowest attachments alphabetically
sub rule_pick_chain_from_valid
{
    my( $graph, $chains, @trees ) = @_;

    my $chosen_chains =
        pick_chain_with_lowest_attachments_alphabetically( $graph, $chains, @trees );

    my( $chain ) = sort compare_arrays @$chosen_chains;
    return $chain;
}

# Subroutine selects chain that has the lowest attachments by alpabetical naming
sub pick_chain_with_lowest_attachments_alphabetically
{
    my( $graph, $chains, @trees ) = @_;

    my $trees_copy = clone \@trees;
    my $index = 0;
    my @locant_placing;
    foreach my $tree (@$trees_copy){
        my %structure = %{clone $tree};

        # Reference to parental chain is removed from the boxes
        foreach my $key (keys %structure) {
            shift @{$structure{$key}};
        }

        # Beginning of the structure is found. Then all chains that belongs to
        # the current tree are selected
        my @first = grep{ $structure{$_}->[0] == 0 } keys %structure;
        my @chains_in_the_tree = grep { $_->[0] == $first[0] || $_->[-1] == $first[0] } @$chains;

        # Placings of the locants found for each chain
        for my $chain (@chains_in_the_tree) {
            push @locant_placing,
                 [$index, $chain->[0], $chain->[-1],
                    [find_locant_placing(
                        $graph,
                        \@{$chain},
                        $tree
                    )]
                 ];
        }
        $index++;
    }

    my @vertices = $graph->vertices;

    # Naming of all attachments found
    my @attachments;
    for (my $i = 0; $i < scalar(@locant_placing); $i++) {
        my @curr_chain = grep { $_->[0] == $locant_placing[$i][1] && $_->[-1] == $locant_placing[$i][2] } @$chains;
        my @attachments_only;
        for my $locant (@{$locant_placing[$i][3]}) {
            my( $vertex ) = grep {$_->{number} == $curr_chain[0][$locant-1]} @vertices;
            my @curr_neighbours = $graph->neighbours( $vertex );
            
            for my $neighbour (@curr_neighbours) {
                next if any { $neighbour->{number} eq $_ } @{$curr_chain[0]};
                my $graph_copy = $graph->copy;
                $graph_copy->delete_edge( $vertex, $neighbour );
                my $attachment_name = get_chain( $graph_copy, $neighbour );
                my $prefix = ($attachment_name =~ /^\(/) ? 'yl)' : 'yl';
                push @attachments_only, $attachment_name . $prefix;
            }
        }
        # Replacing systematic IUPAC attachment names with their preferrable
        # ones.
        for my $att_name (@attachments_only) {
            if( exists $preferrable_names{$att_name} ) {
                $att_name = $preferrable_names{$att_name};
            }
        }
        push @attachments, [clone( $curr_chain[0] ), \@attachments_only];
    }
    my @sorted_attachments = sort sort_attachments @attachments;
    my $correct_attach = $sorted_attachments[0][1];

    # All chains that have the same - alpabetically lowest attachments selected
    my @correct_chains_all = grep { join( ',', @{$_->[1]} ) eq join( ',', @$correct_attach ) } @attachments;
    my @correct_chains = map { $_->[0] } @correct_chains_all;
    return \@correct_chains;
}

# Returns array that contains numbers of vertices that are in main chain
sub save_main_chain_vertices_in_array
{
    my( $curr_vertex, $all_vertices, $structure ) = @_;

    if( $structure->{$curr_vertex}[0] == $curr_vertex ) {
        return $all_vertices;
    } else {
        push @$all_vertices, $structure->{$curr_vertex}[0];
        return save_main_chain_vertices_in_array(
                $structure->{$curr_vertex}->[0],
                $all_vertices,
                $structure );
    }
}

# Returns array that contains lengths of all side chains
sub find_lengths_of_side_chains
{
    my( $graph, $curr_vertex, $main_chain_vertices, $side_chain_lengths, $structure, $atoms_left ) = @_;
    return sort @$side_chain_lengths unless $atoms_left;

    my @vertices = $graph->vertices;
    my( $vertex ) = grep { $_->{number} == $curr_vertex } @vertices;
    my @curr_neighbours = $graph->neighbours( $vertex );
    if( scalar @curr_neighbours == 1 ) {
        $graph->delete_vertex( $vertex );
        $atoms_left--;
        find_lengths_of_side_chains(
            $graph,
            $curr_neighbours[0]->{number},
            $main_chain_vertices,
            $side_chain_lengths,
            $structure,
            $atoms_left
        );
    } else {
        my @side_chain_neighbours;
        my $next_chain_vertex;

        # Find all neighours of the chain that does not exist in main chain and
        # the next chain to be analyzed
        foreach my $neigh ( @curr_neighbours ) {
            if( any { $neigh->{number} eq $_ } @$main_chain_vertices ) {
                $next_chain_vertex = $neigh;
            } else {
                push @side_chain_neighbours, $neigh;
            }
        }
        $graph->delete_vertex( $vertex );
        $atoms_left--;

        # For each side chain neighbour, find their chain lengths
        foreach my $neighbour ( @side_chain_neighbours ) {
            push @$side_chain_lengths,
                 BFS_calculate_chain_length( $graph, $neighbour );
        }

        find_lengths_of_side_chains(
            $graph,
            $next_chain_vertex->{number},
            $main_chain_vertices,
            $side_chain_lengths,
            $structure,
            $atoms_left
        );
    }
}

# TODO: Try to merge all subroutines

# Find placings of all locants in the chain
sub find_locant_placing
{
    my( $graph, $main_chain, $structure ) = @_;

    # Code is destructive, need to make a copy before execution:
    $graph = $graph->copy;

    my @vertices = $graph->vertices;
    my @places_of_locants;
    my $vertex_number = scalar @$main_chain;

    for my $curr_vertex ( reverse @$main_chain ) {
        my( $vertex ) = grep { $_->{number} == $curr_vertex } @vertices;
        my @curr_neighbours = $graph->neighbours( $vertex );
        return @places_of_locants unless scalar @curr_neighbours;

        $graph->delete_vertex( $vertex );
        if( scalar @curr_neighbours > 1 ) {
            foreach my $neigh (@curr_neighbours) {
                next if any { $neigh->{number} eq $_ } @$main_chain;
                push @places_of_locants, $vertex_number;
            }
        }
        $vertex_number--;
    }
}

# Returns number of side chains
sub find_number_of_side_chains
{
    my( $graph, $main_chain, $structure ) = @_;

    # Code is destructive, need to make a copy before execution:
    $graph = $graph->copy;

    my @vertices = $graph->vertices;
    my $number_of_side_chains = 0;

    for my $curr_vertex ( reverse @$main_chain ) {
        my( $vertex ) = grep { $_->{number} == $curr_vertex } @vertices;
        my @curr_neighbours = $graph->neighbours( $vertex );
        return $number_of_side_chains unless scalar @curr_neighbours;

        $graph->delete_vertex( $vertex );
        if( scalar @curr_neighbours > 1 ) {
            foreach my $neigh (@curr_neighbours) {
                next if any { $neigh->{number} eq $_ } @$main_chain;
                $number_of_side_chains++;
            }
        }
    }
}

sub find_number_of_branched_side_chains
{
    my( $graph, $main_chain, $structure ) = @_;

    # Code is destructive, need to make a copy before execution:
    $graph = $graph->copy;

    my @vertices = $graph->vertices;
    my $number_of_branched_side_chains = 0;

    for my $curr_vertex ( reverse @$main_chain ) {
        my( $vertex ) = grep {$_->{number} == $curr_vertex} @vertices;
        my @curr_neighbours = $graph->neighbours( $vertex );
        return $number_of_branched_side_chains unless scalar @curr_neighbours;

        $graph->delete_vertex( $vertex );
        if( scalar @curr_neighbours > 1 ) {
            foreach my $neigh (@curr_neighbours) {
                next if any { $neigh->{number} eq $_ } @$main_chain;
                $number_of_branched_side_chains +=
                    BFS_is_chain_branched( $graph, $neigh );
            }
        }
    }
}

# Finds center (or two centers) of a tree graph
sub graph_center
{
    my( $graph ) = @_;

    $graph = $graph->copy;
    my $nvertices = scalar $graph->vertices;
    while( $graph->vertices > 2 ) {
        $graph->delete_vertices( grep { $graph->degree( $_ ) == 1 }
                                      $graph->vertices );
        my $nvertices_now = scalar $graph->vertices;
        if( $nvertices_now == $nvertices ) {
            # Safeguard for cycles and/or isolated vertices
            die 'cannot find center of cyclic or isolated graphs';
        }
        $nvertices = $nvertices_now;
    }
    return $graph->vertices;
}

# Finds longest paths in a tree graph. The subroutine finds all longest
# paths originating at graph center(s) and produces all their combinations.
sub graph_longest_paths
{
    my( $graph ) = @_;

    my @centers = graph_center( $graph );
    my @longest_paths;
    if( @centers == 1 ) {
        # Single-centered graph
        # Removing the center from longest path parts, to be added later
        my @longest_path_parts = map { [ @{$_}[1..-1] ] }
                                     graph_longest_paths_from_vertex( $graph, $centers[0] );
        for my $i (0..$#longest_path_parts) {
            for my $j (0..$#longest_path_parts) {
                next if $i == $j;
                push @longest_paths, [ reverse( @{$longest_path_parts[$i]} ),
                                       $centers[0],
                                       @{$longest_path_parts[$j]} ];
            }
        }
    } else {
        # Double-centered graph
        $graph = $graph->copy;
        $graph->delete_edge( @centers );
        my @longest_path_parts1 = graph_longest_paths_from_vertex( $graph, $centers[0] );
        my @longest_path_parts2 = graph_longest_paths_from_vertex( $graph, $centers[1] );
        for my $i (0..$#longest_path_parts1) {
            for my $j (0..$#longest_path_parts2) {
                push @longest_paths, [ reverse( @{$longest_path_parts1[$i]} ),
                                                @{$longest_path_parts2[$j]} ];
            }
        }
    }

    return @longest_paths;
}

sub graph_longest_paths_from_vertex
{
    my( $graph, $vertex ) = @_;

    my %from   = ( $vertex => undef );
    my %length = ( $vertex => 0 );
    my $bfs = Graph::Traversal::BFS->new(
        $graph,
        tree_edge =>
            sub {
                my( $u, $v ) = @_;
                ( $u, $v ) = ( $v, $u ) if exists $from{$v};
                $from{$v} = $u;
                $length{$v} = $length{$u} + 1;
            },
        start => $vertex,
    );
    $bfs->bfs;

    my @furthest_leaves;
    my $furthest_distance = 0;
    for my $vertex ( $graph->vertices ) {
        next unless exists $length{$vertex}; # May happen in disconnected graphs
        if(      $length{$vertex} < $furthest_distance ) {
            next;
        } elsif( $length{$vertex} == $furthest_distance ) {
            push @furthest_leaves, $vertex;
        } else {
            @furthest_leaves = ( $vertex );
            $furthest_distance = $length{$vertex};
        }
    }

    # Backtrack starting from the furthest leaves to collect the longest
    # paths. In the returned result path, starting vertex is the first.
    my @longest_paths;
    for my $vertex ( @furthest_leaves ) {
        my @path;
        while( $vertex ) {
            push @path, $vertex;
            $vertex = $from{$vertex};
        }
        push @longest_paths, [ reverse @path ];
    }

    return @longest_paths;
}

# Sorts locant placings from lowest to biggest
sub compare_locant_placings {
    my @first  = @{$a->[3]};
    my @second = @{$b->[3]};
    my @index = (0..scalar @first-1);

    foreach( @index ) {
        return $first[$_] <=> $second[$_] if $first[$_] <=> $second[$_];
    }
    return 0;
}

# Sorts side chain legths from lowest to biggest
sub compare_side_chain_lengths {
    my @first  = @{$a->[3]};
    my @second = @{$b->[3]};
    my @index = (0..scalar @first-1);

    foreach( @index ) {
        return $first[$_] <=> $second[$_] if $first[$_] <=> $second[$_];
    }
    return 0;
}

# Sorts given names only based on alphabetical part of the name
sub compare_only_aphabetical {
    my $a_alpha = $a;
    my $b_alpha = $b;
    $a_alpha =~ s/[^a-zA-Z]+//g;
    $b_alpha =~ s/[^a-zA-Z]+//g;

    if( $a_alpha eq 'tertbutyl' ) {
        $a_alpha = 'butyl';
    } elsif ($b_alpha eq 'tertbutyl') {
        $b_alpha = 'butyl';
    }
    return $a_alpha cmp $b_alpha;
}

# Sorts given names only based on alphabetical part of the name
sub sort_attachments {
    my @first = @{@{$a}[1]};
    my @second = @{@{$b}[1]};
    my @index = (0..scalar @first-1);

    foreach( @index ){
        my $first_alpha  = $first[$_];
        my $second_alpha = $second[$_];

        $first_alpha  =~ s/[^a-zA-Z]+//g;
        $second_alpha =~ s/[^a-zA-Z]+//g;

        if( $first_alpha eq 'tertbutyl') {
            $first_alpha = 'butyl'
        } elsif ($second_alpha eq 'tertbutyl') {
            $second_alpha = 'butyl'
        }
        if ($first_alpha lt $second_alpha) {return 1}
        elsif ($first_alpha gt $second_alpha) {return -1}
    }
    return 0;
}

# Sorts arrays from lowest to biggest by values
sub compare_arrays {
    my @first = @{$a};
    my @second = @{$b};
    my @index = (0..scalar @first-1);

    foreach( @index ){
        if ($first[$_] > $second[$_]) {return 1}
        elsif ($first[$_] < $second[$_]) {return -1}
    }
    return 0;
}

# According to https://en.wikipedia.org/wiki/IUPAC_numerical_multiplier
sub IUPAC_numerical_multiplier
{
    my( $N, $is_middle ) = @_;

    my $ones      = $N % 10;
    my $tens      = int( $N /   10 ) % 10;
    my $hundreds  = int( $N /  100 ) % 10;
    my $thousands = int( $N / 1000 ) % 10;

    my @prefix = ( '', 'hen', 'di', 'tri', 'tetra', 'penta', 'hexa', 'hepta', 'octa', 'nona' );

    return 'heni' if $N == 1 && $is_middle;
    return 'mono' if $N == 1;
    return 'do'   if $N == 2 && $is_middle;
    return 'tria' if $N == 3 && $is_middle;

    if( $N < 10 ) {
        my $value = $prefix[$ones];
        $value =~ s/a$// unless $is_middle;
        return $value;
    }

    return 'dec'    . ($is_middle ? 'a' : '') if $N == 10;
    return 'undec'  . ($is_middle ? 'a' : '') if $N == 11;
    return 'dodec'  . ($is_middle ? 'a' : '') if $N == 12;
    return 'tridec' . ($is_middle ? 'a' : '') if $N == 13;
    return IUPAC_numerical_multiplier( $ones, 1 ) . 'dec' . ($is_middle ? 'a' : '') if $N < 20;
    return 'icos'   . ($is_middle ? 'a' : '') if $N == 20;
    return 'henicos'. ($is_middle ? 'a' : '') if $N == 21;
    return 'docos'  . ($is_middle ? 'a' : '') if $N == 22;
    return 'tricos' . ($is_middle ? 'a' : '') if $N == 23;
    return IUPAC_numerical_multiplier( $ones, 1 ) . 'cos' . ($is_middle ? 'a' : '') if $N < 30;

    if( $N < 100 ) {
        if( $ones == 1 ) {
            return 'hen' . IUPAC_numerical_multiplier( $tens, 1 ) . 'cont' . ($is_middle ? 'a' : '');
        } elsif ( $ones == 3 ) {
            return 'tri' . IUPAC_numerical_multiplier( $tens, 1 ) . 'cont' . ($is_middle ? 'a' : '');
        } else {
            return IUPAC_numerical_multiplier( $ones, 1 ) .
                   IUPAC_numerical_multiplier( $tens, 1 ) .
                   'cont' . ($is_middle ? 'a' : '');
        }
    }

    return 'trihect' if $N == 103;
    return IUPAC_numerical_multiplier( int( $tens . $ones ), 1 ) . 'hect' if $N < 200;
    return IUPAC_numerical_multiplier( int( $tens . $ones ), 1 ) . $prefix[$hundreds] . 'ct' . ($is_middle ? 'a' : '') if $N <  1000;
    return IUPAC_numerical_multiplier( int( $hundreds . $tens . $ones ), 1 ) . 'kili'                     if $N <  2000;
    return IUPAC_numerical_multiplier( int( $hundreds . $tens . $ones ), 1 ) . $prefix[$thousands] . 'li' if $N < 10000;
    die "cannot generate IUPAC numerical multiplier for $N\n";
}

sub alkane_chain_name
{
    my( $N ) = @_;

    my @names = qw( ? meth eth prop but );

    return $names[$N] if $N < @names;
    return IUPAC_numerical_multiplier( $N );
}

1;
