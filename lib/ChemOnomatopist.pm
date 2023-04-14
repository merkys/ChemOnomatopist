package ChemOnomatopist;

use strict;
use warnings;

# ABSTRACT: Give molecule a name
# VERSION

use ChemOnomatopist::Chain;
use ChemOnomatopist::ChainHalf;
use ChemOnomatopist::Group;
use ChemOnomatopist::Group::Carbonyl;
use ChemOnomatopist::Group::Carboxyl;
use ChemOnomatopist::Group::Hydroxy;
use ChemOnomatopist::Util::Graph qw(
    BFS_calculate_chain_length
    BFS_is_chain_branched
    graph_center
    graph_longest_paths_from_vertex
    tree_branch_positions
    tree_number_of_branches
);
use Chemistry::OpenSMILES::Writer qw( write_SMILES );
use Clone qw( clone );
use Graph::Nauty qw( canonical_order );
use Graph::Traversal::BFS;
use Graph::Undirected;
use List::Util qw( all any max sum0 uniq );
use Scalar::Util qw( blessed );
use Set::Object qw( set );

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
    my( $what, $use_new_method ) = @_;

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

    # Check for unsupported elements.
    if( any { !is_element( $_, 'C' ) && !is_element( $_, 'H' ) }
            $graph->vertices ) {
        die "cannot handle atoms other than C and H now\n";
    }

    my $order;
    if( $use_new_method ) {
        $order = [ map { $_->{number} } select_main_chain_new( $graph->copy ) ];
    } else {
        ( $order ) = select_main_chain( $graph->copy );
    }
    return get_chain_2( $graph->copy,
                        $order,
                        { choose_direction => 1 } ) . 'ane';
}

# get_sidechain_name() receives a graph and a position to start the chain in it.
# From that position it finds the longest chain and returns the constructed name.
sub get_sidechain_name
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
        unshift @chain, $end;

        my $min = $end;
        for my $neighbour ($graph->neighbours( $end )) {
            next if !exists $order{$neighbour};
            next if $order{$neighbour} >= $order{$min};
            $min = $neighbour;
        }
        last if $min eq $end;
        $end = $min;
    }
    $graph->delete_path( @chain );

    # Examine the attachments to the main chain: delete the edges
    # connecting them to the main chain, at the same time giving them
    # names according to their lengths via calls to get_sidechain_name()
    my %attachments;
    for my $i (0..$#chain) {
        my $atom = $chain[$i];
        for my $neighbour ($graph->neighbours( $atom )) {
            $graph->delete_edge( $atom, $neighbour );
            next if is_element( $neighbour, 'H' );

            my $attachment_name = get_sidechain_name( $graph, $neighbour ) . 'yl';
            $attachment_name = bracket( $attachment_name ) if $attachment_name =~ /^[0-9]/;
            push @{$attachments{$attachment_name}}, $i;
        }
    }

    # Collecting names of all the attachments
    my $name = '';
    for my $attachment_name (sort { $a cmp $b } keys %attachments) {
        my $number = IUPAC_numerical_multiplier( scalar @{$attachments{$attachment_name}} );
        $number = '' if $number eq 'mono';
        $number .= 'a' unless $number =~ /^(|\?|.*i)$/;
        $name = $name ? $name . '-' : $name;
        $name .= join( ',', map { $_ + 1 } @{$attachments{$attachment_name}} ) .
                 '-' . $number . $attachment_name;
    }

    return $name . alkane_chain_name( scalar @chain );
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

    # Disconnect the main chain: this way every main chain atom remains
    # connected only to the side chains.
    $graph->delete_path( @chain );

    # Examine the attachments to the main chain: delete the edges
    # connecting them to the main chain, at the same time giving them
    # names according to their lengths via calls to get_sidechain_name()
    my %attachments;
    for my $i (0..$#chain) {
        my $atom = $chain[$i];
        for my $neighbour ($graph->neighbours( $atom )) {
            $graph->delete_edge( $atom, $neighbour );
            next if is_element( $neighbour, 'H' );

            my $attachment_name = get_sidechain_name( $graph, $neighbour ) . 'yl';
            $attachment_name = bracket( $attachment_name ) if $attachment_name =~ /^[0-9]/;
            push @{$attachments{$attachment_name}}, $i;
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
        next unless @{$attachments{$preferrable_names{$att_name}}} > 1;
        $attachments{'(' . $preferrable_names{$att_name} . ')'} =
                        $attachments{$preferrable_names{$att_name}};
        delete $attachments{$preferrable_names{$att_name}};
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
            $number = '' if $number eq 'mono';
            $number .= 'a' unless $number =~ /^(|\?|.*i)$/;
        }
        $name .= join( ',', map { $_ + 1 } @{$attachments{$attachment_name}} ) .
                 '-' . $number . $attachment_name;
    }

    return $name . alkane_chain_name( scalar @chain );
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

    return;
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
    my ( $graph, $start ) = @_;

    my $carbon_graph = $graph->copy;
    $carbon_graph->delete_vertices( grep { !is_element( $_, 'C') } $carbon_graph->vertices );

    my $is_any_visited;
    my %lengths;
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

    # Make a carbon-only graph
    my $carbon_graph = $graph->copy;
    $carbon_graph->delete_vertices( grep { !is_element( $_, 'C' ) } $carbon_graph->vertices );

    # Going through every vertice in "farthest" array and creating tree-like structures
    my @all_trees;
    for (my $i = 0; $i < scalar @farthest; $i++) {
        my %tree = ( $farthest[$i] => [ $farthest[$i], 0 ] );

        my( $vertex ) = grep { $_->{number} eq $farthest[$i] } $carbon_graph->vertices;
        # Start creation of the tree from all the starting vertices
        push @all_trees, \%{create_tree( $carbon_graph->copy, $vertex, \%tree )};
    }

    my $trees;

    # Extracts arrays of all longest chains with numbers of vertices (in order)
    # from each tree-like structure
    my @main_chains = prepare_paths( @all_trees );

    # From all possible main chains in tree-like structures, subroutine returns
    # the ones that has the greatest number of side chains. Also returns only the
    # trees that still have some possible main chains after selection
    ( $trees, @main_chains ) = rule_greatest_number_of_side_chains(
                                    $carbon_graph->copy,
                                    \@main_chains,
                                    @all_trees,
                               );
    return @{$main_chains[0]} if scalar @{$main_chains[0]} == 1;

    # If more than one chain is left, second rule is applied.
    # From all main chains left, subroutine selects all that have the
    # lowest numbered locants. Also the trees that have possible main
    # chains returned
    ( $trees, @main_chains ) = rule_lowest_numbered_locants(
                                    $carbon_graph->copy,
                                    @main_chains,
                                    @$trees,
                               );
    return @{$main_chains[0]} if scalar @{$main_chains[0]} == 1;

    # If more than one chain is left, third rule is applied.
    # From all main chains left, subroutine selects all that have the
    # most carbons in side chains. Also the trees that have possible main
    # chains returned
    ( $trees, @main_chains ) = rule_most_carbon_in_side_chains(
                                    $carbon_graph->copy,
                                    @main_chains,
                                    @$trees,
                               );
    return @{$main_chains[0]} if scalar @{$main_chains[0]} == 1;

    # If more than one chain is left, fourth rule is applied.
    # From all main chains left, subroutine selects all that have
    # the least branched side chains. Also the trees that have
    # possible main chains returned
    ( $trees, @main_chains ) = rule_least_branched_side_chains(
                                    $carbon_graph->copy,
                                    @main_chains,
                                    @$trees,
                               );
    return @{$main_chains[0]} if scalar @{$main_chains[0]} == 1;

    # If more than one chain is left, program states that there are
    # few chains that are identical by all the rules and selects
    # one from all that are left.
    # One main chain is picked from all that are left
    my $main_chain = rule_pick_chain_from_valid(
                        $carbon_graph->copy,
                        @main_chains,
                        @$trees,
                     );
    return $main_chain;
}

# Selects the main chain by evaluating its parts
sub select_main_chain_new
{
    my( $tree ) = @_;

    # Remove non-carbon atoms
    $tree = $tree->copy;
    $tree->delete_vertices( grep { !is_element( $_, 'C' ) } $tree->vertices );

    # Here the candidate halves for the longest (and "best") path are placed in @path_parts.
    # Each of candidate halves start with center atom.
    my @center = graph_center( $tree );
    my @path_parts;
    if( @center == 1 ) {
        # Longest path has odd length
        for my $path ( graph_longest_paths_from_vertex( $tree, $center[0] ) ) {
            push @path_parts,
                 ChemOnomatopist::ChainHalf->new( $tree, scalar @center, @$path );
        }
    } else {
        # Longest path has even length
        # Graph copy is needed as it has to be modified to find the longest paths from both centers.
        my $copy = $tree->copy;
        $copy->delete_edge( @center );
        for my $vertex ( @center ) {
            push @path_parts,
                 map { ChemOnomatopist::ChainHalf->new( $tree, scalar @center, @$_ ) }
                     graph_longest_paths_from_vertex( $copy, $vertex );
        }
    }

    for my $rule ( sub { return @_ },
                   \&rule_greatest_number_of_side_chains_new ) {
        my @path_parts_now = $rule->( @path_parts );

        # CHECK: Can a rule cause disappearance of parts?
        # FIXME: Methane will have a single part
        next if @path_parts_now < 2;

        @path_parts = @path_parts_now; # Narrow down the selection

        # If two chains cannot be chosen now, pass on to the next rule
        next unless @path_parts == 2;
        next unless $path_parts[0] ne $path_parts[1];

        # FIXME: There should be a nicer way to do this...
        # FIXME: rule_lowest_numbered_locants_new() may return more than one chain, need to sort them by attachment names
        my @chains = rule_lowest_numbered_locants_new( @path_parts );
        return $chains[0]->vertices;
    }

    # At this point we are left with a set of longest chains all having the same number of side chains

    my @chains = rule_lowest_numbered_locants_new( @path_parts );

    for my $rule ( sub { return @_ },
                   \&rule_most_carbon_in_side_chains_new,
                   \&rule_least_branched_side_chains_new,
                   \&pick_chain_with_lowest_attachments_alphabetically_new ) {
        my @chains_now = $rule->( @chains );

        # CHECK: Can a rule cause disappearance of all chains?
        next unless @chains_now;

        @chains = @chains_now; # Narrow down the selection

        # If a single chain cannot be chosen now, pass on to the next rule
        next unless @chains == 1;

        return $chains[0]->vertices;
    }

    # TODO: Handle the case when none of the rules select proper chains
    return ();
}

sub rule_greatest_number_of_side_chains_new
{
    my( @path_parts ) = @_;

    my @sorted_values = sort uniq map { $_->number_of_branches } @path_parts;

    my @path_parts_now;
    my $seen_groups = set();
    for my $value (reverse @sorted_values) {
        last if $seen_groups->size >= 2;
        push @path_parts_now,
             grep { $_->number_of_branches == $value &&
                    !$seen_groups->has( $_->group ) } @path_parts;
        $seen_groups = set( map { $_->group } @path_parts_now );
    }

    return @path_parts_now;
}

sub rule_lowest_numbered_locants_new
{
    my( @path_parts ) = @_;

    my @chains;
    my $min_value;

    for my $part1 (@path_parts) {
        for my $part2 (@path_parts) {
            next if $part1->group eq $part2->group;
            my $chain = ChemOnomatopist::Chain->new( $part1, $part2 );
            my $value = $chain->locant_positions;
            if( @chains ) {
                if( $value < $min_value ) {
                    @chains = ( $chain );
                    $min_value = $value;
                } elsif( $value == $min_value ) {
                    push @chains, $chain;
                }
            } else {
                push @chains, $chain;
                $min_value = $value;
            }
        }
    }

    return @chains;
}

sub rule_most_carbon_in_side_chains_new
{
    my( @chains ) = @_;

    my( $max_value ) = reverse sort uniq map { $_->number_of_carbons } @chains;
    return grep { $_->number_of_carbons == $max_value } @chains;
}

sub rule_least_branched_side_chains_new
{
    my( @chains ) = @_;

    my( $min_value ) = sort uniq map { $_->number_of_branches_in_sidechains } @chains;
    return grep { $_->number_of_branches_in_sidechains == $min_value } @chains;
}

# FIXME: This rule is dumb now: it just returns the first chain it gets
sub pick_chain_with_lowest_attachments_alphabetically_new
{
    return $_[0];
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
            push @all_chains, save_main_chain_vertices_in_array(
                                $ending,
                                [$ending],
                                $tree
                              );
        }
    }

    # Adds reverted chains if they are not present yet as the longest chains
    for my $chain (@all_chains) {
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
        my( $first ) = grep { $structure{$_}->[0] == 0 } keys %structure;
        my @structure_chains = grep { $_->[0] == $first || $_->[-1] == $first } @$chains;

        # Structure with index of the tree, beginning and ending of the chain,
        # lengths of side chains of the chain created for each tree
        for my $chain (@structure_chains) {
            push @side_chain_lengths,
                 [$index, $chain->[0], $chain->[-1],
                    [find_lengths_of_side_chains(
                        $graph->copy,
                        $chain->[-1],
                        \@{$chain},
                        [],
                        $tree,
                        scalar @{$chain}
                    )]
                 ];
        }
        $index++;
    }

    # All chains that have the highest number of carbons in side chains are selected and returned
    my @sorted_final = sort compare_locant_placings @side_chain_lengths;

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

        # Beginning of the structure is found. Then all chains that belong to
        # the current tree are selected
        my( $first ) = grep { $structure{$_}->[0] == 0 } keys %structure;
        my @chains_in_the_tree = grep { $_->[0] == $first || $_->[-1] == $first } @$chains;

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

    my( $chain ) = sort compare_arrays
                        pick_chain_with_lowest_attachments_alphabetically( $graph, $chains, @trees );
    return $chain;
}

# Subroutine selects chain that has the lowest attachments by alpabetical naming
sub pick_chain_with_lowest_attachments_alphabetically
{
    my( $graph, $chains, @trees ) = @_;

    # Locant placements are found for all trees
    my @attachments;
    for my $i (0..$#trees) {
        my %structure = %{clone $trees[$i]};

        # Reference to parental chain is removed from the boxes
        for my $key (keys %structure) {
            shift @{$structure{$key}};
        }

        # Beginning of the structure is found. Then all chains that belongs to
        # the current tree are selected
        my( $first ) = grep { $structure{$_}->[0] == 0 } keys %structure;
        my @chains_in_the_tree = grep { $_->[0] == $first || $_->[-1] == $first } @$chains;

        # Placings of the locants found for each chain
        for my $chain (@chains_in_the_tree) {
            my @attachments_only;
            for my $locant (uniq find_locant_placing( $graph, $chain )) {
                my( $vertex ) = grep { $_->{number} == $chain->[$locant-1] } $graph->vertices;

                # Cycle through non-mainchain neighbours:
                for my $neighbour ($graph->neighbours( $vertex )) {
                    next if any { $neighbour->{number} eq $_ } @$chain;

                    # Find the name for a sidechain
                    my $graph_copy = $graph->copy;
                    $graph_copy->delete_edge( $vertex, $neighbour );
                    my $attachment_name = get_sidechain_name( $graph_copy, $neighbour ) . 'yl';
                    $attachment_name = "($attachment_name)" if $attachment_name =~ /^[0-9]/;

                    # Replace systematic IUPAC attachment names with their preferrable ones
                    if( exists $preferrable_names{$attachment_name} ) {
                        $attachment_name = $preferrable_names{$attachment_name};
                    }

                    push @attachments_only, $attachment_name;
                }
            }

            push @attachments, [clone( $chain ), \@attachments_only];
        }
    }

    # All chains that have the same - alpabetically lowest attachments selected
    my( $best ) = sort { cmp_attachments( $a->[1], $b->[1] ) } @attachments;
    my @correct_chains_all = grep { join( ',', @{$_->[1]} ) eq
                                    join( ',', @{$best->[1]} ) }
                                  @attachments;
    return map { $_->[0] } @correct_chains_all;
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

# Find placings of all locants in the chain.
# Returns an array of indices denoting locant placements.
# This code treats vertices with degrees larger than 2 as having sidechain attachments.
sub find_locant_placing
{
    my( $graph, $main_chain ) = @_;

    # Indices received instead of vertices, transform them.
    # This later on should be removed.
    if( @$main_chain && !ref $main_chain->[0] ) {
        my %vertices_by_id = map { ( $_->{number} => $_ ) } $graph->vertices;
        $main_chain = [ map { $vertices_by_id{$_} } @$main_chain ];
    }

    # Visit all attachments and memorize their attachment positions
    my @locants;
    for my $i (0..$#$main_chain) {
        my $vertex = $main_chain->[$i];
        next unless $graph->degree( $vertex ) > 2;
        push @locants, ( $i ) x ( $graph->degree( $vertex ) - 2 );
    }

    return @locants;
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

# Sorts locant placings from lowest to largest
# This had code identical to compare_side_chain_lengths(), thus calls to the latter have been redirected here.
sub compare_locant_placings {
    my @A = @$a;
    my @B = @$b;

    if( @A >= 4 && ref $A[3] ) {
        # This is the "old" data structure
        @A = @{$A[3]};
        @B = @{$B[3]};
    }

    for (0..$#A) {
        return $A[$_] <=> $B[$_] if $A[$_] <=> $B[$_];
    }

    return 0;
}

# Sorts given names only based on alphabetical part of the name
sub compare_only_aphabetical {
    my $a_alpha = $a;
    my $b_alpha = $b;
    $a_alpha =~ s/[^a-zA-Z]+//g;
    $b_alpha =~ s/[^a-zA-Z]+//g;

    # FIXME: Not sure how to sort 'butyl' and 'tertbutyl', but stable
    # order is important.
    unless( grep( { $_ eq 'butyl'     } ( $a_alpha, $b_alpha ) ) &&
            grep( { $_ eq 'tertbutyl' } ( $a_alpha, $b_alpha ) ) ) {
        $a_alpha = 'butyl' if $a_alpha eq 'tertbutyl';
        $b_alpha = 'butyl' if $b_alpha eq 'tertbutyl';
    }

    return $a_alpha cmp $b_alpha;
}

# Sorts given names only based on alphabetical part of the name
sub cmp_attachments
{
    my( $a, $b ) = @_;
    my @a = @{$a};
    my @b = @{$b};

    for (0..$#a) {
        my $a_alpha = $a[$_];
        my $b_alpha = $b[$_];

        $a_alpha =~ s/[^a-zA-Z]+//g;
        $b_alpha =~ s/[^a-zA-Z]+//g;

        $a_alpha = 'butyl' if $a_alpha eq 'tertbutyl';
        $b_alpha = 'butyl' if $b_alpha eq 'tertbutyl';

        return $b_alpha cmp $a_alpha if $b_alpha cmp $a_alpha;
    }

    return 0;
}

# Sorts arrays from lowest to biggest by values
sub compare_arrays {
    my @first  = @$a;
    my @second = @$b;
    my @index  = (0..scalar @first-1);

    foreach( @index ){
        return $first[$_] <=> $second[$_] if $first[$_] <=> $second[$_];
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

sub bracket
{
    my( $name ) = @_;
    return $name =~ /\(/ ? "[$name]" : "($name)";
}

1;
