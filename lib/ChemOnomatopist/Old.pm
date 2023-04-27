package ChemOnomatopist::Old;

use strict;
use warnings;

# ABSTRACT: Give molecule a name
# VERSION

use parent ChemOnomatopist::;

use ChemOnomatopist::Util::Graph qw(
    BFS_calculate_chain_length
    BFS_is_chain_branched
);
use Clone qw( clone );
use Graph::Traversal::BFS;
use Graph::Undirected;
use List::Util qw( all any max min sum0 uniq );
use Scalar::Util qw( blessed );

# There must be a nicer way to handle calls to parent...
sub AUTOLOAD {
    our $AUTOLOAD;
    my $call = $AUTOLOAD;
    $call =~ s/.*:://;
    return if $call eq 'DESTROY';
    return ChemOnomatopist->can( $call )->( @_ );
}

no warnings 'recursion';

our @numbers = @ChemOnomatopist::numbers;
our @numberskis = @ChemOnomatopist::numberskis;

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
    if( $graph->edges != $graph->vertices - 1 ) {
        # If it is not a tree, than the graph has cycles, and we have to
        # do our best to recognise them. To make it easier, hydrogen atoms
        # are removed here for now.
        $graph->delete_vertices( grep { $_->{symbol} eq 'H' } $graph->vertices );
        if( $graph->edges != $graph->vertices ) {
            die "cannot handle cycles with attachments for now\n";
        }
        if( any { uc $_->{symbol} ne 'C' } $graph->vertices ) {
            die "cannot handle heterocycles for now\n";
        }
        if( all { $_->{symbol} eq 'C' } $graph->vertices ) {
            # Cycloalkane detected
            return 'cyclo' . alkane_chain_name( scalar $graph->vertices ) . 'ane';
        }
        if( ( all { $_->{symbol} eq 'c' } $graph->vertices ) &&
            ( scalar $graph->vertices ) =~ /^(4|6|8|10|12|14|16)$/ ) {
            # Annulene detected
            return 'cyclo' . $numbers[scalar $graph->vertices] .
                   $numbers[scalar $graph->vertices / 2] . 'ene';
        }
        # No other types of graphs with cycles can be processed for now
        die "only limited set of homocycles is supported for now\n";
    }

    # Check for unsupported elements.
    if( any { !is_element( $_, 'C' ) && !is_element( $_, 'H' ) }
            $graph->vertices ) {
        die "cannot handle atoms other than C and H now\n";
    }

    my( $order ) = select_main_chain( $graph->copy );
    my @chain;
    for my $curr_vertex (@$order) {
        my( $vertex ) = grep { $_->{number} == $curr_vertex } $graph->vertices;
        push @chain, $vertex;
    }
    return get_mainchain_name( $graph->copy, \@chain );
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

    my( $chain ) = sort { cmp_arrays( $a, $b ) }
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
                    my $attachment_name = get_sidechain_name( $graph_copy, $neighbour );
                    $attachment_name = "($attachment_name)" if $attachment_name =~ /^[0-9]/;

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

1;
