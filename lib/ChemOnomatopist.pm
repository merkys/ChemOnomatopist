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
use Scalar::Util qw( blessed );
use Clone qw( clone );

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
        my $bfs = Graph::Traversal::BFS->new(
            $carbon_graph,
            pre => sub { if(!$is_any_visited) {
                            $lengths{$_[0]->{number}} = 0; $is_any_visited = 1  } },
            tree_edge => sub { if ( !defined $lengths{$_[0]->{number}} ) {
                ( $_[0], $_[1] ) = ( $_[1], $_[0] );
            }
            $lengths{$_[1]->{number}} = $lengths{$_[0]->{number}} + 1},
            start => $start);
        my @order = $bfs->bfs;
        return @order;
    }
    else{
        my $bfs = Graph::Traversal::BFS->new( $carbon_graph );
        my @order = $bfs->bfs;
        return @order;
    }
}

sub BFS_calculate_chain_length
{
    my $length = 0;
    my ( $graph, $start ) = @_;
    my $bfs = Graph::Traversal::BFS->new( $graph, start => $start);
    my @order = $bfs->bfs;
    return scalar @order;
}

sub BFS_is_chain_branched
{
    my $branched = 0;
    my ( $graph, $start ) = @_;

    my $graph_copy = $graph->copy;

    my $bfs = Graph::Traversal::BFS->new(
        $graph,
        pre => sub {
            my @neighbours = $graph_copy->neighbours( $_[0] );
            if ( scalar @neighbours > 1 ) { $branched = 1; }
            $graph_copy->delete_vertex($_[0]);
        },
        start => $start);
    my @order = $bfs->bfs;
    return $branched;
}

# Returns tree like structures for all the longest paths
sub create_structure
{
    my @all_trees;
    my ( $graph ) = @_;
    my @order = BFS_order_carbons_only($graph);

    my $start = $order[-1];

    my @second_order = BFS_order_carbons_only($graph, $start);

    my $end = $second_order[-1];

    # Finding all farthest vertices
    my @farthest = grep { $lengths{$_} eq $lengths{$end->{number}} } keys %lengths;
    push(@farthest, $start->{number});

    if (scalar @farthest == 2) {
        return @order;
    }

    for (my $i = 0; $i < scalar(@farthest); $i++) {
        my %tree;
        my @value_array;
        push(@value_array, $farthest[$i]);
        push(@value_array, 0);
        $tree{$farthest[$i]} = [@value_array];

        my $carbon_graph = $graph->copy;
        $carbon_graph->delete_vertices( grep {!is_element( $_, 'C') } $carbon_graph->vertices );

        my @vertice = grep { $_->{number} eq $farthest[$i] } $carbon_graph->vertices;

        # Start creation of the tree from all the starting vertices
        push(@all_trees, \%{create_tree($carbon_graph, $vertice[0], \%tree)});
    }
    my @main_chains = prepare_paths(@all_trees);
    my $trees;
    my $carbon_graph = $graph->copy;
    $carbon_graph->delete_vertices(
        grep {!is_element( $_, 'C') } $carbon_graph->vertices
    );
    ($trees, @main_chains) = rule_greatest_number_of_side_chains($carbon_graph, \@main_chains, @all_trees);

    if (scalar @{$main_chains[0]} != 1){
        my @trr = @{$trees};
        my $carbon_graph = $graph->copy;
        $carbon_graph->delete_vertices(
            grep {!is_element( $_, 'C') } $carbon_graph->vertices
        );
        ($trees, @main_chains) = rule_lowest_numbered_locants(
                                    $carbon_graph, @main_chains, @trr
                                 );
        if (scalar @{$main_chains[0]} != 1){
            my @trr = @{$trees};
            my $carbon_graph = $graph->copy;
            $carbon_graph->delete_vertices(
                grep {!is_element( $_, 'C') } $carbon_graph->vertices
            );
            ($trees, @main_chains) = rule_most_carbon_in_side_chains(
                                        $carbon_graph, @main_chains, @trr
                                    );
            if (scalar @{$main_chains[0]} != 1){
                my @trr = @{$trees};
                my $carbon_graph = $graph->copy;
                $carbon_graph->delete_vertices(
                    grep {!is_element( $_, 'C') } $carbon_graph->vertices
                );
                ($trees, @main_chains) = rule_least_branched_side_chains($carbon_graph, @main_chains, @trr);
                if(scalar @{$main_chains[0]} != 1){
                    #pick one from all equal chains
                }
            }
            else{
                #return created order
            }
        }
        else{
            #return created order
        }
    }
    else{
        #return created order
    }
}
# Creating tree like structure for all the longest paths in molecule
sub create_tree
{
    my ( $graph, $atom, $tree ) = @_;

    my @neighbours = $graph->neighbours( $atom );

    my @array = @ { $tree->{$atom->{number}}};
    my @new_array;

    push(@new_array, $atom->{number});

    shift @array;

    foreach my $arr (@array){
        push(@new_array, $arr + 1);
    }

    if (scalar @neighbours == 1) {
        unless (exists $tree->{$neighbours[0]->{number}}){
            $tree->{$neighbours[0]->{number}} = [@new_array];
            $graph->delete_vertex($atom);
            create_tree($graph, $neighbours[0], $tree);
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
    return $tree;
}

sub prepare_paths
{
    my ( @trees ) = @_;

    my $trees_copy = clone(\@trees);

    my @all_chains;

    foreach my $tree (@{$trees_copy}){

        my %structure = %{$tree};

        my @sorted = sort {
                            $structure{$a}->[1] <=> $structure{$b}->[1]
                          } keys %structure;

        my $last = $sorted[-1];
        my @chain_ending = grep{ $structure{$_}->[1] == $structure{$last}->[1] }
                                                                keys %structure;

        foreach my $ending (@chain_ending){
            my @vertex_array = ($ending);
            my %structure2 = %{$tree};
            push (@all_chains,
                save_main_chain_vertices_in_array(
                    $ending, \@vertex_array, \%structure2
                )
            )
        }
    }
    return @all_chains;
}

# Tries to find the chain which has the greatest number of side chains
sub rule_greatest_number_of_side_chains
{
    my ( $graph, $chains, @trees ) = @_;

    my $trees_copy = clone(\@trees);
    my $index = 0;
    my @number_of_side_chains = ();

    foreach my $tree (@{$trees_copy}){

        my %structure = %{clone $tree};
        my %structure2 = %{$tree};

        foreach my $key (keys %structure)
        {
            shift @ { $structure{$key} };
        }

        my @first = grep{ $structure{$_}->[0] == 0 } keys %structure;
        my @chains_in_the_tree = grep {@{$_}[0] == $first[0]} @{$chains};

        for my $chain (@chains_in_the_tree){
            my $graph_copy = $graph->copy;
            push ( @number_of_side_chains,
                [$index, @{$chain}[0], @{$chain}[-1],
                    find_number_of_side_chains(
                        $graph_copy,
                        \@{$chain},
                        \%structure2
                    )
                ]
            )
        }
        $index += 1;
    }

    my @sorted_paths = sort {
                                $a->[3] <=> $b->[3]
                            } @number_of_side_chains;

    my $path_length = $sorted_paths[-1][3];
    my @longest_paths = grep {$_->[3] == $path_length} @number_of_side_chains;
    my %seen;
    my @uniq_longest_paths = grep { !$seen{$_->[0]}++ } @longest_paths;
    my @result = @trees[map {$_->[0]} @uniq_longest_paths];
    my @eligible_chains;

    for my $chain (@{$chains}){
        if (grep {$_->[1] == $chain->[0] and $_->[2] == $chain->[-1]} @longest_paths){
            push(@eligible_chains, $chain);
        }
    }
    return \@result, \@eligible_chains;
}

# Tries to find the chain which has the lowest-numbered locants
sub rule_lowest_numbered_locants
{
    my ( $graph, $chains, @trees ) = @_;

    my $trees_copy = clone(\@trees);
    my $index = 0;
    my @locant_placing = ();

    foreach my $tree (@{$trees_copy}){

        my %structure = %{clone $tree};
        my %structure2 = %{$tree};

        foreach my $key (keys %structure)
        {
            shift @ { $structure{$key} };
        }

        my @first = grep{ $structure{$_}->[0] == 0 } keys %structure;
        my @chains_in_the_tree = grep {@{$_}[0] == $first[0]} @{$chains};

        for my $chain (@chains_in_the_tree){
            my $graph_copy = $graph->copy;
            push ( @locant_placing,
                [$index, @{$chain}[0], @{$chain}[-1],
                    [find_locant_placing(
                        $graph_copy,
                        \@{$chain},
                        \%structure2
                    )]
                ]
            )
        }
        $index += 1;
    }

    my @sorted_paths = sort compare_locant_placings reverse @locant_placing;

    my $lowest_locants = $sorted_paths[0][3];
    my @lowest_locants_paths = grep {
                    join("", ( @{$_->[3]})) == join("", (@{$lowest_locants}))
                                    } @locant_placing;
    my %seen;
    my @uniq_lowest_locants_paths = grep { !$seen{$_->[0]}++ } @lowest_locants_paths;
    my @result = @trees[map {$_->[0]} @uniq_lowest_locants_paths];

    my @eligible_chains;

    for my $chain (@{$chains}){
        if (grep {$_->[1] == $chain->[0] and $_->[2] == $chain->[-1]} @lowest_locants_paths){
            push(@eligible_chains, $chain);
        }
    }
    return \@result, \@eligible_chains;
}

# Tries to find chain that has the greatest number of carbon atoms in the smaller
# side chains
sub rule_most_carbon_in_side_chains
{
    my ( $graph, $chains, @trees ) = @_;

    my $trees_copy = clone(\@trees);
    my @side_chain_lengths = ();
    my $index = 0;

    foreach my $tree (@{$trees_copy}){

        my %structure = %{ clone $tree};
        my %structure2 = %{$tree};
        my @all_vertices = keys %structure;

        foreach my $key (keys %structure)
        {
            shift @ { $structure{$key} };
        }

        my @sorted = sort {
                            @{$structure{$a}} <=> @{$structure{$b}}
                           or
                            $structure{$a}->[0] cmp $structure{$b}->[0]
                          } keys %structure;
        my $last = $sorted[-1];

        my @first = grep{ $structure{$_}->[0] == 0 } keys %structure;

        my @structure_chains = grep{@{$_}[0] == $first[0]} @{$chains};

        for my $chain (@structure_chains){
            my $graph_copy = $graph->copy;
            my @side_chain_length = ();
            push ( @side_chain_lengths,
                [$index, @{$chain}[0], @{$chain}[-1],
                    [find_lengths_of_side_chains(
                        $graph_copy,
                        @{$chain}[-1],
                        \@{$chain},
                        \@side_chain_length,
                        \%structure2
                    )]
                ]
            )
        }
        $index += 1;
    }

    my @sorted_final = sort compare_side_chain_lengths @side_chain_lengths;

    my $last = $sorted_final[-1][3];
    my @greatest_no_of_side_chains_paths = grep {
                    join("", @{@{$_}[3]}) eq join("", @{$last})
                                    } @sorted_final;
    my @eligible_chains;

    for my $chain (@{$chains}){
        if (grep {$_->[1] == $chain->[0] and $_->[2] == $chain->[-1]}
            @greatest_no_of_side_chains_paths
        ){
            push(@eligible_chains, $chain);
        }
    }
    my %seen;
    my @uniq_side_chain_paths =
                grep { !$seen{$_->[0]}++ } @greatest_no_of_side_chains_paths;
    my @result = @trees[map {$_->[0]} @uniq_side_chain_paths];

    return \@result, \@eligible_chains;
}

# Tries to find chain that have the least branched side chains
sub rule_least_branched_side_chains
{
    my ( $graph, $chains, @trees ) = @_;

    my $trees_copy = clone(\@trees);
    my $index = 0;
    my @number_of_branched_side_chains = ();

    foreach my $tree (@{$trees_copy}){

        my %structure = %{clone $tree};
        my %structure2 = %{$tree};

        foreach my $key (keys %structure)
        {
            shift @ { $structure{$key} };
        }

        my @first = grep{ $structure{$_}->[0] == 0 } keys %structure;
        my @chains_in_the_tree = grep {@{$_}[0] == $first[0]} @{$chains};

        for my $chain (@chains_in_the_tree){
            my $graph_copy = $graph->copy;
            push ( @number_of_branched_side_chains,
                [$index, @{$chain}[0], @{$chain}[-1],
                    find_number_of_branched_side_chains(
                        $graph_copy,
                        \@{$chain},
                        \%structure2
                    )
                ]
            )
        }
        $index += 1;
    }
    my @sorted_paths = sort {
                                $a->[3] <=> $b->[3]
                            } @number_of_branched_side_chains;

    my $path_length = $sorted_paths[-1][3];
    my @longest_paths = grep {$_->[3] == $path_length} @number_of_branched_side_chains;
    my %seen;
    my @uniq_longest_paths = grep { !$seen{$_->[0]}++ } @longest_paths;
    my @result = @trees[map {$_->[0]} @uniq_longest_paths];
    my @eligible_chains;

    for my $chain (@{$chains}){
        if (grep {$_->[1] == $chain->[0] and $_->[2] == $chain->[-1]} @longest_paths){
            push(@eligible_chains, $chain);
        }
    }
    return \@result, \@eligible_chains;
}

# Returns array that contains numbers of vertices that are in side chains
sub remove_main_chain_vertices_from_array
{
    my ( $curr_vertex, $all_vertices, $structure ) = @_;

    if ($structure->{$curr_vertex}->[0] == $curr_vertex) {
        my @other_vertices = grep { $_ != $curr_vertex } @{ $all_vertices};
        return @other_vertices;
    }
    else {
        my @other_vertices = grep { $_ != $curr_vertex } @{ $all_vertices};
        remove_main_chain_vertices_from_array(
        $structure->{$curr_vertex}->[0], \@other_vertices , $structure);
    }
}

# Returns array that contains numbers of vertices that are in main chain
sub save_main_chain_vertices_in_array
{
    my ( $curr_vertex, $all_vertices, $structure ) = @_;

    if ($structure->{$curr_vertex}->[0] == $curr_vertex){
        return $all_vertices;
    }
    else {
        push (@{$all_vertices}, $structure->{$curr_vertex}->[0]);
        save_main_chain_vertices_in_array(
        $structure->{$curr_vertex}->[0], $all_vertices , $structure);
    }
}

# Returns array that contains lengths of all side chains
sub find_lengths_of_side_chains
{
    my ( $graph, $curr_vertex, $all_vertices, $vertex_array, $structure ) = @_;
    if ($structure->{$curr_vertex}->[0] == $curr_vertex){
        return sort @{$vertex_array};
    }
    else {
        my @vertices = $graph->vertices();
        my @vertex = grep {$_->{number} == $curr_vertex} @vertices;
        my @curr_neighbours = $graph->neighbours( $vertex[0] );
        if (scalar @curr_neighbours == 1 ){
            $graph->delete_vertex($vertex[0]);
            find_lengths_of_side_chains(
                $graph,
                $curr_neighbours[0]->{number},
                $all_vertices,
                $vertex_array,
                $structure
            );
        }
        else {
            my @side_chain_neighbours = ();
            my $next_chain_vertex;
            foreach my $neigh (@curr_neighbours) {
                if (grep { $neigh->{number} eq $_ } @{$all_vertices}) {
                   $next_chain_vertex = $neigh;
                }
                else {
                    push (@side_chain_neighbours, $neigh);
                }
            }
            $graph->delete_vertex($vertex[0]);

            foreach my $neighbour (@side_chain_neighbours) {
                push(@{$vertex_array},
                        BFS_calculate_chain_length($graph, $neighbour));
            }

            find_lengths_of_side_chains(
                $graph,
                $next_chain_vertex->{number},
                $all_vertices,
                $vertex_array,
                $structure
            );
        }
    }
}

#TODO: Try to merge all subroutines

sub find_locant_placing
{
    my ( $graph, $main_chain, $structure ) = @_;

    my @vertices = $graph->vertices();
    my @places_of_locants = ();
    my @reverted_main_chain = reverse @{$main_chain};
    my $vertex_number = scalar @reverted_main_chain;

    for my $curr_vertex ( @reverted_main_chain ){
        my @vertex = grep {$_->{number} == $curr_vertex} @vertices;
        my @curr_neighbours = $graph->neighbours( $vertex[0] );
        if ( scalar @curr_neighbours == 0 ){
            return @places_of_locants;
        }
        elsif ( scalar @curr_neighbours == 1 ){
            $graph->delete_vertex($vertex[0]);
        }
        else{
            $graph->delete_vertex($vertex[0]);
            foreach my $neigh (@curr_neighbours) {
                if (grep { $neigh->{number} eq $_ } @{$main_chain}) {
                    next;
                }
                else{
                    push(@places_of_locants, $vertex_number);
                }
            }
        }
        $vertex_number -= 1;
    }
}

sub find_number_of_side_chains
{
    my ( $graph, $main_chain, $structure ) = @_;

    my @vertices = $graph->vertices();
    my $number_of_side_chains = 0;
    my @reverted_main_chain = reverse @{$main_chain};

    for my $curr_vertex ( @reverted_main_chain ){
        my @vertex = grep {$_->{number} == $curr_vertex} @vertices;
        my @curr_neighbours = $graph->neighbours( $vertex[0] );
        if ( scalar @curr_neighbours == 0 ){
            return $number_of_side_chains;
        }
        elsif ( scalar @curr_neighbours == 1 ){
            $graph->delete_vertex($vertex[0]);
        }
        else{
            $graph->delete_vertex($vertex[0]);
            foreach my $neigh (@curr_neighbours) {
                if (grep { $neigh->{number} eq $_ } @{$main_chain}) {
                    next;
                }
                else{
                    $number_of_side_chains += 1;
                }
            }
        }
    }
}

sub find_number_of_branched_side_chains
{
    my ( $graph, $main_chain, $structure ) = @_;

    my @vertices = $graph->vertices();
    my $number_of_branched_side_chains = 0;
    my @reverted_main_chain = reverse @{$main_chain};

    for my $curr_vertex ( @reverted_main_chain ){
        my @vertex = grep {$_->{number} == $curr_vertex} @vertices;
        my @curr_neighbours = $graph->neighbours( $vertex[0] );
        if ( scalar @curr_neighbours == 0 ){
            return $number_of_branched_side_chains;
        }
        if ( scalar @curr_neighbours == 1 ){
            $graph->delete_vertex($vertex[0]);
        }
        else{
            $graph->delete_vertex($vertex[0]);
            foreach my $neigh (@curr_neighbours) {
                if (grep { $neigh->{number} eq $_ } @{$main_chain}) {
                    next;
                }
                else{
                    $number_of_branched_side_chains +=
                                    BFS_is_chain_branched($graph, $neigh);
                }
            }
        }
    }
}

sub compare_locant_placings{
    my @first = @{$a->[3]};
    my @second = @{$b->[3]};
    my @index = (0..scalar @first-1);

    foreach(@index){
        if ($first[$_] > $second[$_]) {return 1}
        elsif ($first[$_] < $second[$_]) {return -1}
    }
    {return 0}
}

sub compare_side_chain_lengths{
    my @first = @{@{$a}[3]};
    my @second = @{@{$b}[3]};
    my @index = (0..scalar @first-1);

    foreach(@index){
        if ($first[$_] > $second[$_]) {return 1}
        elsif ($first[$_] < $second[$_]) {return -1}
    }
    {return 0}
}

1;
