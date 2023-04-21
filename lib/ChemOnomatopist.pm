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
use Graph::Traversal::BFS;
use Graph::Undirected;
use List::Util qw( all any max min sum0 uniq );
use Scalar::Util qw( blessed );
use Set::Object qw( set );

no warnings 'recursion';

our @numbers = ( '?', '', 'di', 'tri', 'tetra', 'penta',
                 'hexa', 'hepta', 'octa', 'nona', 'deca',
                 'undeca', 'dodeca', 'trideca', 'tetradeca',
                 'pentadeca', 'hexadeca', 'heptadeca', 'octadeca', 'nonadeca',
                 'icosa', 'henicosa', 'docosa', 'tricosa',
                 'tetracosa', 'pentacosa', 'hexacosa',
                 'heptacosa', 'octacosa', 'nonacosa', 'triaconta',
                 'hentriaconta', 'dotriaconta', 'tritriaconta', 'tetratriaconta',
                 'pentatriaconta', 'hexatriaconta', 'heptatriaconta',
                 'octatriaconta', 'nonatriaconta', 'tetraconta' );

our @numberskis = ( '?', '', 'bis', 'tris', 'tetrakis', 'pentakis',
                    'hexakis', 'heptakis', 'octakis', 'nonakis', 'decakis' );

our %preferrable_names = (
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

    my $order = [ map { $_->{number} } select_main_chain( $graph->copy ) ];
    return get_mainchain_name( $graph->copy, $order ) . 'ane';
}

# get_sidechain_name() receives a graph and a position to start the chain in it.
# From that position it finds the longest chain and returns the constructed name.
sub get_sidechain_name
{
    my( $graph, $start, $options ) = @_;

    $options = {} unless $options;

    # FIXME: If the chain branches at $start, it should instead be named separately
    # and then marked as attached to the parent chain on this atom.

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

sub get_mainchain_name
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

        if( @chain > 3 ) {
            # As propane cannot have other attachment site but 2, the site is usually omitted
            $name .= join( ',', map { $_ + 1 } @{$attachments{$attachment_name}} ) . '-';
        }
        $name .= $number . $attachment_name;
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

# BFS is performed for the given graph after all vertices that are not carbons
# removed
sub BFS_order_carbons_only
{
    my( $graph, $start ) = @_;

    $graph = $graph->copy;
    $graph->delete_vertices( grep { !is_element( $_, 'C' ) } $graph->vertices );

    my $bfs;
    if( $start ) {
        $bfs = Graph::Traversal::BFS->new( $graph, start => $start );
    } else {
        $bfs = Graph::Traversal::BFS->new( $graph );
    }
    return $bfs->bfs;
}

# Selects the main chain by evaluating its parts
sub select_main_chain
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
                 ChemOnomatopist::ChainHalf->new( $tree, undef, @$path );
        }
    } else {
        # Longest path has even length
        # Graph copy without center edge is required by graph_longest_paths_from_vertex()
        my $copy = $tree->copy;
        $copy->delete_edge( @center );
        for my $vertex ( @center ) {
            push @path_parts,
                 map { ChemOnomatopist::ChainHalf->new( $tree, (grep { $_ ne $vertex } @center), @$_ ) }
                     graph_longest_paths_from_vertex( $copy, $vertex );
        }
    }

    return $path_parts[0]->vertices if @path_parts == 1; # methane

    # Generate all possible chains.
    # FIXME: This needs optimisation.
    my @chains;
    for my $part1 (@path_parts) {
        for my $part2 (@path_parts) {
            next if $part1->group eq $part2->group;
            push @chains, ChemOnomatopist::Chain->new( $part1, $part2 );
        }
    }

    for my $rule ( sub { return @_ },
                   \&rule_greatest_number_of_side_chains, # After this rule we are left with a set of longest chains all having the same number of side chains
                   \&rule_lowest_numbered_locants,
                   \&rule_most_carbon_in_side_chains,
                   \&rule_least_branched_side_chains,
                   \&pick_chain_with_lowest_attachments_alphabetically ) {
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

sub rule_greatest_number_of_side_chains
{
    my( @chains ) = @_;

    my( $max_value ) = sort { $b <=> $a }
                       uniq map { $_->number_of_branches }
                                @chains;
    return grep { $_->number_of_branches == $max_value } @chains;
}

sub rule_lowest_numbered_locants
{
    my( @chains ) = @_;

    my( $max_value ) = sort { cmp_arrays( [ $a->branch_positions ],
                                          [ $b->branch_positions ] ) }
                            @chains;
    return grep { !cmp_arrays( [ $_->branch_positions ],
                               [ $max_value->branch_positions ] ) }
                @chains;
}

sub rule_most_carbon_in_side_chains
{
    my( @chains ) = @_;

    my( $max_value ) = sort { $b <=> $a }
                       uniq map { $_->number_of_carbons }
                                @chains;
    return grep { $_->number_of_carbons == $max_value } @chains;
}

sub rule_least_branched_side_chains
{
    my( @chains ) = @_;

    my( $min_value ) = sort uniq map { $_->number_of_branches_in_sidechains } @chains;
    return grep { $_->number_of_branches_in_sidechains == $min_value } @chains;
}

sub pick_chain_with_lowest_attachments_alphabetically
{
    my( @chains ) = @_;

    my @locant_names = map { [ $_->locant_names ] } @chains;
    my @sorted = reverse sort { cmp_attachments( $locant_names[$a], $locant_names[$b] ) }
                              0..$#locant_names;
    return $chains[$sorted[0]];
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

# Sorts arrays from lowest to biggest by values
sub cmp_arrays
{
    my( $a, $b ) = @_;
    my @first  = @$a;
    my @second = @$b;
    my @index  = (0..scalar @first-1);

    foreach( @index ) {
        return $first[$_] <=> $second[$_] if $first[$_] <=> $second[$_];
    }

    return 0;
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

        my @A = ref $a_alpha ? sort @$a_alpha : ( $a_alpha );
        my @B = ref $b_alpha ? sort @$b_alpha : ( $b_alpha );

        for (0..min( scalar( @A ), scalar( @B ) )-1) {
            my $a_alpha = $A[$_];
            my $b_alpha = $B[$_];

            $a_alpha =~ s/[^a-zA-Z]+//g;
            $b_alpha =~ s/[^a-zA-Z]+//g;

            $a_alpha = 'butyl' if $a_alpha eq 'tertbutyl';
            $b_alpha = 'butyl' if $b_alpha eq 'tertbutyl';

            return $b_alpha cmp $a_alpha if $b_alpha cmp $a_alpha;
        }
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
