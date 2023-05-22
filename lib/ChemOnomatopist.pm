package ChemOnomatopist;

use strict;
use warnings;

# ABSTRACT: Give molecule a name
# VERSION

use ChemOnomatopist::Chain;
use ChemOnomatopist::Chain::Circular;
use ChemOnomatopist::Chain::VertexArray;
use ChemOnomatopist::ChainHalf;
use ChemOnomatopist::Elements qw( %elements );
use ChemOnomatopist::Group;
use ChemOnomatopist::Group::Aldehyde;
use ChemOnomatopist::Group::Amino;
use ChemOnomatopist::Group::Carbonyl;
use ChemOnomatopist::Group::Carboxyl;
use ChemOnomatopist::Group::Ester;
use ChemOnomatopist::Group::Hydroperoxide;
use ChemOnomatopist::Group::Hydroxy;
use ChemOnomatopist::Group::Imino;
use ChemOnomatopist::Name;
use ChemOnomatopist::Util qw( copy );
use ChemOnomatopist::Util::Graph qw(
    BFS_calculate_chain_length
    BFS_is_chain_branched
    graph_center
    graph_cycle_core
    graph_has_cycle
    graph_longest_paths_from_vertex
    graph_path_between_vertices
    tree_branch_positions
    tree_number_of_branches
);
use Graph::Traversal::DFS;
use Graph::Undirected;
use List::Util qw( all any max min sum0 uniq );
use Scalar::Util qw( blessed );
use Set::Object qw( set );

no warnings 'recursion';

sub wjoin(@);

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

    find_groups( $graph );

    my $main_chain = select_mainchain( $graph );
    return get_mainchain_name( $graph, $main_chain );
}

# get_sidechain_name() receives a graph and a position to start the chain in it.
# From that position it finds the longest chain and returns the constructed name.
sub get_sidechain_name
{
    my( $graph, $start, $options ) = @_;

    $options = {} unless $options;

    # Groups that cannot be included in the chain do not matter
    my $branches_at_start = grep { !blessed $_ || $_->is_carbon }
                                 $graph->neighbours( $start );

    my @chain = select_sidechain( $graph, $start );

    # TODO: Handle the case when none of the rules select proper chains
    die "could not select a chain\n" unless @chain;

    # TODO: Bond orders are not handled yet
    for (0..$#chain-1) {
        next unless $graph->has_edge_attributes( $chain[$_], $chain[$_ + 1] );
        die "cannot handle such compounds for now\n";
    }

    # Handle non-carbon substituents
    if( @chain == 1 && !is_element( $chain[0], 'C' ) ) {
        if( blessed $chain[0] ) {
            return ChemOnomatopist::Name->new( $chain[0]->prefix );
        } elsif( exists $elements{$chain[0]->{symbol}} ) {
            my $element = $elements{$chain[0]->{symbol}}->{prefix};
            $element =~ s/a$/o/; # TODO: Is this a general rule? BBv2 seems silent.
            return ChemOnomatopist::Name->new( $element );
        }
    }

    $graph->delete_path( @chain );

    # Examine the attachments to the main chain: delete the edges
    # connecting them to the main chain, at the same time giving them
    # names according to their lengths via calls to get_sidechain_name()
    my %attachments;
    my %attachment_objects;
    for my $i (0..$#chain) {
        my $atom = $chain[$i];
        for my $neighbour ($graph->neighbours( $atom )) {
            $graph->delete_edge( $atom, $neighbour );

            my $attachment_name = get_sidechain_name( $graph, $neighbour );
            push @{$attachments{$attachment_name}}, $i;
            $attachment_objects{$attachment_name} = $attachment_name;
        }
    }

    # Collecting names of all the attachments
    my $name = ChemOnomatopist::Name->new;
    for my $attachment_name (sort { $a cmp $b } keys %attachments) {
        my $attachment = $attachment_objects{$attachment_name};

        $name->append_locants( map { $_ + 1 } @{$attachments{$attachment_name}} ) if @chain > 1;

        if( @{$attachments{$attachment_name}} > 1 ) {
            my $number = IUPAC_numerical_multiplier( scalar @{$attachments{$attachment_name}} );
            $number .= 'a' unless $number =~ /^(|\?|.*i)$/;
            $name->append_multiplier( $number );

            # FIXME: More rules from BBv2 P-16.3.4 should be added
            if( $attachment->has_substituent_locant || # BBv2 P-16.3.4 (a)
                $attachment->starts_with_multiplier || # BBv2 P-16.3.4 (c)
                $attachment =~ /^dec/ ||               # BBv2 P-16.3.4 (d)
                $attachment =~ /^[0-9]/ ) {
                $attachment->bracket;
            }
        }

        $name .= $attachment;
    }
    $name .= alkane_chain_name( scalar @chain );

    if( $branches_at_start > 1 ) {
        my( $branch_point ) = grep { $chain[$_] == $start } 0..$#chain;
        $name .= 'an';
        $name->append_substituent_locant( $branch_point + 1 );
    }

    $name .= 'yl';
    $name->bracket if $name =~ /hydroxymethyl$/; # FIXME: Ugly fix

    return $name;
}

sub get_mainchain_name
{
    my( $graph, $chain, $options ) = @_;

    my @vertices = $graph->vertices;
    my @chain = blessed $chain ? $chain->vertices : @$chain;
    my $most_senior_group = most_senior_group( $graph );

    # Disconnect the main chain: this way every main chain atom remains
    # connected only to the side chains.
    $graph = copy $graph;
    if( blessed $chain && $chain->isa( ChemOnomatopist::Chain::Circular:: ) ) {
        $graph->delete_cycle( @chain );
    } else {
        # TODO: Bond orders are not handled yet
        for (0..$#chain-1) {
            next unless $graph->has_edge_attributes( $chain[$_], $chain[$_ + 1] );
            die "cannot handle such compounds for now\n";
        }
        $graph->delete_path( @chain );
    }

    my %attachments;
    my %heteroatoms;
    my %attachment_objects;
    # Examine the main chain
    # This is skipped for cycles as they should be capable of naming themselves
    if( !blessed $chain || !$chain->isa( ChemOnomatopist::Chain::Circular:: ) ) {
        for my $i (0..$#chain) {
            my $atom = $chain[$i];
            if( blessed $atom ) {
                next if $most_senior_group && $atom->isa( $most_senior_group );
                push @{$attachments{$atom->prefix}}, $i;
                $attachment_objects{$atom->prefix} = ChemOnomatopist::Name->new( $atom->prefix );
            } elsif( !is_element( $atom, 'C' ) &&
                     exists $atom->{symbol} &&
                     exists $elements{$atom->{symbol}} ) {
                push @{$heteroatoms{$atom->{symbol}}}, $i;
            }
        }
    }

    # Examine the attachments to the main chain: delete the edges
    # connecting them to the main chain, at the same time giving them
    # names according to their lengths via calls to get_sidechain_name()
    my @senior_group_attachments;
    for my $i (0..$#chain) {
        my $atom = $chain[$i];
        for my $neighbour ($graph->neighbours( $atom )) {
            $graph->delete_edge( $atom, $neighbour );

            if( $most_senior_group && blessed $neighbour && $neighbour->isa( $most_senior_group ) ) {
                push @senior_group_attachments, $i;
            } else {
                my $attachment_name = get_sidechain_name( $graph, $neighbour );
                push @{$attachments{$attachment_name}}, $i;
                $attachment_objects{$attachment_name} = $attachment_name;
            }
        }
    }

    # Collecting names of all the attachments
    my @order = sort { cmp_only_aphabetical( $a, $b ) || $a cmp $b } keys %attachments;
    my $name = ChemOnomatopist::Name->new;
    for my $i (0..$#order) {
        my $attachment_name = $order[$i];
        my $attachment = $attachment_objects{$attachment_name};

        $name->append_locants( map { $_ + 1 } @{$attachments{$attachment_name}} );

        # FIXME: More rules from BBv2 P-16.3.4 should be added
        if( $attachment !~ /^[\(\[\{]/ &&
            ( $attachment->starts_with_multiplier || # BBv2 P-16.3.4 (c)
              $attachment =~ /^[0-9]/ ) ) {
              $attachment->bracket;
        }

        if( @{$attachments{$attachment_name}} > 1 ) {
            my $number;
            if( $attachment =~ /^[\(\[\{]/ ) {
                $number = IUPAC_complex_numerical_multiplier( scalar @{$attachments{$attachment_name}} );
            } else {
                $number = IUPAC_numerical_multiplier( scalar @{$attachments{$attachment_name}} );
                $number .= 'a' unless $number =~ /^(|\?|.*i)$/;
            }
            $name .= $number;

            # BBv2 P-16.3.4 (a)
            if( $attachment !~ /^[\(\[\{]/ &&
                ( $attachment =~ /^dec/ || # BBv2 P-16.3.4 (d)
                  $attachment->has_substituent_locant ) ) {
                $attachment->bracket;
            }
        }

        $name .= $attachment;
    }

    # Collecting names of all heteroatoms
    for my $element (sort { $elements{$a}->{seniority} <=> $elements{$b}->{seniority} }
                          keys %heteroatoms) {

        $name->append_locants( map { $_ + 1 } @{$heteroatoms{$element}} ) if @chain > 1;

        if( @{$heteroatoms{$element}} > 1 ) {
            my $number = IUPAC_numerical_multiplier( scalar @{$heteroatoms{$element}} );
            $number .= 'a' unless $number =~ /^(|\?|.*i)$/;
            $name .= $number;
        }

        $name->append_element( $elements{$element}->{prefix} );
    }

    if( blessed $chain && $chain->can( 'name' ) ) {
        $name .= $chain->name;
    } else {
        $name .= alkane_chain_name( scalar @chain ) . 'ane';
    }

    if( $most_senior_group ) {
        if( $most_senior_group->is_carbon ) {
            # Most senior group is carbon, thus it is in the chain as well
            my @senior_group_positions =
                grep { blessed $chain[$_] && $chain[$_]->isa( $most_senior_group ) }
                     0..$#chain;
            @senior_group_attachments = sort { $a <=> $b } @senior_group_attachments,
                                                           @senior_group_positions;
        }

        my $number = IUPAC_numerical_multiplier( scalar @senior_group_attachments );
        $number = '' if $number eq 'mono';
        $number .= 'a' unless $number =~ /^(|\?|.*i)$/;

        # Terminal locants are not cited for 1 or 2 senior group attachments according to BBv2 P-14.3.4.1
        if( @senior_group_attachments && @chain > 1 &&
            (!$most_senior_group->is_carbon || @senior_group_attachments > 2) ) {
            $name->append_locants( map { $_ + 1 } @senior_group_attachments );
        }
        $name = wjoin $name,
                      $number,
                      ( @senior_group_attachments > 2 ? $most_senior_group->multisuffix : $most_senior_group->suffix );
    }

    return $name;
}

sub find_groups
{
    my( $graph ) = @_;

    for my $atom ($graph->vertices) {
        my @neighbours = $graph->neighbours( $atom );
        my @C = grep { is_element( $_, 'C' ) } @neighbours;
        my @H = grep { is_element( $_, 'H' ) } @neighbours;
        my @O = grep { is_element( $_, 'O' ) } @neighbours;

        # Detecting amino
        if( is_element( $atom, 'N' ) && @neighbours == 3 && @C == 1 && @H == 2 ) {
            my $amino = ChemOnomatopist::Group::Amino->new( @C );
            $graph->add_edge( @C, $amino );
            $graph->delete_vertices( $atom, @H );
        } elsif( is_element( $atom, 'N' ) && @neighbours == 2 && @C == 1 && @H == 1 ) {
            # Detecting imino
            # FIXME: Check also for double bond
            my $imino = ChemOnomatopist::Group::Imino->new( @C );
            $graph->add_edge( @C, $imino );
            $graph->delete_vertices( $atom, @H );
        }

        # Detecting carbonyl
        if( is_element( $atom, 'O' ) && @neighbours == 1 && @C == 1 ) {
            my $carbonyl = ChemOnomatopist::Group::Carbonyl->new( @C );
            $graph->add_edge( @C, $carbonyl );
            $graph->delete_vertices( $atom );
        # Detecting hydroxy
        } elsif( is_element( $atom, 'O' ) && @neighbours == 2 && @C == 1 && @H == 1 ) {
            my $hydroxy = ChemOnomatopist::Group::Hydroxy->new( @C );
            $graph->add_edge( @C, $hydroxy );
            $graph->delete_vertices( $atom, @H );
        # Detecting hydroperoxide
        } elsif( is_element( $atom, 'O' ) && @neighbours == 2 && @H == 1 && @O == 1 ) {
            my @C = grep { is_element( $_, 'C' ) } $graph->neighbours( @O );
            if( @C == 1 ) {
                my $hydroperoxide = ChemOnomatopist::Group::Hydroperoxide->new( @C );
                $graph->add_edge( @C, $hydroperoxide );
                $graph->delete_vertices( $atom, @H, @O );
            }
        }
    }

    # Second pass is needed to build on top of these trivial groups
    for my $atom ($graph->vertices) {
        my @neighbours = $graph->neighbours( $atom );
        my @groups = grep { blessed $_ && $_->isa( ChemOnomatopist::Group:: ) }
                          @neighbours;
        my @H = grep { is_element( $_, 'H' ) } @neighbours;

        # Detecting aldehyde
        if( is_element( $atom, 'C' ) && @groups == 1 && @H == 1 &&
            $groups[0]->isa( ChemOnomatopist::Group::Carbonyl:: ) ) {
            my $aldehyde = ChemOnomatopist::Group::Aldehyde->new( $atom );
            $graph->delete_vertices( @groups, @H ); # FIXME: Be careful!
            $graph->add_edges( map { $aldehyde, $_ } $graph->neighbours( $atom ) );
            $graph->delete_vertex( $atom );
        # Detecting carboxyl
        } elsif( is_element( $atom, 'C' ) &&
            (any { $_->isa( ChemOnomatopist::Group::Carbonyl:: ) } @groups) &&
            (any { $_->isa( ChemOnomatopist::Group::Hydroxy::  ) } @groups) ) {
            my $carboxyl = ChemOnomatopist::Group::Carboxyl->new( $atom );
            $graph->delete_vertices( @groups ); # FIXME: Be careful!
            $graph->add_edges( map { $carboxyl, $_ } $graph->neighbours( $atom ) );
            $graph->delete_vertex( $atom );
        }
    }

    # Hydrogen atoms are no longer important
    $graph->delete_vertices( grep { is_element( $_, 'H' ) } $graph->vertices );

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
            } elsif( $element eq 'H' ) {
                return '';
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

# Selects the main chain by evaluating its parts
sub select_mainchain
{
    my( $graph ) = @_;

    # Find the most senior group, undefined if alkane
    my $most_senior_group = most_senior_group( $graph );

    my @chains;
    if( graph_has_cycle( $graph ) ) {
        # If it is not a tree, than the graph has cycles, and we have to
        # do our best to recognise them. To make it easier, hydrogen atoms
        # are removed here for now.

        my $core = graph_cycle_core( $graph );
        if( any { $core->degree( $_ ) > 2 } $core->vertices ) {
            die "cannot handle cyclic compounds other than monocycles\n";
        }

        # FIXME: For now we generate all possible traversals of the same cycle.
        #        This is not optimal, some caching could be introduced.
        my @vertices = Graph::Traversal::DFS->new( $core )->dfs;
        for (0..$#vertices) {
            push @chains,
                 ChemOnomatopist::Chain::Circular->new( $graph, @vertices );
            push @vertices, shift @vertices;
        }
        @vertices = reverse @vertices;
        for (0..$#vertices) {
            push @chains,
                 ChemOnomatopist::Chain::Circular->new( $graph, @vertices );
            push @vertices, shift @vertices;
        }
    } elsif( $most_senior_group ) {
        # TODO: Select a chain containing most of the senior groups
        my @groups = grep { blessed( $_ ) && $_->isa( $most_senior_group ) } $graph->vertices;
        my @carbons = uniq map { $_->C } @groups; # FIXME: Carbons with the most attachments should be preferred

        if( @carbons == 1 ) {
            # As the starting position is known, it is enough to take the "side chain"
            # containing this particular carbon:
            my @vertices = select_sidechain( $graph, @carbons );
            push @chains, ChemOnomatopist::Chain::VertexArray->new( $graph, @vertices ),
                          ChemOnomatopist::Chain::VertexArray->new( $graph, reverse @vertices );
        } else {
            my @paths;
            my $max_value;
            for my $i (0..$#carbons) {
                for my $j (($i+1)..$#carbons) {
                    my @path = graph_path_between_vertices( $graph, $carbons[$i], $carbons[$j] );
                    my $value = (set( @carbons ) * set( @path ))->size;
                    if(      !defined $max_value || $max_value < $value ) {
                        @paths = ( \@path );
                        $max_value = $value;
                    } elsif( $max_value == $value ) {
                        push @paths, \@path;
                    }
                }
            }

            # Construct all chains having all possible extensions to both sides of the selected path
            my %longest_paths;
            for my $path (@paths) {
                my $copy = copy $graph;
                $copy->delete_path( @$path );
                $copy->delete_vertices( grep { blessed $_ && !is_element( $_, 'C' ) } $copy->vertices );

                my $A = shift @$path;
                my $B = pop @$path;

                if( !exists $longest_paths{$A} ) {
                    $longest_paths{$A} = [ graph_longest_paths_from_vertex( $copy, $A ) ];
                }
                if( !exists $longest_paths{$B} ) {
                    $longest_paths{$B} = [ graph_longest_paths_from_vertex( $copy, $B ) ];
                }

                for my $i (0..$#{$longest_paths{$A}}) {
                    for my $j (0..$#{$longest_paths{$B}}) {
                        push @chains,
                             ChemOnomatopist::Chain::VertexArray->new( $graph,
                                                                       reverse( @{$longest_paths{$A}->[$i]} ),
                                                                       @$path,
                                                                       @{$longest_paths{$B}->[$j]} ),
                             ChemOnomatopist::Chain::VertexArray->new( $graph,
                                                                       reverse( @{$longest_paths{$B}->[$j]} ),
                                                                       reverse( @$path ),
                                                                       @{$longest_paths{$A}->[$j]} );
                    }
                }
            }

            @chains = rule_most_groups( $most_senior_group, @chains );
        }
    } else {
        # Here the candidate halves for the longest (and "best") path are placed in @path_parts.
        # Each of candidate halves start with center atom.
        my @center = graph_center( $graph );
        my @path_parts;
        if( @center == 1 ) {
            # Longest path has odd length
            for my $path ( graph_longest_paths_from_vertex( $graph, $center[0] ) ) {
                push @path_parts,
                     ChemOnomatopist::ChainHalf->new( $graph, undef, @$path );
            }
        } else {
            # Longest path has even length
            # Graph copy without center edge is required by graph_longest_paths_from_vertex()
            my $copy = copy $graph;
            $copy->delete_edge( @center );
            for my $vertex ( @center ) {
                push @path_parts,
                     map { ChemOnomatopist::ChainHalf->new( $graph, (grep { $_ ne $vertex } @center), @$_ ) }
                         graph_longest_paths_from_vertex( $copy, $vertex );
            }
        }

        return @path_parts if @path_parts == 1; # methane

        # Generate all possible chains.
        # FIXME: This needs optimisation.
        for my $part1 (@path_parts) {
            for my $part2 (@path_parts) {
                next if $part1->group eq $part2->group;
                push @chains, ChemOnomatopist::Chain->new( $part1, $part2 );
            }
        }
    }

    my $chain = filter_chains( @chains );
    my @vertices = $chain->vertices;

    # If there is at least one of carbon-based senior group attachment,
    # it means both ends are already senior, prompting to follow the
    # exception of three or more carbon-based groups.
    if( $most_senior_group &&
        $most_senior_group->is_carbon &&
        $chain->number_of_groups( $most_senior_group ) ) {
        shift @vertices;
        pop @vertices;
        $chain = ChemOnomatopist::Chain::VertexArray->new( $graph, @vertices );
    }

    return $chain;
}

# Selects the best side chain
sub select_sidechain
{
    my( $graph, $start ) = @_;

    # Do this for non-carbons for now in order to represent attachments
    return $start unless is_element( $start, 'C' );

    # Cleaning the graph from the heteroatom leaves
    my $C_graph = copy $graph;
    $C_graph->delete_vertices( grep { blessed $_ && !is_element( $_, 'C' ) } $C_graph->vertices );
    while( my @leaves = grep { !is_element( $_, 'C' ) && $C_graph->degree( $_ ) == 1 } $C_graph->vertices ) {
        $C_graph->delete_vertices( @leaves );
    }

    my @path_parts;
    for my $neighbour ($C_graph->neighbours( $start )) {
        my $graph_copy = copy $C_graph;
        $graph_copy->delete_edge( $start, $neighbour );
        for my $path ( graph_longest_paths_from_vertex( $graph_copy, $neighbour ) ) {
            push @path_parts,
                 ChemOnomatopist::ChainHalf->new( $graph, undef, $start, @$path );
        }
    }

    my @chains;
    if(      $C_graph->degree( $start ) > 1 ) {
        # FIXME: Deduplicate: copied from select_mainchain()
        # Generate all possible chains.
        # FIXME: This needs optimisation.
        for my $part1 (@path_parts) {
            for my $part2 (@path_parts) {
                next if $part1->group eq $part2->group;
                push @chains, ChemOnomatopist::Chain->new( $part1, $part2 );
            }
        }
    } elsif( $C_graph->degree( $start ) == 1 ) {
        @chains = map { ChemOnomatopist::Chain::VertexArray->new( $graph, $_->vertices ) } @path_parts;
    } else {
        return ( $start );
    }

    # From BBv2 P-29.2
    my $rule_lowest_free_valence = sub {
        my( @chains ) = @_;

        my @chains_now;
        my $lowest_locant;

        for my $chain (@chains) {
            my @vertices = $chain->vertices;
            my( $locant ) = grep { $vertices[$_] == $start } 0..$#vertices;
            if( @chains_now ) {
                if( $lowest_locant > $locant ) {
                    @chains_now = ( $chain );
                    $lowest_locant = $locant;
                } elsif( $lowest_locant == $locant ) {
                    push @chains_now, $chain;
                }
            } else {
                @chains_now = ( $chain );
                $lowest_locant = $locant;
            }
        }

        return @chains_now;
    };

    for my $rule ( sub { return @_ },
                   \&rule_longest_chains,
                   \&rule_greatest_number_of_side_chains, # After this rule we are left with a set of longest chains all having the same number of side chains
                   $rule_lowest_free_valence,
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
        last;
    }

    # TODO: Handle the case when none of the rules select proper chains
    return ();
}

# TODO: Should reflect the order described in BBv2 P-44
sub filter_chains
{
    my( @chains ) = @_;

    for my $rule ( sub { return @_ },
                   # P-44.1.1: Maximum number of substituents of principal characteristic group.
                   #           This is not needed as select_mainchain() returns such chains.

                   # TODO: P-44.1.2: Concerns rings

                   # P-44.3.1: Maximum number of heteroatoms of any kind
                   \&rule_most_heteroatoms,
                   # P-44.3.2: Maximum number of skeletal atoms
                   \&rule_longest_chains,
                   # P-44.3.3: Maximum number of the most senior heteroatom
                   \&rule_most_senior_heteroatoms,

                   # TODO: P-44.4.1.1: Maximum number of multiple bonds
                   # TODO: P-44.4.1.2: Maximum number of double bonds
                   # TODO: P-44.4.1.3: Nonstandard bonding numbers
                   # TODO: P-44.4.1.4: Concerns rings
                   # P-44.4.1.5: Lowest locants for heteroatoms in skeletal chain
                   \&rule_lowest_numbered_heteroatoms,
                   # TODO: P-44.4.1.6: Lowest locants for heteroatoms in skeletal chain according to heteroatom seniority
                   # TODO: P-44.4.1.7: Concerns rings
                   # P-44.4.1.8: Lowest locants for suffix groups
                   \&rule_lowest_numbered_senior_groups,
                   # TODO: P-44.4.1.9: Concerns rings
                   # TODO: P-44.4.1.10: Lowest locants for prefixes/suffixes expressing degrees of hydrogenation
                   # TODO: P-44.4.1.11: Concerns isotopes
                   # TODO: P-44.4.1.12: Concerns stereogenic centers

                   # TODO: P-45.1: Multiplication of identical senior parent structures

                   # P-45.2.1: Maximum number of prefix substituents
                   #           FIXME: This includes suffix substituents now
                   \&rule_greatest_number_of_side_chains,
                   # P-45.2.2: Lowest locants for prefix substituents
                   #           FIXME: This includes suffix substituents now
                   \&rule_lowest_numbered_locants,
                   # TODO: P-45.2.3: Lowest locants for prefix substituents in their order of citation in the name
                   # TODO: P-45.3: Nonstandard bond numbers
                   # TODO: P-45.4: Concerns isotopes
                   # TODO: P-45.5: Alphanumerical order of names (maybe covered by P-45.2.3 already?)
                   # TODO: P-45.6: Concerns stereochemistry

                   # TODO: Put these in correct order:
                   \&rule_most_carbon_in_side_chains,
                   \&rule_least_branched_side_chains,
                   \&pick_chain_with_lowest_attachments_alphabetically ) {
        my @chains_now = $rule->( @chains );

        if( 0 ) {
            require Sub::Identify;
            print STDERR '>>> ', Sub::Identify::sub_name( $rule ), "\n";
        }

        # CHECK: Can a rule cause disappearance of all chains?
        next unless @chains_now;

        @chains = @chains_now; # Narrow down the selection

        # If a single chain cannot be chosen now, pass on to the next rule
        return shift @chains if @chains == 1;
    }

    # TODO: Handle the case when none of the rules select proper chains
    return ();
}

sub rule_most_groups
{
    my( $class, @chains ) = @_;

    my( $max_value ) = sort { $b <=> $a }
                       map { $_->number_of_groups( $class ) } @chains;
    return grep { $_->number_of_groups( $class ) == $max_value } @chains;
}

sub rule_lowest_numbered_senior_groups
{
    my( @chains ) = @_;
    my( $max_value ) = sort { cmp_arrays( [ $a->most_senior_group_positions ],
                                          [ $b->most_senior_group_positions ] ) }
                            @chains;
    return grep { !cmp_arrays( [ $_->most_senior_group_positions ],
                               [ $max_value->most_senior_group_positions ] ) }
                @chains;
}

# This rule is employed only if longest chains are not already preselected
sub rule_longest_chains
{
    my( @chains ) = @_;

    my( $max_value ) = sort { $b <=> $a }
                       uniq map { $_->length } @chains;
    return grep { $_->length == $max_value } @chains;
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

sub rule_most_heteroatoms
{
    my( @chains ) = @_;

    my( $max_value ) = sort { $b <=> $a }
                       map { $_->number_of_heteroatoms } @chains;
    return grep { $_->number_of_heteroatoms == $max_value } @chains;
}

sub rule_most_senior_heteroatoms
{
    my( @chains ) = @_;

    my( $max_value ) = sort { cmp_heteroatoms( $a, $b ) }
                       map  { [ $_->heteroatoms ] } @chains;
    return grep { !cmp_heteroatoms( [ $_->heteroatoms ], $max_value ) } @chains;
}

sub rule_lowest_numbered_heteroatoms
{
    my( @chains ) = @_;

    my( $max_value ) = sort { cmp_arrays( [ $a->heteroatom_positions ],
                                          [ $b->heteroatom_positions ] ) }
                            @chains;
    return grep { !cmp_arrays( [ $_->heteroatom_positions ],
                               [ $max_value->heteroatom_positions ] ) }
                @chains;
}

sub pick_chain_with_lowest_attachments_alphabetically
{
    my( @chains ) = @_;

    my @locant_names = map { [ $_->locant_names ] } @chains;
    my @sorted = sort { cmp_attachments( $locant_names[$a], $locant_names[$b] ) }
                      0..$#locant_names;
    return $chains[$sorted[0]];
}

sub most_senior_group
{
    my( @vertices ) = @_;

    if( @vertices == 1 && blessed $vertices[0] && $vertices[0]->isa( Graph::Undirected:: ) ) {
        # Graph given instead of an array of vertices
        @vertices = $vertices[0]->vertices;
    }

    # Find the most senior group, undefined if alkane
    for my $group (@ChemOnomatopist::Group::order) {
        next unless any { blessed $_ && $_->isa( $group ) } @vertices;
        return $group;
    }

    return;
}

# Given two lists of heteroatoms, return the one with the most senior ones
sub cmp_heteroatoms
{
    my( $A, $B ) = @_;

    my( %A, %B );
    for (@$A) { $A{$_}++ }
    for (@$B) { $B{$_}++ }

    my @elements = sort { $elements{$a}->{seniority} <=> $elements{$b}->{seniority} }
                   grep { $elements{$_}->{seniority} >= 5 } # O and after
                        keys %elements;
    for (@elements) {
        $A{$_} = 0 unless $A{$_};
        $B{$_} = 0 unless $B{$_};

        return $B{$_} <=> $A{$_} if $B{$_} <=> $A{$_};
    }

    return 0;
}

# Sorts given names only based on alphabetical part of the name.
# tert compounds are ordered according to BBv2 P-14.5 which says:
# "[t]he preferred order for alphanumerical order is: nonitalic Roman letters > italic letters > Greek letters."
sub cmp_only_aphabetical
{
    my( $a, $b ) = @_;

    $a =~ s/[^a-zA-Z]+//g;
    $b =~ s/[^a-zA-Z]+//g;

    my $a_has_tert = $a =~ s/^tert(butyl)$/$1/;
    my $b_has_tert = $b =~ s/^tert(butyl)$/$1/;

    return $a_has_tert <=> $b_has_tert if $a_has_tert <=> $b_has_tert;
    return $a cmp $b;
}

# Sorts arrays from lowest to biggest by values
sub cmp_arrays
{
    my( $a, $b ) = @_;

    for (0..min( scalar( @$a ), scalar( @$b ) )-1) {
        return $a->[$_] <=> $b->[$_] if $a->[$_] <=> $b->[$_];
    }

    return @$a <=> @$b;
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

        my @A = ref $a_alpha eq 'ARRAY' ? sort @$a_alpha : ( $a_alpha );
        my @B = ref $b_alpha eq 'ARRAY' ? sort @$b_alpha : ( $b_alpha );

        for (0..min( scalar( @A ), scalar( @B ) )-1) {
            my $a_alpha = $A[$_];
            my $b_alpha = $B[$_];

            $a_alpha =~ s/[^a-zA-Z]+//g;
            $b_alpha =~ s/[^a-zA-Z]+//g;

            my $a_has_tert = $a_alpha =~ s/^tert(butyl)$/$1/;
            my $b_has_tert = $b_alpha =~ s/^tert(butyl)$/$1/;

            return $a_has_tert <=> $b_has_tert if $a_has_tert <=> $b_has_tert;
            return $a_alpha cmp $b_alpha if $a_alpha cmp $b_alpha;
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

    if( $N < 1000 ) {
        my $prefix = int( $tens . $ones ) == 1 ? 'hen' : IUPAC_numerical_multiplier( int( $tens . $ones ), 1 );
        return $prefix . 'hect' if $N < 200;
        return $prefix . $prefix[$hundreds] . 'ct' . ($is_middle ? 'a' : '');
    }

    return IUPAC_numerical_multiplier( int( $hundreds . $tens . $ones ), 1 ) . 'kili'                     if $N <  2000;
    return IUPAC_numerical_multiplier( int( $hundreds . $tens . $ones ), 1 ) . $prefix[$thousands] . 'li' if $N < 10000;
    die "cannot generate IUPAC numerical multiplier for $N\n";
}

sub IUPAC_complex_numerical_multiplier
{
    my( $N ) = @_;

    my @multipliers = ( undef, '', 'bis', 'tris' );
    return $multipliers[$N] if $N < @multipliers;
    return IUPAC_numerical_multiplier( $N, 1 ) . 'kis';
}

sub alkane_chain_name
{
    my( $N ) = @_;

    my @names = qw( ? meth eth prop but );

    return $names[$N] if $N < @names;
    return IUPAC_numerical_multiplier( $N );
}

sub wjoin(@)
{
    my( @parts ) = grep { $_ ne '' } @_;

    for (0..(@parts-2)) {
        next if $parts[$_] eq 'di' || $parts[$_] eq 'tri';
        $parts[$_] =~ s/[aeiou](-[0-9,]+-|)$/$1/ if $parts[$_ + 1] =~ /^[aeiou]/;
    }

    return join '', @parts;
}

1;
