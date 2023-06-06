package ChemOnomatopist;

use strict;
use warnings;

# ABSTRACT: Give molecule a name
# VERSION

use ChemOnomatopist::Chain;
use ChemOnomatopist::Chain::Circular;
use ChemOnomatopist::Chain::FromHalves;
use ChemOnomatopist::ChainHalf;
use ChemOnomatopist::Elements qw( %elements );
use ChemOnomatopist::Group;
use ChemOnomatopist::Group::Aldehyde;
use ChemOnomatopist::Group::Amino;
use ChemOnomatopist::Group::Carboxyl;
use ChemOnomatopist::Group::Cyanide;
use ChemOnomatopist::Group::Ester;
use ChemOnomatopist::Group::Hydroperoxide;
use ChemOnomatopist::Group::Hydroxy;
use ChemOnomatopist::Group::Imino;
use ChemOnomatopist::Group::Ketone;
use ChemOnomatopist::Group::Monocycle;
use ChemOnomatopist::Group::Monospiro;
use ChemOnomatopist::Group::Nitro;
use ChemOnomatopist::Group::Nitroso;
use ChemOnomatopist::Group::XO3;
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
use Chemistry::OpenSMILES qw( is_double_bond );
use Graph::Traversal::DFS;
use Graph::Undirected;
use List::Util qw( all any max min sum0 uniq );
use Scalar::Util qw( blessed );
use Set::Object qw( set );

no warnings 'recursion';

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

    # TODO: Extend to other subclasses of ChemOnomatopist::Group::
    if( blessed $start &&
        ( $start->isa( ChemOnomatopist::Group::Monocycle:: ) ||
          $start->isa( ChemOnomatopist::Group::Monospiro:: ) ) ) {
        return ChemOnomatopist::Name->new( $start->prefix );
    }

    $options = {} unless $options;

    # Groups that cannot be included in the chain do not matter
    my $branches_at_start = grep { !blessed $_ || $_->is_carbon }
                                 $graph->neighbours( $start );

    my $chain = select_sidechain( $graph, $start );
    my @chain = blessed $chain && $chain->can( 'vertices' ) ? $chain->vertices : $chain;

    # TODO: Handle the case when none of the rules select proper chains
    die "could not select a chain\n" unless @chain;

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

    $graph = copy $graph;
    $graph->delete_path( @chain );

    # Examine the attachments to the main chain: delete the edges
    # connecting them to the main chain, at the same time giving them
    # names according to their lengths via calls to get_sidechain_name()
    my %attachments;
    my %attachment_objects;
    for my $i (0..$#chain) {
        my $atom = $chain[$i];
        for my $neighbour ($graph->neighbours( $atom )) {
            my $bond = '-';
            if(         $graph->has_edge_attribute( $atom, $neighbour, 'bond' ) ) {
                $bond = $graph->get_edge_attribute( $atom, $neighbour, 'bond' );
            }
            $graph->delete_edge( $atom, $neighbour );

            my $attachment_name = get_sidechain_name( $graph, $neighbour );
            $attachment_name .= 'idene' if $bond eq '=';
            push @{$attachments{$attachment_name}}, $i;
            $attachment_objects{$attachment_name} = $attachment_name;
        }
    }

    # Collecting names of all the attachments
    my $name = ChemOnomatopist::Name->new;
    for my $attachment_name (sort { $a cmp $b } keys %attachments) {
        my $attachment = $attachment_objects{$attachment_name};

        if( @chain > 1 &&
            ( !blessed $chain ||
              !$chain->can( 'max_valence' ) ||
              scalar keys %attachments > 1 ||
              @{$attachments{$attachment_name}} != $chain->max_valence - 1 ) ) {
            if( blessed $chain && $chain->isa( ChemOnomatopist::Chain:: ) ) {
                $name->append_locants( $chain->locants( @{$attachments{$attachment_name}} ) );
            } else {
                $name->append_locants( map { $_ + 1 } @{$attachments{$attachment_name}} );
            }
        }

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
    $name .= unbranched_chain_name( blessed $chain && $chain->can( 'vertices' ) ? $chain : \@chain );
    $name->{name} =~ s/(an)?e$//; # FIXME: Dirty

    if( $branches_at_start > 1 ) {
        my( $branch_point ) = grep { $chain[$_] == $start } 0..$#chain;
        $name .= 'an' unless $name->{name} =~ /-en$/; # FIXME: Dirty
        $name->append_substituent_locant( $branch_point + 1 );
    }

    $name .= 'yl';
    $name->bracket if $name =~ /hydroxymethyl$/; # FIXME: Dirty

    return $name;
}

sub get_mainchain_name
{
    my( $graph, $chain, $options ) = @_;

    my @vertices = $graph->vertices;
    my @chain = $chain->vertices;
    my @groups = most_senior_groups( $graph->vertices );
    my $most_senior_group = blessed $groups[0] if @groups;

    # Disconnect the main chain: this way every main chain atom remains
    # connected only to the side chains.
    if( $chain->isa( ChemOnomatopist::Chain:: ) ) {
        $graph = copy $chain->graph;
        $graph->delete_edges( map { @$_ } $graph->subgraph( \@chain )->edges );
    } else {
        $graph = copy $graph;
        $graph->delete_path( @chain );
    }

    my %attachments;
    my %heteroatoms;
    my %attachment_objects;
    # Examine the main chain
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

    # Examine the attachments to the main chain: delete the edges
    # connecting them to the main chain, at the same time giving them
    # names according to their lengths via calls to get_sidechain_name()
    my @senior_group_attachments;
    for my $i (0..$#chain) {
        my $atom = $chain[$i];
        for my $neighbour ($graph->neighbours( $atom )) {
            my $bond = '-';
            if(         $graph->has_edge_attribute( $atom, $neighbour, 'bond' ) ) {
                $bond = $graph->get_edge_attribute( $atom, $neighbour, 'bond' );
            }
            $graph->delete_edge( $atom, $neighbour );

            if( $most_senior_group && blessed $neighbour && $neighbour->isa( $most_senior_group ) ) {
                push @senior_group_attachments, $i;
            } else {
                my $attachment_name = get_sidechain_name( $graph, $neighbour );
                $attachment_name .= 'idene' if $bond eq '=';

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

        # Locants are not important in single-substituted homogeneous cycles
        if( @order > 1 || @{$attachments{$attachment_name}} > 1 || $chain->needs_substituent_locants ) {
            $name->append_locants( $chain->locants( @{$attachments{$attachment_name}} ) );
        }

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

        if( $chain->needs_heteroatom_locants ) {
            $name->append_locants( $chain->locants( @{$heteroatoms{$element}} ) );
        }

        if( $chain->needs_heteroatom_names ) {
            if( @{$heteroatoms{$element}} > 1 ) {
                my $number = IUPAC_numerical_multiplier( scalar @{$heteroatoms{$element}} );
                $number .= 'a' unless $number =~ /^(|\?|.*i)$/;
                $name .= $number;
            }

            $name->append_element( $elements{$element}->{prefix} );
        }
    }

    if( $chain->isa( ChemOnomatopist::Group:: ) ) {
        $name .= $chain->suffix;
    } else {
        $name .= unbranched_chain_name( $chain );
    }

    if( $most_senior_group && !$most_senior_group->isa( ChemOnomatopist::Chain:: ) ) {
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
        if( $chain->needs_suffix_locant ) {
            $name->append_locants( $chain->locants( @senior_group_attachments ) );
        }
        $name->append_multiplier( $number );
        $name->append_suffix( @senior_group_attachments > 2 ? $groups[0]->multisuffix : $groups[0]->suffix );
    }

    $name =~ s/benzen-1-ol$/phenol/;
    $name = 'benzoic acid' if $name eq 'benzenoic acid';
    $name = 'toluene'      if $name eq 'methylbenzene';
    $name =~ s/^(\d,\d-)dimethylbenzene$/$1xylene/;

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

        # N-based groups
        if( is_element( $atom, 'N' ) && @neighbours == 3 && @C == 1 && @H == 2 ) {
            # Detecting amino
            my $amino = ChemOnomatopist::Group::Amino->new( @C );
            $graph->add_edge( @C, $amino );
            $graph->delete_vertices( $atom, @H );
        } elsif( is_element( $atom, 'N' ) && @neighbours == 2 && @C == 1 && @H == 1 &&
                 is_double_bond( $graph, $atom, @C ) ) {
            # Detecting imino
            my $imino = ChemOnomatopist::Group::Imino->new( @C );
            $graph->add_edge( @C, $imino );
            $graph->delete_vertices( $atom, @H );
        } elsif( is_element( $atom, 'N' ) && @C == 1 && @O == 2 && $atom->{charge} && $atom->{charge} == 1 &&
                 (any {  is_double_bond( $graph, $atom, $_ ) } @O) &&
                 (any { !is_double_bond( $graph, $atom, $_ ) && $_->{charge} && $_->{charge} == -1 } @O) ) {
            # Detecting nitro
            my $nitro = ChemOnomatopist::Group::Nitro->new( @C );
            $graph->add_edge( @C, $nitro );
            $graph->delete_vertices( $atom, @O );
        } elsif( is_element( $atom, 'N' ) && @neighbours == 1 && @C == 1 &&
                 is_triple_bond( $graph, $atom, @C ) ) {
            # Detecting cyanide
            my( $C ) = grep { is_element( $_, 'C' ) } $graph->neighbours( @C );
            my $cyanide = ChemOnomatopist::Group::Cyanide->new( $C );
            $graph->add_edge( $C, $cyanide );
            $graph->delete_vertices( $atom, @C );
        }

        # O-based groups
        if( is_element( $atom, 'O' ) && @neighbours == 2 && @C == 1 && @H == 1 ) {
            # Detecting hydroxy
            my $hydroxy = ChemOnomatopist::Group::Hydroxy->new( @C );
            $graph->add_edge( @C, $hydroxy );
            $graph->delete_vertices( $atom, @H );
        } elsif( is_element( $atom, 'O' ) && @neighbours == 2 && @H == 1 && @O == 1 ) {
            # Detecting hydroperoxide
            my @C = grep { is_element( $_, 'C' ) } $graph->neighbours( @O );
            if( @C == 1 ) {
                my $hydroperoxide = ChemOnomatopist::Group::Hydroperoxide->new( @C );
                $graph->add_edge( @C, $hydroperoxide );
                $graph->delete_vertices( $atom, @H, @O );
            }
        }

        # Ketones and their chalcogen analogues
        if( @neighbours == 1 && @C == 1 && is_double_bond( $graph, $atom, @C ) &&
            any { is_element( $atom, $_ ) } ( 'O', 'S', 'Se', 'Te' ) ) {
            my $ketone = ChemOnomatopist::Group::Ketone->new( @C, $atom );
            $graph->add_edge( @C, $ketone );
            $graph->delete_vertices( $atom );
        }

        # Nitroso and its analogues
        if( @neighbours == 2 && @C == 1 && @O == 1 && is_double_bond( $graph, $atom, @O ) &&
            any { is_element( $atom, $_ ) } qw( Br Cl F I N ) ) {
            my $nitroso = ChemOnomatopist::Group::Nitroso->new( @C, $atom );
            $graph->add_edge( @C, $nitroso );
            $graph->delete_vertices( $atom, @O );
        }

        # XO3
        if( @neighbours == 4 && @C == 1 && @O == 3 && (all { is_double_bond( $graph, $atom, $_ ) } @O) &&
            any { is_element( $atom, $_ ) } qw( Br Cl F I ) ) {
            my $XO3 = ChemOnomatopist::Group::XO3->new( @C, $atom );
            $graph->add_edge( @C, $XO3 );
            $graph->delete_vertices( $atom, @O );
        }
    }

    # Due to the issue in Graph, bridges() returns strings instead of real objects.
    # Graph issue: https://github.com/graphviz-perl/Graph/issues/29
    my %vertices_by_name = map { $_ => $_ } $graph->vertices;
    my $cut_vertices = set( map { $vertices_by_name{$_} } map { @$_ } $graph->bridges );

    # Second pass is needed to build on top of these trivial groups
    for my $atom ($graph->vertices) {
        my @neighbours = $graph->neighbours( $atom );
        my @groups = grep { blessed $_ && $_->isa( ChemOnomatopist::Group:: ) }
                          @neighbours;
        my @C = grep { is_element( $_, 'C' ) } @neighbours;
        my @H = grep { is_element( $_, 'H' ) } @neighbours;
        my @O = grep { is_element( $_, 'O' ) } @neighbours;

        if( is_element( $atom, 'C' ) && @groups == 1 && @H == 1 &&
            $groups[0]->isa( ChemOnomatopist::Group::Ketone:: ) ) {
            # Detecting aldehyde
            my $aldehyde = ChemOnomatopist::Group::Aldehyde->new( $atom );
            $graph->delete_vertices( @groups, @H ); # FIXME: Be careful!
            $graph->add_edges( map { $aldehyde, $_ } $graph->neighbours( $atom ) );
            $graph->delete_vertex( $atom );
        } elsif( is_element( $atom, 'C' ) && @C == 1 &&
            (any { $_->isa( ChemOnomatopist::Group::Ketone:: ) } @groups) &&
            (any { $_->isa( ChemOnomatopist::Group::Hydroxy::  ) } @groups) ) {
            # Detecting carboxyl
            my $carboxyl = ChemOnomatopist::Group::Carboxyl->new( @C );
            $graph->delete_vertices( $atom, @groups ); # FIXME: Be careful!
            $graph->add_edges( $carboxyl, @C );
        } elsif( is_element( $atom, 'C' ) && @groups == 1 && @C == 1 && @O == 2 &&
                 $groups[0]->isa( ChemOnomatopist::Group::Ketone:: ) &&
                 all { $cut_vertices->has( $_ ) } @O ) {
            # Detecting esters
            # Both oxygens have to be cut vertices to avoid one being in a ring
            $graph->delete_vertices( @groups );
            my( $hydroxylic ) = grep { $_ != $atom } map { $graph->neighbours( $_ ) } grep { !blessed $_ } @O;
            my( $acid ) = @C;
            my $ester = ChemOnomatopist::Group::Ester->new( $hydroxylic, $acid );
            $graph->add_edges( $ester, $hydroxylic );
            $graph->add_edges( $ester, $acid );
            $graph->delete_vertices( $atom, @O );
        }
    }

    # Hydrogen atoms are no longer important.
    # They are demoted to hydrogen counts for future reference.
    for my $H (grep { is_element( $_, 'H' ) } $graph->vertices) {
        my( $parent ) = $graph->neighbours( $H );
        $parent->{hcount} = 0 unless exists $parent->{hcount};
        $parent->{hcount}++;
        $graph->delete_vertex( $H );
    }

    # Detecting monocyclic compounds
    if( graph_has_cycle( $graph ) ) {
        # If it is not a tree, then the graph has cycles, and we have to do our best to recognise them.

        my $core = graph_cycle_core( $graph );
        my %vertices_by_degree;
        for my $vertex ($core->vertices) {
            my $degree = $core->degree( $vertex );
            $vertices_by_degree{$degree} = [] unless $vertices_by_degree{$degree};
            push @{$vertices_by_degree{$degree}}, $vertex;
        }

        # The cycle object is given the original graph to retain the original atom-atom relations
        my $compound;
        if(      join( ',', sort keys %vertices_by_degree ) eq '2' ) {
            $compound = ChemOnomatopist::Group::Monocycle->new( copy $graph, Graph::Traversal::DFS->new( $core )->dfs );
        } elsif( join( ',', sort keys %vertices_by_degree ) eq '2,4' && @{$vertices_by_degree{4}} == 1 ) {
            # BBv2 P-24.2.1 Monospiro alicyclic ring systems
            $compound = ChemOnomatopist::Group::Monospiro->new( copy $graph, $core->vertices );
        } elsif( join( ',', sort keys %vertices_by_degree ) eq '2,3' && @{$vertices_by_degree{3}} == 2 &&
                 $core->has_edge( @{$vertices_by_degree{3}} ) && $core->is_edge_connected ) {
            # Ortho-fused as defined in BBv2 P-25.3.1.1.1
            # Edge connected means that it has no bridges
            $compound = ChemOnomatopist::Group::Bicycle->new( copy $graph, $core->vertices );
        } else {
            die "cannot handle cyclic compounds other than monocycles and monospiro\n";
        }

        for my $neighbour ( map { $graph->neighbours( $_ ) } $core->vertices ) {
            $graph->add_edge( $compound, $neighbour );

            # Reattach groups
            if( blessed $neighbour &&
                $neighbour->isa( ChemOnomatopist::Group:: ) &&
                $core->has_vertex( $neighbour->C ) ) {
                $neighbour->{C} = $compound;
            }
        }
        $graph->delete_vertices( $core->vertices );
        $graph->delete_edge( $compound, $compound ); # May have been added, must be removed
    }

    return;
}

# Check if an object or Perl hash is of certain chemical element
sub is_element
{
    my( $atom, $element ) = @_;
    return unless ref $atom;

    $element = ucfirst $element;

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
           ucfirst $atom->{symbol} eq $element;
}

sub is_triple_bond
{
    my( $graph, $A, $B ) = @_;
    return $graph->has_edge_attribute( $A, $B, 'bond' ) &&
           $graph->get_edge_attribute( $A, $B, 'bond' ) eq '#';
}

# Given a graph, selects the main chain.
# The returned chain is an object of ChemOnomatopist::ChainHalf or its subclasses.
sub select_mainchain
{
    my( $graph ) = @_;

    # Find the most senior group, undefined if alkane
    my @groups = most_senior_groups( $graph->vertices );
    my $most_senior_group = blessed $groups[0] if @groups;

    my @chains;
    if( @groups ) {
        # TODO: Select a chain containing most of the senior groups
        my @carbons = uniq map { $_->C } @groups; # FIXME: Carbons with the most attachments should be preferred

        if( $most_senior_group->isa( ChemOnomatopist::Chain:: ) ) {
            return shift @groups;
        } elsif( @carbons == 1 ) {
            if( blessed $carbons[0] && @groups == 1 &&
                ( $carbons[0]->isa( ChemOnomatopist::Group::Monocycle:: ) ||
                  $carbons[0]->isa( ChemOnomatopist::Group::Monospiro:: ) ) ) {
                # For senior attachments to cycles
                push @chains, @carbons;
            } else {
                # As the starting position is known, it is enough to take the "side chain"
                # containing this particular carbon:
                my $chain = select_sidechain( $graph, @carbons );
                my @vertices = blessed $chain && $chain->can( 'vertices' ) ? $chain->vertices : $chain;
                push @chains, ChemOnomatopist::Chain->new( $graph, @vertices ),
                              ChemOnomatopist::Chain->new( $graph, reverse @vertices );
            }
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
                             ChemOnomatopist::Chain->new( $graph,
                                                          reverse( @{$longest_paths{$A}->[$i]} ),
                                                          @$path,
                                                          @{$longest_paths{$B}->[$j]} ),
                             ChemOnomatopist::Chain->new( $graph,
                                                          reverse( @{$longest_paths{$B}->[$j]} ),
                                                          reverse( @$path ),
                                                          @{$longest_paths{$A}->[$j]} );
                    }
                }
            }
            die "cannot determine the parent structure\n" unless @chains;

            @chains = rule_most_groups( $most_senior_group, @chains );
        }
    } else {
        # Here the candidate halves for the longest (and "best") path are placed in @path_parts.
        # Each of candidate halves start with center atom.
        my $subgraph = copy $graph;
        $subgraph->delete_vertices( grep { blessed $_ &&
                                           $_->isa( ChemOnomatopist::Group:: ) &&
                                           !$_->is_part_of_chain } $subgraph->vertices );
        my @center = graph_center( $subgraph );
        my @path_parts;
        if( @center == 1 ) {
            # Longest path has odd length
            for my $path ( graph_longest_paths_from_vertex( $subgraph, $center[0] ) ) {
                push @path_parts,
                     ChemOnomatopist::ChainHalf->new( $graph, undef, @$path );
            }
        } else {
            # Longest path has even length
            # Graph copy without center edge is required by graph_longest_paths_from_vertex()
            my $copy = copy $subgraph;
            $copy->delete_edge( @center );
            for my $vertex ( @center ) {
                push @path_parts,
                     map { ChemOnomatopist::ChainHalf->new( $graph, (grep { $_ ne $vertex } @center), @$_ ) }
                         graph_longest_paths_from_vertex( $copy, $vertex );
            }
        }

        return shift @path_parts if @path_parts == 1; # methane

        # Generate all possible chains.
        # FIXME: This needs optimisation.
        for my $part1 (@path_parts) {
            for my $part2 (@path_parts) {
                next if $part1->group eq $part2->group;
                push @chains, ChemOnomatopist::Chain::FromHalves->new( $part1, $part2 );
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
        !$chain->isa( ChemOnomatopist::Group:: ) &&
         $chain->number_of_groups( $most_senior_group ) ) {

        shift @vertices;
        pop @vertices;
        $chain = ChemOnomatopist::Chain->new( $graph, @vertices );
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
                push @chains, ChemOnomatopist::Chain::FromHalves->new( $part1, $part2 );
            }
        }
    } elsif( $C_graph->degree( $start ) == 1 ) {
        @chains = map { ChemOnomatopist::Chain->new( $graph, $_->vertices ) } @path_parts;
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

        return shift @chains;
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

                   # P-44.4.1.1: Maximum number of multiple bonds
                   \&rule_most_multiple_bonds,
                   # P-44.4.1.2: Maximum number of double bonds
                   \&rule_most_double_bonds,
                   # TODO: P-44.4.1.3: Nonstandard bonding numbers
                   # TODO: P-44.4.1.4: Concerns rings
                   # P-44.4.1.5: Lowest locants for heteroatoms in skeletal chain
                   \&rule_lowest_numbered_heteroatoms,
                   # P-44.4.1.6: Lowest locants for heteroatoms in skeletal chain according to heteroatom seniority
                   \&rule_lowest_numbered_most_senior_heteroatoms,
                   # TODO: P-44.4.1.7: Concerns rings
                   # P-44.4.1.8: Lowest locants for suffix groups
                   \&rule_lowest_numbered_senior_groups,
                   # TODO: P-44.4.1.9: Concerns rings
                   # TODO: P-44.4.1.10: Lowest locants for prefixes/suffixes expressing degrees of hydrogenation
                   #                    This is not fully implemented now
                   \&rule_lowest_numbered_multiple_bonds,
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

sub rule_lowest_numbered_multiple_bonds
{
    my @chains = @_;
    my( $max_value ) = sort { cmp_arrays( $a, $b ) }
                       map  {  [ $_->multiple_bond_positions ] } @chains;
    return grep { !cmp_arrays( [ $_->multiple_bond_positions ],
                               $max_value ) }
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

    my( $max_value ) = sort { cmp_heteroatom_counts( $a, $b ) }
                       map  { [ $_->heteroatoms ] } @chains;
    return grep { !cmp_heteroatom_counts( [ $_->heteroatoms ], $max_value ) } @chains;
}

sub rule_most_multiple_bonds
{
    my( @chains ) = @_;

    my( $max_value ) = sort { $b <=> $a }
                       map { $_->number_of_multiple_bonds } @chains;
    return grep { $_->number_of_multiple_bonds == $max_value } @chains;
}

sub rule_most_double_bonds
{
    my( @chains ) = @_;

    my( $max_value ) = sort { $b <=> $a }
                       map { $_->number_of_double_bonds } @chains;
    return grep { $_->number_of_double_bonds == $max_value } @chains;
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

# Chains given for this rule have the same number of heteroatoms at the same positions.
sub rule_lowest_numbered_most_senior_heteroatoms
{
    my( @chains ) = @_;

    my( $max_value ) = sort { cmp_heteroatom_seniority( [ $a->heteroatoms ],
                                                        [ $b->heteroatoms ] ) }
                            @chains;

    return grep { !cmp_heteroatom_seniority( [ $_->heteroatoms ],
                                             [ $max_value->heteroatoms ] ) }
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

sub most_senior_groups
{
    my( @vertices ) = @_;

    if( @vertices == 1 && blessed $vertices[0] && $vertices[0]->isa( Graph::Undirected:: ) ) {
        # Graph given instead of an array of vertices
        @vertices = $vertices[0]->vertices;
    }

    my @groups = grep { blessed $_ && $_->isa( ChemOnomatopist::Group:: ) && !$_->is_prefix_only } @vertices;
    return unless @groups;

    my( $most_senior_group ) = sort { ChemOnomatopist::Group::cmp( $a, $b ) } @groups;
    return grep { $_->isa( blessed $most_senior_group ) } @groups;
}

# Given two lists of heteroatoms, return the one with the most senior ones
sub cmp_heteroatom_counts
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

sub cmp_heteroatom_seniority
{
    my( $A, $B ) = @_;

    for (0..$#$A) {
        next unless $elements{$A->[$_]}->{seniority} <=> $elements{$B->[$_]}->{seniority};
        return      $elements{$A->[$_]}->{seniority} <=> $elements{$B->[$_]}->{seniority};
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

    if( $N < 10 ) {
        my $value = $prefix[$ones];
        $value =~ s/a$// unless $is_middle;
        return $value;
    }

    return 'dec'    . ($is_middle ? 'a' : '') if $N == 10;
    return 'undec'  . ($is_middle ? 'a' : '') if $N == 11;
    return IUPAC_numerical_multiplier( $ones, 1 ) . 'dec' . ($is_middle ? 'a' : '') if $N < 20;
    return 'icos'   . ($is_middle ? 'a' : '') if $N == 20;
    return IUPAC_numerical_multiplier( $ones, 1 ) . 'cos' . ($is_middle ? 'a' : '') if $N < 30;

    if( $N < 100 ) {
        return ($ones == 1 ? $prefix[$ones] : IUPAC_numerical_multiplier( $ones, 1 )) .
               IUPAC_numerical_multiplier( $tens, 1 ) .
               ($tens == 3 ? 'a' : '') .
               'cont' . ($is_middle ? 'a' : '');
    }

    if( $N < 1000 ) {
        my $prefix = int( $tens . $ones ) == 1 ? $prefix[$ones] : IUPAC_numerical_multiplier( int( $tens . $ones ), 1 );
        $prefix[1] = 'he';
        return $prefix . $prefix[$hundreds] . 'ct' . ($is_middle ? 'a' : '');
    }

    if( $N < 10000 ) {
        $prefix[0] = 'ki';
        return IUPAC_numerical_multiplier( int( $hundreds . $tens . $ones ), 1 ) .
               $prefix[$thousands] . 'li';
    }

    die "cannot generate IUPAC numerical multiplier for $N\n";
}

sub IUPAC_complex_numerical_multiplier
{
    my( $N ) = @_;

    my @multipliers = ( undef, '', 'bis', 'tris' );
    return $multipliers[$N] if $N < @multipliers;
    return IUPAC_numerical_multiplier( $N, 1 ) . 'kis';
}

sub alkane_chain_name($)
{
    my( $N ) = @_;

    my @names = qw( ? meth eth prop but );

    return $names[$N] if $N < @names;
    return IUPAC_numerical_multiplier( $N );
}

sub unbranched_chain_name($)
{
    my( $chain ) = @_;

    my @chain = blessed $chain ? $chain->vertices : @$chain;

    my $name = ChemOnomatopist::Name->new;
    $name->append_stem( alkane_chain_name scalar @chain );

    my( @double, @triple );
    if( blessed $chain ) {
        my @bonds = $chain->bonds;
        @double = grep { $bonds[$_] eq '=' } 0..$#bonds;
        @triple = grep { $bonds[$_] eq '#' } 0..$#bonds;
        if( @double ) {
            if( ( !$chain->isa( ChemOnomatopist::Chain::Circular:: ) && @chain > 2 )||
                @double > 1 || @triple ) {
                $name->append_locants( map { $_ + 1 } @double );
            }
            $name->append_multiplier( IUPAC_numerical_multiplier( scalar @double ) ) if @double > 1;
            $name .= 'en';
        }
        if( @triple ) {
            if( ( !$chain->isa( ChemOnomatopist::Chain::Circular:: ) && @chain > 2 ) ||
                @triple > 1 || @double ) {
                $name->append_locants( map { $_ + 1 } @triple );
            }
            $name->append_multiplier( IUPAC_numerical_multiplier( scalar @triple ) ) if @triple > 1;
            $name .= 'yn';
        }
    }
    $name .= 'an' unless @double || @triple;
    $name .= 'e';
    return $name;
}

1;
