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
use ChemOnomatopist::Group::AcylHalide;
use ChemOnomatopist::Group::Aldehyde;
use ChemOnomatopist::Group::Amide;
use ChemOnomatopist::Group::Amine;
use ChemOnomatopist::Group::Carboxyl;
use ChemOnomatopist::Group::Cyanide;
use ChemOnomatopist::Group::Ester;
use ChemOnomatopist::Group::Guanidine;
use ChemOnomatopist::Group::Hydrazine;
use ChemOnomatopist::Group::Hydroperoxide;
use ChemOnomatopist::Group::Hydroxy;
use ChemOnomatopist::Group::Imino;
use ChemOnomatopist::Group::Ketone;
use ChemOnomatopist::Group::Monocycle;
use ChemOnomatopist::Group::Monospiro;
use ChemOnomatopist::Group::Nitro;
use ChemOnomatopist::Group::Nitroso;
use ChemOnomatopist::Group::Polyacene;
use ChemOnomatopist::Group::Polyaphene;
use ChemOnomatopist::Group::SulfinicAcid;
use ChemOnomatopist::Group::Sulfinyl;
use ChemOnomatopist::Group::SulfonicAcid;
use ChemOnomatopist::Group::Sulfonyl;
use ChemOnomatopist::Group::XO3;
use ChemOnomatopist::Group::Xanthene;
use ChemOnomatopist::MolecularGraph;
use ChemOnomatopist::Name;
use ChemOnomatopist::Name::Part::AlkaneANSuffix;
use ChemOnomatopist::Name::Part::Stem;
use ChemOnomatopist::Util qw( copy );
use ChemOnomatopist::Util::Graph qw(
    BFS_calculate_chain_length
    BFS_is_chain_branched
    cyclic_components
    graph_center
    graph_cycle_core
    graph_cycles
    graph_has_cycle
    graph_longest_paths_from_vertex
    graph_path_between_vertices
    graph_replace
    graph_without_edge_attributes
    subgraph
    tree_branch_positions
    tree_number_of_branches
);
use ChemOnomatopist::Util::SMILES qw( cycle_SMILES_explicit );
use Chemistry::OpenSMILES qw(
    %bond_symbol_to_order
    is_double_bond
    is_ring_atom
    is_single_bond
    is_triple_bond
);
use Graph::Nauty qw( are_isomorphic );
use Graph::SSSR;
use Graph::Traversal::DFS;
use Graph::Undirected;
use List::Util qw( all any max min sum0 uniq );
use Scalar::Util qw( blessed );
use Set::Object qw( set );

no warnings 'recursion';

our $CAUTIOUS = '';
our $DEBUG = '';

sub get_name
{
    my( $what ) = @_;

    # Detect the type of the input data
    my( $graph );
    if( blessed $what && $what->isa( Graph::Undirected:: ) ) {
        $graph = ChemOnomatopist::MolecularGraph->new( $what );
    } else {
        # Assume SMILES string
        die "cannot handle stereochemistry now\n" if $what =~ /[\\\/@]/;

        require Chemistry::OpenSMILES::Parser;
        my $parser = Chemistry::OpenSMILES::Parser->new;
        my @graphs = map { ChemOnomatopist::MolecularGraph->new( $_ ) }
                         $parser->parse( $what );
        die "separate molecular entities are not handled yet\n" if @graphs > 1;
        $graph = shift @graphs;
    }
    die "nothing supplied for get_name()\n" unless $graph;

    if( any { $_ ne 'H' && !exists $elements{$_} }
        map { ucfirst $_->{symbol} } $graph->vertices ) {
        die "unknown elements detected\n";
    }

    find_groups( $graph );

    my $main_chain = select_mainchain( $graph );
    return get_mainchain_name( $graph, $main_chain );
}

# get_sidechain_name() receives a graph and a position to start the chain in it.
# From that position it finds the longest chain and returns the constructed name.
sub get_sidechain_name
{
    my( $graph, $parent, $start ) = @_;

    # Record the type of parent bond
    my $parent_bond = '-' if $parent;
    if(                $graph->has_edge_attribute( $parent, $start, 'bond' ) ) {
        $parent_bond = $graph->get_edge_attribute( $parent, $start, 'bond' );
    }

    # Groups that cannot be included in the chain do not matter
    my $branches_at_start = grep { !blessed $_ || $_->is_carbon }
                            grep { !$parent || $_ != $parent }
                                 $graph->neighbours( $start );

    my $chain;
    if( blessed $start && $start->isa( ChemOnomatopist::Chain:: ) ) {
        $chain = $start;
        $chain->parent( $parent ) if $parent;
    } elsif( $graph->groups( $start ) ) {
        ( $chain ) = $graph->groups( $start );
        $chain->parent( $parent ) if $parent;
    } else {
        $chain = select_sidechain( $graph, $parent, $start );
    }
    my @chain = $chain->vertices;

    # Handle non-carbon substituents
    if( @chain == 1 && $graph->degree( @chain ) == 0 + defined $parent && !blessed $chain[0] &&
        !is_element( $chain[0], 'C' ) && exists $elements{$chain[0]->{symbol}} ) {
        my $element = $elements{$chain[0]->{symbol}}->{prefix};
        $element =~ s/a$/o/; # TODO: Is this a general rule? BBv2 seems silent.
        return ChemOnomatopist::Name::Part::Element->new( $element )->to_name;
    }

    # Examine the attachments to the main chain: delete the edges
    # connecting them to the main chain, at the same time giving them
    # names according to their lengths via calls to get_sidechain_name()
    my %attachments;
    my %attachment_objects;
    my %heteroatoms;
    my %isotopes;
    for my $i (0..$#chain) {
        my $atom = $chain[$i];

        if( !blessed $atom && !is_element( $atom, 'C' ) &&
            $atom->{symbol} && exists $elements{$atom->{symbol}} ) {
            push @{$heteroatoms{$atom->{symbol}}}, $i;
        }

        if( !blessed $atom && exists $atom->{isotope} ) {
            push @{$isotopes{$atom->{isotope} . element( $atom )}}, $i;
        }

        if( exists $atom->{h_isotope} && grep { defined $_ } @{$atom->{h_isotope}} ) {
            for (grep { defined $_ } @{$atom->{h_isotope}}) {
                push @{$isotopes{$_ . 'H'}}, $i;
            }
        }

        for my $neighbour ($graph->neighbours( $atom )) {
            next if any { $neighbour == $_ } @chain; # Skip atoms from this chain
            next if $parent && $neighbour == $parent;
            my $attachment_name = get_sidechain_name( $graph, $atom, $neighbour );
            push @{$attachments{$attachment_name}}, $i;
            $attachment_objects{$attachment_name} = $attachment_name;
        }
    }

    # Collecting names of all the attachments
    my $name = ChemOnomatopist::Name->new;
    for my $attachment_name (sort { $a cmp $b } keys %attachments) {
        my $attachment = $attachment_objects{$attachment_name};

        if( $chain->needs_substituent_locants ) {
            $name->append_locants( $chain->locants( @{$attachments{$attachment_name}} ) );
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
        } else {
            if( $chain->needs_substituent_locants &&
                !$attachment->is_enclosed &&
                (!$attachment->is_simple || $attachment->starts_with_locant) ) {
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

    # Attaching isotopes
    my $isotopes = '';
    for my $isotope (sort { cmp_isotopes( $a, $b ) } keys %isotopes) {
        if( $chain->needs_substituent_locants ) { # FIXME: Is this right?
            $isotopes .= join( ',', $chain->locants( @{$isotopes{$isotope}} ) ) . '-';
        }

        $isotopes .= $isotope;
        $isotopes .= scalar @{$isotopes{$isotope}} if @{$isotopes{$isotope}} > 1 || $isotope =~ /H$/;
    }
    $name .= "($isotopes)" if $isotopes ne '';

    if( $chain->isa( ChemOnomatopist::Group:: ) ) {
        my $prefix = $chain->prefix( $parent );
        # All groups are most likely stems
        $prefix = ChemOnomatopist::Name::Part::Stem->new( $prefix )->to_name unless blessed $prefix;
        $name .= $prefix;
    } elsif( @chain == 1 && blessed $chain[0] ) {
        my $prefix = $chain[0]->prefix( $parent );
        # All group-containing chains are most likely stems
        $prefix = ChemOnomatopist::Name::Part::Stem->new( $prefix )->to_name unless blessed $prefix;
        if( $chain[0]->isa( ChemOnomatopist::Group::Sulfinyl:: ) ) { # BBv2 P-63.6
            $name->[-1] =~ s/yl$/ane/;
        }
        $name .= $prefix;
    } else {
        if( $chain->isa( ChemOnomatopist::Chain::Ether:: ) && $name->has_locant ) {
            $name->bracket;
        }
        $name .= $chain->prefix( $parent );
        pop @$name if $name->[-1] eq 'e'; # FIXME: Dirty
        pop @$name if $name->[-1] eq 'an';

        if( $branches_at_start > 1 ) {
            my( $branch_point ) = grep { $chain[$_] == $start } 0..$#chain;
            if( $branch_point || !$chain->is_saturated ) {
                # According to BBv2 P-29.2 (1)
                $name .= 'an' unless $name =~ /-(di|tri)?en$/; # FIXME: Dirty
                $name->append_substituent_locant( $branch_point + 1 );
            }
        }

        $name .= 'yl' unless $name =~ /y$/;
        $name->bracket if $name =~ /hydroxymethyl$/; # FIXME: Dirty
    }

    if( @chain > 1 || !blessed $chain[0] ) { # Groups converted to chains can add suffixes themselves
        $name .= 'idene' if $parent_bond && $parent_bond eq '=';
        $name .= 'idyne' if $parent_bond && $parent_bond eq '#';
    }

    $name = ChemOnomatopist::Name->new( 'benzyl' ) if $name eq 'phenylmethyl';

    return $name;
}

sub get_mainchain_name
{
    my( $graph, $chain ) = @_;

    my @vertices = $graph->vertices;
    my @chain = $chain->vertices;
    my @groups = most_senior_groups( $graph );
    my $most_senior_group = blessed $groups[0] if @groups;

    # Collect the heteroatoms and isotopes in the chain
    my %heteroatoms;
    my %isotopes;
    for my $i (0..$#chain) {
        my $atom = $chain[$i];
        next if blessed $atom;

        if( !is_element( $atom, 'C' ) &&
            defined element( $atom ) &&
            exists $elements{element( $atom )} ) {
            push @{$heteroatoms{element( $atom )}}, $i;
        }

        if( exists $atom->{isotope} ) {
            push @{$isotopes{$atom->{isotope} . element( $atom )}}, $i;
        }

        if( exists $atom->{h_isotope} && grep { defined $_ } @{$atom->{h_isotope}} ) {
            for (grep { defined $_ } @{$atom->{h_isotope}}) {
                push @{$isotopes{$_ . 'H'}}, $i;
            }
        }
    }

    $graph = copy $chain->graph;

    # Collect the substituents
    my %attachments;
    my %attachment_objects;
    my @senior_group_attachments;
    for my $i (0..$#chain) {
        my $atom = $chain[$i];
        for my $neighbour ($graph->neighbours( $atom )) {
            next if any { $_ == $neighbour } @chain; # Skip atoms from this chain
            my( $group ) = $graph->groups( $neighbour );
            if( grep { $_ == $neighbour } @groups ) {
                push @senior_group_attachments, $i;
            } else {
                my $attachment_name = get_sidechain_name( $graph, $atom, $group ? $group : $neighbour );
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

        if( $chain->needs_substituent_locants ) {
            $name->append_locants( $chain->locants( @{$attachments{$attachment_name}} ) );
        }

        # FIXME: More rules from BBv2 P-16.3.4 and P-16.5.1 should be added
        if( !$attachment->is_enclosed &&
            ( $attachment->starts_with_multiplier || # BBv2 P-16.3.4 (c)
              $attachment =~ /^[0-9]/ ) ) {
              $attachment->bracket;
        }

        if( @{$attachments{$attachment_name}} > 1 ) {
            my $number;
            if( $attachment->is_enclosed ) {
                $number = IUPAC_complex_numerical_multiplier( scalar @{$attachments{$attachment_name}} );
            } else {
                $number = IUPAC_numerical_multiplier( scalar @{$attachments{$attachment_name}} );
                $number .= 'a' unless $number =~ /^(|\?|.*i)$/;
            }
            $name .= $number;

            # BBv2 P-16.3.4 (a)
            if( !$attachment->is_enclosed &&
                ( $attachment =~ /^dec/ || # BBv2 P-16.3.4 (d)
                  $attachment->has_substituent_locant ) ) {
                $attachment->bracket;
            }
        } else {
            # This is an attempt to implement rules from P-16.5.1.
            # However, they are quite vague, thus there is not much of guarantee the following code is correct.
            if( !$attachment->is_enclosed &&
                ($attachment->has_locant || !$attachment->is_simple) &&
                $chain->needs_substituent_locants &&
                $attachment ne 'tert-butyl' ) {
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

    # Attaching isotopes
    my $isotopes = '';
    for my $isotope (sort { cmp_isotopes( $a, $b ) } keys %isotopes) {
        if( $chain->needs_substituent_locants ) { # FIXME: Is this right?
            $isotopes .= join( ',', $chain->locants( @{$isotopes{$isotope}} ) ) . '-';
        }

        $isotopes .= $isotope;
        $isotopes .= scalar @{$isotopes{$isotope}} if @{$isotopes{$isotope}} > 1 || $isotope =~ /H$/;
    }
    $name .= "($isotopes)" if $isotopes ne '';

    $name .= $chain->suffix;

    if( $most_senior_group && !$most_senior_group->isa( ChemOnomatopist::Chain:: ) ) {
        if( $groups[0]->is_carbon ) {
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
        if( $chain->isa( ChemOnomatopist::Group::Monocycle:: ) ||
            $chain->isa( ChemOnomatopist::Group::Monospiro:: ) ) {
            $name->append_suffix( $groups[0]->suffix_if_cycle_substituent );
        } elsif( @senior_group_attachments > 2 ) {
            $name->append_suffix( $groups[0]->multisuffix );
        } else {
            $name->append_suffix( $groups[0]->suffix );
        }
    }

    $name =~ s/benzen(-1-)?ol$/phenol/;
    $name = 'anisole'      if $name eq 'methoxybenzene';         # BBv2 P-63.2.4.1
    $name = 'benzoic acid' if $name eq 'benzenecarboxylic acid'; # BBv2 P-65.1.1.1
    $name = 'toluene'      if $name eq 'methylbenzene';
    $name =~ s/^(\d,\d-)dimethylbenzene$/$1xylene/;
    $name = 'formic acid' if $name eq 'methanoic acid';
    $name =~ s/ethanoic acid$/acetic acid/; # BBv2 P-65.1.1.1
    $name =~ s/phenylamine$/aniline/;

    return $name;
}

sub find_groups
{
    my( $graph ) = @_;

    # First pass is to detect guanidine and hydrazine
    for my $atom ($graph->vertices) {
        my @neighbours = $graph->neighbours( $atom );
        my @H = grep { is_element( $_, 'H' ) } @neighbours;
        my @N = grep { is_element( $_, 'N' ) } @neighbours;

        if( is_element( $atom, 'C' ) && @neighbours == 3 && @N == 3 &&
            !is_ring_atom( $graph, $atom, -1 ) ) {
            # Detecting guanidine
            my $guanidine = ChemOnomatopist::Group::Guanidine->new( copy $graph, $atom );
            graph_replace_all( $graph, $guanidine, $atom, @N );
        }

        if( is_element( $atom, 'N' ) && @neighbours == 3 && @N == 1 && @H == 2 &&
            all { !is_ring_atom( $graph, $_, -1 ) } ( $atom, @N ) ) {
            # Detecting hydrazine
            my $hydrazine = ChemOnomatopist::Group::Hydrazine->new( $graph, $atom, @N );
            $graph->add_group( $hydrazine );
        }
    }

    for my $atom ($graph->vertices) {
        next if $graph->groups( $atom );

        my @neighbours = $graph->neighbours( $atom );
        my @C  = grep { is_element( $_, 'C' ) } @neighbours;
        my @H  = grep { is_element( $_, 'H' ) } @neighbours;
        my @N  = grep { is_element( $_, 'N' ) } @neighbours;
        my @O  = grep { is_element( $_, 'O' ) } @neighbours;
        my @S  = grep { is_element( $_, 'S' ) } @neighbours;
        my @Se = grep { is_element( $_, 'Se' ) } @neighbours;
        my @Te = grep { is_element( $_, 'Te' ) } @neighbours;

        # N-based groups
        if( is_element( $atom, 'N' ) && @C == 1 && @O == 2 && $atom->{charge} && $atom->{charge} == 1 &&
                 (any {  is_double_bond( $graph, $atom, $_ ) } @O) &&
                 (any { !is_double_bond( $graph, $atom, $_ ) && $_->{charge} && $_->{charge} == -1 } @O) ) {
            # Detecting nitro
            my $nitro = ChemOnomatopist::Group::Nitro->new( @C );
            graph_replace_all( $graph, $nitro, $atom, @O );
        } elsif( is_element( $atom, 'N' ) && @neighbours == 3 && !is_ring_atom( $graph, $atom, -1 ) ) {
            # Detecting amines
            my $amine = ChemOnomatopist::Group::Amine->new( @C );
            graph_replace_all( $graph, $amine, $atom, @H );
        } elsif( is_element( $atom, 'N' ) && @neighbours == 2 && @C == 1 && @H == 1 &&
                 is_double_bond( $graph, $atom, @C ) ) {
            # Detecting imino
            my $imino = ChemOnomatopist::Group::Imino->new( @C );
            graph_replace_all( $graph, $imino, $atom, @H );
        } elsif( is_element( $atom, 'N' ) && @neighbours == 1 && @C == 1 &&
                 $graph->degree( @C ) >= 2 &&
                 is_triple_bond( $graph, $atom, @C ) ) {
            # Detecting cyanide
            my( $C ) = grep { $_ != $atom } $graph->neighbours( @C );
            my $cyanide = ChemOnomatopist::Group::Cyanide->new( $C );
            graph_replace_all( $graph, $cyanide, $atom, @C );
        }

        # Hydroxy groups and their chalcogen analogues
        if( @neighbours == 2 && ( @C || @N || @O || @S || @Se || @Te ) && @H == 1 &&
            any { is_element( $atom, $_ ) } qw( O S Se Te ) ) {
            my $hydroxy = ChemOnomatopist::Group::Hydroxy->new( @C, @N, @O, @S, @Se, @Te, $atom );
            graph_replace_all( $graph, $hydroxy, $atom, @H );
        }

        # Ketones and their chalcogen analogues
        if( @neighbours == 1 && @C == 1 && is_double_bond( $graph, $atom, @C ) &&
            any { is_element( $atom, $_ ) } qw( O S Se Te ) ) {
            my $ketone = ChemOnomatopist::Group::Ketone->new( @C, $atom );
            graph_replace_all( $graph, $ketone, $atom );
        }

        # Nitroso and its analogues
        if( @neighbours == 2 && @C == 1 && @O == 1 && is_double_bond( $graph, $atom, @O ) &&
            any { is_element( $atom, $_ ) } qw( Br Cl F I N ) ) {
            my $nitroso = ChemOnomatopist::Group::Nitroso->new( @C, $atom );
            graph_replace_all( $graph, $nitroso, $atom, @O );
        }

        # XO3
        if( @neighbours == 4 && @C == 1 && @O == 3 && (all { is_double_bond( $graph, $atom, $_ ) } @O) &&
            any { is_element( $atom, $_ ) } qw( Br Cl F I ) ) {
            my $XO3 = ChemOnomatopist::Group::XO3->new( @C, $atom );
            graph_replace_all( $graph, $XO3, $atom, @O );
        }

        # Sulfinyl group and its analogues
        if( @neighbours == 3 && @O == 1 && is_double_bond( $graph, $atom, @O ) &&
            !is_ring_atom( $graph, $atom ) && any { is_element( $atom, $_ ) } qw( S Se Te ) ) {
            my $sulfinyl = ChemOnomatopist::Group::Sulfinyl->new( $atom );
            graph_replace_all( $graph, $sulfinyl, $atom, @O );
        }

        # Sulfonyl group and its analogues
        if( @neighbours == 4 && @O == 2 && (all { is_double_bond( $graph, $atom, $_ ) } @O) &&
            !is_ring_atom( $graph, $atom ) && any { is_element( $atom, $_ ) } qw( S Se Te ) ) {
            my $sulfonyl = ChemOnomatopist::Group::Sulfonyl->new( $atom );
            graph_replace_all( $graph, $sulfonyl, $atom, @O );
        }
    }

    if( any { !blessed $_ && exists $_->{charge} } $graph->vertices ) {
        die "cannot handle charges for now\n";
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
            $groups[0]->isa( ChemOnomatopist::Group::Ketone:: ) &&
            all { $graph->degree( $_ ) == 1 } @H ) {
            # Detecting aldehyde
            my $aldehyde = ChemOnomatopist::Group::Aldehyde->new( $atom, @groups );
            $graph->delete_vertices( @groups, @H );
            $graph->add_edges( map { $aldehyde, $_ } $graph->neighbours( $atom ) );
            $graph->delete_vertex( $atom );
        } elsif( is_element( $atom, 'C' ) && @neighbours == 3 && @groups >= 1 && @O == 2 &&
            (any { $_->isa( ChemOnomatopist::Group::Ketone:: ) } @groups) &&
            ( (any { $_->isa( ChemOnomatopist::Group::Hydroxy:: ) } @groups) ||
              (all { $graph->degree( $_ ) == 1 } @O) ) ) {
            # Detecting carboxyl
            my( $parent ) = grep { !is_element( $_, 'O' ) && !is_element( $_, 'H' ) } @neighbours;
            my $carboxyl = ChemOnomatopist::Group::Carboxyl->new( $parent );
            $graph->delete_vertices( $atom, @O );
            $graph->add_vertex( $carboxyl );
            $graph->add_edges( $carboxyl, $parent ) if $parent;
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
        } elsif( is_element( $atom, 'C' ) && @groups == 2 &&
                 (any { $_->isa( ChemOnomatopist::Group::Amine:: ) } @groups) &&
                 (any { $_->isa( ChemOnomatopist::Group::Ketone:: ) } @groups) ) {
            # Detecting amides
            my $amide = ChemOnomatopist::Group::Amide->new( $atom );
            my( $amine ) = grep { $_->isa( ChemOnomatopist::Group::Amine:: ) } @groups;
            my( $ketone ) = grep { $_->isa( ChemOnomatopist::Group::Ketone:: ) } @groups;
            $graph->delete_vertices( $ketone );
            graph_replace_all( $graph, $amide, $amine );
        } elsif( is_element( $atom, 'C' ) && @neighbours == 3 && @C == 1 &&
                 @groups == 1 && $groups[0]->isa( ChemOnomatopist::Group::Ketone:: ) && is_element( @groups, 'O' ) &&
                 element(   grep { !blessed $_ && !is_element( $_, 'C' ) } @neighbours ) =~ /^(F|Cl|Br|I)$/ ) {
            my( $halide ) = grep { !blessed $_ && !is_element( $_, 'C' ) } @neighbours;
            my $acyl_halide = ChemOnomatopist::Group::AcylHalide->new( @C, $halide );
            graph_replace_all( $graph, $acyl_halide, $atom, grep { blessed $_ || !is_element( $_, 'C' ) } @neighbours );
        }

        if( !blessed $atom && is_element( $atom, 'N' ) &&
            @neighbours - @H >= 2 && !is_ring_atom( $graph, $atom, -1 ) &&
            !$graph->groups( $atom ) ) {
            die "cannot process secondary and tertiary amines yet\n";
        }

        # Detecting sulfinic acids
        if( !blessed $atom && is_element( $atom, 'S' ) && @neighbours == 3 && @C == 1 && @O == 2 &&
            @groups == 1 && $groups[0]->isa( ChemOnomatopist::Group::Hydroxy:: ) &&
            is_double_bond( $graph, $atom, grep { !blessed $_ } @O ) ) {
            my $acid = ChemOnomatopist::Group::SulfinicAcid->new( @C );
            graph_replace_all( $graph, $acid, $atom, @O );
        }
        # Detecting sulfonic acids
        if( !blessed $atom && is_element( $atom, 'S' ) && @neighbours == 4 && @C == 1 && @O == 3 &&
            @groups == 1 && $groups[0]->isa( ChemOnomatopist::Group::Hydroxy:: ) &&
            all { is_double_bond( $graph, $atom, $_ ) } grep { !blessed $_ } @O ) {
            my $acid = ChemOnomatopist::Group::SulfonicAcid->new( @C );
            graph_replace_all( $graph, $acid, $atom, @O );
        }

        if( !blessed $atom && @C == 1 && @groups == 1 &&
            ( is_element( $atom, 'O' ) || is_element( $atom, 'S' ) || is_element( $atom, 'Se' ) || is_element( $atom, 'Te' ) ) && 
            $groups[0]->isa( ChemOnomatopist::Group::Hydroxy:: ) ) {
            my $hydroperoxide = ChemOnomatopist::Group::Hydroperoxide->new( @C, $atom, @groups );
            graph_replace_all( $graph, $hydroperoxide, $atom, @groups );
        }
    }

    # Hydrogen atoms are no longer important.
    # They are demoted to hydrogen counts for future reference.
    for my $H (grep { is_element( $_, 'H' ) } $graph->vertices) {
        my( $parent ) = $graph->neighbours( $H );
        $parent->{hcount} = 0 unless exists $parent->{hcount};
        $parent->{hcount}++;
        push @{$parent->{h_isotope}}, $H->{isotope};
        $graph->delete_vertex( $H );
    }

    # Detecting cyclic compounds
    my @ring_systems = cyclic_components( $graph );

    # Aromatising mancudes - experimental
    for my $core (@ring_systems) {
        my %valences;
        for my $edge ($core->edges) {
            my $order = 1;
            if( $core->has_edge_attribute( @$edge, 'bond' ) ) {
                $order = $bond_symbol_to_order{$core->get_edge_attribute( @$edge, 'bond' )};
            }
            $valences{$edge->[0]} += $order;
            $valences{$edge->[1]} += $order;
        }

        my %uniq = map { join( '', sort @$_ ) => $_ } Graph::SSSR::get_SSSR( $core, 8 );
        my @aromatic;
        for my $cycle (values %uniq) {
            next if any { $core->degree( $_ ) > 3 } @$cycle;

            my @v2 = grep { $core->degree( $_ ) == 2 } @$cycle;
            my @v3 = grep { $core->degree( $_ ) == 3 } @$cycle;

            my @nonaromatic_atoms = grep { $valences{$_} < 3 } @v2;
            if( @nonaromatic_atoms ) {
                next if @nonaromatic_atoms > 1;
                next unless element( $nonaromatic_atoms[0] ) =~ /^[NOP]$/;
            }

            push @aromatic, [ Graph::Traversal::DFS->new( subgraph( $core, @$cycle ) )->dfs ];
        }
        for my $cycle (@aromatic) {
            my @vertices = @$cycle;
            for (0..$#vertices) {
                $graph->set_edge_attribute( $vertices[$_],
                                            $vertices[($_ + 1) % @vertices],
                                            'bond',
                                            ':' );
                if( $vertices[$_]->{symbol} =~ /^(Se|As|[BCNOPS])$/ ) {
                    $vertices[$_]->{symbol} = lcfirst $vertices[$_]->{symbol};
                }
            }
        }
    }

    for my $core (@ring_systems) {
        my %vertices_by_degree;
        for my $vertex ($core->vertices) {
            my $degree = $core->degree( $vertex );
            $vertices_by_degree{$degree} = [] unless $vertices_by_degree{$degree};
            push @{$vertices_by_degree{$degree}}, $vertex;
        }

        my $compound;
        if(      join( ',', sort keys %vertices_by_degree ) eq '2' ) {
            # Monocycles
            $compound = ChemOnomatopist::Group::Monocycle->new( $graph, Graph::Traversal::DFS->new( $core )->dfs );
        } elsif( ChemOnomatopist::Group::Monospiro->has_form( $core ) ) {
            # BBv2 P-24.2.1 Monospiro alicyclic ring systems
            $compound = ChemOnomatopist::Group::Monospiro->new( $graph, $core->vertices );
        } elsif( ChemOnomatopist::Group::Bicycle->has_form( $core ) ) {
            # Ortho-fused as defined in BBv2 P-25.3.1.1.1
            $compound = ChemOnomatopist::Group::Bicycle->new( $graph, $core->vertices );
        } elsif( join( ',', sort keys %vertices_by_degree ) eq '2,3' ) {
            # Fused ring systems of three or more rings
            # Graph::SSSR v0.1.0 does not know how to return unique rings
            my %uniq = map { join( '', sort @$_ ) => $_ } Graph::SSSR::get_SSSR( $core, 8 );
            my @cycles = map { ChemOnomatopist::Group::Monocycle->new( copy $graph, Graph::Traversal::DFS->new( $_ )->dfs ) }
                         map { subgraph( $core, @$_ ) }
                             values %uniq;
            if( (grep {  $_->is_benzene }  @cycles) == 2 &&
                (grep { !$_->is_benzene }  @cycles) == 1 &&
                (all  {  $_->length == 6 } @cycles) &&
                (any  { !$_->is_hydrocarbon } @cycles) ) {
                $compound = ChemOnomatopist::Group::Xanthene->new( $graph, @cycles );
            } elsif( @cycles >= 4 &&
                     (all { $_->length == 6 && $_->is_hydrocarbon } @cycles) &&
                     are_isomorphic( graph_without_edge_attributes( $core ),
                                     ChemOnomatopist::Group::Polyacene->ideal_graph( scalar $core->vertices ),
                                     sub { return 'C' } ) ) {
                $compound = ChemOnomatopist::Group::Polyacene->new( $graph, @cycles );
            } elsif( @cycles >= 4 &&
                     (all { $_->length == 6 && $_->is_hydrocarbon } @cycles) &&
                     are_isomorphic( graph_without_edge_attributes( $core ),
                                     ChemOnomatopist::Group::Polyaphene->ideal_graph( scalar $core->vertices ),
                                     sub { return 'C' } ) ) {
                $compound = ChemOnomatopist::Group::Polyaphene->new( $graph, @cycles );
            } else {
                die "cannot handle complicated cyclic compounds\n";
            }
        } else {
            die "cannot handle complicated cyclic compounds\n";
        }
        $graph->add_group( $compound );
    }

    return;
}

# Derive the chemical element of atom or group representation
# TODO: Replace object methods is_...
sub element
{
    my( $atom_or_group ) = @_;
    return undef unless ref $atom_or_group;

    if( !blessed $atom_or_group ) {
        die "unknown value '$atom_or_group' given for element()\n" unless ref $atom_or_group eq 'HASH';
        return ucfirst $atom_or_group->{symbol};
    }

    if( $atom_or_group->isa( 'Chemistry::Atom' ) ) { # PerlMol Atom
        return $atom_or_group->symbol;
    }

    return $atom_or_group->element;
}

# Check if an object or Perl hash is of certain chemical element
sub is_element
{
    my( $atom, $element ) = @_;
    return unless ref $atom;

    $element = ucfirst $element;

    if( blessed $atom ) {
        return element( $atom ) && element( $atom ) eq $element;
    }

    return ref $atom eq 'HASH' &&
           exists $atom->{symbol} &&
           ucfirst $atom->{symbol} eq $element;
}

# Given a graph, selects the main chain.
# The returned chain is an object of ChemOnomatopist::Chain or its subclasses.
sub select_mainchain
{
    my( $graph ) = @_;

    # Find the most senior group, undefined if alkane
    my @groups = most_senior_groups( $graph );
    my $most_senior_group = blessed $groups[0] if @groups;

    my @chains;
    if( @groups ) {
        # Select a chain containing most of the senior groups
        # FIXME: Parents with the most attachments should be preferred
        my @parents = uniq grep { defined $_ } map { $_->C } @groups;

        # Prefer circular structures
        if( @parents > 1 && (grep { blessed $_ && $_->isa( ChemOnomatopist::Chain:: ) } @parents) == 1 ) {
            @parents =       grep { blessed $_ && $_->isa( ChemOnomatopist::Chain:: ) } @parents;
        } elsif( @parents && uniq( map { $graph->groups( $_ ) } @parents ) == 1 ) {
            @parents = $graph->groups( $parents[0] );
        }

        if( $most_senior_group->isa( ChemOnomatopist::Chain:: ) ) {
            @chains = map { $_->can( 'candidates' ) ? $_->candidates : $_ } @groups;
        } elsif( @parents == 1 ) {
            if( blessed $parents[0] && $parents[0]->can( 'candidates' ) ) {
                @chains = $parents[0]->candidates;
            } else {
                # As the starting position is known, it is enough to take the "sidechain"
                # containing this particular parent:
                my $chain = select_sidechain( $graph, undef, @parents );
                my @vertices = $chain->can( 'vertices' ) ? $chain->vertices : $chain;
                push @chains, ChemOnomatopist::Chain->new( $graph, undef, @vertices ),
                              ChemOnomatopist::Chain->new( $graph, undef, reverse @vertices );
            }
        } elsif( @parents ) {
            my $copy = $graph->copy;
            $copy->delete_vertices( map { $_->vertices } $copy->groups );
            my @paths;
            my $max_value;
            for my $i (0..$#parents) {
                for my $j (($i+1)..$#parents) {
                    my @path = graph_path_between_vertices( $copy, $parents[$i], $parents[$j] );
                    next unless @path;
                    my $value = (set( @parents ) * set( @path ))->size;
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
                $copy->delete_vertices( grep { !is_element( $_, 'C' ) && !$_->is_part_of_chain }
                                        grep { blessed $_ }
                                             $copy->vertices );

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
                                                          undef,
                                                          reverse( @{$longest_paths{$A}->[$i]} ),
                                                          @$path,
                                                          @{$longest_paths{$B}->[$j]} ),
                             ChemOnomatopist::Chain->new( $graph,
                                                          undef,
                                                          reverse( @{$longest_paths{$B}->[$j]} ),
                                                          reverse( @$path ),
                                                          @{$longest_paths{$A}->[$j]} );
                    }
                }
            }
            die "cannot determine the parent structure\n" unless @chains;

            @chains = rule_most_groups( $most_senior_group, @chains );
        } elsif( @groups ) {
            # Attempt to build chains from functional groups
            @chains = map { ChemOnomatopist::Chain->new( $graph, undef, $_ ) } @groups;
        } else {
            die "cannot determine the parent structure\n";
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

    # This is needed to detect ethers.
    # However, it clears the cache of chains, thus is quite suboptimal.
    if( blessed $chain eq ChemOnomatopist::Chain::FromHalves:: ) {
        $chain = ChemOnomatopist::Chain->new( $graph, undef, @vertices );
    }

    # Replace the original chain with the selected candidate
    if( $chain->isa( ChemOnomatopist::Group:: ) && $chain->candidate_for ) {
        graph_replace_all( $graph, $chain, $chain->candidate_for );
    }

    # If there is at least one of carbon-based senior group attachment,
    # it means both ends are already senior, prompting to follow the
    # exception of three or more carbon-based groups.
    if( $most_senior_group && $groups[0]->is_carbon &&
        !$chain->isa( ChemOnomatopist::Group:: ) &&
         $chain->number_of_groups( $most_senior_group ) ) {

        shift @vertices;
        pop @vertices;
        $chain = ChemOnomatopist::Chain->new( $graph, undef, @vertices );
    }

    return $chain;
}

# Selects the best side chain
sub select_sidechain
{
    my( $graph, $parent, $start ) = @_;

    # Do this for non-carbons for now in order to represent attachments
    if( !is_element( $start, 'C' ) && $graph->degree( $start ) == 0 + defined $parent ) {
        return ChemOnomatopist::Chain->new( $graph, $parent, $start );
    }

    # Chalcogen analogues of ethers
    if( $graph->degree( $start ) == 1 + defined $parent && grep { element( $start ) && element( $start ) eq $_ } qw( S Se Te ) ) {
        return ChemOnomatopist::Chain->new( $graph, $parent, $start );
    }

    my $C_graph = copy $graph;
    $C_graph->delete_edge( $start, $parent ) if $parent;
    # Delete formed chains and non-carbon leaves
    # FIXME: Some other chains should as well be excluded
    $C_graph->delete_vertices( grep { $_ != $start && !is_element( $_, 'C' ) && $C_graph->degree( $_ ) == 1 } $C_graph->vertices );
    $C_graph->delete_vertices( grep { $_ != $start }
                               map  { $_->vertices }
                                    $C_graph->groups );

    return ChemOnomatopist::Chain->new( $graph, $parent, $start ) unless $C_graph->degree( $start );

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
    if(      $C_graph->degree( $start ) > 1 && (!blessed $start || !$start->is_terminal) ) {
        # FIXME: Deduplicate: copied from select_mainchain()
        # Generate all possible chains.
        # FIXME: This needs optimisation.
        for my $part1 (@path_parts) {
            for my $part2 (@path_parts) {
                next if $part1->group eq $part2->group;
                push @chains, ChemOnomatopist::Chain::FromHalves->new( $part1, $part2 );
            }
        }
    } else {
        @chains = map { ChemOnomatopist::Chain->new( $graph, $parent, $_->vertices ) }
                      @path_parts;
    }

    die "cannot select a sidechain\n" unless @chains;

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
                   \&rule_most_multiple_bonds,
                   \&rule_most_double_bonds,
                   \&rule_lowest_numbered_multiple_bonds,
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

    die "cannot select a sidechain\n";
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
                   # P-44.3.3: Maximum number of the most senior skeletal heteroatom
                   \&rule_greatest_number_of_most_senior_heteroatoms,

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

        if( $DEBUG ) {
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

sub graph_replace_all
{
    my( $graph, $new, @old ) = @_;

    # Replace in the main graph
    graph_replace( $graph, $new, @old );

    my $old = set( @old );
    for my $vertex ($graph->vertices) {
        next if $vertex == $new;
        next if $old->has( $vertex );
        next unless blessed $vertex;

        # Update parents
        if( $vertex->isa( ChemOnomatopist::Group:: ) && $vertex->C && $old->has( $vertex->C ) ) {
            $vertex->{C} = $new;
        }

        # Update internal graphs
        if( $vertex->isa( ChemOnomatopist::Group::Bicycle:: ) ||
            $vertex->isa( ChemOnomatopist::Group::Guanidine:: ) ||
            $vertex->isa( ChemOnomatopist::Group::Monocycle:: ) ||
            $vertex->isa( ChemOnomatopist::Group::Monospiro:: ) ||
            $vertex->isa( ChemOnomatopist::Group::Xanthene:: ) ) {
            graph_replace( $vertex->graph, $new, @old );
        }
    }

    return $graph;
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

sub rule_greatest_number_of_most_senior_heteroatoms
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

    my @groups;
    if( @vertices == 1 && blessed $vertices[0] && $vertices[0]->isa( Graph::Undirected:: ) ) {
        # Graph given instead of an array of vertices
        @groups = $vertices[0]->groups if $vertices[0]->isa( ChemOnomatopist::MolecularGraph:: );
        @vertices = $vertices[0]->vertices;
    }

    push @groups,
         grep { blessed $_ && $_->isa( ChemOnomatopist::Group:: ) && !$_->is_prefix_only } @vertices;
    return unless @groups;

    my( $most_senior_group ) = sort { ChemOnomatopist::Group::cmp( $a, $b ) } @groups;
    return grep { !ChemOnomatopist::Group::cmp( $_, $most_senior_group ) } @groups;
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

# BBv2 P-82.2.1
sub cmp_isotopes
{
    my @mass;
    my @symbol;
    for (@_) {
        if( /^(\d+)(\D+)$/ ) {
            push @mass, $1;
            push @symbol, $2;
        }
    }
    return $symbol[0] cmp $symbol[1] || $mass[0] cmp $mass[1];
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

    die "alkane chain of zero length detected\n" unless $N;

    my @names = qw( ? meth eth prop but );

    return $names[$N] if $N < @names;
    return IUPAC_numerical_multiplier( $N );
}

sub unbranched_chain_name($)
{
    my( $chain ) = @_;

    my @chain = $chain->vertices;

    my $name = ChemOnomatopist::Name->new;
    if( $chain->length == 1 && !blessed $chain[0] && !is_element( @chain, 'C' ) ) {
        $name .= 'ne'; # Leaving element prefix appending to the caller
        return $name;
    }

    my @bonds = $chain->bonds;
    my @double = grep { $bonds[$_] eq '=' } 0..$#bonds;
    my @triple = grep { $bonds[$_] eq '#' } 0..$#bonds;

    # BBv2 P-63.2.2.2
    if( $chain->parent && (all { !blessed $_ } @chain) && is_element( $chain[0], 'O' ) &&
        !@double && !@triple && all { is_element( $_, 'C' ) } @chain[1..$#chain] ) {
        $name->append_stem( alkane_chain_name( $chain->length - 1 ) );
        $name .= 'oxy';
        return $name;
    }

    if( $chain->isa( ChemOnomatopist::Chain::Amide:: ) ||
        $chain->isa( ChemOnomatopist::Chain::Amine:: ) ) {
        $name->append_stem( alkane_chain_name scalar grep { !blessed $_ } $chain->vertices );
    } elsif( (any { is_element( $_, 'C' ) } @chain) ||
        scalar( uniq map { element $_ } @chain ) > 1 ) {
        $name->append_stem( alkane_chain_name $chain->length );
    }

    if( @double ) {
        $name .= 'a' if @double >= 2; # BBv2 P-16.8.2
        if( $chain->needs_multiple_bond_locants || @double > 1 || @triple ) {
            $name->append_locants( $chain->bond_locants( @double ) );
        }
        if( @double > 1 ) {
            my $multiplier = IUPAC_numerical_multiplier scalar @double;
            $multiplier .= 'a' unless $multiplier =~ /i$/; # BBv2 P-31.1.1.2
            $name->append_multiplier( $multiplier );
        }
        $name .= 'en';
    }
    if( @triple ) {
        $name .= 'a' if @triple >= 2 && !@double; # BBv2 P-16.8.2
        if( $chain->needs_multiple_bond_locants || @triple > 1 || @double ) {
            $name->append_locants( $chain->bond_locants( @triple ) );
        }
        if( @triple > 1 ) {
            my $multiplier = IUPAC_numerical_multiplier scalar @triple;
            $multiplier .= 'a' unless $multiplier =~ /i$/; # BBv2 P-31.1.1.2
            $name->append_multiplier( $multiplier );
        }
        $name .= 'yn';
    }

    $name .= ChemOnomatopist::Name::Part::AlkaneANSuffix->new( 'an' ) unless @double || @triple;
    $name .= 'e';
    return $name;
}

1;
