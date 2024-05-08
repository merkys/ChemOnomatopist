package ChemOnomatopist::Chain;

# ABSTRACT: Chain of atoms
# VERSION

use strict;
use warnings;

use ChemOnomatopist;
use ChemOnomatopist::Chain::Amide;
use ChemOnomatopist::Chain::Amine;
use ChemOnomatopist::Chain::Ether;
use ChemOnomatopist::Charge;
use ChemOnomatopist::Elements qw( %elements );
use ChemOnomatopist::Group::Carboxyl;
use ChemOnomatopist::Group::Ether;
use ChemOnomatopist::Isotope;
use ChemOnomatopist::Name::Part::AlkaneANSuffix;
use ChemOnomatopist::Name::Part::Isotope;
use ChemOnomatopist::Util::SMILES qw( path_SMILES );
use Graph::Traversal::DFS;
use List::Util qw( all any first none sum0 uniq );
use Scalar::Util qw( blessed );
use Set::Object qw( set );

sub vertices();

sub new
{
    my( $class, $graph, $parent, @vertices ) = @_;

    # TODO: For now only chains with a single oxygen atom are detected as ethers.
    # First of all, ChemOnomatopist::Chain::Ether is not very clever.
    # Second, BBv2 P-63.2.4.1 does not draw a clear line between ether naming and skeletal replacement nomenclature.
    my $self;
    if( (grep { blessed $_ && $_->isa( ChemOnomatopist::Group::Ether:: ) } @vertices) == 1 ) {
        $self = ChemOnomatopist::Chain::Ether->new( $graph, $parent, @vertices );
    } elsif( blessed $vertices[0] && $vertices[0]->isa( ChemOnomatopist::Group::Amine:: ) &&
        none { blessed $_ && $_->isa( ChemOnomatopist::Group::Amine:: ) } @vertices[1..$#vertices] ) {
        my $amine = shift @vertices;
        my $chain = bless { vertices => \@vertices, graph => $graph }, $class;
        $self = ChemOnomatopist::Chain::Amine->new( $graph, $chain, $amine );
    } else {
        $self = { vertices => \@vertices, graph => $graph, cache => {} };
        $self->{parent} = $parent if $parent;
        $self = bless $self, $class;
    }
    return $self;
}

# Accessors

sub graph()
{
    my( $self ) = @_;
    return $self->{graph};
}

sub parent(;$)
{
    my( $self, $parent ) = @_;
    my $old_parent = exists $self->{parent} ? $self->{parent} : undef;
    $self->{parent} = $parent if $parent; # TODO: Maybe invalidate the related cache
    return $old_parent;
}

sub parent_locant()
{
    my( $self ) = @_;
    my $parent = $self->parent;
    return $parent unless $parent;

    my $graph = $self->graph;
    my @vertices = $self->vertices;
    return first { $graph->has_edge( $vertices[$_], $parent ) } 0..$#vertices;
}

sub substituents()
{
    my( $self ) = @_;

    my $vertices = set( $self->vertices );

    my @substituents;
    for my $vertex ($self->vertices) {
        for my $neighbour ($self->graph->neighbours( $vertex )) {
            next if $vertices->has( $neighbour );
            push @substituents, $neighbour;
        }
    }

    return @substituents;
}

sub vertices()
{
    my( $self ) = @_;
    return @{$self->{vertices}};
}

# Properties

sub backbone_SMILES()
{
    my( $self ) = @_;
    return path_SMILES( $self->graph, $self->vertices );
}

sub bonds()
{
    my( $self ) = @_;

    return @{$self->{bonds}} if $self->{bonds};

    my $graph = $self->graph;
    my @vertices = $self->vertices;

    my @bonds;
    for (0..$self->length-2) {
        if( $graph->has_edge_attribute( @vertices[$_ .. $_+1], 'bond' ) ) {
            push @bonds, $graph->get_edge_attribute( @vertices[$_ .. $_+1], 'bond' );
        } else {
            push @bonds, '-';
        }
    }

    $self->{bonds} = \@bonds;
    return @bonds;
}

sub branch_positions()
{
    my( $self ) = @_;

    return @{$self->{branch_positions}} if $self->{branch_positions};

    my $graph = $self->_disconnected_chain_graph;
    my @vertices = $self->vertices;

    my @branch_positions =
        map { ( $_ ) x $graph->degree( $vertices[$_] ) }
            grep { $graph->degree( $vertices[$_] ) }
                 0..$#vertices;

    $self->{branch_positions} = \@branch_positions;
    return @branch_positions;
}

sub group_positions
{
    my( $self, $class ) = @_;

    return @{$self->{group_positions}{$class}} if $self->{group_positions}{$class};

    my $graph = $self->_disconnected_chain_graph;
    my @vertices = $self->vertices;

    my @group_positions;
    for (0..$#vertices) {
        my $groups = grep { blessed $_ && $_->isa( $class ) }
                          $graph->neighbours( $vertices[$_] );
        next unless $groups;
        push @group_positions, ( $_ ) x $groups;
    }

    $self->{group_positions}{$class} = \@group_positions;
    return @group_positions;
}

sub heteroatom_positions()
{
    my( $self ) = @_;

    return @{$self->{heteroatom_positions}} if $self->{heteroatom_positions};

    my @vertices = $self->vertices;
    my @heteroatom_positions;
    for (0..$#vertices) {
        next if blessed $vertices[$_] && !$vertices[$_]->isa( ChemOnomatopist::Group::Ether:: );
        next if ChemOnomatopist::is_element( $vertices[$_], 'C' );
        push @heteroatom_positions, $_;
    }

    $self->{heteroatom_positions} = \@heteroatom_positions;
    return @heteroatom_positions;
}

# CHECKME: Can chains have indicated hydrogens?
sub indicated_hydrogens() { my @hydrogen_positions; return @hydrogen_positions }

sub nonstandard_valence_positions()
{
    my( $self ) = @_;

    return @{$self->{nonstandard_valence_positions}} if $self->{nonstandard_valence_positions};

    my @vertices = $self->vertices;
    my @nonstandard_valence_positions;
    for (0..$#vertices) {
        next if blessed $vertices[$_];
        next if ChemOnomatopist::element( $vertices[$_] ) eq 'C';
        next if ChemOnomatopist::element( $vertices[$_] ) eq 'N';
        next unless exists $vertices[$_]->{valence};
        push @nonstandard_valence_positions, $_;
    }

    $self->{nonstandard_valence_positions} = \@nonstandard_valence_positions;
    return @nonstandard_valence_positions;
}

sub is_hydrocarbon()
{
    my( $self ) = @_;
    return $self->number_of_heteroatoms == 0;
}

sub is_saturated()
{
    my( $self ) = @_;
    return all { $_ eq '-' } $self->bonds;
}

sub is_substituted()
{
    my( $self ) = @_;
    return $self->number_of_branches > 0;
}

# Returns maximum number of substitutable locations
sub max_valence()
{
    my( $self ) = @_;

    my $max_valence = 0;
    for my $vertex ($self->vertices) {
        my $element = ChemOnomatopist::element( $vertex );
        next unless $element;
        next if !exists $elements{$element};
        next if !exists $elements{$element}->{standard_bonding_number};
        $max_valence += $elements{$element}->{standard_bonding_number};
    }

    for my $bond ($self->bonds) {
        $max_valence -= 2 if $bond eq '-';
        $max_valence -= 3 if $bond eq ':'; # Aromatic bond is 1.5
        $max_valence -= 4 if $bond eq '=';
        $max_valence -= 6 if $bond eq '#';
        $max_valence -= 8 if $bond eq '$';
    }

    return $max_valence;
}

sub most_senior_groups()
{
    my( $self ) = @_;
    return ChemOnomatopist::most_senior_groups( $self->vertices, $self->substituents );
}

sub most_senior_group_positions()
{
    my( $self ) = @_;

    return @{$self->{most_senior_group_positions}} if $self->{most_senior_group_positions};

    my @groups = ChemOnomatopist::most_senior_groups( $self->vertices, $self->substituents );
    $self->{most_senior_group_positions} = [];
    return () unless @groups;

    my $groups = set( @groups );
    my @vertices = $self->vertices;
    my $vertices = set( @vertices );
    my @positions;
    for (0..$#vertices) {
        my $vertex = $vertices[$_];
        push @positions, $_ if $groups->has( $vertex );
        for my $neighbour ($self->{graph}->neighbours( $vertex )) {
            next if $vertices->has( $neighbour );
            push @positions, $_ if $groups->has( $neighbour );
        }
    }

    $self->{most_senior_group_positions} = \@positions;
    return @positions;
}

sub multiple_bond_positions()
{
    my( $self ) = @_;
    my @bonds = $self->bonds;
    return grep { $bonds[$_] =~ /^[=#\$]$/ } 0..$#bonds;
}

sub double_bond_positions()
{
    my( $self ) = @_;
    my @bonds = $self->bonds;
    return grep { $bonds[$_] eq '=' } 0..$#bonds;
}

sub needs_multiple_bond_locants()
{
    my( $self ) = @_;
    return $self->length > 2;
}

sub needs_multiple_bond_suffix()
{
    my( $self ) = @_;
    my( $first ) = $self->vertices;
    return 1 unless blessed $first;
    return $first->needs_multiple_bond_suffix;
}

sub needs_heteroatom_locants()
{
    my( $self ) = @_;
    return '' if $self->length == 1;

    my @vertices = $self->vertices;
    # Check if this is -oxy substituent
    if( $self->parent && !$self->number_of_branches && $self->number_of_heteroatoms == 1 &&
        $vertices[0]->{symbol} eq 'O' ) {
        return '';
    }

    if(      scalar( uniq $self->heteroatoms ) == 1 ) {
        return $self->number_of_heteroatoms != $self->length;
    } elsif( scalar( uniq $self->heteroatoms ) >  1 ) {
        return 1;
    }
}

sub needs_heteroatom_names()
{
    my( $self ) = @_;

    my @vertices = $self->vertices;
    # Check if this is -oxy substituent
    if( $self->parent && !$self->number_of_branches && $self->number_of_heteroatoms == 1 &&
        $vertices[0]->{symbol} eq 'O' ) {
        return '';
    }

    # Chalcogen analogues of ethers
    if( @vertices == 1 && $self->parent && grep { ChemOnomatopist::element( @vertices ) eq $_ } qw( S Se Te ) ) {
        return '';
    }

    return 1;
}

sub needs_suffix_locant()
{
    my( $self ) = @_;

    # BBv2 P-14.3.4.2 (a): mononuclear parent hydrides do not need locants
    return '' if $self->length == 1;

    return 1 if $self->isotopes;

    # BBv2 P-14.3.4.2 (b): monosubstituted homogeneous chains consisting of only two identical atoms do not need locants
    return '' if $self->length == 2 && $self->number_of_branches == 1;

    my @most_senior_groups = $self->most_senior_groups;
    return '' unless @most_senior_groups;
    return 1 if $self->number_of_heteroatoms; # P-15.4.3.2.3: Characteristic groups cited as suffixes are given locants
    return 1 if @most_senior_groups > 2;
    return 1 if !ChemOnomatopist::element( $most_senior_groups[0] );
    return 1 if  ChemOnomatopist::element( $most_senior_groups[0] ) ne 'C';
    return '';
}

sub needs_substituent_locants()
{
    my( $self ) = @_;
    return '' if $self->length == 1;

    return 1 if $self->number_of_isotopes;

    # FIXME: Make sure the substituents are of the same kind
    return '' if scalar( $self->substituents ) == $self->max_valence;

    # Ad-hoc fix for acetic acid substituents
    if( $self->length == 2 &&
        any { blessed $_ && $_->isa( ChemOnomatopist::Group::Carboxyl:: ) }
            $self->vertices ) {
        return '';
    }

    return 1;
}

sub needs_charge_locants()  { &needs_substituent_locants }
sub needs_isotope_locants() { &needs_substituent_locants }

sub charges()
{
    my( $self ) = @_;
    return @{$self->{charges}} if $self->{charges};

    my @vertices = $self->vertices;
    my @charges;
    for my $i (0..$#vertices) {
        next if blessed $vertices[$i];
        next unless $vertices[$i]->{charge};
        push @charges, ChemOnomatopist::Charge->new( $vertices[$i]->{charge},
                                                     $i,
                                                     $self->locants( $i ) );
    }

    $self->{charges} = \@charges;
    return @charges;
}

sub heteroatoms()
{
    my( $self ) = @_;
    my @vertices = $self->vertices;
    return map { ChemOnomatopist::element( $vertices[$_] ) }
               $self->heteroatom_positions;
}

sub nonstandard_valences()
{
    my( $self ) = @_;
    my @vertices = $self->vertices;
    return map { $vertices[$_]->{valence} }
               $self->nonstandard_valence_positions;
}

sub isotopes()
{
    my( $self ) = @_;
    return @{$self->{isotopes}} if $self->{isotopes};

    my @vertices = $self->vertices;
    my @isotopes;
    for my $i (0..$#vertices) {
        next if blessed $vertices[$i];

        if( exists $vertices[$i]->{isotope} ) {
            push @isotopes, ChemOnomatopist::Isotope->new( ChemOnomatopist::element( $vertices[$i] ),
                                                           $vertices[$i]->{isotope},
                                                           $i,
                                                           $self->locants( $i ) );
        }
        if( exists $vertices[$i]->{h_isotope} ) {
            for (@{$vertices[$i]->{h_isotope}}) {
                next unless defined $_;
                push @isotopes, ChemOnomatopist::Isotope->new( 'H', $_, $i, $self->locants( $i ) );
            }
        }
    }

    $self->{isotopes} = \@isotopes;
    return @isotopes;
}

sub length()
{
    my( $self ) = @_;
    return scalar $self->vertices;
}

sub locant_names()
{
    my( $self ) = @_;

    return @{$self->{locant_names}} if $self->{locant_names};

    my $graph = $self->_disconnected_chain_graph->copy;

    my @locants;
    for my $vertex ($self->vertices) {
        my @current_locants;
        for my $neighbour ($graph->neighbours( $vertex )) {
            $graph->delete_edge( $vertex, $neighbour );
            next if $self->parent && $self->parent == $neighbour;
            if( blessed $neighbour ) {
                push @current_locants, $neighbour->prefix;
            } else {
                push @current_locants,
                     ChemOnomatopist::get_sidechain_name( $graph, $vertex, $neighbour );
            }
        }
        push @locants, \@current_locants;
    }

    $self->{locant_names} = \@locants;
    return @locants;
}

sub locants(@)
{
    my $self = shift;
    return map { $_ + 1 } @_;
}

sub bond_locants(@)
{
    my $self = shift;
    return map { $_ + 1 } @_;
}

sub number_of_branches_in_sidechains()
{
    my( $self ) = @_;

    return $self->{number_of_branches_in_sidechains} if exists $self->{number_of_branches_in_sidechains};

    my $graph = $self->_disconnected_chain_graph->copy;
    my @vertex_neighbours = map { $graph->neighbours( $_ ) } $self->vertices;
    $graph->delete_vertices( $self->vertices );

    # has_vertex() is used to filter out neighbours within the chain
    my $number = sum0 map { $_ > 2 ? $_ - 2 : 0 }
                          map  { $graph->degree( $_ ) }
                          map  { Graph::Traversal::DFS->new( $graph, start => $_ )->dfs }
                          grep { $graph->has_vertex( $_ ) }
                              @vertex_neighbours;

    $self->{number_of_branches_in_sidechains} = $number;
    return $number;
}

sub number_of_carbons()
{
    my( $self ) = @_;

    return $self->{number_of_carbons} if exists $self->{number_of_carbons};

    my $graph = $self->_disconnected_chain_graph;

    my $C = grep { ChemOnomatopist::is_element( $_, 'C' ) }
            map  { Graph::Traversal::DFS->new( $graph, start => $_ )->dfs }
                 $self->vertices;

    # Since main chain carbons are included in the count, they have to be subtracted.
    $C -= $self->length;

    $self->{number_of_carbons} = $C;
    return $C;
}

sub number_of_charges()
{
    my( $self ) = @_;
    return scalar $self->charges;
}

sub number_of_branches()
{
    my( $self ) = @_;
    return scalar $self->branch_positions;
}

sub number_of_double_bonds()
{
    my( $self ) = @_;
    return scalar grep { $_ eq '=' } $self->bonds;
}

sub number_of_groups
{
    my( $self, $class ) = @_;
    return scalar $self->group_positions( $class );
}

sub number_of_heteroatoms()
{
    my( $self ) = @_;
    return scalar $self->heteroatom_positions;
}

sub number_of_indicated_hydrogens()
{
    my( $self ) = @_;
    return scalar $self->indicated_hydrogens;
}

sub number_of_isotopes()
{
    my( $self ) = @_;
    return scalar $self->isotopes;
}

sub number_of_multiple_bonds()
{
    my( $self ) = @_;
    return scalar grep { $_ =~ /^[=#\$]$/ } $self->bonds;
}

sub number_of_nonstandard_valence_positions()
{
    my( $self ) = @_;
    return scalar $self->nonstandard_valence_positions;
}

sub charge_part()
{
    my( $self ) = @_;
    return ChemOnomatopist::Name->new unless $self->charges;

    my @negative = grep { $_->charge < 0 } $self->charges;
    my @positive = grep { $_->charge > 0 } $self->charges;

    my $name = ChemOnomatopist::Name->new;

    if( @positive ) {
        @positive = map { ( $_->locant ) x abs $_->charge } @positive;
        $name->append_locants( @positive ) if $self->needs_charge_locants;
        if( @positive > 1 ) {
            $name->append_multiplier( ChemOnomatopist::IUPAC_complex_numerical_multiplier( scalar @positive ) );
            $name .= '(';
        }
        $name .= @negative ? 'ium' : 'ylium';
        $name .= ')' if @positive > 1;
    }

    if( @negative ) {
        @negative = map { ( $_->locant ) x abs $_->charge } @negative;
        $name->append_locants( @negative ) if $self->needs_charge_locants;
        if( @negative > 1 ) {
            $name->append_multiplier( ChemOnomatopist::IUPAC_numerical_multiplier( scalar @negative ) );
        }
        $name .= 'ide';
    }

    return $name;
}

sub indicated_hydrogens_part() { ChemOnomatopist::Name->new }

sub isotope_part()
{
    my( $self ) = @_;

    my @isotopes = sort { $a->element cmp $b->element ||
                          $a->mass_number <=> $b->mass_number }
                        $self->isotopes;
    return '' unless @isotopes;

    my @order;
    my %freq;
    for my $isotope (@isotopes) {
        my $key = $isotope->mass_number . $isotope->element;
        if( !$freq{$key} ) {
            $freq{$key} = [];
            push @order, $key;
        }
        push @{$freq{$key}}, $isotope;
    }

    my @vertices = $self->vertices;
    my $isotopes = '';
    for my $key (@order) {
        $isotopes .= ',' if $isotopes;
        $isotopes .= join ',', map { $_->locant } @{$freq{$key}} if $self->needs_isotope_locants;
        $isotopes .= '-' if $self->needs_isotope_locants;
        $isotopes .= $key;
        if( @{$freq{$key}} > 1 ||
            ( $key =~ /H$/ && $vertices[$freq{$key}->[0]{index}]->{hcount} > 1 ) ) {
            $isotopes .= scalar @{$freq{$key}};
        }
    }

    return ChemOnomatopist::Name::Part::Isotope->new( "($isotopes)" );
}

sub prefix()
{
    my( $self ) = @_;

    if( $self->length == 1 ) {
        my $vertex = $self->{vertices}[0];
        return $vertex->prefix if blessed $vertex;

        # Chalcogen analogues of ethers
        my $element = ChemOnomatopist::element( $vertex );
        return ChemOnomatopist::Name->new( 'sulfan' ) if $element eq 'S';
        return ChemOnomatopist::Name->new( 'selan'  ) if $element eq 'Se';
        return ChemOnomatopist::Name->new( 'tellan' ) if $element eq 'Te';
    }

    my $name = $self->suffix;
    $name->pop_e;
    pop @$name if $name->ends_with_alkane_an_suffix;
    return $name . 'yl';
}

sub suffix()
{
    my( $self ) = @_;

    my @chain = $self->vertices;

    my $name = ChemOnomatopist::Name->new;
    if( $self->length == 1 && !blessed $chain[0] && ChemOnomatopist::element( @chain ) ne 'C' ) {
        return $name . 'ne'; # Leaving element prefix appending to the caller
    }

    # CHECKME: Not sure if calling prefix() is correct
    return $chain[0]->prefix if $self->length == 1 && blessed $chain[0];

    my @bonds = $self->bonds;
    my @double = grep { $bonds[$_] eq '=' } 0..$#bonds;
    my @triple = grep { $bonds[$_] eq '#' } 0..$#bonds;

    # BBv2 P-63.2.2.2
    if( $self->parent && @chain && (all { !blessed $_ } @chain) && ChemOnomatopist::element( @chain ) eq 'O' &&
        !@double && !@triple && all { ChemOnomatopist::element( $_ ) eq 'C' } @chain[1..$#chain] ) {
        $name->append_stem( ChemOnomatopist::alkane_chain_name( $self->length - 1 ) );
        $name .= 'oxy';
        return $name;
    }

    if( $self->isa( ChemOnomatopist::Chain::Amide:: ) ||
        $self->isa( ChemOnomatopist::Chain::Amine:: ) ) {
        $name->append_stem( ChemOnomatopist::alkane_chain_name( scalar grep { !blessed $_ } $self->vertices ) );
    } elsif( (any { ChemOnomatopist::is_element( $_, 'C' ) } @chain) ||
        scalar( uniq map { ChemOnomatopist::element( $_ ) } @chain ) > 1 ) {
        $name->append_stem( ChemOnomatopist::alkane_chain_name( $self->length ) );
    }

    if( @double ) {
        $name .= 'a' if @double >= 2; # BBv2 P-16.8.2
        if( $self->needs_multiple_bond_locants ) {
            $name->append_locants( $self->bond_locants( @double ) );
        }
        if( @double > 1 ) {
            my $multiplier = ChemOnomatopist::IUPAC_numerical_multiplier( scalar @double );
            $multiplier .= 'a' unless $multiplier =~ /i$/; # BBv2 P-31.1.1.2
            $name->append_multiplier( $multiplier );
        }
        $name .= 'en';
    }
    if( @triple ) {
        $name .= 'a' if @triple >= 2 && !@double; # BBv2 P-16.8.2
        if( $self->needs_multiple_bond_locants ) {
            $name->append_locants( $self->bond_locants( @triple ) );
        }
        if( @triple > 1 ) {
            my $multiplier = ChemOnomatopist::IUPAC_numerical_multiplier( scalar @triple );
            $multiplier .= 'a' unless $multiplier =~ /i$/; # BBv2 P-31.1.1.2
            $name->append_multiplier( $multiplier );
        }
        $name .= 'yn';
    }

    $name .= ChemOnomatopist::Name::Part::AlkaneANSuffix->new( 'an' ) unless @double || @triple;
    $name .= 'e';
    return $name;
}

sub vertex_ids
{
    my $self = shift;
    my %ids = map { $self->{vertices}[$_] => $_ } 0..$self->length-1;
    return map { exists $ids{$_} ? $ids{$_} : undef } @_;
}

sub _cmp_instances
{
    my( $A, $B ) = @_;
    # For now, just compare the sizes of cycles.
    # TODO: Proper ordering should be implemented as per BBv2 P-25.8.1
    return $B->length <=> $A->length;
}

sub _disconnected_chain_graph()
{
    my( $self ) = @_;

    my $graph = $self->graph->copy;
    my $vertices = set( $self->vertices );
    $graph->delete_edges( map { @$_ } grep { $vertices->has( @$_ ) } $graph->edges );

    return $graph;
}

1;
