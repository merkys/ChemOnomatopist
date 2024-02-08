package ChemOnomatopist::Chain;

# ABSTRACT: Chain of atoms
# VERSION

use strict;
use warnings;

use ChemOnomatopist;
use ChemOnomatopist::Chain::Ether;
use ChemOnomatopist::Elements qw( %elements );
use ChemOnomatopist::Group::Carboxyl;
use ChemOnomatopist::Group::Ether;
use ChemOnomatopist::Util::SMILES qw( path_SMILES );
use Graph::Traversal::DFS;
use List::Util qw( all any sum0 uniq );
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
    } elsif( blessed $vertices[0] && $vertices[0]->isa( ChemOnomatopist::Group::Amine:: ) ) {
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
        next if ChemOnomatopist::is_element( $vertices[$_], 'C' );
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
    if( @vertices == 1 && grep { ChemOnomatopist::element( @vertices ) eq $_ } qw( S Se Te ) ) {
        return '';
    }

    return 1;
}

sub needs_suffix_locant()
{
    my( $self ) = @_;

    # BBv2 P-14.3.4.2 (a): mononuclear parent hydrides do not need locants
    return '' if $self->length == 1;

    # BBv2 P-14.3.4.2 (b): monosubstituted homogeneous chains consisting of only two identical atoms do not need locants
    return '' if $self->length == 2 && $self->number_of_branches == 1;

    my @most_senior_groups = $self->most_senior_groups;
    return '' unless @most_senior_groups;
    return 1 if $self->number_of_heteroatoms; # P-15.4.3.2.3: Characteristic groups cited as suffixes are given locants
    return 1 if !$most_senior_groups[0]->is_carbon || @most_senior_groups > 2;
    return '';
}

sub needs_substituent_locants()
{
    my( $self ) = @_;
    return '' if $self->length == 1;

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

sub number_of_multiple_bonds()
{
    my( $self ) = @_;
    return scalar grep { $_ =~ /^[=#\$]$/ } $self->bonds;
}

sub prefix()
{
    my( $self ) = @_;

    # Chalcogen analogues of ethers
    if( $self->length == 1 ) {
        return ChemOnomatopist::Name->new( 'sulfan' ) if ChemOnomatopist::is_element( $self->vertices, 'S' );
        return ChemOnomatopist::Name->new( 'selan'  ) if ChemOnomatopist::is_element( $self->vertices, 'Se' );
        return ChemOnomatopist::Name->new( 'tellan' ) if ChemOnomatopist::is_element( $self->vertices, 'Te' );
    }

    return ChemOnomatopist::unbranched_chain_name( $self ); # FIXME: Add proper suffix
}

sub suffix()
{
    my( $self ) = @_;
    return ChemOnomatopist::unbranched_chain_name( $self );
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
