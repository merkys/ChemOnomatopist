package ChemOnomatopist::Group;

use strict;
use warnings;

use ChemOnomatopist::Chain::Circular;
use ChemOnomatopist::Group::AcylHalide;
use ChemOnomatopist::Group::Aldehyde;
use ChemOnomatopist::Group::Amide;
use ChemOnomatopist::Group::Amine;
use ChemOnomatopist::Group::Carbonitrile;
use ChemOnomatopist::Group::Carboxyl;
use ChemOnomatopist::Group::Cyanide;
use ChemOnomatopist::Group::Ester;
use ChemOnomatopist::Group::Ether;
use ChemOnomatopist::Group::Guanidine;
use ChemOnomatopist::Group::Hydrazide;
use ChemOnomatopist::Group::Hydrazine;
use ChemOnomatopist::Group::Hydroperoxide;
use ChemOnomatopist::Group::Hydroxy;
use ChemOnomatopist::Group::Imino;
use ChemOnomatopist::Group::Ketone;
use ChemOnomatopist::Group::SulfinicAcid;
use ChemOnomatopist::Group::SulfonicAcid;
use ChemOnomatopist::Group::Urea;

use List::Util qw( any );
use Scalar::Util qw( blessed );

# ABSTRACT: Chemical group
# VERSION

# Order from BBv2 P-41
our @order = (
    # Radicals
    # Radical anions
    # Radical cations
    # Anions
    # Zwitterions
    # Cations
    # Acids (BBv2 P-42)
    ChemOnomatopist::Group::Carboxyl::,
    ChemOnomatopist::Group::SulfonicAcid::,
    ChemOnomatopist::Group::SulfinicAcid::,
    ChemOnomatopist::Group::AcylHalide::, # FIXME: Is this correct?
    # Anhydrides
    ChemOnomatopist::Group::Ester::,
    # Acid halides and pseudohalides
    # Amides
    ChemOnomatopist::Group::Urea::,
    ChemOnomatopist::Group::Amide::,
    ChemOnomatopist::Group::Guanidine::,
    # Hydrazides
    ChemOnomatopist::Group::Hydrazide::,
    # Imides
    # 14. Nitriles
    ChemOnomatopist::Group::Carbonitrile::,
    ChemOnomatopist::Group::Cyanide::,
    ChemOnomatopist::Group::Aldehyde::,
    ChemOnomatopist::Group::Ketone::,
    ChemOnomatopist::Group::Hydroxy::,
    ChemOnomatopist::Group::Hydroperoxide::,
    ChemOnomatopist::Group::Amine::,
    ChemOnomatopist::Group::Imino::,

    # TODO: Some are omitted

    # TODO: Classes denoted by the senior atom in heterane nomenclature should go here
    # 21. Nitrogen compounds
    ChemOnomatopist::Group::Hydrazine::,
    # 41. Ethers, then sulfides, sulfoxides, sulfones; then selenides, selenoxides, etc.
    ChemOnomatopist::Group::Ether::,
);

sub new
{
    my( $class, $element ) = @_;
    return bless { element => $element }, $class;
}

sub element() { return $_[0]->{element} }

sub is_carbon()
{
    my( $self ) = @_;
    return $self->element && $self->element eq 'C';
}

sub is_nitrogen()
{
    my( $self ) = @_;
    return $self->element && $self->element eq 'N';
}

sub is_oxygen()
{
    my( $self ) = @_;
    return $self->element && $self->element eq 'O';
}

sub is_part_of_chain() { return '' }

# Certain groups can only be expressed as prefixes
sub is_prefix_only() { return '' }

# Certain groups can only be terminal in chains
sub is_terminal() { return '' }

sub needs_heteroatom_locants { return 1 }
sub needs_heteroatom_names { return 1 }
sub needs_multiple_bond_suffix { return 1 }

sub prefix() { return '' }
sub suffix() { return $_[0]->is_prefix_only ? undef : '' }
sub multisuffix() { return $_[0]->suffix }
sub suffix_if_cycle_substituent() { return $_[0]->suffix }

sub candidate_for()
{
    my( $self ) = @_;
    return undef unless exists $self->{candidate_for};
    return $self->{candidate_for};
}

sub rule_greatest_number_of_most_senior_heteroatoms
{
    my( @chains ) = @_;

    # This order is taken from BBv2 P-41 and is different from order in %elements
    my @element_order = qw( N P As Sb Bi Si Ge Sn Pb B Al Ga In Tl O S Se Te );
    my %element_order = map { $element_order[$_] => $_ } 0..$#element_order;

    my( $most_senior ) = sort { $element_order{$a} <=> $element_order{$b} }
                         grep { exists $element_order{$_} }
                         map  { $_->heteroatoms } @chains;
    return @chains unless $most_senior;

    my( $max_value ) = reverse sort map { scalar( grep { $_ eq $most_senior } $_->heteroatoms ) }
                                        @chains;
    return grep { scalar( grep { $_ eq $most_senior } $_->heteroatoms ) == $max_value }
                @chains;
}

# Compare seniority of two objects
sub cmp
{
    my( $A, $B ) = @_;

    my( $A_pos ) = grep { $A->isa( $order[$_] ) } 0..$#order;
    my( $B_pos ) = grep { $B->isa( $order[$_] ) } 0..$#order;

    # Clear distinction exists
    if( defined $A_pos && defined $B_pos && $A_pos <=> $B_pos ) {
        return $A_pos <=> $B_pos;
    }

    # Any of the objects is in the priority list
    if( defined $A_pos ^ defined $B_pos ) {
        return defined $B_pos <=> defined $A_pos;
    }

    # Same class; class should know how to compare
    if( blessed $A eq blessed $B ) {
        return $A->_cmp_instances( $B );
    }

    # BBv2 P-41
    # First, the chain with the most senior atom wins
    # FIXME: Select just by seniority, not by number
    my @chains = rule_greatest_number_of_most_senior_heteroatoms( $A, $B );
    return ($chains[0] == $B) * 2 - 1 if @chains;

    # Second, the order is heterocycles, polyheteroatom, heteroatom
    if( $A->isa( ChemOnomatopist::Chain::Circular:: ) + 0 ^
        $B->isa( ChemOnomatopist::Chain::Circular:: ) + 0 ) {
        return $B->isa( ChemOnomatopist::Chain::Circular:: ) <=>
               $A->isa( ChemOnomatopist::Chain::Circular:: );
    }

    # TODO: The remaining rules from P-41

    die "cannot compare\n";
}

# Two instances of the same group are thought to be of the same seniority
sub _cmp_instances { return 0 }

1;
