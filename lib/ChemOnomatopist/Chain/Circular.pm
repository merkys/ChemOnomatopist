package ChemOnomatopist::Chain::Circular;

# ABSTRACT: Chain whose first and last members are connected
# VERSION

use strict;
use warnings;

use ChemOnomatopist;
use ChemOnomatopist::Chain; # FIXME: Not sure why it is needed
use ChemOnomatopist::Chain::Monocycle;
use ChemOnomatopist::Elements qw( %elements );
use ChemOnomatopist::Name::Part::NondetachablePrefix;
use ChemOnomatopist::Util::SMILES qw( cycle_SMILES );
use Chemistry::OpenSMILES qw( is_single_bond );
use List::Util qw( all any uniq );
use Scalar::Util qw( blessed );
use Set::Object qw( set );

use parent ChemOnomatopist::Chain::;

sub new
{
    my( $class, $graph, @vertices ) = @_;
    return bless { graph => $graph, vertices => \@vertices }, $class;
}

# Selecting the candidate with the lowest alphabetical order
sub backbone_SMILES()
{
    my( $self ) = @_;

    my @vertices = $self->vertices;
    my @candidates;
    for (0..$#vertices) {
        push @candidates, cycle_SMILES( $self->graph, @vertices );
        push @vertices, shift @vertices;
    }
    @vertices = reverse @vertices;
    for (0..$#vertices) {
        push @candidates, cycle_SMILES( $self->graph, @vertices );
        push @vertices, shift @vertices;
    }

    my( $SMILES ) = sort @candidates;
    return $SMILES;
}

sub is_Hantzsch_Widman()
{
    my( $self ) = @_;
    return $self->length >= 3 && $self->length <= 10 && # "no more than ten ring members"
           $self->number_of_heteroatoms &&              # "containing one or more heteroatoms"
           any { $_->{symbol} =~ /^[cC]$/ } $self->vertices;
}

# FIXME: What to do with furan and others?
sub is_aromatic()
{
    my( $self ) = @_;
    return 1 if all { $_ eq ':' } $self->bonds;
    return '';
}

sub is_benzene()
{
    my( $self ) = @_;
    return $self->is_aromatic && $self->is_homogeneous && $self->length == 6;
}

sub is_heterocycle()
{
    my( $self ) = @_;
    return $self->number_of_heteroatoms > 0;
}

sub is_homogeneous()
{
    my( $self ) = @_;

    if( any { blessed $_ } $self->vertices ) {
        die "cannot process cycles with groups as their members\n";
    }

    my @elements = map { $_->{symbol} } $self->vertices;
    my @bonds = $self->bonds;

    my( $element ) = @elements;
    my( $bond ) = @bonds;

    return '' if any { $_ ne $element } @elements;
    return  1 if $self->is_aromatic;

    return '' if any { $_ ne $bond } @bonds;
    return  1;
}

sub is_monoreplaced()
{
    my( $self ) = @_;
    return $self->number_of_heteroatoms == 1;
}

my %five_membered_aromatic_single_heteroatom = (
    N => 'pyrrole',
    O => 'furan',
    S => 'thiophene',
);

sub name()
{
    my( $self ) = @_;

    my $graph = $self->graph;
    my $SMILES = $self->backbone_SMILES;
    my %names = %ChemOnomatopist::Chain::Monocycle::names;

    # Check the preserved names
    if( $self->length == 5 && $self->number_of_double_bonds &&
        exists $five_membered_aromatic_single_heteroatom{join( '', $self->heteroatoms )} ) {
        return ChemOnomatopist::Name->new( $five_membered_aromatic_single_heteroatom{join( '', $self->heteroatoms )} );
    }
    return ChemOnomatopist::Name->new( $names{$SMILES} ) if exists $names{$SMILES};

    # Check for aromatic notation
    if( $SMILES =~ /:/ ) {
        $SMILES =~ s/([a-z]):/'' . uc( $1 )/ge;
        $SMILES =~ s/\[([a-z]{1,2})\]:/'[' . uc( $1 ) . ']'/ge;
        for my $SMILES_for_name (keys %names) {
            next unless $SMILES_for_name =~ /=/;
            my $name = $names{$SMILES_for_name};
            $SMILES_for_name =~ s/=//g;
            return ChemOnomatopist::Name->new( $name ) if $SMILES eq $SMILES_for_name;
        }
    }

    # Check for annulenes
    if( $self->is_hydrocarbon && $self->is_aromatic &&
        $self->length =~ /^(4|6|8|10|12|14|16)$/ ) {
        return ChemOnomatopist::Name->new(
                    'cyclo' .
                    ChemOnomatopist::IUPAC_numerical_multiplier( $self->length, 1 ) .
                    ChemOnomatopist::IUPAC_numerical_multiplier( $self->length / 2, 1 ) . 'ene' );
    }

    # Check for cycloalkanes
    if( $self->is_hydrocarbon ) {
        my $name = ChemOnomatopist::Name::Part::NondetachablePrefix->new( 'cyclo' )->to_name;
        $name .= $self->SUPER::suffix;
        return $name;
    }

    if( $self->is_Hantzsch_Widman ) {
        # Hantzsch-Widman names (BBv2 P-22.2.2.1)

        # Collect the types of heteroatoms and their attachment positions
        my %heteroatoms;
        my @vertices = $self->vertices;
        for my $i (0..$#vertices) {
            next if ChemOnomatopist::is_element( $vertices[$i], 'C' );
            my $symbol = ucfirst $vertices[$i]->{symbol};
            $heteroatoms{$symbol} = [] unless $heteroatoms{$symbol};
            push @{$heteroatoms{$symbol}}, $i;
        }

        my $least_senior_element;
        my @heteroatom_locants;
        for my $element (sort { $elements{$a}->{seniority} <=> $elements{$b}->{seniority} }
                              keys %heteroatoms) {
            push @heteroatom_locants, @{$heteroatoms{$element}};
            $least_senior_element = $element;
        }

        my $name = ChemOnomatopist::Name->new;
        unless(  @heteroatom_locants == 1 ||
                (scalar keys %heteroatoms == 1 && @heteroatom_locants == $self->length - 1) ) {
            # Locants are omitted according to BBv2 P-22.2.2.1.7
            $name->append_locants( map { $_ + 1 } @heteroatom_locants );
        }

        for my $element (sort { $elements{$a}->{seniority} <=> $elements{$b}->{seniority} }
                              keys %heteroatoms) {
            if( @{$heteroatoms{$element}} > 1 ) {
                $name->append_multiplier( ChemOnomatopist::IUPAC_numerical_multiplier( scalar @{$heteroatoms{$element}} ) );
            }
            $name->append_element( exists $elements{$element}->{HantzschWidman}
                                        ? $elements{$element}->{HantzschWidman}
                                        : $elements{$element}->{prefix} );
        }
        $name->[-1] =~ s/a$//;

        if(      $self->length <= 5 ) {
            my @stems = ( 'ir', 'et', 'ol' );
            $name .= ChemOnomatopist::Name::Part::Stem->new( $stems[$self->length - 3] );
            if( $self->is_saturated ) {
                $name .= $heteroatoms{N} ? 'idine' : 'ane';
            } elsif( $self->length == 3 ) {
                $name .= $heteroatoms{N} ? 'ene' : 'ine';
            } else {
                $name .= 'e';
            }
            return $name;
        } elsif( $self->length == 6 ) {
            if(      ($elements{$least_senior_element}->{seniority} >= 5 &&
                      $elements{$least_senior_element}->{seniority} <= 8) || $least_senior_element eq 'Bi' ) {
                $name .= $self->is_saturated ? 'ane' : 'ine';
            } elsif( ($elements{$least_senior_element}->{seniority} >= 16 &&
                      $elements{$least_senior_element}->{seniority} <= 19) || $least_senior_element eq 'N' ) {
                $name .= $self->is_saturated ? 'inane' : 'ine';
            } else {
                $name .= $self->is_saturated ? 'inane' : 'inine';
            }
            return $name;
        } elsif( $self->length >= 7 ) {
            my @stems = ( 'ep', 'oc', 'on', 'ec' );
            $name .= ChemOnomatopist::Name::Part::Stem->new( $stems[$self->length - 7] );
            $name .= $self->is_saturated ? 'ane' : 'ine';
            return $name;
        }
    }

    my $name = ChemOnomatopist::Name->new( 'cyclo' );
    $name .= $self->SUPER::suffix;
    return $name;
}

sub needs_indicated_hydrogens() { '' }

sub needs_multiple_bond_locants()
{
    my( $self ) = @_;
    return 1 if $self->number_of_charges;
    return 1 if $self->number_of_branches || $self->parent;
    return 1 if scalar( uniq map { $_->{symbol} } $self->vertices ) > 1;
    return $self->number_of_multiple_bonds > 1 && $self->number_of_multiple_bonds < $self->length;
}

sub needs_substituent_locants()
{
    my( $self ) = @_;
    # BBv2 P-14.3.4.2 (c): monosubstituted homogeneous cycles do not need locants
    return '' if $self->is_homogeneous && $self->number_of_branches == 1;

    if( $self->is_homogeneous && $self->number_of_branches >= $self->max_valence - 1 ) {
        # If there is only one kind of substituents, locants are not needed
        return '' if uniq( map { "$_" } map { @$_ } $self->locant_names ) == 1;
    }
    return 1;
}

sub needs_charge_locants()
{
    my( $self ) = @_;
    return 1 if $self->number_of_charges > 1;
    return 1 if $self->number_of_branches;
    return !&is_homogeneous;
}

sub needs_isotope_locants()
{
    my( $self ) = @_;
    return 1 if $self->number_of_isotopes > 1;
    return 1 if $self->number_of_branches;
    return !&is_homogeneous;
}

sub needs_suffix_locant()
{
    my( $self ) = @_;
    return $self->needs_substituent_locants;
}

sub bonds()
{
    my( $self ) = @_;
    my @bonds = $self->SUPER::bonds;

    my $graph = $self->graph;
    my @vertices = $self->vertices;
    if( $graph->has_edge_attribute( $vertices[0], $vertices[-1], 'bond' ) ) {
        push @bonds, $graph->get_edge_attribute( $vertices[0], $vertices[-1], 'bond' );
    } else {
        push @bonds, '-';
    }
    return @bonds;
}

sub indicated_hydrogens()
{
    my( $self ) = @_;

    my @positions;
    my @vertices = $self->vertices;
    my $graph = $self->graph;
    for my $i (0..$#vertices) {
        # Rough interpretation of BBv2 P-14.7.1 and P-22.2.2.1.4
        next unless $vertices[$i]->{hcount};
        if( ChemOnomatopist::element( $vertices[$i] ) eq 'C' ) {
            if(  $vertices[$i]->{hcount} == 2 ||
                ($graph->degree( $vertices[$i] ) == 3 &&
                 $vertices[$i]->{hcount} == 1 ) ) {
                push @positions, $i;
            }
        } elsif( ChemOnomatopist::element( $vertices[$i] ) eq 'N' &&
            $vertices[$i]->{hcount} && $vertices[$i]->{hcount} == 1 ) {
            push @positions, $i;
        }
    }

    return @positions;
}

# Implemented according to BBv2 P-25.3.3.1.1
sub locants(@)
{
    my $self = shift;
    my @vertices = $self->vertices;
    my $graph = $self->graph->subgraph( @vertices );

    return map { $_ + 1 } @_ if $graph->vertices == $graph->edges;

    my %locant_map;
    my $pos = 0;
    my $letter = 'a';
    for my $i (0..$#vertices) {
        if( $graph->degree( $vertices[$i] ) == 2 ||
            !ChemOnomatopist::element( $vertices[$i] ) ||
             ChemOnomatopist::element( $vertices[$i] ) ne 'C' ) {
            $pos++;
            $locant_map{$i} = $pos;
            $letter = 'a';
        } else {
            $locant_map{$i} = $pos . $letter;
            $letter++;
        }
    }

    return map { $locant_map{$_} } @_;
}

# In aromatic systems Kekule bonds have to be calculated, otherwise seniority rules may fail.
sub number_of_double_bonds()
{
    my( $self ) = @_;
    return $self->SUPER::number_of_double_bonds unless $self->is_aromatic;
    return int( $self->length / 2 );
}

sub number_of_multiple_bonds()
{
    my( $self ) = @_;
    return $self->SUPER::number_of_multiple_bonds unless $self->is_aromatic;
    return $self->number_of_double_bonds;
}

sub number_of_rings()
{
    my( $self ) = @_;
    return scalar $self->cycles if $self->can( 'cycles' );
    return 1;
}

sub indicated_hydrogens_part()
{
    my( $self ) = @_;
    my $part = ChemOnomatopist::Name->new;
    return $part unless $self->number_of_indicated_hydrogens;

    if( $self->needs_indicated_hydrogens ) {
        my @indicated_hydrogens = $self->indicated_hydrogens;
        my $single_H = shift @indicated_hydrogens if @indicated_hydrogens % 2;
        if( @indicated_hydrogens ) {
            if( $self->number_of_indicated_hydrogens < $self->length ) {
                $part->append_locants( $self->locants( @indicated_hydrogens ) );
            }
            $part .= ChemOnomatopist::IUPAC_numerical_multiplier( scalar @indicated_hydrogens );
            $part->[-1] .= 'a' unless $part =~ /[ai]$/;
            $part .= 'hydro';
        }
        if( defined $single_H ) {
            $part->append_locants( map { $_ . 'H' } $self->locants( $single_H ) );
        }
    }

    return $part;
}

sub _cmp_instances
{
    my( $A, $B ) = @_;

    # BBv2 P-44.2.1 (a)
    if( $A->is_heterocycle <=> $B->is_heterocycle ) {
        return $B->is_heterocycle <=> $A->is_heterocycle;
    }

    # BBv2 P-44.2.1 (b)
    if( set( $A->heteroatoms )->has( 'N' ) <=> set( $B->heteroatoms )->has( 'N' ) ) {
        return set( $B->heteroatoms )->has( 'N' ) <=> set( $A->heteroatoms )->has( 'N' );
    }

    # BBv2 P-44.2.1 (e)
    return $B->length <=> $A->length if $A->length <=> $B->length;

    # BBv2 P-44.2.1 (f)
    return scalar( $B->heteroatoms ) <=> scalar( $A->heteroatoms );
}

1;
