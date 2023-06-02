package ChemOnomatopist::Chain::Circular;

use strict;
use warnings;

# ABSTRACT: Circular chain
# VERSION

use ChemOnomatopist;
use ChemOnomatopist::Chain; # FIXME: Not sure why it is needed
use ChemOnomatopist::Elements qw( %elements );
use ChemOnomatopist::Group::Monocycle;
use ChemOnomatopist::Util::SMILES qw( cycle_SMILES );
use List::Util qw( all any );

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
        push @candidates, cycle_SMILES( $self->{graph}, @vertices );
        push @vertices, shift @vertices;
    }
    @vertices = reverse @vertices;
    for (0..$#vertices) {
        push @candidates, cycle_SMILES( $self->{graph}, @vertices );
        push @vertices, shift @vertices;
    }

    my( $SMILES ) = sort @candidates;
    return $SMILES;
}

sub is_aromatic()
{
    my( $self ) = @_;
    return 1 if join( '', $self->bonds ) =~ /^((-=)+|(=-)+)$/; # FIXME: Simple aromaticity detection
    return 1 if all { $_ eq ':' } $self->bonds;
    return '';
}

sub is_homogeneous()
{
    my( $self ) = @_;

    my @elements = map { $_->{symbol} } $self->vertices;
    my @bonds = $self->bonds;

    my( $element ) = @elements;
    my( $bond ) = @bonds;

    return '' if any { $_ ne $element } @elements;
    return 1 if join( '', @bonds ) =~ /^((-=)+|(=-)+)$/; # FIXME: Simple aromaticity detection

    return '' if any { $_ ne $bond } @bonds;
    return 1;
}

# Returns undef for cycles not needing special treatment.
# Those can be handled by ChemOnomatopist::get_mainchain_name().
sub name()
{
    my( $self ) = @_;

    my $graph = $self->graph;
    my $SMILES = $self->backbone_SMILES;
    my %names = %ChemOnomatopist::Group::Monocycle::names;

    # Check the preserved names
    return $names{$SMILES} if exists $names{$SMILES};

    # Check for aromatic notation
    if( $SMILES =~ /:/ ) {
        $SMILES =~ s/([a-z]):/'' . uc( $1 )/ge;
        $SMILES =~ s/\[([a-z]{1,2})\]:/'[' . uc( $1 ) . ']'/ge;
        for my $SMILES_for_name (keys %names) {
            next unless $SMILES_for_name =~ /=/;
            my $name = $names{$SMILES_for_name};
            $SMILES_for_name =~ s/=//g;
            return $name if $SMILES eq $SMILES_for_name;
        }
    }

    # Check for annulenes
    if( $self->is_hydrocarbon && $self->is_aromatic &&
        $self->length =~ /^(4|6|8|10|12|14|16)$/ ) {
        return 'cyclo' .
               ChemOnomatopist::IUPAC_numerical_multiplier( $self->length, 1 ) .
               ChemOnomatopist::IUPAC_numerical_multiplier( $self->length / 2, 1 ) . 'ene';
    }

    # Check for cycloalkanes
    if( $self->is_hydrocarbon ) {
        return 'cyclo' . ChemOnomatopist::unbranched_chain_name( $self );
    }

    if( $self->length >= 3 && $self->length <= 10 &&
        any { $_->{symbol} =~ /^[cC]$/ } $self->vertices ) {
        # Hantzsch-Widman names (BBv2 P-22.2.2.1)

        # Select the best numbering for heteroatoms
        my @chains;
        my @vertices = $self->vertices;
        my $cycle = $self->{graph}->subgraph( \@vertices ); # TODO: Add attributes
        for (0..$#vertices) {
            push @chains,
                 ChemOnomatopist::Chain::Circular->new( $cycle, @vertices );
            push @vertices, shift @vertices;
        }
        @vertices = reverse @vertices;
        for (0..$#vertices) {
            push @chains,
                 ChemOnomatopist::Chain::Circular->new( $cycle, @vertices );
            push @vertices, shift @vertices;
        }
        # CHECKME: Additional rules from ChemOnomatopist::filter_chains() might still be needed
        my( $chain ) = sort { _cmp( $a, $b ) } @chains;
        @vertices = $chain->vertices;

        # Collect the types of heteroatoms and their attachment positions
        my %heteroatoms;
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
            $name->append_element( exists $elements{$element}->{HantzschWidman} ? $elements{$element}->{HantzschWidman} : $elements{$element}->{prefix} );
        }
        $name->{name} =~ s/a$//;

        if(      $self->length <= 5 ) {
            my @stems = ( 'ir', 'et', 'ol' );
            $name .= $stems[$self->length - 3];
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
                return $name . $self->is_saturated ? 'ane' : 'ine';
            } elsif( ($elements{$least_senior_element}->{seniority} >= 16 &&
                      $elements{$least_senior_element}->{seniority} <= 19) || $least_senior_element eq 'N' ) {
                return $name . $self->is_saturated ? 'inane' : 'ine';
            } else {
                return $name . $self->is_saturated ? 'inane' : 'inine';
            }
        } elsif( $self->length >= 7 ) {
            my @stems = ( 'ep', 'oc', 'on', 'ec' );
            $name .= $stems[$self->length - 7];
            $name .= $self->is_saturated ? 'ane' : 'ine';
            return $name;
        }
    }

    return 'cyclo' . ChemOnomatopist::unbranched_chain_name( $self );
}

# sub branch_positions() # TODO: Maybe need to add 1 to all returned positions?

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

sub _disconnected_chain_graph()
{
    my( $self ) = @_;

    my $graph = $self->graph->copy; # FIXME: Call our copy() to preserve bonds
    $graph->delete_cycle( $self->vertices );

    return $graph;
}

sub _cmp
{
    my( $A, $B ) = @_;

    my @A_heteroatoms = $A->heteroatoms;
    my @A_positions = $A->heteroatom_positions;

    @A_positions = map { $A_positions[$_] }
                       sort { $elements{$A_heteroatoms[$a]}->{seniority} <=>
                              $elements{$A_heteroatoms[$b]}->{seniority} } 0..$#A_positions;

    my @B_heteroatoms = $B->heteroatoms;
    my @B_positions = $B->heteroatom_positions;

    @B_positions = map { $B_positions[$_] }
                       sort { $elements{$B_heteroatoms[$a]}->{seniority} <=>
                              $elements{$B_heteroatoms[$b]}->{seniority} } 0..$#B_positions;

    return ChemOnomatopist::cmp_arrays( \@A_positions, \@B_positions );
}

1;
