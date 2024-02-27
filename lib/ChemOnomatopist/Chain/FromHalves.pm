package ChemOnomatopist::Chain::FromHalves;

use strict;
use warnings;

# ABSTRACT: Chain formed by two halves
# VERSION

use Chemistry::OpenSMILES qw( is_single_bond );
use Clone qw( clone );
use List::Util qw( all sum sum0 );

use parent ChemOnomatopist::Chain::;

sub new
{
    my( $class, @halves ) = @_;
    # FIXME: This check might be too weak
    die "halves must belong to the same graph\n" if $halves[0]->graph ne $halves[1]->graph;
    return bless { halves => \@halves }, $class;
}

sub graph()
{
    my( $self ) = @_;
    return $self->{halves}[0]->graph;
}

sub halves()
{
    my( $self ) = @_;
    return @{$self->{halves}};
}

sub branch_positions()
{
    my( $self ) = @_;
    my @half0_positions = $self->{halves}[0]->branch_positions;
    my @half1_positions = $self->{halves}[1]->branch_positions;
    # If path parts start at the same atom, its attachments get duplicated
    @half1_positions = grep { $_ } @half1_positions unless $self->{halves}[0]{other_center};
    return ( map { $self->{halves}[0]->length - $_ - 1 }                                 reverse @half0_positions ),
           ( map { $self->{halves}[1]->length + $_ - !defined $self->{halves}[0]{other_center} } @half1_positions );
}

# FIXME: Somewhy this fails '2,8-dioxa-4,5-dithia-11-selenadodecane' test in t/16_heteroatoms.t
#~ sub heteroatom_positions()
#~ {
    #~ my( $self ) = @_;
    #~ my @half0_positions = $self->{halves}[0]->heteroatom_positions;
    #~ my @half1_positions = $self->{halves}[1]->heteroatom_positions;
    #~ # If path parts start at the same atom, its attachments get duplicated
    #~ @half1_positions = grep { $_ } @half1_positions unless $self->{halves}[0]{other_center};
    #~ return ( map { $self->{halves}[0]->length - $_ - 1 }                                 reverse @half0_positions ),
           #~ ( map { $self->{halves}[1]->length + $_ - !defined $self->{halves}[0]{other_center} } @half1_positions );
#~ }

sub most_senior_group_positions()
{
    my( $self ) = @_;
    my @half0_positions = $self->{halves}[0]->most_senior_group_positions;
    my @half1_positions = $self->{halves}[1]->most_senior_group_positions;
    # If path parts start at the same atom, its attachments get duplicated
    @half1_positions = grep { $_ } @half1_positions unless $self->{halves}[0]{other_center};
    return ( map { $self->{halves}[0]->length - $_ - 1 }                                 reverse @half0_positions ),
           ( map { $self->{halves}[1]->length + $_ - !defined $self->{halves}[0]{other_center} } @half1_positions );
}

sub bonds()
{
    my( $self ) = @_;
    my @bonds = reverse $self->{halves}[0]->bonds;

    if( $self->{halves}[0]->{other_center} ) {
        my @centers = map { $_->{other_center} } $self->halves;
        if( $self->graph->has_edge_attribute( @centers, 'bond' ) ) {
            push @bonds, $self->graph->get_edge_attribute( @centers, 'bond' );
        } else {
            push @bonds, '-';
        }
    }

    push @bonds, $self->{halves}[1]->bonds;
    return @bonds;
}

sub locant_names()
{
    my( $self ) = @_;
    return reverse( $self->{halves}[0]->locant_names ),
           $self->{halves}[1]->locant_names;
}

sub isotopes()
{
    my( $self ) = @_;
    # Cloning is needed in order not to affect the original arrays
    my @half0_isotopes = map { clone $_ } $self->{halves}[0]->isotopes;
    my @half1_isotopes = map { clone $_ } $self->{halves}[1]->isotopes;
    @half1_isotopes = grep { $_->{index} } @half1_isotopes unless $self->{halves}[0]{other_center};
    for (@half0_isotopes) {
        $_->{index} = $self->{halves}[0]->length - $_->{index} - 1;
        ( $_->{locant} ) = $self->locants( $_->{index} );
    }
    for (@half1_isotopes) {
        $_->{index} = $self->{halves}[1]->length + $_ - !defined $self->{halves}[0]{other_center};
        ( $_->{locant} ) = $self->locants( $_->{index} );
    }
    return reverse( @half0_isotopes ), @half1_isotopes;
}

# Not sure why this has to be overriden
sub number_of_branches()
{
    my( $self ) = @_;
    return int sum0 map { $_->number_of_branches } $self->halves;
}

sub vertices()
{
    my( $self ) = @_;
    my @A = $self->{halves}[0]->vertices;
    my @B = $self->{halves}[1]->vertices;
    # If there is only one center atom, it appears in both chains
    shift @B unless $self->{halves}[0]->{other_center};
    my @vertices = ( reverse( @A ), @B ); # Otherwise scalar is returned sometimes
    return @vertices;
}

1;
