package ChemOnomatopist::Chain::FromHalves;

# ABSTRACT: Chain formed from two halves
# VERSION

use strict;
use warnings;

use parent ChemOnomatopist::Chain::;

use Chemistry::OpenSMILES qw( is_single_bond );
use Clone qw( clone );
use List::Util qw( all sum sum0 );

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

sub shares_start() { !defined $_[0]->{halves}[0]{other_center} }

sub branch_positions()
{
    my( $self ) = @_;
    my @half0_positions = $self->{halves}[0]->branch_positions;
    my @half1_positions = $self->{halves}[1]->branch_positions;
    # If path parts start at the same atom, its attachments get duplicated
    @half1_positions = grep { $_ } @half1_positions if $self->shares_start;
    return ( map { $self->{halves}[0]->length - $_ - 1 }           reverse @half0_positions ),
           ( map { $self->{halves}[1]->length + $_ - $self->shares_start } @half1_positions );
}

sub most_senior_group_positions()
{
    my( $self ) = @_;
    my @half0_positions = $self->{halves}[0]->most_senior_group_positions;
    my @half1_positions = $self->{halves}[1]->most_senior_group_positions;
    # If path parts start at the same atom, its attachments get duplicated
    @half1_positions = grep { $_ } @half1_positions if $self->shares_start;
    return ( map { $self->{halves}[0]->length - $_ - 1 }           reverse @half0_positions ),
           ( map { $self->{halves}[1]->length + $_ - $self->shares_start } @half1_positions );
}

sub bonds()
{
    my( $self ) = @_;
    my @bonds = reverse $self->{halves}[0]->bonds;

    if( !$self->shares_start ) {
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
    @half1_isotopes = grep { $_->{index} } @half1_isotopes if $self->shares_start;
    for (@half0_isotopes) {
        $_->{index} = $self->{halves}[0]->length - $_->{index} - 1;
        ( $_->{locant} ) = $self->locants( $_->{index} );
    }
    for (@half1_isotopes) {
        $_->{index} = $self->{halves}[1]->length + $_->{index} - $self->shares_start;
        ( $_->{locant} ) = $self->locants( $_->{index} );
    }
    return my @isotopes = ( reverse( @half0_isotopes ), @half1_isotopes );
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
    shift @B if $self->shares_start;
    my @vertices = ( reverse( @A ), @B ); # Otherwise scalar is returned sometimes
    return @vertices;
}

1;
