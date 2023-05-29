package ChemOnomatopist::Chain;

use strict;
use warnings;

# ABSTRACT: Longest chain in a compound
# VERSION

use List::Util qw( sum0 );

sub AUTOLOAD {
    our $AUTOLOAD;
    my $call = $AUTOLOAD;
    $call =~ s/.*:://;
    return if $call eq 'DESTROY';
    if( $call =~ /^number_/ ) {
        return int sum0 map { $_->can( $call )->( $_ ) } @{$_[0]->{halves}};
    } else {
        warn "$call called"; # This may be a source of future problems, it is better to throw a warning here
        return;
    }
}

sub new
{
    my( $class, @halves ) = @_;
    return bless { halves => \@halves }, $class;
}

sub length()
{
    my( $self ) = @_;
    my( $A, $B ) = @{$self->{halves}};
    return $A->length + $B->length - !$A->{other_center};
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

sub heteroatom_positions()
{
    my( $self ) = @_;
    my @half0_positions = $self->{halves}[0]->heteroatom_positions;
    my @half1_positions = $self->{halves}[1]->heteroatom_positions;
    # If path parts start at the same atom, its attachments get duplicated
    @half1_positions = grep { $_ } @half1_positions unless $self->{halves}[0]{other_center};
    return ( map { $self->{halves}[0]->length - $_ - 1 }                                 reverse @half0_positions ),
           ( map { $self->{halves}[1]->length + $_ - !defined $self->{halves}[0]{other_center} } @half1_positions );
}

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
        my $graph = $self->{halves}[0]->{graph};
        my @centers = map { $_->{other_center} } @{$self->{halves}};
        if( $graph->has_edge_attribute( @centers, 'bond' ) ) {
            push @bonds, $graph->get_edge_attribute( @centers, 'bond' );
        } else {
            push @bonds, '-';
        }
    }

    push @bonds, $self->{halves}[1]->bonds;
    return @bonds;
}

sub heteroatoms()
{
    my( $self ) = @_;
    my @vertices = $self->vertices;
    return map { $vertices[$_]->{symbol} } $self->heteroatom_positions;
}

sub locant_names()
{
    my( $self ) = @_;
    return reverse( $self->{halves}[0]->locant_names ),
           $self->{halves}[1]->locant_names;
}

sub multiple_bond_positions()
{
    my( $self ) = @_;
    my @bonds = $self->bonds;
    return grep { $bonds[$_] =~ /^[=#\$]$/ } 0..$#bonds;
}

sub vertices()
{
    my( $self ) = @_;
    my @A = $self->{halves}[0]->vertices;
    my @B = $self->{halves}[1]->vertices;
    # If there is only one center atom, it appears in both chains
    shift @B unless $self->{halves}[0]->{other_center};
    return reverse( @A ), @B;
}

sub number_of_double_bonds()
{
    my( $self ) = @_;
    return scalar grep { $_ eq '=' } $self->bonds;
}

sub number_of_multiple_bonds()
{
    my( $self ) = @_;
    return scalar grep { $_ =~ /^[=#\$]$/ } $self->bonds;
}

# Generators

sub reversed()
{
    my( $self ) = @_;

    return ChemOnomatopist::Chain->new( reverse @{$self->{halves}} );
}

1;
