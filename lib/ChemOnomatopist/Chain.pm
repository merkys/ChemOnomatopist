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
    return $A->length + $B->length - !defined $A->{other_center};
}

sub branch_positions()
{
    my( $self ) = @_;
    my @half0_positions = $self->{halves}[0]->branch_positions;
    my @half1_positions = $self->{halves}[1]->branch_positions;
    # If longest path has odd length, the center atom appears in all chains
    @half1_positions = grep { $_ } @half1_positions if $self->length % 2;
    return ( map { $self->{halves}[0]->length - $_ - 1 }         reverse @half0_positions ),
           ( map { $self->{halves}[1]->length + $_ - $self->length % 2 } @half1_positions );
}

sub locant_names()
{
    my( $self ) = @_;
    return reverse( $self->{halves}[0]->locant_names ),
           $self->{halves}[1]->locant_names;
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

1;
