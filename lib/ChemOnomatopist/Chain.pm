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
        return sum0 map { $_->can( $call )->( $_ ) } @{$_[0]->{halves}};
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
    return $A->length + $B->length - ($A->number_of_centers == 1);
}

sub branch_positions()
{
    my( $self ) = @_;
    my @half0_positions = $self->{halves}[0]->branch_positions;
    my @half1_positions = $self->{halves}[1]->branch_positions;
    # If longest path has odd length, the center atom appears in all chains
    @half1_positions = grep @half1_positions if $self->length % 2;
    return ( map { $self->{halves}[0]->length - $_ - 1 } reverse @half0_positions ),
           ( map { $self->{halves}[1]->length + $_ }             @half1_positions );
}

sub locant_positions()
{
    my( $self ) = @_;
    return $self->{halves}[0]->locant_positions_backward +
           $self->{halves}[1]->locant_positions_forward;
}

sub vertices()
{
    my( $self ) = @_;
    my @half1 = $self->{halves}[0]->vertices;
    my @half2 = $self->{halves}[1]->vertices;
    # If longest path has odd length, the center atom appears in all chains
    shift @half1 if $half1[0] eq $half2[0];
    return reverse( @half1 ), @half2;
}

1;
