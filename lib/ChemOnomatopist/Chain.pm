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
        return sum0 map { $_->can( $call )->() } @{$_[0]->{halves}};
    } else {
        return;
    }
}

sub new
{
    my( $class, @halves ) = @_;
    return bless { halves => \@halves }, $class;
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
