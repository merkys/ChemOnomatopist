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
    my( $class, $is_center_single, @halves ) = @_;
    return bless { halves => \@halves, is_center_single => $is_center_single }, $class;
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
    my @vertices = reverse $self->{halves}[0]->vertices;
    # If longest path has odd length, the center atom appears in all chains
    pop @vertices if $self->{is_center_single};
    push @vertices, $self->{halves}[1]->vertices;
    return @vertices;
}

1;
