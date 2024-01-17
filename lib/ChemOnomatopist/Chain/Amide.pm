package ChemOnomatopist::Chain::Amide;

# ABSTRACT: Amide chain
# VERSION

use strict;
use warnings;

use ChemOnomatopist;
use ChemOnomatopist::Name;
use List::Util qw( first );
use Scalar::Util qw( blessed );

sub AUTOLOAD {
    our $AUTOLOAD;
    my $call = $AUTOLOAD;
    $call =~ s/.*:://;
    return if $call eq 'DESTROY';
    my $self = shift;
    return $self->{chain}->can( $call )->( $self->{chain}, @_ );
}

sub new
{
    my( $class, $graph, $chain, $amide ) = @_;
    return bless { graph => $graph, chain => $chain, amide => $amide };
}

sub vertices()
{
    my $self = shift;
    my @vertices = ( $self->{amide}, $self->{chain}->vertices );
    return @vertices;
}

sub locants(@)
{
    my $self = shift;
    return map { $_ ? $_ : 'N' } @_;
}

sub needs_substituent_locants() { 1 }

1;
