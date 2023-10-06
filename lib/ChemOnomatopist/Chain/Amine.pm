package ChemOnomatopist::Chain::Amine;

use strict;
use warnings;

# ABSTRACT: Amine chain
# VERSION

use ChemOnomatopist;
use ChemOnomatopist::Name;
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
    my( $class, $graph, $chain, $amine ) = @_;
    return bless { graph => $graph, chain => $chain, amine => $amine };
}

sub vertices()
{
    my $self = shift;
    my @vertices = ( $self->{amine}, $self->{chain}->vertices );
    return @vertices;
}

sub locants(@)
{
    my $self = shift;
    return map { $_ ? $_ : 'N' } @_;
}

sub needs_substituent_locants() { return 1 }

1;
