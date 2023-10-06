package ChemOnomatopist::Chain::Amine;

use strict;
use warnings;

# ABSTRACT: Amine chain
# VERSION

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

sub suffix()
{
    my( $self ) = @_;
    my $suffix = $self->{chain}->suffix;
    return $suffix unless $self->{chain}->needs_suffix_locant;
    return $suffix if $self->{chain}->length == 2; # Ad-hoc fix for ethanamines

    my $neighbour = first { $self->graph->has_edge( $self->{amine}, $_ ) }
                          $self->{chain}->vertices;

    return $suffix->append_locants( $self->{chain}->locants( $self->vertex_ids( $neighbour ) ) );
}

1;
