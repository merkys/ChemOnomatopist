package ChemOnomatopist::Chain::Amine;

# ABSTRACT: Amine chain
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
    return map { $_ ? $self->{chain}->locants( $_ - 1 ) : 'N' } @_;
}

sub is_hydrocarbon() { '' }

sub needs_substituent_locants()
{
    my( $self ) = @_;
    return $self->{chain}->length > 0;
}

sub prefix()
{
    my( $self ) = @_;
    my $prefix = $self->{chain}->prefix;
    pop @$prefix if $prefix->[-1] eq 'e'; # FIXME: Dirty
    pop @$prefix if $prefix->[-1] eq 'an';
    $prefix .= 'yl' unless $prefix =~ /[oy]$/;
    $prefix->append_stem( 'amino' );
    return $prefix;
}

sub suffix()
{
    my( $self ) = @_;
    my $suffix = $self->{chain}->suffix;
    return $suffix unless $self->{chain}->needs_suffix_locant;
    return $suffix if $self->{chain}->length == 2; # Ad-hoc fix for ethanamines

    my $neighbour = first { $self->graph->has_edge( $self->{amine}, $_ ) }
                          $self->{chain}->vertices;
    die "cannot perceive connectivity in an amino chain\n" unless defined $neighbour;

    return $suffix->append_locants( $self->{chain}->locants( $self->vertex_ids( $neighbour ) ) );
}

1;
