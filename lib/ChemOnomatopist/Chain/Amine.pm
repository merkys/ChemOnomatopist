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
    return bless { graph => $graph, chain => $chain, amine => $amine }, $class;
}

sub isa
{
    my( $self, $class ) = @_;
    return 1 if $class eq ChemOnomatopist::Chain::Amine::;
    return $self->{chain}->isa( $class );
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

# Need to adjust positions by 1 to accommodate the amino group as the first vertex
sub heteroatom_positions
{
    my( $self ) = @_;
    return map { $_ + 1 } $self->{chain}->heteroatom_positions;
}

sub is_hydrocarbon() { '' }

sub needs_substituent_locants() { $_[0]->{chain}->length > 0 }

sub prefix()
{
    my( $self ) = @_;
    return $self->{amine}->prefix unless $self->length;
    return $self->{chain}->prefix->append_stem( 'amino' );
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
