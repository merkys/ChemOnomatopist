package ChemOnomatopist::Chain::Carboxamide;

# ABSTRACT: Carboxamide chain
# VERSION

use strict;
use warnings;

use parent ChemOnomatopist::Chain::;

use ChemOnomatopist::Name;
use List::Util qw( first );

sub new
{
    my( $class, $graph, $amide, $C, $chain ) = @_;
    $chain->parent( $C );
    return bless { graph => $graph,
                   chain => $chain,
                   vertices => [ $amide, $C, $chain->vertices ] }, $class;
}

sub needs_heteroatom_locants() { '' }
sub needs_heteroatom_names() { '' }

sub locants(@)
{
    my $self = shift;
    return map { $_ > 1 ? $_ - 1 : $_ ? '?' : 'N' } @_;
}

# FIXME: This is a source of possible failures
sub prefix() { ChemOnomatopist::Name->new( 'benzamido' ) }

sub suffix()
{
    my( $self ) = @_;
    return 'benz' if $self->{chain}->is_benzene;

    my $suffix = $self->{chain}->suffix;
    if( !$self->{chain}->isa( ChemOnomatopist::Chain::Monocycle:: ) ||
         $self->{chain}->needs_substituent_locants ) {
        my @vertices = $self->{chain}->vertices;
        my $locant = first { $self->graph->has_edge( $self->{vertices}[1], $vertices[$_] ) }
                           0..$#vertices;
        $suffix->append_locants( $locant + 1 );
    }
    $suffix .= 'carbo';
    $suffix .= 'x' if $self->{vertices}[0]{ketone}->element eq 'O';
    return $suffix;
}

1;
