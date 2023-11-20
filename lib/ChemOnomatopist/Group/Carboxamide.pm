package ChemOnomatopist::Group::Carboxamide;

use strict;
use warnings;

# ABSTRACT: Carboxamide chain
# VERSION

use parent ChemOnomatopist::Group::, ChemOnomatopist::Chain::;

use ChemOnomatopist::Chain::Monocycle;
use List::Util qw( first );

sub new
{
    my( $class, $graph, $amide, $chain ) = @_;
    $chain->parent( $amide->{vertices}[2] );
    return bless { graph => $graph,
                   chain => $chain,
                   vertices => [ $amide->vertices, $chain->vertices ] }, $class;
}

sub needs_heteroatom_locants() { return '' }
sub needs_heteroatom_names() { return '' }

sub locants(@)
{
    my $self = shift;
    return map { $_ > 2 ? $_ - 2 : $_ ? '?' : 'N' } @_;
}

# FIXME: This is a source of possible failures
sub prefix() { return 'benzamido' }

sub suffix()
{
    my( $self ) = @_;
    return 'benzamide' if $self->{chain}->is_benzene;

    my $suffix = $self->{chain}->suffix;
    if( !$self->{chain}->isa( ChemOnomatopist::Chain::Monocycle:: ) ||
         $self->{chain}->needs_substituent_locants ) {
        my @vertices = $self->{chain}->vertices;
        my $locant = first { $self->graph->has_edge( $self->{vertices}[2], $vertices[$_] ) }
                           0..$#vertices;
        $suffix->append_locants( $locant + 1 );
    }
    $suffix .= 'carboxamide';
    return $suffix;
}

1;
