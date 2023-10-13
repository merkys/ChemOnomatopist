package ChemOnomatopist::Chain::Benzamide;

use strict;
use warnings;

# ABSTRACT: Benzamide chain
# VERSION

use parent ChemOnomatopist::Chain::;

use List::Util qw( first );

sub new
{
    my( $class, $graph, $amide, $C, $benzene ) = @_;
    $benzene->parent( $C );
    return bless { graph => $graph,
                   benzene => $benzene,
                   vertices => [ $amide, $C, $benzene->vertices ] }, $class;
}

sub needs_heteroatom_locants() { return '' }
sub needs_heteroatom_names() { return '' }

sub locants(@)
{
    my $self = shift;
    return map { $_ > 1 ? $_ - 1 : $_ ? '?' : 'N' } @_;
}

sub prefix() { return 'benzamido' }

sub suffix()
{
    my( $self ) = @_;
    return 'benz' if $self->{benzene}->is_benzene;

    my $suffix = $self->{benzene}->suffix;
    my @vertices = $self->{benzene}->vertices;
    my $locant = first { $self->graph->has_edge( $self->{vertices}[1], $vertices[$_] ) }
                       0..$#vertices;
    $suffix->append_locants( $locant + 1 );
    $suffix .= 'carbox';
    return $suffix;
}

1;
