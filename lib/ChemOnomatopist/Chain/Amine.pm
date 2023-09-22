package ChemOnomatopist::Chain::Amine;

use strict;
use warnings;

# ABSTRACT: Amine chain
# VERSION

use parent ChemOnomatopist::Chain::;

use ChemOnomatopist;
use ChemOnomatopist::Name;
use Scalar::Util qw( blessed );

sub new
{
    my( $class, $graph, $parent, @vertices ) = @_;
    my $self = { vertices => \@vertices, graph => $graph };
    $self->{parent} = $parent if $parent;
    return bless $self, $class;
}

sub locants(@)
{
    my $self = shift;
    return map { $_ ? $_ : 'N' } @_;
}

sub suffix()
{
    my( $self ) = @_;
    return '' if $self->length == 1;
    return $self->SUPER::suffix;
}

1;
