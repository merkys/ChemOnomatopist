package ChemOnomatopist::Chain::Amide;

use strict;
use warnings;

# ABSTRACT: Amide chain
# VERSION

use parent ChemOnomatopist::Chain::;

use ChemOnomatopist;
use ChemOnomatopist::Group::Amide;
use ChemOnomatopist::Name;
use List::Util qw( all );
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

sub bond_locants(@)
{
    my $self = shift;
    return @_;
}

sub suffix()
{
    my( $self ) = @_;
    return '' if $self->length == 1;

    my $name = $self->SUPER::suffix;
    my @vertices = $self->vertices;
    if( all { blessed $_ && $_->isa( ChemOnomatopist::Group::Amide:: ) }
            ( $vertices[0], $vertices[-1] ) ) {
        $name->append_multiplier( 'di' );
    }
    return $name;
}

1;
