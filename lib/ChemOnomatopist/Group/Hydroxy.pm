package ChemOnomatopist::Group::Hydroxy;

use strict;
use warnings;

# ABSTRACT: Hydroxy group
# VERSION

use parent ChemOnomatopist::Group::;

# From BBv2 P-63.1.5
my %prefixes = ( O => 'hydroxy', S => 'sulfanyl', Se => 'selanyl', Te => 'tellanyl' );
my %suffixes = ( O => 'ol', S => 'thiol', Se => 'selenol', Te => 'tellurol' );

sub prefix
{
    my( $self ) = @_;
    return $prefixes{$self->element};
}

sub suffix
{
    my( $self ) = @_;
    return $suffixes{$self->element};
}

sub _cmp_instances
{
    my( $A, $B ) = @_;
    return $A->element cmp $B->element
}

1;
