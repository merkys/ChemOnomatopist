package ChemOnomatopist::Isotope;

# ABSTRACT: Isotope
# VERSION

use strict;
use warnings;

use ChemOnomatopist::Util;

sub new
{
    my( $class, $element, $mass_number, $index, $locant ) = @_;
    return bless { element => $element,
                   mass_number => $mass_number,
                   index => $index,
                   locant => $locant }, $class;
}

sub atomic_number()
{
    my( $self ) = @_;
    return ChemOnomatopist::Util::atomic_number( $self->element );
}

sub element()     { $_[0]->{element} }
sub index()       { $_[0]->{index} }
sub locant()      { $_[0]->{locant} }
sub mass_number() { $_[0]->{mass_number} }

1;
