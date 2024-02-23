package ChemOnomatopist::Isotope;

# ABSTRACT: Isotope
# VERSION

use strict;
use warnings;

sub new
{
    my( $class, $element, $atomic_number, $locant ) = @_;
    return bless { element => $element,
                   atomic_number => $atomic_number,
                   locant => $locant }, $class;
}

sub element()       { $_[0]->{element} }
sub atomic_number() { $_[0]->{atomic_number} }
sub locant()        { $_[0]->{locant} }

1;
