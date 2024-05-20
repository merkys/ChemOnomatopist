package ChemOnomatopist::Isotope;

# ABSTRACT: Isotope
# VERSION

use strict;
use warnings;

use Chemistry::Isotope qw( isotope_abundance );

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
    my $abundance = isotope_abundance( $self->element );
    my( $most_abundant ) = sort { $abundance->{$b} <=> $abundance->{$a} }
                                keys %$abundance;
    return $most_abundant + 0;
}

sub element()       { $_[0]->{element} }
sub index()         { $_[0]->{index} }
sub locant()        { $_[0]->{locant} }
sub mass_number()   { $_[0]->{mass_number} }

1;
