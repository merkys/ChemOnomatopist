package ChemOnomatopist::Charge;

# ABSTRACT: Charged atom
# VERSION

use strict;
use warnings;

sub new
{
    my( $class, $charge, $index, $locant ) = @_;
    return bless { charge => $charge,
                   index => $index,
                   locant => $locant }, $class;
}

sub charge() { $_[0]->{charge} }
sub index()  { $_[0]->{index} }
sub locant() { $_[0]->{locant} }

1;
