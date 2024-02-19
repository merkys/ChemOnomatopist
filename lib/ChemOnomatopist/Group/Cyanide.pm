package ChemOnomatopist::Group::Cyanide;

# ABSTRACT: Cyanide group
# VERSION

use strict;
use warnings;

use parent ChemOnomatopist::Group::;

sub prefix() { 'cyano' }
sub suffix() { 'nitrile' }
sub multisuffix() { 'carbonitrile' }

1;
