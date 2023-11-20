package ChemOnomatopist::Group::Cyanide;

use strict;
use warnings;

# ABSTRACT: Cyanide group
# VERSION

use parent ChemOnomatopist::Group::;

sub prefix { return 'cyano' }
sub suffix { return 'nitrile' }

1;
