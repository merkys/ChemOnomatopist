package ChemOnomatopist::Group::Carbonitrile;

# ABSTRACT: Carbonitrile group
# VERSION

use strict;
use warnings;

use parent ChemOnomatopist::Group::;

sub prefix { return 'cyano' }
sub suffix { return 'carbonitrile' }

1;
