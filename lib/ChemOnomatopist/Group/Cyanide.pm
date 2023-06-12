package ChemOnomatopist::Group::Cyanide;

use strict;
use warnings;

# ABSTRACT: Cyanide group
# VERSION

use parent ChemOnomatopist::Group::;

sub prefix { return 'cyano' }
sub suffix { return 'carbonitrile' }

sub is_part_of_chain() { return '' }

1;
