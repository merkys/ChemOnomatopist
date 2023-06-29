package ChemOnomatopist::Group::Amide;

use strict;
use warnings;

# ABSTRACT: Amide group
# VERSION

use parent ChemOnomatopist::Group::;

sub element() { return 'C' }

sub is_part_of_chain() { return 1 }

sub prefix { return 'amido' } # FIXME: Not sure if really
sub suffix { return 'amide' }

1;
