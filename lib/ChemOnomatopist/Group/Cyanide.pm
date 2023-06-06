package ChemOnomatopist::Group::Cyanide;

use strict;
use warnings;

# ABSTRACT: Cyanide group
# VERSION

use parent ChemOnomatopist::Group::;

sub prefix { return 'cyano' }
sub suffix { return 'cyanide' }

sub is_part_of_chain() { return '' }
sub is_prefix_only() { return 1 } # FIXME: Possibly can act as prefix

1;
