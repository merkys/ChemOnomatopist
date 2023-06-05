package ChemOnomatopist::Group::Nitroso;

use strict;
use warnings;

# ABSTRACT: Nitroso group
# VERSION

use parent ChemOnomatopist::Group::;

sub prefix { return 'nitroso' }
sub suffix { return undef } # Cannot act as suffix

sub is_prefix_only() { return 1 }

1;
