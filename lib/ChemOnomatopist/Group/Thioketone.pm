package ChemOnomatopist::Group::Thioketone;

use strict;
use warnings;

# ABSTRACT: Thioketone group
# VERSION

use parent ChemOnomatopist::Group::;

sub prefix { return 'sulfanylidene' };
sub suffix { return 'thione' };

1;
