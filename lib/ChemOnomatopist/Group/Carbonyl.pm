package ChemOnomatopist::Group::Carbonyl;

use strict;
use warnings;

# ABSTRACT: Carbonyl group
# VERSION

use parent ChemOnomatopist::Group::;

sub is_oxygen { return 1 }

sub prefix { return 'oxo' };
sub suffix { return 'one' };

1;
