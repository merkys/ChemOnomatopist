package ChemOnomatopist::Group::Amide;

use strict;
use warnings;

# ABSTRACT: Amide group
# VERSION

use parent ChemOnomatopist::Group::;

sub element() { return 'N' }

sub is_terminal() { return 1 }

sub prefix { return 'amido' }
sub suffix { return 'amide' }

1;
