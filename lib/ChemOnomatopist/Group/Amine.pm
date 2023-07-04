package ChemOnomatopist::Group::Amine;

use strict;
use warnings;

# ABSTRACT: Amino group
# VERSION

use parent ChemOnomatopist::Group::;

sub element() { return 'N' }

sub prefix { return 'amino' }
sub suffix { return 'amine' }

1;