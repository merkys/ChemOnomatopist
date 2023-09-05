package ChemOnomatopist::Group::SulfonicAcid;

use strict;
use warnings;

# ABSTRACT: Sulfonic acid group
# VERSION

use parent ChemOnomatopist::Group::;

sub element() { return 'S' }

# From BBv2 P-65.3.0
sub prefix() { return 'sulfo' }
sub suffix() { return 'sulfonic acid' }

1;
