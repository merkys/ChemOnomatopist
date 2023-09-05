package ChemOnomatopist::Group::SulfinicAcid;

use strict;
use warnings;

# ABSTRACT: Sulfinic acid group
# VERSION

use parent ChemOnomatopist::Group::;

sub element() { return 'S' }

# From BBv2 P-65.3.0
sub prefix() { return 'sulfino' }
sub suffix() { return 'sulfinic acid' }

1;
