package ChemOnomatopist::Group::SulfinicAcid;

# ABSTRACT: Sulfinic acid group
# VERSION

use strict;
use warnings;

use parent ChemOnomatopist::Group::;

# From BBv2 P-65.3.0
sub prefix() { 'sulfino' }
sub suffix() { 'sulfinic acid' }

1;
