package ChemOnomatopist::Group::SulfinicAcid;

use strict;
use warnings;

# ABSTRACT: Sulfinic acid group
# VERSION

use parent ChemOnomatopist::Group::;

sub element() { return 'S' }

sub prefix() { return 'sulfino' }
sub suffix() { return 'sulfinic acid' }

1;
