package ChemOnomatopist::Group::Amino;

use strict;
use warnings;

# ABSTRACT: Amino group
# VERSION

use parent ChemOnomatopist::Group::;

sub prefix { return 'amino' }
sub suffix { return 'amine' }

1;
