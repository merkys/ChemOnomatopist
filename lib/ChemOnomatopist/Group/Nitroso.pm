package ChemOnomatopist::Group::Nitroso;

use strict;
use warnings;

# ABSTRACT: Nitroso group
# VERSION

use parent ChemOnomatopist::Group::;

sub prefix { return 'nitroso' }
sub suffix { return '' } # CHECKME: Is it specified anywhere?

1;
