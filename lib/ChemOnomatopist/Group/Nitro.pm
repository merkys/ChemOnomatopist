package ChemOnomatopist::Group::Nitro;

use strict;
use warnings;

# ABSTRACT: Nitro group
# VERSION

use parent ChemOnomatopist::Group::;

sub prefix { return 'nitro' }
sub suffix { return '' } # CHECKME: Is it specified anywhere?

1;
