package ChemOnomatopist::Group::Cyanide;

use strict;
use warnings;

# ABSTRACT: Cyanide group
# VERSION

use parent ChemOnomatopist::Group::;

# FIXME: When part of chain, suffix is 'nitrile'; when an attachment - 'carbonitrile'.

sub prefix { return 'cyano' }
sub suffix { return 'nitrile' }

1;
