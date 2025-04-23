package ChemOnomatopist::Name::Part::Locants;

use strict;
use warnings;

# ABSTRACT: Locants of a chemical name
# VERSION

use parent ChemOnomatopist::Name::Part::;

sub is_numeric() { $_[0]->{value} =~ /^\]?\d+(,\d+)*[-\]]?$/ }

1;
