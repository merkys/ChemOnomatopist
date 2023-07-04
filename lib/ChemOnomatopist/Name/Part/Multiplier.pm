package ChemOnomatopist::Name::Part::Multiplier;

use strict;
use warnings;

# ABSTRACT: Multiplier of a chemical name
# VERSION

use parent ChemOnomatopist::Name::Part::;

use overload '""'  => sub { return $_[0]->{value} };

1;
