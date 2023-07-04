package ChemOnomatopist::Name::Part::Stem;

use strict;
use warnings;

# ABSTRACT: Stem of a chemical name
# VERSION

use parent ChemOnomatopist::Name::Part::;

use overload '""'  => sub { return $_[0]->{value} };

1;
