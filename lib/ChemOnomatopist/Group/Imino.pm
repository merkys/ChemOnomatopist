package ChemOnomatopist::Group::Imino;

use strict;
use warnings;

# ABSTRACT: Imino group
# VERSION

use parent ChemOnomatopist::Group::;

sub prefix { return 'imino' }
sub suffix { return 'imine' }

1;
