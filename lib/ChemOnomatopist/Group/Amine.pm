package ChemOnomatopist::Group::Amine;

use strict;
use warnings;

# ABSTRACT: Amino group
# VERSION

use parent ChemOnomatopist::Group::;

sub element() { return 'N' }

sub is_part_of_chain() { return 1 }
sub is_terminal() { return 1 }

sub prefix { return 'amino' }
sub suffix { return 'amine' }

1;
