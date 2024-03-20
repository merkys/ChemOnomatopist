package ChemOnomatopist::Group::Amine;

# ABSTRACT: Amino group
# VERSION

use strict;
use warnings;

use parent ChemOnomatopist::Group::;

sub element() { 'N' }

sub is_terminal() { 1 }

sub prefix() { 'amino' }
sub suffix() { 'amine' }

1;
