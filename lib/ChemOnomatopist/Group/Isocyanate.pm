package ChemOnomatopist::Group::Isocyanate;

# ABSTRACT: Isocyanate group
# VERSION

use strict;
use warnings;

use parent ChemOnomatopist::Group::;

use ChemOnomatopist::Name;

sub is_prefix_only() { 1 }

sub prefix() { ChemOnomatopist::Name->new( 'isocyanato' ) }
sub suffix() { ChemOnomatopist::Name->new( 'isocyanate' ) }

1;
