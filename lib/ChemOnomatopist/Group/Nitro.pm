package ChemOnomatopist::Group::Nitro;

# ABSTRACT: Nitro group
# VERSION

use strict;
use warnings;

use parent ChemOnomatopist::Group::;

use ChemOnomatopist::Name;

sub element() { 'N' }

sub is_prefix_only() { 1 }

sub prefix { ChemOnomatopist::Name->new( 'nitro' ) }

1;
