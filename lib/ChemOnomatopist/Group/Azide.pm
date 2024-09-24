package ChemOnomatopist::Group::Azide;

# ABSTRACT: Azide group
# VERSION

use strict;
use warnings;

use parent ChemOnomatopist::Group::;

sub element() { 'N' }

sub is_prefix_only() { 1 }

sub prefix { ChemOnomatopist::Name->new( 'azido' ) }

1;
