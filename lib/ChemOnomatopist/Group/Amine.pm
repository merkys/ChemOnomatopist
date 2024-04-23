package ChemOnomatopist::Group::Amine;

# ABSTRACT: Amino group
# VERSION

use strict;
use warnings;

use parent ChemOnomatopist::Group::;

use ChemOnomatopist::Name;

sub element() { 'N' }

sub is_terminal() { 1 }

sub prefix() { ChemOnomatopist::Name->new( 'amino' ) }
sub suffix() { ChemOnomatopist::Name->new( 'amine' ) }

1;
