package ChemOnomatopist::Group::SulfonicAcid;

# ABSTRACT: Sulfonic acid group
# VERSION

use strict;
use warnings;

use parent ChemOnomatopist::Group::;

use ChemOnomatopist::Name;

sub element() { 'S' }

# From BBv2 P-65.3.0
sub prefix() { ChemOnomatopist::Name->new( 'sulfo' ) }
sub suffix() { ChemOnomatopist::Name->new( 'sulfonic acid' ) }

1;
