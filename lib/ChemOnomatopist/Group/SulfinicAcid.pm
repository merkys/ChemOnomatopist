package ChemOnomatopist::Group::SulfinicAcid;

# ABSTRACT: Sulfinic acid group
# VERSION

use strict;
use warnings;

use ChemOnomatopist::Group::SulfonicAcid;
use parent ChemOnomatopist::Group::SulfonicAcid::;

use ChemOnomatopist::Name;

# From BBv2 P-65.3.0
sub element_prefix() { ChemOnomatopist::Name->new( 'sulfino' ) }

1;
