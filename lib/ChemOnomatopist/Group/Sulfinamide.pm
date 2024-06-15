package ChemOnomatopist::Group::Sulfinamide;

# ABSTRACT: Sulfinamide group or its Se/Te equivalent
# VERSION

use strict;
use warnings;

use parent ChemOnomatopist::Group::;

use ChemOnomatopist::Name;

my %prefixes = (
    S  => 'sulfinamide',
    Se => 'seleninamide',
    Te => 'tellurinamide',
);

sub prefix { ChemOnomatopist::Name->new( 'sulfinamido' ) } # FIXME: May be incorrect
sub suffix { ChemOnomatopist::Name->new( $prefixes{$_[0]->element} ) }

1;
