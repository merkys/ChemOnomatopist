package ChemOnomatopist::Group::Isocyanate;

# ABSTRACT: Isocyanate group
# VERSION

use strict;
use warnings;

use parent ChemOnomatopist::Group::;

use ChemOnomatopist::Name;

# Should be prefix-only as per P-61.8
sub is_prefix_only() { 1 }

my %infix = (
    O  => '',
    S  => 'thio',
    Se => 'seleno',
    Te => 'telluro',
);

sub prefix() { ChemOnomatopist::Name->new( 'iso' . $infix{$_[0]->element} . 'cyanato' ) }
sub suffix() { ChemOnomatopist::Name->new( 'iso' . $infix{$_[0]->element} . 'cyanate' ) }

1;
