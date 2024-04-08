package ChemOnomatopist::Group::Peroxide;

# ABSTRACT: Peroxide group
# VERSION

use strict;
use warnings;

use parent ChemOnomatopist::Group::;

use ChemOnomatopist::Name;

sub new
{
    my( $class, @atoms ) = @_;
    return bless { atoms => \@atoms }, $class;
}

sub element() { 'O' }

sub prefix() { 'peroxy' }
sub suffix() { 'peroxolate' }

1;
