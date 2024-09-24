package ChemOnomatopist::Group::Imine;

# ABSTRACT: Imine group
# VERSION

use strict;
use warnings;

use parent ChemOnomatopist::Group::;

use ChemOnomatopist::Name;

sub new
{
    my( $class, $atom ) = @_;
    return bless { atom => $atom }, $class;
}

sub element() { 'N' }
sub charge() { ChemOnomatopist::charge( $_[0]->{atom} ) }
sub is_terminal() { 1 }

sub needs_multiple_bond_suffix() { '' }

sub prefix() { ChemOnomatopist::Name->new( 'imino' ) }
sub suffix() { ChemOnomatopist::Name->new( 'imine' ) }

1;
