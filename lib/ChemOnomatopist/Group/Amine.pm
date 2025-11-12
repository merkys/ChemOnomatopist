package ChemOnomatopist::Group::Amine;

# ABSTRACT: Amino group
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

sub prefix() { ChemOnomatopist::Name->new( 'amino' ) }
sub suffix() { ChemOnomatopist::Name->new( $_[0]->charge && $_[0]->charge == -1 ? 'aminide' : 'amine' ) }

1;
