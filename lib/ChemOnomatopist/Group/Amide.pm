package ChemOnomatopist::Group::Amide;

use strict;
use warnings;

# ABSTRACT: Amide group
# VERSION

use parent ChemOnomatopist::Group::;

sub new
{
    my( $class, $parent, $ketone ) = @_;
    return bless { parent => $parent, ketone => $ketone }, $class;
}

sub element() { return 'N' }

sub is_terminal() { return 1 }

sub prefix { return 'amido' }
sub suffix { return 'amide' }

1;
