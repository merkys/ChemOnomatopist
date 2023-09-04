package ChemOnomatopist::Group::Hydroperoxide;

use strict;
use warnings;

# ABSTRACT: Hydroperoxide group
# VERSION

use parent ChemOnomatopist::Group::;

sub new
{
    my( $class, $carbon, @atoms ) = @_;
    return bless { C => $carbon, atoms => \@atoms }, $class;
}

sub element() { return $_[0]->{atoms}[0]{symbol} }

sub prefix { return 'hydroperoxy' }
sub suffix { return 'peroxol' }

1;
