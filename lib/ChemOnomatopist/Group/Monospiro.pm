package ChemOnomatopist::Group::Monospiro;

use strict;
use warnings;

# ABSTRACT: Monospiro compound
# VERSION

use parent ChemOnomatopist::Group::;

use ChemOnomatopist::Chain::Circular;

sub new
{
    my( $class, $graph, @vertices ) = @_;
    return bless { graph => $graph, vertices => \@vertices }, $class;
}

1;
