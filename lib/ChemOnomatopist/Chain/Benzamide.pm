package ChemOnomatopist::Chain::Benzamide;

use strict;
use warnings;

# ABSTRACT: Benzamide chain
# VERSION

use parent ChemOnomatopist::Chain::;

sub new
{
    my( $class, $graph, @vertices ) = @_;
    return bless { graph => $graph, vertices => \@vertices }, $class;
}

1;
