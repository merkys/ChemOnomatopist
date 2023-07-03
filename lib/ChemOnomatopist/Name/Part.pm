package ChemOnomatopist::Name::Part;

use strict;
use warnings;

# ABSTRACT: Semantic part of a chemical name
# VERSION

sub new
{
    my( $class, $value ) = @_;
    return bless { value => $value }, $class;
}

1;
