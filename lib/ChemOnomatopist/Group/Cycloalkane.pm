package ChemOnomatopist::Group::Cycloalkane;

use strict;
use warnings;

# ABSTRACT: Cycloalkane group
# VERSION

use parent ChemOnomatopist::Group::;

sub new {
    my( $class, $size ) = @_;
    return bless { size => $size }, $class;
}

1;
