package ChemOnomatopist::Comparable::Array::Numeric;

# ABSTRACT: Comparable array of numbers
# VERSION

use strict;
use warnings;

use ChemOnomatopist::Util qw( cmp_arrays );

use overload '<=>' => \&cmp_arrays;

sub new
{
    my $class = shift;
    return bless \@_, $class;
}

1;
