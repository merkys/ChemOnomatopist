package ChemOnomatopist::Group::Carbonyl;

use strict;
use warnings;

# ABSTRACT: Carbonyl group
# VERSION

use parent ChemOnomatopist::Group::;

use Scalar::Util qw( blessed );

sub is_oxygen { return 1 }

sub prefix { return 'oxo' };
sub suffix { return 'one' };

1;
