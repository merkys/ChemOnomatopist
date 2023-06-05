package ChemOnomatopist::Group::Nitro;

use strict;
use warnings;

# ABSTRACT: Nitro group
# VERSION

use parent ChemOnomatopist::Group::;

sub prefix { return 'nitro' }
sub suffix { return undef } # Cannot act as suffix

sub is_prefix_only() { return 1 }

1;
