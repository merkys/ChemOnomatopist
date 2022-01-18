package ChemOnomatopist::Group;

use strict;
use warnings;

# ABSTRACT: Chemical group
# VERSION

sub new
{
    my( $class, $self ) = @_;
    return bless $self, $class;
}

# Neither of these by default
sub is_carbon { return '' }
sub is_oxygen { return '' }

1;
