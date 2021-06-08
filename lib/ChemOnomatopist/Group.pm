package ChemOnomatopist::Group;

use strict;
use warnings;

sub new
{
    my( $class, $self ) = @_;
    return bless $self, $class;
}

# Not carbon by default
sub is_carbon { return '' }

1;
