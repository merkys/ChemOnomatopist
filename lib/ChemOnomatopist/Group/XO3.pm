package ChemOnomatopist::Group::XO3;

use strict;
use warnings;

# ABSTRACT: XO3 group
# VERSION

use parent ChemOnomatopist::Group::;

use ChemOnomatopist::Elements qw( %elements );

sub is_prefix_only() { return 1 }

# Compiled from BBv2 Table 5.1 (P-59.1.9)
sub prefix {
    my( $self ) = @_;
    my $prefix = 'per' . $elements{$self->element}{prefix};
    $prefix =~ s/a$/yl/;
    return $prefix;
}

1;
