package ChemOnomatopist::Group::XO3;

# ABSTRACT: XO3 group
# VERSION

use strict;
use warnings;

use ChemOnomatopist::Elements qw( %elements );
use ChemOnomatopist::Name;

use parent ChemOnomatopist::Group::;

sub is_prefix_only() { 1 }

# Compiled from BBv2 Table 5.1 (P-59.1.9)
sub prefix
{
    my( $self ) = @_;
    my $prefix = 'per' . $elements{$self->element}{prefix};
    $prefix =~ s/a$/yl/;
    return ChemOnomatopist::Name->new( $prefix );
}

1;
