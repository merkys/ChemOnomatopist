package ChemOnomatopist::Group::Nitroso;

use strict;
use warnings;

# ABSTRACT: Nitroso group or its analogue
# VERSION

use parent ChemOnomatopist::Group::;

use ChemOnomatopist::Elements qw( %elements );

sub new
{
    my( $class, $carbon, $atom ) = @_;
    return bless { C => $carbon, atom => $atom }, $class;
}

sub is_prefix_only() { return 1 }

# Compiled from BBv2 Table 5.1 (P-59.1.9)
sub prefix {
    my( $self ) = @_;
    return 'nitroso' if $self->{atom}{symbol} eq 'N';
    my $prefix = $elements{$self->{atom}{symbol}}{prefix};
    $prefix =~ s/a$/osyl/;
    return $prefix;
}

1;
