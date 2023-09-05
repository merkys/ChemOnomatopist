package ChemOnomatopist::Group::AcylHalide;

use strict;
use warnings;

# ABSTRACT: Acyl halide group
# VERSION

use ChemOnomatopist::Elements qw( %elements );

use parent ChemOnomatopist::Group::;

sub new
{
    my( $class, $carbon, $halide ) = @_;
    return bless { C => $carbon, halide => $halide }, $class;
}

sub element() { return 'C' }

# FIXME
sub prefix() { return 'sulfino' }

sub suffix()
{
    my( $self ) = @_;
    my $name = 'oyl ' . $elements{$self->{halide}{symbol}}->{prefix};
    $name =~ s/a$/ide/;
    return $name;
}

1;
