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

sub _cmp_instances
{
    my( $A, $B ) = @_;
    return $A->{halide}{symbol} cmp $B->{halide}{symbol};
}

1;
