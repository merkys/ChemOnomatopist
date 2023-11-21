package ChemOnomatopist::Group::AcylHalide;

use strict;
use warnings;

# ABSTRACT: Acyl halide group
# VERSION

use ChemOnomatopist::Elements qw( %elements );
use List::Util qw( any );
use Scalar::Util qw( blessed );

use parent ChemOnomatopist::Group::;

sub new
{
    my( $class, $halide ) = @_;
    return bless { halide => $halide }, $class;
}

sub element() { return 'C' }

sub prefix()
{
    my( $self ) = @_;
    my $name = 'carbono';
    $name .= blessed $self->{halide} ? $self->{halide}->prefix : $elements{$self->{halide}{symbol}}->{prefix};
    $name =~ s/a$/idoyl/;
    return $name;
}

sub suffix()
{
    my( $self ) = @_;
    my $name = 'oyl ';
    $name .= blessed $self->{halide} ? $self->{halide}->prefix : $elements{$self->{halide}{symbol}}->{prefix};
    $name =~ s/[ao]$/ide/;
    return $name;
}

sub _cmp_instances
{
    my( $A, $B ) = @_;
    return if any { blessed $_ } ( $A, $B ); # FIXME: Need proper sorting
    return $A->{halide}{symbol} cmp $B->{halide}{symbol};
}

1;
