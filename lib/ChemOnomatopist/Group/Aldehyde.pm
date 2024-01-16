package ChemOnomatopist::Group::Aldehyde;

use strict;
use warnings;

# ABSTRACT: Aldehyde group
# VERSION

use parent ChemOnomatopist::Group::;

sub new
{
    my( $class, $ketone ) = @_;
    return bless { ketone => $ketone }, $class;
}

sub element() { 'C' }

sub is_part_of_chain() { 1 }

sub prefix { 'formyl' }

sub suffix()
{
    my( $self ) = @_;
    my $name = $self->{ketone}->suffix;
    $name =~ s/one$/al/;
    return $name;
}

sub multisuffix { 'carbaldehyde' }
sub suffix_if_cycle_substituent { 'carbaldehyde' }

1;
