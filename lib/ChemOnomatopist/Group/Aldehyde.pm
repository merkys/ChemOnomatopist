package ChemOnomatopist::Group::Aldehyde;

use strict;
use warnings;

# ABSTRACT: Aldehyde group
# VERSION

use parent ChemOnomatopist::Group::;

sub new
{
    my( $class, $carbon, $ketone ) = @_;
    return bless { C => $carbon, ketone => $ketone }, $class;
}

sub element { return 'C' }

sub is_part_of_chain() { return 1 }

sub prefix { return 'formyl' }

sub suffix()
{
    my( $self ) = @_;
    my $name = $self->{ketone}->suffix;
    $name =~ s/one$/al/;
    return $name;
}

sub multisuffix { return 'carbaldehyde' }
sub suffix_if_cycle_substituent { return 'carbaldehyde' }

1;
