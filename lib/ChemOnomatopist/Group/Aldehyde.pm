package ChemOnomatopist::Group::Aldehyde;

# ABSTRACT: Aldehyde group
# VERSION

use strict;
use warnings;

use parent ChemOnomatopist::Group::;

use ChemOnomatopist::Name;

sub new
{
    my( $class, $ketone ) = @_;
    return bless { ketone => $ketone }, $class;
}

sub element() { 'C' }

sub is_part_of_chain() { 1 }

sub prefix { ChemOnomatopist::Name->new( 'formyl' ) }

sub suffix()
{
    my( $self ) = @_;
    my $name = $self->{ketone}->suffix;
    $name =~ s/one$/al/;
    $name = 'selenal' if $name eq 'selal';
    return $name;
}

sub multisuffix { 'carbaldehyde' }
sub suffix_if_cycle_substituent { 'carbaldehyde' }

1;
