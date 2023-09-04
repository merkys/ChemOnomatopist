package ChemOnomatopist::Group::Hydroperoxide;

use strict;
use warnings;

# ABSTRACT: Hydroperoxide group
# VERSION

use List::Util qw( all );

use parent ChemOnomatopist::Group::;

sub new
{
    my( $class, $carbon, @atoms ) = @_;
    return bless { C => $carbon, atoms => \@atoms }, $class;
}

sub element() { return $_[0]->{atoms}[0]{symbol} }

sub prefix { return 'hydroperoxy' }

sub suffix
{
    my( $self ) = @_;
    return 'peroxol' if all { $_->{symbol} eq 'O' } @{$self->{atoms}};
    return '-' . join( '', map { $_->{symbol} } @{$self->{atoms}} ) . '-peroxol'; # FIXME: Missing element names
}

1;
