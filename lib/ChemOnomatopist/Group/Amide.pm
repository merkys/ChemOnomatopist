package ChemOnomatopist::Group::Amide;

use strict;
use warnings;

# ABSTRACT: Amide group
# VERSION

use parent ChemOnomatopist::Group::;

sub new
{
    my( $class, $parent, $ketone ) = @_;
    return bless { parent => $parent, ketone => $ketone }, $class;
}

sub element() { return 'N' }

sub is_terminal() { return 1 }

sub prefix { return 'amido' }

my %infix = (
    S => 'thio',
    Se => 'seleno',
    Te => 'telluro',
);

sub suffix
{
    my( $self ) = @_;
    my $suffix = '';
    if( exists $infix{$self->{ketone}->element} ) {
        $suffix .= $infix{$self->{ketone}->element};
    }
    return $suffix . 'amide';
}

1;
