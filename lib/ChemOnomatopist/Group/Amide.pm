package ChemOnomatopist::Group::Amide;

# ABSTRACT: Amide group
# VERSION

use strict;
use warnings;

use parent ChemOnomatopist::Group::;

use ChemOnomatopist::Name;

sub new
{
    my( $class, $parent, $ketone ) = @_;
    return bless { parent => $parent, ketone => $ketone }, $class;
}

sub element() { 'N' }

sub is_terminal() { 1 }

sub prefix { ChemOnomatopist::Name->new( 'amido' ) }

my %infix = (
    S => 'thio',
    Se => 'seleno',
    Te => 'telluro',
);

sub suffix
{
    my( $self ) = @_;
    my $suffix = ChemOnomatopist::Name->new;
    if( $self->{ketone} && exists $infix{$self->{ketone}->element} ) {
        $suffix .= $infix{$self->{ketone}->element};
    }
    return $suffix . 'amide';
}

1;
