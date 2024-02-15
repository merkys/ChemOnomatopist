package ChemOnomatopist::Group::Carbaldehyde;

# ABSTRACT: Carbaldehyde group
# VERSION

use strict;
use warnings;

use parent ChemOnomatopist::Group::;

use ChemOnomatopist::Elements qw( %elements );
use ChemOnomatopist::Name::Part::Stem;

sub new
{
    my( $class, $aldehyde ) = @_;
    return bless { ketone => $aldehyde->{ketone} }, $class;
}

sub element() { 'C' }

sub prefix()
{
    my( $self ) = @_;
    return 'formyl' if $self->{ketone}->element eq 'O';

    my $name = ChemOnomatopist::Name::Part::Stem->new( 'methane' )->to_name;
    my $element = $elements{$self->{ketone}->element}->{prefix};
    $element =~ s/a$/oyl/;
    $name .= ChemOnomatopist::Name::Part::Stem->new( $element );
    return $name;
}

my %suffixes = ( O => '', S => 'othi', Se => 'oselen', Te => 'otellan' );

sub suffix() { 'carb' . $suffixes{$_[0]->{ketone}->element} . 'aldehyde' }

1;
