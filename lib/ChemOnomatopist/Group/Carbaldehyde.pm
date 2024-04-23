package ChemOnomatopist::Group::Carbaldehyde;

# ABSTRACT: Carbaldehyde group
# VERSION

use strict;
use warnings;

use parent ChemOnomatopist::Group::;

use ChemOnomatopist::Elements qw( %elements );
use ChemOnomatopist::Name;

sub new
{
    my( $class, $aldehyde ) = @_;
    return bless { ketone => $aldehyde->{ketone} }, $class;
}

sub element() { 'C' }

sub prefix()
{
    my( $self ) = @_;
    return ChemOnomatopist::Name->new( 'formyl' ) if $self->{ketone}->element eq 'O';

    my $name = ChemOnomatopist::Name->new( 'methane' );
    my $element = $elements{$self->{ketone}->element}->{prefix};
    $element =~ s/a$/oyl/;
    return $name->append_stem( $element );
}

my %suffixes = ( O => '', S => 'othi', Se => 'oselen', Te => 'otellan' );

sub suffix() { 'carb' . $suffixes{$_[0]->{ketone}->element} . 'aldehyde' }

1;
