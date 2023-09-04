package ChemOnomatopist::Group::Hydroperoxide;

use strict;
use warnings;

# ABSTRACT: Hydroperoxide group
# VERSION

use ChemOnomatopist::Elements qw( %elements );
use List::Util qw( all );

use parent ChemOnomatopist::Group::;

sub new
{
    my( $class, $carbon, @atoms ) = @_;
    return bless { C => $carbon, atoms => \@atoms }, $class;
}

sub element() { return $_[0]->{atoms}[0]{symbol} }

sub prefix
{
    my( $self ) = @_;

    my @elements = map { ChemOnomatopist::element( $_ ) } @{$self->{atoms}};

    return 'hydroperoxy' if all { $_ eq 'O' } @elements;
    return 'disulfanyl'  if all { $_ eq 'S' } @elements;

    my $name = '';
    for my $element (reverse @elements) { # FIXME: Incomplete
        if( $element eq 'O' ) {
            $name .= 'hydr' unless $name;
            $name .= 'oxy';
        }
        $name .= 'sulfanyl' if $element eq 'S';
        $name .= 'selanyl'  if $element eq 'Se';
        $name .= 'tellanyl' if $element eq 'Te';
    }
    return $name;
}

sub suffix
{
    my( $self ) = @_;

    my @elements = map { ChemOnomatopist::element( $_ ) } @{$self->{atoms}};

    return 'peroxol' if all { $_ eq 'O' } @elements;

    my $name;
    if( $elements[0] eq $elements[1] ) {
        $name = 'di' . $elements{$elements[0]}->{prefix};
        $name =~ s/a$/o/;
    } else {
        $name = '-' . join( '', @elements ) . '-' .
                join '', sort map { s/a$/o/; $_ } map { $elements{$_}->{prefix} } grep { $_ ne 'O' } @elements;
    }
    return $name . 'peroxol';
}

1;
