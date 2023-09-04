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

sub prefix { return 'hydroperoxy' }

sub suffix
{
    my( $self ) = @_;
    return 'peroxol' if all { ChemOnomatopist::element( $_ ) eq 'O' } @{$self->{atoms}};

    my @elements = map { ChemOnomatopist::element( $_ ) } @{$self->{atoms}};
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
