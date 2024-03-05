package ChemOnomatopist::Group::Hydroperoxide;

# ABSTRACT: Hydroperoxide group
# VERSION

use strict;
use warnings;

use parent ChemOnomatopist::Group::;

use ChemOnomatopist::Elements qw( %elements );
use ChemOnomatopist::Name;
use List::Util qw( all );

sub new
{
    my( $class, @atoms ) = @_;
    return bless { atoms => \@atoms }, $class;
}

sub element() { $_[0]->{atoms}[0]{symbol} }

sub prefix()
{
    my( $self ) = @_;

    my @elements = map { ChemOnomatopist::element( $_ ) } @{$self->{atoms}};

    return 'hydroperoxy' if all { $_ eq 'O' } @elements;
    if( all { $_ eq 'S' } @elements ) {
        my $name = ChemOnomatopist::Name->new;
        $name->append_multiplier( 'di' );
        $name->append_stem( 'sulfanyl' );
        return $name;
    }

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

sub suffix()
{
    my( $self ) = @_;

    my @elements = map { ChemOnomatopist::element( $_ ) } @{$self->{atoms}};

    my $name = ChemOnomatopist::Name->new;
    if( all { $_ eq 'O' } @elements ) {
        $name .= 'peroxol';
    } else {
        if( $elements[0] eq $elements[1] ) {
            $name->append_multiplier( 'di' );
            my $element_prefix = $elements{$elements[0]}->{prefix};
            $element_prefix =~ s/a$/o/;
            $name .= $element_prefix;
        } else {
            $name = '-' . join( '', @elements ) . '-' .
                    join '', sort map { s/a$/o/; $_ } map { $elements{$_}->{prefix} } grep { $_ ne 'O' } @elements;
        }
        $name .= 'peroxol';
    }
    return $name;
}

1;
