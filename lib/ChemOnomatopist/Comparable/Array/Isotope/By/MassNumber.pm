package ChemOnomatopist::Comparable::Array::Isotope::By::MassNumber;

# ABSTRACT: Comparable array of isotopes
# VERSION

use strict;
use warnings;

use overload '<=>' => \&cmp;

use ChemOnomatopist::Util qw( array_frequencies );
use List::Util qw( uniq );

sub new
{
    my $class = shift;
    return bless \@_, $class;
}

sub cmp
{
    my( $A, $B ) = @_;

    my %A_freq = array_frequencies map { $_->mass_number } @$A;
    my %B_freq = array_frequencies map { $_->mass_number } @$B;

    my @keys = (keys %A_freq, keys %B_freq);
    for (reverse sort uniq @keys) {
        return  1 if !exists $A_freq{$_};
        return -1 if !exists $B_freq{$_};
        return $A_freq{$_} <=> $B_freq{$_} if $A_freq{$_} <=> $B_freq{$_};
    }

    return 0;
}

1;
