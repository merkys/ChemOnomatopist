package ChemOnomatopist::Isotope;

# ABSTRACT: Isotope
# VERSION

use strict;
use warnings;

use ChemOnomatopist::Util qw(
    array_frequencies
    cmp_arrays
);
use Chemistry::Isotope qw( isotope_abundance );
use List::Util qw( uniq );

sub new
{
    my( $class, $element, $mass_number, $index, $locant ) = @_;
    return bless { element => $element,
                   mass_number => $mass_number,
                   index => $index,
                   locant => $locant }, $class;
}

sub atomic_number()
{
    my( $self ) = @_;
    my $abundance = isotope_abundance( $self->element );
    my( $most_abundant ) = sort { $abundance->{$b} <=> $abundance->{$a} }
                                keys %$abundance;
    return $most_abundant + 0;
}

sub element()       { $_[0]->{element} }
sub index()         { $_[0]->{index} }
sub locant()        { $_[0]->{locant} }
sub mass_number()   { $_[0]->{mass_number} }

sub cmp_isotope_lists_by_greater_number_of_nuclides_of_higher_atomic_number
{
    my( $A, $B ) = @_;

    my %A_atomic_number_freq = array_frequencies map { $_->atomic_number } @$A;
    my %B_atomic_number_freq = array_frequencies map { $_->atomic_number } @$B;
    my @keys = (keys %A_atomic_number_freq, keys %B_atomic_number_freq);
    for (reverse sort uniq @keys) {
        return  1 if !exists $A_atomic_number_freq{$_};
        return -1 if !exists $B_atomic_number_freq{$_};
        return $A_atomic_number_freq{$_} <=> $B_atomic_number_freq{$_}
            if $A_atomic_number_freq{$_} <=> $B_atomic_number_freq{$_};
    }

    return 0;
}

sub cmp_isotope_lists_by_greater_number_of_nuclides_of_higher_mass_number
{
    my( $A, $B ) = @_;

    my %A_mass_number_freq = array_frequencies map { $_->mass_number } @$A;
    my %B_mass_number_freq = array_frequencies map { $_->mass_number } @$B;
    my @keys = (keys %A_mass_number_freq, keys %B_mass_number_freq);
    for (reverse sort uniq @keys) {
        return  1 if !exists $A_mass_number_freq{$_};
        return -1 if !exists $B_mass_number_freq{$_};
        return $A_mass_number_freq{$_} <=> $B_mass_number_freq{$_}
            if $A_mass_number_freq{$_} <=> $B_mass_number_freq{$_};
    }

    return 0;
}

1;
