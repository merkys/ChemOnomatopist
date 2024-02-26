package ChemOnomatopist::Isotope;

# ABSTRACT: Isotope
# VERSION

use strict;
use warnings;

use ChemOnomatopist::Util qw(
    array_frequencies
    cmp_arrays
);
use List::Util qw( uniq );

sub new
{
    my( $class, $element, $atomic_number, $locant ) = @_;
    return bless { element => $element,
                   atomic_number => $atomic_number,
                   locant => $locant }, $class;
}

sub element()       { $_[0]->{element} }
sub atomic_number() { $_[0]->{atomic_number} }
sub locant()        { $_[0]->{locant} }

# TODO: Employ in ChemOnomatopist
sub cmp($$)
{
    my( $A, $B ) = @_;
    $A->element cmp $B->element || $A->atomic_number <=> $B->atomic_number;
}

sub cmp_isotope_lists($$)
{
    my( $A, $B ) = @_;

    # BBv3 P-44.4.1.11.1: Senior set is larger
    return @$B <=> @$A if @$B <=> @$A;

    my %A_atomic_number_freq = array_frequencies map { $_->atomic_number } @$A;
    my %B_atomic_number_freq = array_frequencies map { $_->atomic_number } @$B;

    # BBv3 P-44.4.1.11.2: Senior set has greater number of nuclides of higher atomic number
    # CHECKME: P-44.4.1.11.3 seems to be covered by this as well?
    my @keys = (keys %A_atomic_number_freq, keys %B_atomic_number_freq);
    for (reverse sort uniq @keys) {
        return  1 if !exists $A_atomic_number_freq{$_};
        return -1 if !exists $B_atomic_number_freq{$_};
        return $A_atomic_number_freq{$_} <=> $B_atomic_number_freq{$_}
            if $A_atomic_number_freq{$_} <=> $B_atomic_number_freq{$_};
    }

    my $cmp_result = 0;

    # BBv3 P-44.4.1.11.4: Senior set has lower overall locants
    $cmp_result = cmp_arrays( [ sort map { $_->locant } @$A ],
                              [ sort map { $_->locant } @$B ] );
    return $cmp_result if $cmp_result;

    # BBv3 P-44.4.1.11.5: Senior set has lower locants for nuclides of higher atomic number
    $cmp_result = cmp_arrays( [ map { $_->locant } sort { $b->atomic_number <=> $a->atomic_number } @$A ],
                              [ map { $_->locant } sort { $b->atomic_number <=> $a->atomic_number } @$B ] );
    return $cmp_result if $cmp_result;

    return 0;
}

sub collate($@)
{
    my $needs_isotope_locants = shift;
    my @isotopes = sort { ChemOnomatopist::Isotope::cmp( $a, $b ) } @_;
    my @order;
    my %freq;
    for my $isotope (@isotopes) {
        my $key = $isotope->atomic_number . $isotope->element;
        if( !$freq{$key} ) {
            $freq{$key} = [];
            push @order, $key;
        }
        push @{$freq{$key}}, $isotope;
    }
    my $isotopes = '';
    for my $key (@order) {
        $isotopes .= ',' if $isotopes;
        $isotopes .= join ',', map { $_->locant } @{$freq{$key}} if $needs_isotope_locants;
        $isotopes .= '-' if $needs_isotope_locants;
        $isotopes .= $key;
        $isotopes .= scalar @{$freq{$key}} if @{$freq{$key}} > 1 || $key =~ /H$/;
    }
    return "($isotopes)";
}

1;
