package ChemOnomatopist::Util;

# ABSTRACT: Generic utilities
# VERSION

use strict;
use warnings;

use Chemistry::Isotope qw( isotope_abundance );
use Exporter;
use Graph::Undirected;
use List::Util qw( first max min pairs );
use Scalar::Util qw( blessed );

use parent Exporter::;

our @EXPORT_OK = qw(
    all_max
    all_min
    array_frequencies
    atomic_number
    circle_permutations
    cmp_arrays
    copy
    element
    zip
);

sub all_max(&@)
{
    my $code = shift;
    return () unless @_;

    my @values = map { $code->( $_ ) } @_;
    my $max = max @values;

    return map { $_[$_] } grep { $values[$_] == $max } 0..$#_;
}

sub all_min(&@)
{
    my $code = shift;
    return () unless @_;

    my @values = map { $code->( $_ ) } @_;
    my $min = min @values;

    return map { $_[$_] } grep { $values[$_] == $min } 0..$#_;
}

sub array_frequencies(@)
{
    my %frequencies;
    for (@_) {
        $frequencies{$_} = 0 unless exists $frequencies{$_};
        $frequencies{$_}++;
    }
    return %frequencies;
}

sub atomic_number($)
{
    my( $element ) = @_;
    $element = element( $element ) if ref $element;
    return unless $element;

    my $abundance = isotope_abundance( $element );
    my $most_abundant = first { 1 }
                        sort  { $abundance->{$b} <=> $abundance->{$a} }
                              keys %$abundance;
    return $most_abundant + 0;
}

sub circle_permutations(@)
{
    my @permutations;
    for (0..$#_) {
        push @permutations, [ @_ ];
        push @_, shift @_;
    }
    @_ = reverse @_;
    for (0..$#_) {
        push @permutations, [ @_ ];
        push @_, shift @_;
    }
    return @permutations;
}

# Takes two arrays and compares them.
# Comparison is first performed numerically on the corresponding array elements.
# Then, if all corresponding elements are equal, arrays are compared by length. 
sub cmp_arrays($$)
{
    my( $a, $b ) = @_;

    for (0..min( scalar( @$a ), scalar( @$b ) )-1) {
        return $a->[$_] <=> $b->[$_] if $a->[$_] <=> $b->[$_];
    }

    return @$a <=> @$b;
}

sub copy($)
{
    my( $object ) = @_;
    die "can only copy Graph now\n" unless blessed $object && $object->isa( Graph::Undirected:: );

    # Graphs have to be copied with the following code as Graph::copy() loses attributes.
    my $graph = $object;
    my $copy = $graph->copy;
    for my $edge ($graph->edges) {
        next unless $graph->has_edge_attributes( @$edge );
        $copy->set_edge_attributes( @$edge, $graph->get_edge_attributes( @$edge ) );
    }
    return $copy;
}

# Derive the chemical element of atom or group representation
sub element
{
    my( $atom_or_group ) = @_;
    return undef unless ref $atom_or_group;

    if( !blessed $atom_or_group ) {
        die "unknown value '$atom_or_group' given for element()\n" unless ref $atom_or_group eq 'HASH';
        return ucfirst $atom_or_group->{symbol};
    }

    if( $atom_or_group->isa( 'Chemistry::Atom' ) ) { # PerlMol Atom
        return $atom_or_group->symbol;
    }

    return $atom_or_group->element;
}

sub zip(@)
{
    die "odd input to zip\n" if @_ % 2;
    my $N = @_ / 2;
    return map { $_[$_], $_[$_ + $N] } 0..$N-1;
}

1;
