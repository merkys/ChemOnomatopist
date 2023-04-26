#!/usr/bin/perl

use strict;
use warnings;

use ChemOnomatopist;
use Test::More;

my %unsupported = (
    'C=N' => 'cannot handle atoms other than C, H and O now',
    'c1ccccoc1' => 'cannot handle heterocycles for now',
);

plan tests => scalar keys %unsupported;

for my $case (sort keys %unsupported) {
    my $message;
    eval { ChemOnomatopist::get_name( $case ) };
    $message = $@ if $@;
    $message =~ s/\n$// if $message;
    is( $message, $unsupported{$case} );
}
