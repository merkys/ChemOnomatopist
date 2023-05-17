#!/usr/bin/perl

use strict;
use warnings;

use ChemOnomatopist;
use Test::More;

my %unsupported = (
    'CCCC=CCC'  => undef, # FIXME: Should be an error
    'c1ccccoc1' => 'cannot handle complicated monocycles for now',
);

plan tests => scalar keys %unsupported;

for my $case (sort keys %unsupported) {
    my $message;
    eval { ChemOnomatopist::get_name( $case ) };
    $message = $@ if $@;
    $message =~ s/\n$// if $message;
    is( $message, $unsupported{$case} );
}
