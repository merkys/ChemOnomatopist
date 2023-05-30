#!/usr/bin/perl

use strict;
use warnings;

use ChemOnomatopist;
use Test::More;

my %unsupported = (
    # 'N1NNN1'   => 'cannot handle complicated monocycles for now', # FIXME: This is not supported
    'C1=CC=CC2=CC=CC=C12' => 'cannot handle ortho-fused rings for now',
    'C1=CC=CC=2C=CC=3C(=C4C=CC=CC=C4C3)C12' => 'cannot handle cyclic compounds other than monocycles',
    'C1CCCC12CCCCC2' => 'cannot handle monospiro ring systems for now',
);

plan tests => scalar keys %unsupported;

for my $case (sort keys %unsupported) {
    my $message;
    eval { ChemOnomatopist::get_name( $case ) };
    $message = $@ if $@;
    $message =~ s/\n$// if $message;
    is( $message, $unsupported{$case} );
}
