#!/usr/bin/perl

use strict;
use warnings;

use ChemOnomatopist;
use Test::More;

my %unsupported = (
    # 'N1NNN1'   => 'cannot handle complicated monocycles for now', # FIXME: This is not supported
    'C1=CC=C2C(=C1)C=CC3=C2N=CC=C3' => 'cannot handle complicated cyclic compounds', # PubChem 9191
    'C1=CC=CC=2C1=C1C=C3C=C4C=CC=CC4=CC3=CC1=CC2' => 'cannot handle complicated cyclic compounds', # PubChem 67470
    'C1=CC2=C3C=CC(=CC=C4C=CC(=C5C=CC(=CC=C1C=C2)C=C5)C=C4)C=C3' => 'cannot handle complicated cyclic compounds', # PubChem 157100544
);

plan tests => scalar keys %unsupported;

for my $case (sort keys %unsupported) {
    my $message;
    eval { ChemOnomatopist::get_name( $case ) };
    $message = $@ if $@;
    $message =~ s/\n$// if $message;
    is( $message, $unsupported{$case} );
}
