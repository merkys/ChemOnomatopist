#!/usr/bin/perl

use strict;
use warnings;

use ChemOnomatopist;
use Test::More;

my %unsupported = (
    # 'N1NNN1'   => 'cannot handle complicated monocycles for now', # FIXME: This is not supported
    'C1=CC=CC=2C=CC=3C(=C4C=CC=CC=C4C3)C12' => 'cannot handle cyclic compounds other than monocycles and monospiro',
    'CC1=C2C(=CC=C1)C(C(CCS2)C(=O)OC)O' => 'cannot determine the parent structure', # PubChem 54384155
);

plan tests => scalar keys %unsupported;

for my $case (sort keys %unsupported) {
    my $message;
    eval { ChemOnomatopist::get_name( $case ) };
    $message = $@ if $@;
    $message =~ s/\n$// if $message;
    is( $message, $unsupported{$case} );
}
