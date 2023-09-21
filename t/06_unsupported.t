#!/usr/bin/perl

use strict;
use warnings;

use ChemOnomatopist;
use Test::More;

my %unsupported = (
    # 'N1NNN1'   => 'cannot handle complicated monocycles for now', # FIXME: This is not supported
    'C1=CC=CC=2C1=C1C=C3C=C4C=CC=CC4=CC3=CC1=CC2' => 'cannot handle complicated cyclic compounds', # PubChem 67470
    # 'CC1=C2C(=CC=C1)C(C(CCS2)C(=O)OC)O' => 'cannot determine the parent structure', # PubChem 54384155
    # 'COC1=NN=C(C=C1C(=O)O)C2=CC=CC=N2' => 'unknown locant in multicyclic compound', # PubChem 117127049 # This works now, differs only in brackets and is flaky
);

plan tests => scalar keys %unsupported;

for my $case (sort keys %unsupported) {
    my $message;
    eval { ChemOnomatopist::get_name( $case ) };
    $message = $@ if $@;
    $message =~ s/\n$// if $message;
    is( $message, $unsupported{$case} );
}
