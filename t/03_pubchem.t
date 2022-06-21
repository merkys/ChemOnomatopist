#!/usr/bin/perl

use strict;
use warnings;

use ChemOnomatopist;
use Test::More;

if( !$ENV{EXTENDED_TESTING} ) {
    plan skip_all => "Skip \$ENV{EXTENDED_TESTING} is not set\n";
}

open( my $inp, '<', 't/PubChemData' ) or die;

my @cases;
while (my $line = <$inp>) {
    my @fields = split /\s+/, $line;
    push @cases, { id => $fields[0], iupac => $fields[1], smiles => $fields[2] };
}
close $inp;

plan tests => scalar @cases;

for my $case (@cases) {
    is( ChemOnomatopist::get_name( $case->{smiles} ),
        $case->{iupac},
        'ID ' . $case->{id} );
}
