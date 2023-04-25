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
    # Minor regularizations for PubChem names:
    $case->{iupac} =~ s/(di|tri|tetra|penta|hepta)(tert-butyl)/$1\($2\)/g;
    $case->{iupac} =~ s/di\((non|heptadec|hentriacont|undec|tridec|docos|icos|tetradec|pentadec)yl\)/di$1yl/g;
    $case->{iupac} =~ s/tetra\(tridecyl\)/tetratridecyl/g;

    is( ChemOnomatopist::get_name( $case->{smiles} ),
        $case->{iupac},
        'ID ' . $case->{id} );
}
