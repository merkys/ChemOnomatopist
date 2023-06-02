#!/usr/bin/perl

use strict;
use warnings;

use ChemOnomatopist;
use Test::More;

my @cases = (
    # From BBv2 P-25.0
    { smiles => 'C1=CC=CC2=CC=CC=C12', iupac => 'naphthalene' },

    # From BBv2 P-25.2.1
    { smiles => 'N1=CN=CC2=NC=CN=C12', iupac => 'pteridine' },
    { smiles => 'N1=CN=C2N=CNC2=C1', iupac => 'purine', AUTHOR => 1 },
    { smiles => 'N1N=CC2=CC=CC=C12', iupac => '1H-indazole' },

    { smiles => 'C1=CC=C2C(C(C=CC2=C1)O)O', iupac => '1,2-dihydronaphthalene-1,2-diol', AUTHOR => 1 }, # PubChem 362

    { smiles => 'C1=CC=CC=CC2=CC=CC=CC=C12', iupac => 'octalene' }, # From BBv2 P-25.1.2.3

    # From BBv2 P-25.2.2.4
    { smiles => 'C1=COC=CC2=C1C=CC=C2', iupac => '3-benzoxepine', AUTHOR => 1 },
    { smiles => 'O1C=CC2=C1C=CC=C2', iupac => '1-benzofuran', AUTHOR => 1 },
    { smiles => 'C=1OC=C2C1C=CC=C2', iupac => '2-benzofuran', AUTHOR => 1 },
    { smiles => 'C1=CC=COC=CC=CC=COC=CC=CC2=C1C=CC=C2', iupac => '5,12-benzodioxacyclooctadecine', AUTHOR => 1 },
);

@cases = grep { !exists $_->{AUTHOR} } @cases unless $ENV{AUTHOR_TESTING};
plan skip_all => 'No available cases' unless @cases;
plan tests => scalar @cases;

for my $case (@cases) {
    eval { is ChemOnomatopist::get_name( $case->{smiles} ), $case->{iupac} };
    $@ =~ s/\n$// if $@;
    fail $case->{smiles} . ": $@" if $@;
}
