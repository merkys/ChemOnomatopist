#!/usr/bin/perl

use strict;
use warnings;

use ChemOnomatopist;
use Test::More;

my @cases = (
    # From BBv2 P-24.2.1
    { smiles => 'C1CCCC12CCCCC2', iupac => 'spiro[4.5]decane' },
    { smiles => 'C1CCCC12CCCC2',  iupac => 'spiro[4.4]nonane' },

    # From BBv2 P-24.2.4.1.1
    { smiles => 'C1CCCC12OCCCC2', iupac => '6-oxaspiro[4.5]decane' },
    { smiles => 'C1CCCC12NCCOC2', iupac => '9-oxa-6-azaspiro[4.5]decane' },
    { smiles => 'C1CCCC12CSCNC2', iupac => '7-thia-9-azaspiro[4.5]decane' },

    # From BBv2 P-31.1.5.1.1
    { smiles => 'C1=CCCCC12CC=CCC2', iupac => 'spiro[5.5]undeca-1,8-diene' },
    { smiles => 'C1CCCC12C=CCCC2', iupac => 'spiro[4.5]dec-6-ene' },

    { smiles => 'C1C=CCC12CCCCC2', iupac => 'spiro[4.5]dec-2-ene' }, # From BBv2 P-93.5.3.3

    { smiles => 'C1CC2(C1)C(C2(C#N)C#N)(C#N)C#N', iupac => 'spiro[2.3]hexane-1,1,2,2-tetracarbonitrile' }, # PubChem 263661
    { smiles => 'CC(C)(C)C1CCC2(C(C1)(CCCO)O)OCCO2', iupac => '8-tert-butyl-6-(3-hydroxypropyl)-1,4-dioxaspiro[4.5]decan-6-ol', AUTHOR => 1 }, # PubChem 496485
    { smiles => 'CCCN1C(=O)C(NC(=O)C12CCNCC2)CC(C)C', iupac => '3-(2-methylpropyl)-1-propyl-1,4,9-triazaspiro[5.5]undecane-2,5-dione' }, # PubChem 9856956
    { smiles => 'C1CCOC2(C1)CC(=O)CCO2', iupac => '1,7-dioxaspiro[5.5]undecan-4-one' }, # PubChem 11217490
    { smiles => 'C1CCC2(C1)CCOC(=O)C2', iupac => '8-oxaspiro[4.5]decan-9-one', AUTHOR => 1 }, # PubChem 12733330
    { smiles => 'O1CCOC12CCC(CC2)C(C(C)C)=O', iupac => '1-(1,4-dioxaspiro[4.5]decan-8-yl)-2-methylpropan-1-one' }, # PubChem 13353114

    # This might not be a real compound, nevertheless, locants should probably not be added to it
    { smiles => 'FC1(C(C(C(C(C12C(C(C(C2(F)F)(F)F)(F)F)(F)F)(F)F)(F)F)(F)F)(F)F)F', iupac => 'octadecafluorospiro[4.5]decane', AUTHOR => 1 },
);

@cases = grep { !exists $_->{AUTHOR} } @cases unless $ENV{AUTHOR_TESTING};
plan skip_all => 'No available cases' unless @cases;
plan tests => scalar @cases;

for my $case (@cases) {
    my $ok;
    eval { $ok = is ChemOnomatopist::get_name( $case->{smiles} ), $case->{iupac}, $case->{smiles} };
    $@ =~ s/\n$// if $@;
    fail $case->{smiles} . ": $@" if $@;
    diag 'test supposed to fail with AUTHOR_TESTING' if $case->{AUTHOR} && $ok;
}
