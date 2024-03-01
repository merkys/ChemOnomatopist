#!/usr/bin/perl

use strict;
use warnings;

use ChemOnomatopist;
use Test::More;

my @cases = (
    # From BBv3 P-12.1
    { smiles => 'C12CCC(CC1)C2', iupac => 'bicyclo[2.2.1]heptane' },
    { smiles => 'CC1(C2CCC1CC2)C', iupac => '7,7-dimethylbicyclo[2.2.1]heptane' },

    # From BBv3 P-23.3.1
    { smiles => 'C12COCC(CC1)C2', iupac => '3-oxabicyclo[3.2.1]octane' },
    { smiles => 'C12[Se]CC(CC1)C2', iupac => '2-selenabicyclo[2.2.1]heptane' },

    { smiles => 'C12COCC(OC1)CO2', iupac => '3,6,8-trioxabicyclo[3.2.2]nonane' }, # From BBv3 P-23.3.2.1
    { smiles => 'C12OCSC(CC1)C2', iupac => '2-oxa-4-thiabicyclo[3.2.1]octane' }, # From BBv3 P-23.3.2.2
    { smiles => 'C12C[SH2]CC(CC1)C2', iupac => '3λ4-thiabicyclo[3.2.1]octane' }, # From BBv3 P-23.6.1
    { smiles => 'C12[AsH3][AsH]C(CC1)C2', iupac => '2λ5,3-diarsabicyclo[2.2.1]heptane', AUTHOR => 'flaky' }, # From BBv3 P-23.6.2

    # From BBv3 P-23.7
    { smiles => 'C12CC3CC(CC(C1)C3)C2', iupac => 'adamantane', AUTHOR => 1 },
    { smiles => 'C12C3C4C5C3C1C5C24', iupac => 'cubane', AUTHOR => 1 },
    { smiles => 'N12CCC(CC1)CC2', iupac => '1-azabicyclo[2.2.2]octane' },
    { smiles => 'C12C3C4C2C4C31', iupac => 'tetracyclo[2.2.0.02,6.03,5]hexane', AUTHOR => 1 },

    # From BBv3 P-31.1.4.1
    { smiles => 'C12C=CCC(CC1)C2', iupac => 'bicyclo[3.2.1]oct-2-ene' },
    { smiles => 'C12C=CC(C=C1)CC2', iupac => 'bicyclo[2.2.2]octa-2,5-diene', AUTHOR => 1 },

    # From BBv3 P-31.1.4.2
    { smiles => 'C12CCCCC2=CC1', iupac => 'bicyclo[4.2.0]oct-6-ene', AUTHOR => 1 },
    { smiles => 'C12CCCCCCC(=CCCCC1)C2', iupac => 'bicyclo[6.5.1]tetradec-8-ene' },
    { smiles => 'C12=CC=CC=C2C1', iupac => 'bicyclo[4.1.0]hepta-1,3,5-triene', AUTHOR => 1 },
    { smiles => 'C12CC3=CCCC(CC(CC=4CCCC(C1)C4)C2)C3', iupac => 'tetracyclo[7.7.1.13,7.111,15]nonadeca-3,11(18)-diene', AUTHOR => 1 },

    # From BBv3 P-31.1.4.3
    { smiles => 'C12C#CCCCCCCCC=CC=CCC(CC=C1)C2', iupac => 'bicyclo[14.3.1]icosa-11,13,18-trien-2-yne' },
    { smiles => 'C12C=CCCCCCCCC#CC(CCC1)C2', iupac => 'bicyclo[11.3.1]heptadec-2-en-11-yne' },
    { smiles => 'C12C#CC=CC=CCCC(=CCC1)C2', iupac => 'bicyclo[8.3.1]tetradeca-4,6,10-trien-2-yne' },

    # From BBv3 P-31.1.4.4
    { smiles => 'C12SCC(C=C1)CC2', iupac => '2-thiabicyclo[2.2.2]oct-5-ene', AUTHOR => 'flaky' },
    { smiles => 'C12OCC(C=C1)C2', iupac => '2-oxabicyclo[2.2.1]hept-5-ene' },
    { smiles => 'C12CNCC(C=C1)CC2', iupac => '3-azabicyclo[3.2.2]non-6-ene' },
);

@cases = grep { !exists $_->{AUTHOR} } @cases unless $ENV{AUTHOR_TESTING};
plan skip_all => 'No available cases' unless @cases;
plan tests => scalar @cases;

for my $case (@cases) {
    my $ok;
    eval { $ok = is ChemOnomatopist::get_name( $case->{smiles} ), $case->{iupac}, $case->{smiles} };
    $@ =~ s/\n$// if $@;
    fail $case->{smiles} . ": $@" if $@;
    if( $case->{AUTHOR} && $ok ) {
        diag 'test supposed to fail with AUTHOR_TESTING' .
             ( $case->{AUTHOR} !~ /^1$/ ? ': ' . $case->{AUTHOR} : '' );
    }
}
