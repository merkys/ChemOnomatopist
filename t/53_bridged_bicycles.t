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
    { smiles => 'C12[Se]CC(CC1)C2', iupac => '2-selenabicyclo[2.2.1]heptane', AUTHOR => 'flaky' },

    { smiles => 'C12COCC(OC1)CO2', iupac => '3,6,8-trioxabicyclo[3.2.2]nonane', AUTHOR => 'flaky' }, # From BBv3 P-23.3.2.1
    { smiles => 'C12OCSC(CC1)C2', iupac => '2-oxa-4-thiabicyclo[3.2.1]octane' }, # From BBv3 P-23.3.2.2
    { smiles => 'C12C[SH2]CC(CC1)C2', iupac => '3λ4-thiabicyclo[3.2.1]octane' }, # From BBv3 P-23.6.1
    { smiles => 'C12[AsH3][AsH]C(CC1)C2', iupac => '2λ5,3-diarsabicyclo[2.2.1]heptane', AUTHOR => 'flaky' }, # From BBv3 P-23.6.2

    # From BBv3 P-23.7
    { smiles => 'C12CC3CC(CC(C1)C3)C2', iupac => 'adamantane', AUTHOR => 1 },
    { smiles => 'C12C3C4C5C3C1C5C24', iupac => 'cubane', AUTHOR => 1 },
    { smiles => 'N12CCC(CC1)CC2', iupac => '1-azabicyclo[2.2.2]octane' },
    { smiles => 'C12C3C4C2C4C31', iupac => 'tetracyclo[2.2.0.02,6.03,5]hexane', AUTHOR => 1 },
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
