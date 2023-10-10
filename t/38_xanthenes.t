#!/usr/bin/perl

use strict;
use warnings;

use ChemOnomatopist;
use Test::More;

my @cases = (
    # From BBv2 P-25.2.1
    { smiles => 'C1=CC=CC2=NC3=CC=CC=C3C=C12', iupac => 'acridine' },
    { smiles => 'C1=CC=CC=2OC3=CC=CC=C3CC12', iupac => '9H-xanthene' },

    # From BBv2 P-25.2.2.2
    { smiles => 'C1=CC=CC=2OC3=CC=CC=C3OC12', iupac => 'oxanthrene' },
    { smiles => 'C1=CC=CC2=NC3=CC=CC=C3N=C12', iupac => 'phenazine' },

    { smiles => 'C1=CC=CC=2OC3=CC=CC=C3C(C12)=O', iupac => '9H-xanthen-9-one' }, # From Wikipedia Xanthone
    { smiles => 'OC1=CC(=CC=2OC3=CC(=CC(=C3C(C12)=O)C)OC)OC', iupac => '1-hydroxy-3,6-dimethoxy-8-methyl-9H-xanthen-9-one' }, # From Wikipedia Lichexanthone
    { smiles => 'C1=CC=CC=2SC3=CC=CC=C3C(C12)=O', iupac => '9H-thioxanthen-9-one' }, # From Wikipedia Thioxanthone
    { smiles => 'CC(C)C1=CC=CC=2C(C3=CC=CC=C3SC12)=O', iupac => '4-(propan-2-yl)-9H-thioxanthen-9-one' }, # From Wikipedia Isopropylthioxanthone

    { smiles => 'ClC1=CC=2OC3=CC(=C(C=C3OC2C=C1Cl)Cl)Cl', iupac => '2,3,7,8-tetrachlorooxanthrene', AUTHOR => 1 }, # From Wikipedia 2,3,7,8-Tetrachlorodibenzodioxin
    { smiles => 'ClC=1C(=C(C(=C2OC=3C(=C(C(=C(C3OC12)Cl)Cl)Cl)Cl)Cl)Cl)Cl', iupac => 'octachlorooxanthrene', AUTHOR => 1 }, # From Wikipedia Octachlorodibenzodioxin
    { smiles => 'ClC1=C(C(=C(C=2OC3=C(C(=C(C=C3OC12)Cl)Cl)Cl)Cl)Cl)Cl', iupac => '1,2,3,4,6,7,8-heptachlorooxanthrene', AUTHOR => 1 }, # From Wikipedia Heptachlorodibenzo-p-dioxin

    { smiles => 'C1CCN(CC1)CCCC2C3=CC=CC=C3SC4=CC=CC=C24', iupac => '1-[3-(9H-thioxanthen-9-yl)propyl]piperidine' }, # PubChem 155569813
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
