#!/usr/bin/perl

use strict;
use warnings;

use ChemOnomatopist;
use Test::More;

my @cases = (
    # BBv3 P-31.2.2
    { smiles => 'C1(=CC=CC2=CC=CC=C12)C1CC2=CC=CC=C2CC1', iupac => '1\',2\',3\',4\'-tetrahydro-1,2\'-binaphthalene', AUTHOR => 1 },
    { smiles => 'N1=CCCCC=C1', iupac => '4,5-dihydro-3H-azepine' },
    { smiles => 'N=1CCCC1', iupac => '3,4-dihydro-2H-pyrrole', AUTHOR => 1 },
    { smiles => 'C1#CC=CC=C1', iupac => '1,2-didehydrobenzene', AUTHOR => 1 },

    # BBv3 P-31.2.3.1
    { smiles => 'C1=CCC=CC1', iupac => 'cyclohexa-1,4-diene' },
    { smiles => 'C1=CCCCC1', iupac => 'cyclohexene' },
    { smiles => 'N1CC=CC=C1', iupac => '1,2-dihydropyridine', AUTHOR => 1 },
    { smiles => 'S1C=CNCCC1', iupac => '4,5,6,7-tetrahydro-1,4-thiazepine' },
    { smiles => 'N1CC=CC=CC1', iupac => '2,7-dihydro-1H-azepine' },
    { smiles => 'P1CCC=C1', iupac => '2,3-dihydro-1H-phosphole', AUTHOR => 1 },

    # BBv3 P-31.2.3.2
    { smiles => 'P1C=CC=C1', iupac => '1H-phosphole' },
    { smiles => 'P1CCCC1', iupac => 'phospholane' },
    { smiles => '[SiH]1=CC=CC=C1', iupac => 'siline' },
    { smiles => '[SiH2]1CCCCC1', iupac => 'silinane' },
    { smiles => 'O1C=CC=C1', iupac => 'furan' },
    { smiles => 'O1CCCC1', iupac => 'oxolane' },
    { smiles => 'N1=CC=CC=C1', iupac => 'pyridine' },
    { smiles => 'N1CCCCC1', iupac => 'piperidine' },
    { smiles => 'O1CC=CC=C1', iupac => '2H-pyran' },
    { smiles => 'O1CCCCC1', iupac => 'oxane' },
    { smiles => 'S1C=CN=CC=C1', iupac => '1,4-thiazepine' },
    { smiles => 'S1CCNCCC1', iupac => '1,4-thiazepane' },

    # BBv3 P-31.2.3.3.1
    { smiles => 'C1CCC2=CC=CC=C12', iupac => '2,3-dihydro-1H-indene' },
    { smiles => 'N1CCC2=CC=CC=C12', iupac => '2,3-dihydro-1H-indole' },
    { smiles => 'C1NCC2=CC=CC=C12', iupac => '2,3-dihydro-1H-isoindole' },
    { smiles => 'O1CCCC2=C1C=CC=C2', iupac => '3,4-dihydro-2H-1-benzopyran' },
    { smiles => 'S1CCCC2=C1C=CC=C2', iupac => '3,4-dihydro-2H-1-benzothiopyran' },
    { smiles => '[Se]1CCCC2=C1C=CC=C2', iupac => '3,4-dihydro-2H-1-benzoselenopyran' },
    { smiles => '[Te]1CCCC2=C1C=CC=C2', iupac => '3,4-dihydro-2H-1-benzotelluropyran' },
    { smiles => 'C1OCCC2=C1C=CC=C2', iupac => '3,4-dihydro-1H-2-benzopyran' },
    { smiles => 'C1SCCC2=C1C=CC=C2', iupac => '3,4-dihydro-1H-2-benzothiopyran' },
    { smiles => 'C1[Se]CCC2=C1C=CC=C2', iupac => '3,4-dihydro-1H-2-benzoselenopyran' },
    { smiles => 'C1[Te]CCC2=C1C=CC=C2', iupac => '3,4-dihydro-1H-2-benzotelluropyran' },

    # BBv3 P-31.2.3.3.2
    { smiles => 'C1C=CCC2=CC=CC=C12', iupac => '1,4-dihydronaphthalene' },
    { smiles => 'C1=CC=CC2=C1C=CCCC2', iupac => '6,7-dihydro-5H-benzo[7]annulene' },
    { smiles => 'C1CCCC2CCCCC12', iupac => 'decahydronaphthalene' },
    { smiles => 'C1CCCC2CC3CCCCC3CC12', iupac => 'tetradecahydroanthracene' },
    { smiles => 'C1CCC2CC3C1C1C4C5C(C(C3CCC2)C1)CCCC(C5)CCC4', iupac => 'octadecahydro-7,14-methano-4,6:8,10-dipropanodicyclohepta[a,d][8]annulene' },
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
