#!/usr/bin/perl

use strict;
use warnings;

use ChemOnomatopist;
use Test::More;

my @cases = (
    { smiles => 'O1C(CCC1)=O', iupac => 'oxolan-2-one' }, # Table 1

    { smiles => 'C=C1C=CC(C=C1C1C(C(N2CCCCC12)=O)C(=O)C1CN(CCC1)C)=O', iupac => '1-(6-methylene-3-oxo-cyclohexa-1,4-dienyl)-2-(1-methyl-piperidine-3-carbonyl)-hexahydro-indolizin-3-one', AUTHOR => 1 }, # Figure 1

    # Table 3
    { smiles => 'CC1=CC=C(C=C1)C1=CC=CC=C1', iupac => '4-methylbiphenyl', AUTHOR => 1 },
    { smiles => 'C1CCCC2C3C4CCCCC4C(C12)C3', iupac => 'tetradecahydro-9,10-methanoanthracene', AUTHOR => 1 },
    { smiles => 'COCCOCCOCCOC', iupac => '2,5,8,11-tetraoxadodecane' },
    { smiles => 'O(CC(=O)O)CC(=O)O', iupac => '2,2\'-oxydiacetic acid', AUTHOR => 1 },

    { smiles => 'CC1CC2(CCC1)CCCCC2', iupac => '2-methylspiro[5.5]undecane', AUTHOR => 1 }, # Figure 4

    # Table 4
    { smiles => 'BrCC(CCC)CCl', iupac => '1-bromo-2-(chloromethyl)pentane' },
    { smiles => 'ClC1C(C1)OC1C(C1)F', iupac => '1-chloro-2-(2-fluorocyclopropoxy)cyclopropane', AUTHOR => 1 },
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
