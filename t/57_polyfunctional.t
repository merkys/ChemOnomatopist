#!/usr/bin/perl

use strict;
use warnings;

use ChemOnomatopist;
use Test::More;

my @cases = (
    # From BBv3 P-63.7
    { smiles => 'CC(CS)(C)OS', iupac => '2-methyl-2-(sulfanyloxy)propane-1-thiol' },
    { smiles => 'OCCC1=C(C=CC=C1)O', iupac => '2-(2-hydroxyethyl)phenol' },
    { smiles => 'OC1=C(C=CC=C1)C(CO)O', iupac => '1-(2-hydroxyphenyl)ethane-1,2-diol', AUTHOR => 1 },
    { smiles => '[SeH]OCCOO', iupac => '2-(selanyloxy)ethane-1-peroxol' },
    { smiles => 'NCC(C)(OO)C', iupac => '1-amino-2-methylpropane-2-peroxol' },
    { smiles => 'CS(=O)(=O)CCO', iupac => '2-(methanesulfonyl)ethan-1-ol', AUTHOR => 1 },
    { smiles => 'O(O)C1C(CCCC1)(O)OOC1C(CCCC1)=O', iupac => '2-[(2-hydroperoxy-1-hydroxycyclohexyl)peroxy]cyclohexan-1-one' },
    { smiles => 'NCCO', iupac => '2-aminoethan-1-ol' },
    { smiles => 'C(C)OOCCOC', iupac => '1-(ethylperoxy)-2-methoxyethane', AUTHOR => 1 },
    { smiles => 'COCCCSC', iupac => '1-methoxy-3-(methylsulfanyl)propane', AUTHOR => 1 },
    { smiles => 'CSSC(=CCCC)SC', iupac => '1-(methyldisulfanyl)-1-(methylsulfanyl)pent-1-ene', AUTHOR => 1 },
    { smiles => 'CO[Si](CCCS)(OC)OC', iupac => '3-(trimethoxysilyl)propane-1-thiol', AUTHOR => 1 },
    { smiles => 'C(C)SC=C(SCCC)SCCC', iupac => '1-{[2-(ethylsulfanyl)-1-(propylsulfanyl)ethen-1-yl]sulfanyl}propane', AUTHOR => 1 },
    { smiles => 'CC(CC)N(C(C)(CC)O)C(C)CC', iupac => '2-[di(butan-2-yl)amino]butan-2-ol', AUTHOR => 1 },
    { smiles => 'N[SiH](C1(CC(CC1)O)[SiH2]CN)C', iupac => '3-[amino(methyl)silyl]-3-[(aminomethyl)silyl]cyclopentan-1-ol', AUTHOR => 1 },
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
