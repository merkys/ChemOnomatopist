#!/usr/bin/perl

use strict;
use warnings;

use ChemOnomatopist;
use Test::More;

my @cases = (
    # From BBv3 P-44.1.1
    { smiles => 'C1(CCCCC1)CCC(=O)O', iupac => '3-cyclohexylpropanoic acid' },
    { smiles => 'C(CC)C=1C=C(C(=O)O)C=CC1', iupac => '3-propylbenzoic acid' },
    { smiles => 'ClCCCCC(CCO)C(C)O', iupac => '3-(4-chlorobutyl)pentane-1,4-diol' },
    { smiles => 'N(N)C(=O)O', iupac => 'hydrazinecarboxylic acid', AUTHOR => 1 },
    { smiles => '[SiH3]CCC(=O)O', iupac => '3-silylpropanoic acid', AUTHOR => 1 },
    { smiles => 'C(C)[SiH2]C(=O)O', iupac => 'ethylsilanecarboxylic acid', AUTHOR => 1 },
    { smiles => 'C(CCC)OCCOCC(SCCSCC(=O)O)(CCSCCSCC(=O)O)COCCOCCCC', iupac => '7,7-bis[(2-butoxyethoxy)methyl]-3,6,10,13-tetrathiapentadecane-1,15-dioic acid', AUTHOR => 1 },
    { smiles => '[SiH2](CC[SiH3])CC[SiH3]', iupac => '[silanediyldi(ethane-2,1-diyl)]bis(silane)', AUTHOR => 1 },
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
