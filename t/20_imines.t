#!/usr/bin/perl

use strict;
use warnings;

use ChemOnomatopist;
use Test::More;

my @cases = (
    # From BBv2 P-62.3.1.1
    { smiles => 'C(CCCCC)=N', iupac => 'hexan-1-imine' },
    { smiles => 'CN=CC', iupac => 'N-methylethanimine' },
    { smiles => 'ClC1=CC=C(C=C1)N=CC1=CC=C(C=C1)Cl', iupac => 'N,1-bis(4-chlorophenyl)methanimine', AUTHOR => 1 },
    { smiles => 'ClC1=CC=C(C=C1)C=NC1=CC=C(N)C=C1', iupac => '4-{[(4-chlorophenyl)methylidene]amino}aniline', AUTHOR => 1 },
    { smiles => 'S1C(CCC1)=N', iupac => 'thiolan-2-imine' },
    { smiles => 'C1C(C=CC2=CC=CC=C12)=N', iupac => 'naphthalen-2(1H)-imine', AUTHOR => 1 },

    # From BBv2 P-62.3.1.2
    { smiles => 'N=C1CC(CC(C1)=N)=O', iupac => '3,5-diiminocyclohexan-1-one' },
    { smiles => 'N=C(CC1CC(CCC1)C(=O)O)C', iupac => '3-(2-iminopropyl)cyclohexane-1-carboxylic acid' },
    { smiles => 'N=C1CCC(N1)=O', iupac => '5-iminopyrrolidin-2-one' },
    { smiles => 'N=C1C=CC(C=C1)=O', iupac => '4-iminocyclohexa-2,5-dien-1-one' },

    # From BBv2 P-62.3.1.3
    { smiles => 'CP=N', iupac => '1-methylphosphanimine', AUTHOR => 1 },
    { smiles => 'C[Si](=NC1=CC=CC=C1)C', iupac => '1,1-dimethyl-N-phenylsilanimine', AUTHOR => 1 },
    { smiles => 'CN=[SiH]CC(=O)OC', iupac => 'methyl [(methylimino)silyl]acetate', AUTHOR => 1 },

    { smiles => 'C(C(=N)C(=O)O)C(=O)O', iupac => '2-iminobutanedioic acid' },
    { smiles => 'CC(C)C=N', iupac => '2-methylpropan-1-imine' },
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
