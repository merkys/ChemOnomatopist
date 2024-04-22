#!/usr/bin/perl

use strict;
use warnings;

use ChemOnomatopist;
use Test::More;

my @cases = (
    # From BBv3 P-73.2.2.1.1
    { smiles => '[CH3+]', iupac => 'methylium' },
    { smiles => 'C1(=CC=CC=C1)[Si+](C1=CC=CC=C1)C1=CC=CC=C1', iupac => 'triphenylsilylium' },
    { smiles => '[CH2+]CC', iupac => 'propylium' },
    { smiles => '[CH+]1CCC1', iupac => 'cyclobutylium' },

    # From BBv3 P-73.2.2.1.2
    { smiles => 'C1(=CC=CC=C1)[S+]', iupac => 'phenylsulfanylium' },
    { smiles => 'CNN=[N+]', iupac => '3-methyltriaz-1-en-1-ylium' },
    { smiles => 'C[Si]([Si+]([Si](C)(C)C)C)(C)C', iupac => 'heptamethyltrisilan-2-ylium' },
    { smiles => 'O1[C+]=CC=C1', iupac => 'furan-2-ylium' },
    { smiles => 'C1CCCC12CC[CH+]CC2', iupac => 'spiro[4.5]decan-8-ylium' },
    { smiles => '[CH2+]C[CH2+]', iupac => 'propane-1,3-bis(ylium)' },
    { smiles => 'CN([N+2])C', iupac => '2,2-dimethylhydrazine-1,1-bis(ylium)' },
    { smiles => 'C[C+2]C', iupac => 'propane-2,2-bis(ylium)' },
    { smiles => '[CH+]1[CH+]C=C1', iupac => 'cyclobut-3-ene-1,2-bis(ylium)' },
    { smiles => '[CH+]1C=CC=C1', iupac => 'cyclopenta-2,4-dien-1-ylium' },
    { smiles => 'O=C1[N+]C(CC1)=O', iupac => '2,5-dioxopyrrolidin-1-ylium' },
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
