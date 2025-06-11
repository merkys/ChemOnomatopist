#!/usr/bin/perl

use strict;
use warnings;

use ChemOnomatopist;
use Test::More;

my @cases = (
    # From BBv3 P-66.1.1.2
    { smiles => 'CS(=O)(=O)N', iupac => 'methanesulfonamide' },
    { smiles => 'CC(CC)S(=O)N', iupac => 'butane-2-sulfinamide' },
    { smiles => 'O1C(=CC=C1)[Se](=O)N', iupac => 'furan-2-seleninamide' },
    { smiles => 'N1(CCCC1)S(=O)(=O)N', iupac => 'pyrrolidine-1-sulfonamide' },

    # From BBv3 P-66.1.1.4.2
    { smiles => 'CNS(=O)C=1C=C2C=CC(=CC2=CC1)C(=O)O', iupac => '6-[(methylamino)sulfinyl]naphthalene-2-carboxylic acid' },

    # Hand-crafted tests, confirmed with OPSIN
    { smiles => 'C(CCCC)S(N)=[Se]', iupac => 'pentane-1-sulfinoselenoamide' },
    { smiles => 'C(CCCC)S(N)=[Te]', iupac => 'pentane-1-sulfinotelluroamide' },
    { smiles => 'C(CCCC)S(=O)(N)=S', iupac => 'pentane-1-sulfonothioamide' },
    { smiles => 'C(CCCC)S(N)(=[Se])=S', iupac => 'pentane-1-sulfonoselenothioamide' },
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
