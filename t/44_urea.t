#!/usr/bin/perl

use strict;
use warnings;

use ChemOnomatopist;
use Test::More;

my @cases = (
    { smiles => 'CNC(=O)N', iupac => 'methylurea' }, # From BBv3 P-14.3.4.3
    { smiles => 'FN(C(N(F)F)=O)F', iupac => 'tetrafluorourea' }, # From BBv3 P-14.3.4.5
    { smiles => 'NC(=S)N', iupac => 'thiourea' }, # From BBv3 P-15.5.3.4.3
    { smiles => 'NC(=O)N', iupac => 'urea' }, # From BBv3 P-34.1.1.5
    { smiles => 'CN(C(=O)N)N=O', iupac => 'N-methyl-N-nitrosourea', AUTHOR => 1 }, # From BBv3 P-61.5.2 # FIXME: detected as hydrazine

    # From BBv3 P-66.1.6.1.1.2
    { smiles => 'CNC(=O)NC', iupac => 'N,N\'-dimethylurea' },
    { smiles => 'CC(C)=NC(=O)N', iupac => 'N-(propan-2-ylidene)urea', AUTHOR => 1 },
    { smiles => 'C(#N)C(CCSC)NC(=O)NC', iupac => 'N-[1-cyano-3-(methylsulfanyl)propyl]-N\'-methylurea', AUTHOR => 1 },

    # From BBv2
    { smiles => 'CC(CC)NC(=[Se])N', iupac => 'N-(butan-2-yl)selenourea', AUTHOR => 1 }, # From BBv3 P-66.1.6.1.3.1
    { smiles => 'C(N)(=N)NC(=O)N', iupac => 'N-carbamimidoylurea', AUTHOR => 1 },

    # From BBv3 P-66.1.6.1.1.3
    { smiles => 'CNC(=O)NC1=C(C2=CC=CC=C2C=C1)C(=O)O', iupac => '2-[(methylcarbamoyl)amino]naphthalene-1-carboxylic acid', AUTHOR => 1 },
    { smiles => 'C(=O)(NC1=CC=C2C=CC(=CC2=C1)S(=O)(=O)O)NC1=CC=C2C=CC(=CC2=C1)S(=O)(=O)O', iupac => '7,7′-[carbonylbis(azanediyl)]di(naphthalene-2-sulfonic acid)', AUTHOR => 1 },
    { smiles => 'C(N)(=O)NC(C1=CC=CC=C1)=O', iupac => 'N-carbamoylbenzamide', AUTHOR => 1 },
    { smiles => 'C(N)(=O)NS(=O)(=O)C1=CC=CC=C1', iupac => 'N-carbamoylbenzenesulfonamide', AUTHOR => 1 },
    { smiles => 'C(N)(=O)NC(CC1=CC=CC=C1)=O', iupac => 'N-carbamoyl-2-phenylacetamide', AUTHOR => 1 },
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
