#!/usr/bin/perl

use strict;
use warnings;

use ChemOnomatopist;
use Test::More;

my @cases = (
    # From BBv2
    { smiles => 'CNC(=O)N', iupac => 'methylurea' },
    { smiles => 'C(#N)C(CCSC)NC(=O)NC', iupac => 'N-[1-cyano-3-(methylsulfanyl)propyl]-N\'-methylurea', AUTHOR => 1 },
    { smiles => 'CC(CC)NC(=[Se])N', iupac => 'N-(butan-2-yl)selenourea', AUTHOR => 1 },
    { smiles => 'C(N)(=N)NC(=O)N', iupac => 'N-carbamimidoylurea', AUTHOR => 1 },
    { smiles => 'CN(C(=O)N)N=O', iupac => 'N-methyl-N-nitrosourea', AUTHOR => 1 },
    { smiles => 'CNC(=O)NC', iupac => 'N,N\'-dimethylurea' },
    { smiles => 'CC(C)=NC(=O)N', iupac => 'N-(propan-2-ylidene)urea', AUTHOR => 1 },
    { smiles => 'FN(C(N(F)F)=O)F', iupac => 'tetrafluorourea' },
    { smiles => 'NC(=S)N', iupac => 'thiourea' },
    { smiles => 'NC(=O)N', iupac => 'urea' },
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
