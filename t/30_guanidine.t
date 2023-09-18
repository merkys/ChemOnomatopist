#!/usr/bin/perl

use strict;
use warnings;

use ChemOnomatopist;
use Test::More;

my @cases = (
    { smiles => 'NC(=N)N', iupac => 'guanidine' },

    # From BBv2 P-66.4.1.2.1.2
    { smiles => 'CN(C(=NC1=CC=CC=C1)N(C)C)C', iupac => "N,N,N',N'-tetramethyl-N''-phenylguanidine" },
    { smiles => 'CN(C(=N)NC)C', iupac => "N,N,N'-trimethylguanidine" },

    { smiles => 'NC(N)=NCCCC(=O)O', iupac => '4-[(diaminomethylidene)amino]butanoic acid' }, # From BBv2 P-66.4.1.2.1.3
    { smiles => 'NC(C(=O)O)CCCNC(N)=N', iupac => '2-amino-5-(carbamimidoylamino)pentanoic acid' }, # From BBv2 P-103.1.1.1, arginine
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
