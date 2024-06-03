#!/usr/bin/perl

use strict;
use warnings;

use ChemOnomatopist;
use Test::More;

my @cases = (
    # From BBv3 P-67.1.1.1
    { smiles => '[AsH2](O)=O', iupac => 'arsinic acid' },
    { smiles => '[AsH2]O', iupac => 'arsinous acid' },
    { smiles => '[AsH](O)(O)=O', iupac => 'arsonic acid' },
    { smiles => '[AsH](O)O', iupac => 'arsonous acid' },
    { smiles => '[As](O)(O)(O)=O', iupac => 'arsoric acid' },
    { smiles => '[As](O)(O)O', iupac => 'arsorous acid' },
    { smiles => '[NH2+](O)[O-]', iupac => 'azinic acid', AUTHOR => 1 },
    { smiles => '[NH+](O)(O)[O-]', iupac => 'azonic acid', AUTHOR => 1 },
    { smiles => 'N(O)O', iupac => 'azonous acid', AUTHOR => 1 },
    { smiles => '[N+](O)(O)(O)[O-]', iupac => 'nitroric acid', AUTHOR => 1 },
    { smiles => 'N(O)(O)O', iupac => 'azorous acid', AUTHOR => 1 },
    { smiles => 'B(O)(O)O', iupac => 'boric acid', AUTHOR => 1 },
    { smiles => 'BO', iupac => 'borinic acid', AUTHOR => 1 },
    { smiles => 'B(O)O', iupac => 'boronic acid', AUTHOR => 1 },
    { smiles => 'Br(=O)(=O)O', iupac => 'bromic acid', AUTHOR => 1 },
    { smiles => 'Br(=O)O', iupac => 'bromous acid', AUTHOR => 1 },
    { smiles => 'Cl(=O)(=O)O', iupac => 'chloric acid', AUTHOR => 1 },
    { smiles => 'Cl(=O)O', iupac => 'chlorous acid', AUTHOR => 1 },
    { smiles => 'BrO', iupac => 'hypobromous acid', AUTHOR => 1 },
    { smiles => 'ClO', iupac => 'hypochlorous acid', AUTHOR => 1 },
    { smiles => 'FO', iupac => 'hypofluorous acid', AUTHOR => 1 },
    { smiles => 'IO', iupac => 'hypoiodous acid', AUTHOR => 1 },
    { smiles => 'I(=O)(=O)O', iupac => 'iodic acid', AUTHOR => 1 },
    { smiles => 'I(=O)O', iupac => 'iodous acid', AUTHOR => 1 },
    { smiles => '[N+](=O)(O)[O-]', iupac => 'nitric acid', AUTHOR => 1 },
    { smiles => 'N(=O)O', iupac => 'nitrous acid', AUTHOR => 1 },
    { smiles => 'Br(=O)(=O)(=O)O', iupac => 'perbromic acid', AUTHOR => 1 },
    { smiles => 'Cl(=O)(=O)(=O)O', iupac => 'perchloric acid', AUTHOR => 1 },
    { smiles => 'F(=O)(=O)(=O)O', iupac => 'perfluoric acid', AUTHOR => 1 },
    { smiles => 'I(=O)(=O)(=O)O', iupac => 'periodic acid', AUTHOR => 1 },
    { smiles => '[PH2](O)=O', iupac => 'phosphinic acid' },
    { smiles => 'PO', iupac => 'phosphinous acid' },
    { smiles => 'P(O)(O)=O', iupac => 'phosphonic acid' },
    { smiles => 'P(O)O', iupac => 'phosphonous acid' },
    { smiles => 'P(O)(O)(O)=O', iupac => 'phosphoric acid' },
    { smiles => 'P(O)(O)O', iupac => 'phosphorous acid' },
    { smiles => '[Se](O)(O)(=O)=O', iupac => 'selenic acid', AUTHOR => 1 },
    { smiles => '[Se](O)(O)=O', iupac => 'selenous acid', AUTHOR => 1 },
    { smiles => '[Si](O)(O)(O)O', iupac => 'silicic acid' },
    { smiles => '[SbH2](O)=O', iupac => 'stibinic acid' },
    { smiles => '[SbH2]O', iupac => 'stibinous acid' },
    { smiles => '[SbH](O)(O)=O', iupac => 'stibonic acid' },
    { smiles => '[SbH](O)O', iupac => 'stibonous acid' },
    { smiles => '[Sb](O)(O)(O)=O', iupac => 'stiboric acid' },
    { smiles => '[Sb](O)(O)O', iupac => 'stiborous acid' },
    { smiles => 'S(O)(O)(=O)=O', iupac => 'sulfuric acid', AUTHOR => 1 },
    { smiles => 'S(O)(O)=O', iupac => 'sulfurous acid', AUTHOR => 1 },
    { smiles => '[Te](O)(O)(=O)=O', iupac => 'telluric acid', AUTHOR => 1 },
    { smiles => '[Te](O)(O)=O', iupac => 'tellurous acid', AUTHOR => 1 },
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
