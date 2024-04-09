#!/usr/bin/perl

use strict;
use warnings;

use ChemOnomatopist;
use Test::More;

my @cases = (
    # From BBv3 P-45.5
    { smiles => 'BrC1=C(C=C(C2=CC=CC=C12)Cl)CCOCC1=C(C2=CC=CC=C2C(=C1)Br)Br', iupac => '1-bromo-4-chloro-2-{2-[(1,4-dibromonaphthalen-2-yl)methoxy]ethyl}naphthalene', AUTHOR => 1 },
    { smiles => 'BrC1=C(NC2=C(C=C(C=C2)Br)Br)C=CC(=C1)Cl', iupac => '2-bromo-4-chloro-N-(2,4-dibromophenyl)aniline', AUTHOR => 1 },
    # { smiles => '', iupac => '13-bromo-14-chloro-2-(3,4-dibromonaphthalen-2-yl)-1(2)-naphthalena-3,5(1,4),7(1)-tribenzenaheptaphane' }, # FIXME: Not parsed by OPSIN
    { smiles => 'FC(C(C)F)C(CCC(=O)O)C(C(C)[N+](=O)[O-])[N+](=O)[O-]', iupac => '4-(1,2-difluoropropyl)-5,6-dinitroheptanoic acid' },
    { smiles => '[81Br]C(C(C(CC(=O)O)C(C)C(C)[81Br])[N+](=O)[O-])C', iupac => '5-(81Br)bromo-3-[3-(81Br)bromobutan-2-yl]-4-nitrohexanoic acid' },
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
