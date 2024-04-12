#!/usr/bin/perl

use strict;
use warnings;

use ChemOnomatopist;
use Test::More;

my @cases = (
    # From BBv3 P-25.1.2.7
    { smiles => 'C1=CC2=CC=CC3=CC=CC1=C23', iupac => 'acenaphthylene' },
    { smiles => 'C1=CC2=CC=CC3=CC4=CC=CC=C4C1=C23', iupac => 'aceanthrylene' },
    { smiles => 'C1=CC=C2C=CC3=CC4=CC=CC=C4C1=C23', iupac => 'acephenanthrylene' },
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
