#!/usr/bin/perl

use strict;
use warnings;

use ChemOnomatopist;
use Test::More;

# Collected by Miglė
my @prefixes = qw(
    ?
    meth
    eth
    prop
    but
    pent
    hex
    hept
    oct
    non
    dec
    undec
    dodec
    tridec
    tetradec
    pentadec
    hexadec
    heptadec
    octadec
    nonadec
    icos
    henicos
    docos
    tricos
    tetracos
    pentacos
    hexacos
    heptacos
    octacos
    nonacos
    triacont
    hentriacont
    dotriacont
    tritriacont
    tetratriacont
    pentatriacont
    hexatriacont
    heptatriacont
    octatriacont
    nonatriacont
    tetracont
    hentetracont
    dotetracont
    tritetracont
    tetratetracont
    pentatetracont
    hexatetracont
    heptatetracont
    octatetracont
    nonatetracont
    pentacont
    henpentacont
    dopentacont
    tripentacont
    tetrapentacont
    pentapentacont
    hexapentacont
    heptapentacont
    octapentacont
    nonapentacont
    hexacont
    henhexacont
    dohexacont
    trihexacont
    tetrahexacont
    pentahexacont
    hexahexacont
    heptahexacont
    octahexacont
    nonahexacont
    heptacont
    henheptacont
    doheptacont
    triheptacont
    tetraheptacont
    pentaheptacont
    hexaheptacont
    heptaheptacont
    octaheptacont
    nonaheptacont
    octacont
    henoctacont
    dooctacont
    trioctacont
    tetraoctacont
    pentaoctacont
    hexaoctacont
    heptaoctacont
    octaoctacont
    nonaoctacont
    nonacont
    hennonacont
    dononacont
    trinonacont
    tetranonacont
    pentanonacont
    hexanonacont
    heptanonacont
    octanonacont
    nonanonacont
    hect
    henhect
    dohect
    trihect
    tetrahect
    pentahect
    hexahect
    heptahect
    octahect
    nonahect
    decahect
    undecahect
    dodecahect
    tridecahect
    tetradecahect
    pentadecahect
    hexadecahect
    heptadecahect
    octadecahect
    nonadecahect
    icosahect
    henicosahect
    docosahect
    tricosahect
    tetracosahect
    pentacosahect
    hexacosahect
    heptacosahect
    octacosahect
    nonacosahect
    triacontahect
    hentriacontahect
    dotriacontahect
    tritriacontahect
    tetratriacontahect
    pentatriacontahect
    hexatriacontahect
    heptatriacontahect
    octatriacontahect
    nonatriacontahect
    tetracontahect
    hentetracontahect
    dotetracontahect
    tritetracontahect
    tetratetracontahect
    pentatetracontahect
    hexatetracontahect
    heptatetracontahect
    octatetracontahect
    nonatetracontahect
    pentacontahect
    henpentacontahect
    dopentacontahect
    tripentacontahect
    tetrapentacontahect
    pentapentacontahect
    hexapentacontahect
    heptapentacontahect
    octapentacontahect
    nonapentacontahect
    hexacontahect
    henhexacontahect
    dohexacontahect
    trihexacontahect
    tetrahexacontahect
    pentahexacontahect
    hexahexacontahect
    heptahexacontahect
    octahexacontahect
    nonahexacontahect
    heptacontahect
    henheptacontahect
    doheptacontahect
    triheptacontahect
    tetraheptacontahect
    pentaheptacontahect
    hexaheptacontahect
    heptaheptacontahect
    octaheptacontahect
    nonaheptacontahect
    octacontahect
    henoctacontahect
    dooctacontahect
    trioctacontahect
    tetraoctacontahect
    pentaoctacontahect
    hexaoctacontahect
    heptaoctacontahect
    octaoctacontahect
    nonaoctacontahect
    nonacontahect
    hennonacontahect
    dononacontahect
    trinonacontahect
    tetranonacontahect
    pentanonacontahect
    hexanonacontahect
    heptanonacontahect
    octanonacontahect
    nonanonacontahect
);

my %cases = (
    # Taken from:
    # https://en.wikipedia.org/w/index.php?title=IUPAC_numerical_multiplier&oldid=1086173027
     241 => 'hentetracontadict',
     411 => 'undecatetract',
     548 => 'octatetracontapentact',
    9267 => 'heptahexacontadictanonali',

    # Taken from BBv2 P-14.2.1.2
     363 => 'trihexacontatrict',
     486 => 'hexaoctacontatetract',
);

plan tests => scalar( @prefixes ) + scalar keys %cases;

for (0..$#prefixes) {
    is( ChemOnomatopist::alkane_chain_name( $_ ),
        $prefixes[$_],
        "Number $_" );
}

for (sort keys %cases) {
    is( ChemOnomatopist::alkane_chain_name( $_ ),
        $cases{$_},
        "Number $_" );
}
