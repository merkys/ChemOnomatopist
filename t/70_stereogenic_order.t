#!/usr/bin/perl

use strict;
use warnings;

use ChemOnomatopist;
use Chemistry::OpenSMILES qw( is_chiral );
use Chemistry::OpenSMILES::Parser;
use List::Util qw( first );
use Test::More;

my @cases = (
    { smiles => 'ClCC[C@H](CC(C)O)O', order => '8,2,4,13' },
);

plan tests => scalar @cases;

my $parser = Chemistry::OpenSMILES::Parser->new;

for my $case (@cases) {
    my( $moiety ) = $parser->parse( $case->{smiles} );
    my $atom = first { is_chiral $_ } $moiety->vertices;
    my @order = sort { ChemOnomatopist::order_by_neighbours( $moiety, $a, $b ) } $moiety->neighbours( $atom );
    is join( ',', map { $_->{number} } @order ), $case->{order};
}
