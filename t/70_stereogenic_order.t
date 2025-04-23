#!/usr/bin/perl

use strict;
use warnings;

use ChemOnomatopist;
use Chemistry::OpenSMILES qw( is_chiral );
use Chemistry::OpenSMILES::Parser;
use List::Util qw( first );
use Test::More;

my @cases = (
    { smiles => 'ClCC[C@H](CC(C)O)O',    order => '8,2,4,13' }, # From BBv3 P-92.1.4.1

    # From BBv3 P-92.1.4.2
    { smiles => 'C[C@H](C=C)O',          order => '4,2,0,8'  },
    { smiles => 'O[C@H](C=O)C',          order => '0,2,4,6'  },

    # From BBv3 P-92.1.4.3
    { smiles => 'C1(CCC[C@]12C=CCC2)=O', order => '0,5,8,3', AUTHOR => 'not implemented yet' },
    { smiles => 'C1(CC1)[C@@H](C(C)C)O', order => '7,0,4,13' },
);

@cases = grep { !exists $_->{AUTHOR} } @cases unless $ENV{AUTHOR_TESTING};
plan skip_all => 'No available cases' unless @cases;
plan tests => scalar @cases;

my $parser = Chemistry::OpenSMILES::Parser->new;

for my $case (@cases) {
    my( $moiety ) = $parser->parse( $case->{smiles} );
    my $atom = first { is_chiral $_ } $moiety->vertices;
    my @order = sort { ChemOnomatopist::order_by_neighbours( $moiety, $atom, $a, $b ) } $moiety->neighbours( $atom );
    my $ok = is join( ',', map { $_->{number} } @order ), $case->{order};
    if( $case->{AUTHOR} && $ok ) {
        diag 'test supposed to fail with AUTHOR_TESTING' .
             ( $case->{AUTHOR} !~ /^1$/ ? ': ' . $case->{AUTHOR} : '' );
    }
}
