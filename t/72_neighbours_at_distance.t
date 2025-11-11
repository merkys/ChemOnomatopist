#!/usr/bin/perl

use strict;
use warnings;

use ChemOnomatopist;
use ChemOnomatopist::Util::Graph qw( neighbours_at_distance );
use Chemistry::OpenSMILES::Parser;
use List::Util qw( first );
use Test::More;
use Set::Object qw( set );

my @cases = (
    { smiles => 'ClCC[C@H](CC(C)O)O', start => 3, distance => 0, neighbours => ';1,11,12;21;5,14,15' },
    { smiles => 'ClCC[C@H](CC(C)O)O', start => 3, distance => 1, neighbours => ';;0,9,10;6,7,16' },
    { smiles => 'ClCC[C@H](CC(C)O)O', start => 3, distance => 2, neighbours => ';;;17,18,19,20' },
    { smiles => 'ClCC[C@H](CC(C)O)O', start => 3, distance => 3, neighbours => ';;;' },
    { smiles => 'ClCC[C@H](CC(C)O)O', start => 3, distance => 4, neighbours => ';;;' },
    { smiles => 'ClCC[C@H](CC(C)O)O', start => 3, distance => 5, neighbours => ';;;' },

    { smiles => 'c1(ccccc1)Cc1ccccc1', start => 6, distance => 0, neighbours => ';;1,5;8,12' },
    { smiles => 'c1(ccccc1)Cc1ccccc1', start => 6, distance => 1, neighbours => ';;2,4,13,17;9,11,20,24' },
    { smiles => 'c1(ccccc1)Cc1ccccc1', start => 6, distance => 2, neighbours => ';;10,10,21,23;3,3,14,16' },
    { smiles => 'c1(ccccc1)Cc1ccccc1', start => 6, distance => 3, neighbours => ';;2,4,15,15;9,11,22,22' },
    { smiles => 'c1(ccccc1)Cc1ccccc1', start => 6, distance => 4, neighbours => ';;1,5,14,16;8,12,21,23' },
);

@cases = grep { !exists $_->{AUTHOR} } @cases unless $ENV{AUTHOR_TESTING};
plan skip_all => 'No available cases' unless @cases;
plan tests => scalar @cases;

my $parser = Chemistry::OpenSMILES::Parser->new;

for my $case (@cases) {
    my( $moiety ) = $parser->parse( $case->{smiles} );
    my $start = first { $_->{number} == $case->{start} } $moiety->vertices;
    my @neighbours;
    for my $atom ($moiety->neighbours( $start )) {
        my @atoms = neighbours_at_distance( $moiety, $atom, $case->{distance}, set( $start ) );
        push @neighbours, join ',', sort { $a <=> $b } map { $_->{number} } @atoms;
    }
    is join( ';', sort @neighbours ), $case->{neighbours};
}
