#!/usr/bin/perl

use strict;
use warnings;

use Chemistry::OpenSMILES::Parser;
use ChemOnomatopist;
use ChemOnomatopist::Chain::Circular;
use Graph::Traversal::DFS;
use Test::More;

my @cases = sort keys %ChemOnomatopist::Chain::Circular::names;

plan tests => scalar @cases;

for my $SMILES (@cases) {
    my $ring = $SMILES;
    $ring =~ s/^(.)/${1}1/;
    $ring .= 1;

    my $parser = Chemistry::OpenSMILES::Parser->new;
    my( $graph ) = $parser->parse( $ring );
    $graph->delete_vertices( grep { $_->{symbol} eq 'H' } $graph->vertices );

    my $cycle = ChemOnomatopist::Chain::Circular->new( $graph, Graph::Traversal::DFS->new( $graph )->dfs );
    is $cycle->backbone_SMILES, $SMILES;
}
