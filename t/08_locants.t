#!/usr/bin/perl

use strict;
use warnings;

use ChemOnomatopist;
use ChemOnomatopist::Chain::FromHalves;
use ChemOnomatopist::ChainHalf;
use ChemOnomatopist::MolecularGraph;
use Test::More;

sub chains($@)
{
    my( $graph, @halves ) = @_;
    @halves = map { ChemOnomatopist::ChainHalf->new( $graph, undef, @$_ ) } @halves;

    my @chains;
    for my $half1 (@halves) {
        for my $half2 (@halves) {
            next if $half1 == $half2;
            push @chains, ChemOnomatopist::Chain::FromHalves->new( $half1, $half2 );
        }
    }

    return @chains;
}

plan tests => 4;

my $graph;
my @paths;

my @atoms = map { { symbol => 'C', number => $_ } } 0..38;

$graph = ChemOnomatopist::MolecularGraph->new;
$graph->add_path( map { $atoms[$_] } ( 1, 0, 11..17 ) );
$graph->add_path( map { $atoms[$_] } ( 1, 22..27 ) );
$graph->add_path( map { $atoms[$_] } ( 1, 32..37 ) );
$graph->add_edge( map { $atoms[$_] } ( 11, 18 ) );
$graph->add_edge( map { $atoms[$_] } ( 25, 28 ) );
$graph->add_edge( map { $atoms[$_] } ( 32, 38 ) );

@paths = ChemOnomatopist::rule_lowest_numbered_locants( chains $graph,
                                                               [ map { $atoms[$_] } ( 0, 11..17 ) ],
                                                               [ map { $atoms[$_] } ( 0, 1, 22..27 ) ],
                                                               [ map { $atoms[$_] } ( 0, 1, 32..37 ) ] );
is scalar( @paths ), 1;
is join( ',', map { $_->{number} } $paths[0]->vertices ), '27,26,25,24,23,22,1,0,11,12,13,14,15,16,17';

my $chain = ChemOnomatopist::select_mainchain( $graph );
is join( ',', map { $_->{number} } $chain->vertices ), '27,26,25,24,23,22,1,0,11,12,13,14,15,16,17';

# Figure 7 from Urbonaitė, 2022.
# In the image, however, one branch is held as having priority over another, while in fact they are equal.
$graph = Graph::Undirected->new;
$graph->add_path( 'A'..'G' );
$graph->add_path( 'B', 'H' );
$graph->add_path( 'D', 'I'..'K' );
$graph->add_path( 'J', 'L' );
$graph->add_path( 'E', 'M' );

@paths = ChemOnomatopist::rule_lowest_numbered_locants( chains $graph,
                                                               [ reverse 'A'..'D' ],
                                                               [ 'D'..'G' ],
                                                               [ 'D', 'I'..'K' ] );
is scalar( @paths ), 2;
