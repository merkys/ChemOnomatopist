#!/usr/bin/perl

use strict;
use warnings;

use ChemOnomatopist;
use Graph;

use Test::More tests => 1;

my $g = Graph->new( undirected => 1 );
$g->add_path( 'A', 'B', 'C', 'D', 'E', 'F', 'G', 'H' );
$g->add_path( 'D', 'DB', 'DC', );
$g->add_path( 'F', 'FB' );
$g->add_path( 'F', 'FC' );

is( ChemOnomatopist::get_name( $g ), '3-(methyl)-3-(methyl)-5-(ethyl)-octane' );
