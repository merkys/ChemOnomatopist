#!/usr/bin/perl

use strict;
use warnings;

use ChemOnomatopist;
use Test::More;

plan tests => 3;

my @sorted;

@sorted = sort { ChemOnomatopist::cmp_attachments( $a, $b ) } ( [ 'butyl' ], [ 'tert-butyl' ] );
is join( ';', map { join ',', @$_ } @sorted ), 'butyl;tert-butyl';

@sorted = sort { ChemOnomatopist::cmp_attachments( $a, $b ) } ( [ 'tricosyl' ], [ 'tert-butyl' ] );
is join( ';', map { join ',', @$_ } @sorted ), 'tricosyl;tert-butyl';

@sorted = sort { ChemOnomatopist::cmp_attachments( $a, $b ) } ( [ 'ethyl' ], [ 'methyl' ] );
is join( ';', map { join ',', @$_ } @sorted ), 'ethyl;methyl';
