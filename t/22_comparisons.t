#!/usr/bin/perl

use strict;
use warnings;

use ChemOnomatopist;
use ChemOnomatopist::Isotope;
use Test::More;

plan tests => 11;

my @sorted;

@sorted = sort { ChemOnomatopist::cmp_attachments( $a, $b ) } ( [ 'butyl' ], [ 'tert-butyl' ] );
is join( ';', map { join ',', @$_ } @sorted ), 'butyl;tert-butyl';

@sorted = sort { ChemOnomatopist::cmp_attachments( $a, $b ) } ( [ 'tricosyl' ], [ 'tert-butyl' ] );
is join( ';', map { join ',', @$_ } @sorted ), 'tricosyl;tert-butyl';

@sorted = sort { ChemOnomatopist::cmp_attachments( $a, $b ) } ( [ 'ethyl' ], [ 'methyl' ] );
is join( ';', map { join ',', @$_ } @sorted ), 'ethyl;methyl';

@sorted = sort { ChemOnomatopist::cmp_heteroatom_seniority( $a, $b ) } ( [ 'N' ], [ 'S' ] );
is join( ';', map { join ',', @$_ } @sorted ), 'S;N';

@sorted = sort { ChemOnomatopist::cmp_heteroatom_seniority( $a, $b ) } ( [ 'O' ], [ 'Si' ] );
is join( ';', map { join ',', @$_ } @sorted ), 'O;Si';

@sorted = sort { ChemOnomatopist::cmp_heteroatom_seniority( $a, $b ) } ( [ 'N' ], [ 'O' ] );
is join( ';', map { join ',', @$_ } @sorted ), 'O;N';

@sorted = sort { ChemOnomatopist::cmp_heteroatom_seniority( $a, $b ) } ( [ 'N' ], [ 'O' ], [ 'P' ] );
is join( ';', map { join ',', @$_ } @sorted ), 'O;N;P';

my @isotopes;

# From BBv3 P-44.4.1.11.1
@isotopes = ( ChemOnomatopist::Isotope->new( 'H',  2, 1 ),
              ChemOnomatopist::Isotope->new( 'C', 14, 1 ) );
is ChemOnomatopist::Isotope::cmp_isotope_lists( [ $isotopes[0], $isotopes[0] ], [ $isotopes[1] ] ), -1;

# From BBv3 P-44.4.1.11.2
@isotopes = ( ChemOnomatopist::Isotope->new( 'C', 14, 1 ),
              ChemOnomatopist::Isotope->new( 'H',  2, 1 ) );
is ChemOnomatopist::Isotope::cmp_isotope_lists( [ $isotopes[0] ], [ $isotopes[1] ] ), -1;

# From BBv3 P-44.4.1.11.3
@isotopes = ( ChemOnomatopist::Isotope->new( 'C', 14, 1 ),
              ChemOnomatopist::Isotope->new( 'C', 13, 1 ) );
is ChemOnomatopist::Isotope::cmp_isotope_lists( [ $isotopes[0] ], [ $isotopes[1] ] ), -1;

# From BBv3 P-44.4.1.11.4
@isotopes = ( ChemOnomatopist::Isotope->new( 'H', 2, 2 ),
              ChemOnomatopist::Isotope->new( 'H', 2, 3 ) );
is ChemOnomatopist::Isotope::cmp_isotope_lists( [ $isotopes[0] ], [ $isotopes[1] ] ), -1;
