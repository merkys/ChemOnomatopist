#!/usr/bin/perl

use strict;
use warnings;

use ChemOnomatopist;
use Test::More;

open (my $f, '<', 't/PubChemData') or die;

my @iupac;
my @smiles;
while (my $line = <$f>) {
  my @elems = split ' ', $line;
  push(@iupac, $elems[1]);
  push(@smiles, $elems[2]);
}

close $f;

plan tests => scalar @iupac;

for my $i (0 .. $#iupac){
    is( ChemOnomatopist::get_name( $smiles[$i] ), $iupac[$i] );
}