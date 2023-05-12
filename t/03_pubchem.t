#!/usr/bin/perl

use strict;
use warnings;

use Chemistry::OpenSMILES::Parser;
use ChemOnomatopist;
use List::Util qw( any );
use Test::More;

if( !$ENV{EXTENDED_TESTING} ) {
    plan skip_all => "Skip \$ENV{EXTENDED_TESTING} is not set\n";
}

open( my $inp, '<', 't/PubChemData' ) or die;

my @cases;
while (<$inp>) {
    my( $id, $iupac, $smiles ) = split /\t/, $_;

    # TODO: The following compounds are not properly named yet
    next if $iupac =~ /edial$/;
    next if $iupac =~ / /;
    next if $iupac =~ /acetyl/;

    next if $smiles =~ /[\[\]=\$\#]/; # TODO: Cannot process these

    my $parser = Chemistry::OpenSMILES::Parser->new;
    my( @graphs ) = $parser->parse( $smiles );

    next if @graphs > 1;
    my $graph = shift @graphs;

    # TODO: Cannot yet process double, triple, ... bonds between carbon atoms
    #~ next if any { $_->[0]{symbol} eq 'C' && $_->[1]{symbol} eq 'C' }
            #~ grep { $graph->has_edge_attributes( @$_ ) }
            #~ $graph->edges;

    #~ next if any { join( '', sort map { $_->{symbol} } $graph->neighbours( $_ ) ) =~ /^(CC|OO)$/ }
            #~ grep { $graph->degree( $_ ) == 2 }
            #~ grep { $_->{symbol} eq 'O' }
            #~ $graph->vertices;

    push @cases, { id => $id, iupac => $iupac, graph => $graph };
}
close $inp;

plan tests => scalar @cases;

for my $case (@cases) {
    # Minor regularizations for PubChem names:
    $case->{iupac} =~ s/(di|tri|tetra|penta|hepta)(tert-butyl)/$1\($2\)/g;
    $case->{iupac} =~ s/di\((non|heptadec|hentriacont|undec|tridec|docos|icos|tetradec|pentadec)yl\)/di$1yl/g;
    $case->{iupac} =~ s/tetra\(tridecyl\)/tetratridecyl/g;

    # The following multiplicative prefixes were simplified by BBv2
    $case->{iupac} =~ s/-tris-/-tri/g;
    $case->{iupac} =~ s/-(tetra|hexa)kis-/-$1/g;

    my $name;
    eval { $name = ChemOnomatopist::get_name( $case->{graph} ); };
    if( $@ ) {
        $@ =~ s/\n$//;
        fail 'ID ' . $case->{id} . " failed: $@";
    } else {
        is $name, $case->{iupac}, 'ID ' . $case->{id};
    }
}
