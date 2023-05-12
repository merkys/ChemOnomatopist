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
while (my $line = <$inp>) {
    my @fields = split /\t/, $line;
    push @cases, { id => $fields[0], iupac => $fields[1], smiles => $fields[2] };
}
close $inp;

plan tests => scalar @cases;

for my $case (@cases) {
    next if $case->{iupac}  =~ /edial$/;
    next if $case->{iupac}  =~ / /;
    next if $case->{iupac}  =~ /acetyl/;

    # Minor regularizations for PubChem names:
    $case->{iupac} =~ s/(di|tri|tetra|penta|hepta)(tert-butyl)/$1\($2\)/g;
    $case->{iupac} =~ s/di\((non|heptadec|hentriacont|undec|tridec|docos|icos|tetradec|pentadec)yl\)/di$1yl/g;
    $case->{iupac} =~ s/tetra\(tridecyl\)/tetratridecyl/g;

    # The following multiplicative prefixes were simplified by BBv2
    $case->{iupac} =~ s/-tris-/-tri/g;
    $case->{iupac} =~ s/-(tetra|hexa)kis-/-$1/g;

    is( ChemOnomatopist::get_name( $case->{smiles} ),
        $case->{iupac},
        'ID ' . $case->{id} );
}

if( 0 ) {
    my $parser = Chemistry::OpenSMILES::Parser->new;
    my( $graph ) = $parser->parse( $case->{smiles} ); # Taking only the first graph

    # Cannot yet process double, triple, ... bonds between carbon atoms
    next if any { $_->[0]{symbol} eq 'C' && $_->[1]{symbol} eq 'C' }
            grep { $graph->has_edge_attributes( @$_ ) }
            $graph->edges;

    next if any { join( '', sort map { $_->{symbol} } $graph->neighbours( $_ ) ) =~ /^(CC|OO)$/ }
            grep { $graph->degree( $_ ) == 2 }
            grep { $_->{symbol} eq 'O' }
            $graph->vertices;

    my $name;
    eval {
        $name = ChemOnomatopist::get_name( $graph );
    };
    if( $@ ) {
        $@ =~ s/\n$//;
        ok 1, 'ID ' . $case->{id} . " skipped: $@";
    } else {
        is $name, $case->{iupac}, 'ID ' . $case->{id};
    }
}
