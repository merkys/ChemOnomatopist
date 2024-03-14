package ChemOnomatopist::Group::Ester;

# ABSTRACT: Ester group
# VERSION

use strict;
use warnings;

use parent ChemOnomatopist::Group::, ChemOnomatopist::Chain::;

sub new
{
    my( $class, $graph, $C, $ketone, $O ) = @_;
    return bless { graph => $graph, vertices => [ $C, $O ] }, $class;
}

sub needs_heteroatom_locants() { '' }
sub needs_heteroatom_names() { '' }

sub needs_substituent_locants() { '' }

sub prefix() { die "cannot use esters as prefixes yet\n" }
sub suffix() { 'anoate' }

1;
