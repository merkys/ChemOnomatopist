package ChemOnomatopist::Group::Guanidine;

use strict;
use warnings;

# ABSTRACT: Guanidine group
# VERSION

use parent ChemOnomatopist::Group::, ChemOnomatopist::Chain::Circular::;

use Chemistry::OpenSMILES qw(
    is_double_bond
    is_single_bond
);

sub new
{
    my( $class, $graph, $atom ) = @_;
    my $self = bless { graph => $graph, vertices => [ $graph->neighbours( $atom ) ] }, $class;
    $graph->delete_vertices( $atom, grep { ChemOnomatopist::is_element( $_, 'H' ) } $graph->vertices );
    return $self;
}

sub needs_heteroatom_locants() { return '' }
sub needs_heteroatom_names() { return '' }

sub prefix { return 'carbamimidoylamino' } # FIXME: Two kinds exist, BBv2 P-66.4.1.2.1.3
sub suffix { return 'guanidine' }

1;
