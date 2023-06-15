package ChemOnomatopist::Group::Guanidine;

use strict;
use warnings;

# ABSTRACT: Guanidine group
# VERSION

use parent ChemOnomatopist::Group::, ChemOnomatopist::Chain::Circular::;

use Chemistry::OpenSMILES qw( is_double_bond );

sub new
{
    my( $class, $graph, $atom ) = @_;
    # Double ! is used here as is_double_bond() returns 1 or undef which is a bug in Chemistry::OpenSMILES?
    my @vertices = sort { !!is_double_bond( $graph, $atom, $a ) <=>
                          !!is_double_bond( $graph, $atom, $b ) }
                        $graph->neighbours( $atom );
    my @orders = map { !!is_double_bond( $graph, $atom, $_ ) } @vertices;
    my $self = bless { graph => $graph, vertices => \@vertices, is_double_bond => \@orders }, $class;
    $graph->delete_vertices( $atom, grep { ChemOnomatopist::is_element( $_, 'H' ) } $graph->vertices );
    return $self;
}

sub needs_heteroatom_locants() { return '' }
sub needs_heteroatom_names() { return '' }
sub needs_substituent_locants() { return 1 } # FIXME: There may be identical substituents, what to do then?

sub locants(@) {
    my $self = shift;
    return map { 'N' . "'" x $_ } @_;
}

sub prefix { return 'carbamimidoylamino' } # FIXME: Two kinds exist, BBv2 P-66.4.1.2.1.3
sub suffix { return 'guanidine' }

1;
