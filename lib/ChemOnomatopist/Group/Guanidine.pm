package ChemOnomatopist::Group::Guanidine;

use strict;
use warnings;

# ABSTRACT: Guanidine group
# VERSION

use parent ChemOnomatopist::Group::;

use Chemistry::OpenSMILES qw(
    is_double_bond
    is_single_bond
);

sub new
{
    my( $class, $graph, $atom ) = @_;

    my @amino = grep { is_single_bond( $graph, $atom, $_ ) } $graph->neighbours( $atom );
    my @imino = grep { is_double_bond( $graph, $atom, $_ ) } $graph->neighbours( $atom );

    my @amino_neighbours = map { [ grep { $_ != $atom } $graph->neighbours( $_ ) ] } @amino;
    my @imino_neighbours = map {   grep { $_ != $atom } $graph->neighbours( $_ )   } @imino;

    my $self = bless { amino_neighbours => \@amino_neighbours,
                       imino_neighbours => \@imino_neighbours,
                       graph => $graph }, $class;

    for (@imino_neighbours, map { @$_ } @amino_neighbours) {
        $graph->add_edge( $self, $_ );
    }

    return $self;
}

sub C() { return $_[0] }

sub prefix { return 'carbamimidoylamino' } # FIXME: Two kinds exist, BBv2 P-66.4.1.2.1.3
sub suffix { return 'guanidine' }

1;
