package ChemOnomatopist::Group::Monocycle::Fused;

use strict;
use warnings;

# ABSTRACT: Monocyclic group which is fused to some other cycle
# VERSION

use parent ChemOnomatopist::Group::Monocycle::;

use ChemOnomatopist;
use List::Util qw( all );

sub new
{
    my( $class, $graph, $system, @vertices ) = @_;

    return bless { graph => $graph, system => \$system, vertices => \@vertices }, $class;
}

sub is_aromatic()
{
    my( $self ) = @_;
    return 1 if $self->SUPER::is_aromatic;
    return $self->{system}->is_aromatic;
}

1;
