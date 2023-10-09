package ChemOnomatopist::Chain::Bicycle::Purine;

use strict;
use warnings;

# ABSTRACT: Purine chain
# VERSION

use parent ChemOnomatopist::Chain::Bicycle::;

sub new
{
    my( $class, $graph, $pyrimidine, $imidazole, @vertices ) = @_;
    return bless { graph => $graph,
                   cycles => [ $pyrimidine, $imidazole ],
                   vertices => \@vertices };
}

sub locants(@)
{
    my $self = shift;
    return map { $_ + 1 } @_;
}

1;
