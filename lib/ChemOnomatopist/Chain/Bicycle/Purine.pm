package ChemOnomatopist::Chain::Bicycle::Purine;

use strict;
use warnings;

# ABSTRACT: Purine chain
# VERSION

use parent ChemOnomatopist::Chain::Bicycle::;

sub new
{
    my( $class, $graph, $pyrimidine, $imidazole ) = @_;

    my @vertices = ( @{$pyrimidine->{vertices}}[1..5],
                     $pyrimidine->{vertices}[0],
                     reverse @{$imidazole->{vertices}}[0..2] );
    return bless { graph => $graph,
                   cycles => [ $pyrimidine, $imidazole ],
                   vertices => \@vertices };
}

sub locants(@)
{
    my $self = shift;
    return map { $_ + 1 } @_;
}

sub suffix() { return ChemOnomatopist::Name->new( 'purine' ) }

1;
