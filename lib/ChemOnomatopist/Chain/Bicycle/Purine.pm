package ChemOnomatopist::Chain::Bicycle::Purine;

# ABSTRACT: Purine chain
# VERSION

use strict;
use warnings;

use parent ChemOnomatopist::Chain::Bicycle::;

sub new
{
    my( $class, $graph, @cycles ) = @_;

    # Checking if the order of cycles is correct
    my( $pyrimidine, $imidazole ) = @cycles;
    if( $pyrimidine->length == 5 ) {
        ( $pyrimidine, $imidazole ) = map { $_->flipped } @cycles
    }

    my @vertices = ( @{$pyrimidine->{vertices}}[1..5],
                     $pyrimidine->{vertices}[0],
                     reverse @{$imidazole->{vertices}}[0..2] );
    return bless { graph => $graph,
                   cycles => [ $pyrimidine, $imidazole ],
                   vertices => \@vertices };
}

sub is_purine() { 1 }

sub locants(@)
{
    my $self = shift;
    return map { $_ + 1 } @_;
}

sub suffix() { ChemOnomatopist::Name->new( 'purine' ) }

1;
