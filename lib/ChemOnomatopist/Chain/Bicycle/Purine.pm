package ChemOnomatopist::Chain::Bicycle::Purine;

# ABSTRACT: Purine chain
# VERSION

use strict;
use warnings;

use List::Util qw( first );

use parent ChemOnomatopist::Chain::Bicycle::;

sub new
{
    my( $class, $graph, @cycles ) = @_;

    my $imidazole  = first { $_->length == 5 } @cycles;
    my $pyrimidine = first { $_->length == 6 } @cycles;
    if( ChemOnomatopist::element( $pyrimidine->{vertices}[0] ) eq 'N' ) {
        $pyrimidine = $pyrimidine->flipped;
    }
    if( $pyrimidine->{vertices}[-1] != $imidazole->{vertices}[-1] ) {
        $imidazole = $imidazole->flipped;
    }

    my @vertices = ( @{$pyrimidine->{vertices}}[1..5],
                     $pyrimidine->{vertices}[0],
                     @{$imidazole->{vertices}}[0..2] );
    return bless { graph => $graph,
                   cycles => \@cycles,
                   vertices => \@vertices }, $class;
}

sub is_purine() { 1 }

sub locants(@)
{
    my $self = shift;
    return map { $_ + 1 } @_;
}

sub suffix() { ChemOnomatopist::Name->new( 'purine' ) }

1;
