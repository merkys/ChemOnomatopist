package ChemOnomatopist::Chain::Benzamide;

use strict;
use warnings;

# ABSTRACT: Benzamide chain
# VERSION

use parent ChemOnomatopist::Chain::;

sub new
{
    my( $class, $graph, $amide, $C, $benzene ) = @_;
    return bless { graph => $graph,
                   benzene => $benzene,
                   vertices => [ $amide, $C, $benzene->vertices ] }, $class;
}

sub candidates
{
    my( $self ) = @_;
    my @candidates = ( $self );

    my @benzene_vertices = $self->{benzene}->vertices;
    push @benzene_vertices, shift @benzene_vertices;
    @benzene_vertices = reverse @benzene_vertices;

    push @candidates,
         ChemOnomatopist::Chain::Benzamide->new( $self->graph,
                                                 @{$self->{vertices}}[0..1],
                                                 ChemOnomatopist::Chain::Monocycle->new( $self->graph,
                                                                                         @benzene_vertices ) );
    $candidates[-1]->{candidate_for} = $self;
    return @candidates;
}

1;
