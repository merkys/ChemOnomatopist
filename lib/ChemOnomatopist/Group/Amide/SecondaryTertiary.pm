package ChemOnomatopist::Group::Amide::SecondaryTertiary;

use strict;
use warnings;

# ABSTRACT: Secondary or tertiary amide
# VERSION

use ChemOnomatopist::Util qw( copy );

use parent ChemOnomatopist::Group::, ChemOnomatopist::Chain::;

sub new
{
    my( $class, $graph, @vertices ) = @_;
    return bless { graph => $graph, vertices => \@vertices }, $class;
}

sub C() { return @_ }

sub candidates()
{
    my( $self ) = @_;

    my @chains;
    for ($self->graph->neighbours( $self )) {
        my $graph = copy $self->graph;
        $graph->delete_edge( $self, $_ );
        my $chain = ChemOnomatopist::select_sidechain( $graph, $self, $_ );
        $chain->{candidate_for} = $self;
        push @chains, ChemOnomatopist::Group::Amide::SecondaryTertiary->new( $self->graph,
                                                                             $self,
                                                                             $chain->vertices );
    }

    return @chains;
}

sub locants(@)
{
    my $self = shift;
    # Locant '1' is the carbon with ketone attachment, it is probably never used
    return map { $_ ? $_ + 1 : 'N' } @_;
}

sub suffix()
{
    my( $self ) = @_;
    my $name = ChemOnomatopist::unbranched_chain_name( $self ); # FIXME: This is incorrect, it does not pick all chain attachments
    $name =~ s/e$//;
    return $name . 'amide';
}

1;
