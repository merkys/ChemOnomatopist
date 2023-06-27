package ChemOnomatopist::Group::Amine::SecondaryTertiary;

use strict;
use warnings;

# ABSTRACT: Secondary or tertiary amine
# VERSION

use List::Util qw( any );

use parent ChemOnomatopist::Group::, ChemOnomatopist::Chain::;

sub new
{
    my( $class, $graph, @vertices ) = @_;
    return bless { graph => $graph, vertices => \@vertices }, $class;
}

sub C()
{
    my( $self ) = @_;
    return $self;
}

sub candidates()
{
    my( $self ) = @_;

    my @chains;
    for ($self->graph->neighbours( $self )) {
        my $chain = ChemOnomatopist::select_sidechain( $self->graph, $self, $_ );
        $chain->{candidate_for} = $self;
        push @chains, ChemOnomatopist::Group::Amine::SecondaryTertiary->new( $self->graph,
                                                                             $self,
                                                                             $chain->vertices );
    }

    return @chains;
}

sub is_nitrogen() { return 1 }

sub locants(@)
{
    my $self = shift;
    return map { $_ ? $_ : 'N' } @_;
}

sub prefix(;$)
{
    my( $self, $parent ) = @_;
    return 'amino' unless $parent;

    my @neighbours = $self->graph->neighbours( $self );
    die "cannot process tertiary amines for now\n" if @neighbours == 3;
    die "cannot process complicated secondary amines for now\n" unless any { $_ == $parent } @neighbours;

    return ChemOnomatopist::get_sidechain_name( $self->graph, $self, grep { $_ != $parent } @neighbours ) . 'amino';
}

sub suffix()
{
    my( $self ) = @_;
    my $remaining_chain = ChemOnomatopist::Chain->new( $self->graph, $self->vertices );
    my $name = ChemOnomatopist::unbranched_chain_name( $remaining_chain );
    $name =~ s/e$//;
    $name .= '-1-' if $self->length > 3;
    return $name . 'amine';
}

1;
