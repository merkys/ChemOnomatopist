package ChemOnomatopist::Group::Amine::SecondaryTertiary;

use strict;
use warnings;

# ABSTRACT: Secondary or tertiary amine
# VERSION

use List::Util qw( any );

use parent ChemOnomatopist::Group::, ChemOnomatopist::Chain::;

sub new
{
    my( $class, $graph, $atom ) = @_;
    return bless { graph => $graph, vertices => [ $atom ] }, $class;
}

sub element() { return 'N' }

sub is_part_of_chain() { return 1 }

sub locants(@)
{
    my $self = shift;
    return map { $_ ? $_ : 'N' } @_;
}

sub prefix(;$)
{
    my( $self, $parent ) = @_;
    return 'amino' unless $parent;

    my @neighbours = $self->graph->neighbours( $self->vertices );
    die "cannot process tertiary amines for now\n" if @neighbours == 3;
    die "cannot process complicated secondary amines for now\n" unless any { $_ == $parent } @neighbours;
    return ChemOnomatopist::Name->new( 'amino' ) if @neighbours == 1;

    my $name = ChemOnomatopist::get_sidechain_name( $self->graph, $self->vertices, grep { $_ != $parent } @neighbours );
    if( $name eq 'phenyl' ) {
        $name = ChemOnomatopist::Name->new( 'anilino' );
    } else {
        $name .= 'amino';
    }
    return $name;
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
