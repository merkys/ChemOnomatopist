package ChemOnomatopist::Group::Hydrazide;

# ABSTRACT: Hydrazide group
# VERSION

use strict;
use warnings;

use ChemOnomatopist::Name;

use parent ChemOnomatopist::Group::, ChemOnomatopist::Chain::;

sub new
{
    my( $class, $graph, $ketone, @vertices ) = @_;
    return bless { graph => $graph, ketone => $ketone, vertices => \@vertices }, $class;
}

sub needs_heteroatom_locants() { '' }
sub needs_heteroatom_names() { '' }
sub needs_substituent_locants()
{
    my( $self ) = @_;
    return $self->number_of_branches > 1 && $self->number_of_branches < $self->max_valence;
}

sub locants(@)
{
    my $self = shift;
    return map { $_ == 0 ? "N'" : $_ == 1 ? 'N' : $_ - 1 } @_;
}

my %suffixes = ( O => '', S => 'thio', Se => 'seleno', Te => 'telluro' );

sub prefix() { ChemOnomatopist::Name->new( $suffixes{$_[0]->{ketone}->element} . 'hydrazidyl' ) }
sub suffix() { ChemOnomatopist::Name->new( $suffixes{$_[0]->{ketone}->element} . 'hydrazide'  ) }

1;
