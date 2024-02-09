package ChemOnomatopist::Group::Urea;

# ABSTRACT: Urea group
# VERSION

use strict;
use warnings;

use ChemOnomatopist::Elements qw( %elements );
use ChemOnomatopist::Group::Amide;
use Scalar::Util qw( blessed );

use parent ChemOnomatopist::Group::, ChemOnomatopist::Chain::;

sub new
{
    my( $class, $graph, @vertices ) = @_;
    my( $ketone_element ) = map  { $_->{ketone}->element }
                            grep { blessed $_ && $_->isa( ChemOnomatopist::Group::Amide:: ) }
                                 @vertices;
    return bless { graph => $graph,
                   vertices => \@vertices,
                   ketone_element => $ketone_element },
                 $class;
}

sub needs_substituent_locants() { $_[0]->number_of_branches > 1 && $_[0]->number_of_branches < 4 }

sub locants(@)
{
    my $self = shift;
    return map { $_ ? 'N' . '\'' x ($_ - 1) : 1 } @_;
}

sub prefix() { 'carbamoylamino' }

sub suffix()
{
    my( $self ) = @_;
    my $name = '';
    if( $self->{ketone_element} ne 'O' ) {
        $name = $elements{$self->{ketone_element}}->{prefix};
        $name =~ s/a$/o/;
    }
    return $name . 'urea';
}

1;
