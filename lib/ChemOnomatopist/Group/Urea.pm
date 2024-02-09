package ChemOnomatopist::Group::Urea;

# ABSTRACT: Urea group
# VERSION

use strict;
use warnings;

use ChemOnomatopist::Elements qw( %elements );
use ChemOnomatopist::Group::Amide;
use List::Util qw( first );
use Scalar::Util qw( blessed );

use parent ChemOnomatopist::Group::, ChemOnomatopist::Chain::;

sub new
{
    my( $class, $graph, @vertices ) = @_;
    my $ketone = first { blessed $_ && $_->isa( ChemOnomatopist::Group::Ketone:: ) } @vertices;
    return bless { graph => $graph,
                   vertices => \@vertices,
                   ketone_element => $ketone->{element} },
                 $class;
}

sub heteroatom_positions() { () }

sub needs_substituent_locants() { $_[0]->number_of_branches > 1 && $_[0]->number_of_branches < 4 }

sub locants(@)
{
    my $self = shift;
    return map { $_ ? 'N' . '\'' x ($_ - 2) : 1 } @_;
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
