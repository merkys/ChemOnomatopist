package ChemOnomatopist::Group::Urea;

# ABSTRACT: Urea group
# VERSION

use strict;
use warnings;

use ChemOnomatopist::Elements qw( %elements );
use ChemOnomatopist::Group::Amide;
use ChemOnomatopist::Name;
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

sub heteroatom_positions() { my @not_applicable }

sub needs_substituent_locants() { $_[0]->number_of_branches > 1 && $_[0]->number_of_branches < 4 }

sub locants(@)
{
    my $self = shift;
    return map { $_ ? 'N' . '\'' x ($_ - 2) : 1 } @_;
}

sub prefix() { ChemOnomatopist::Name->new( 'carbamoylamino' ) }

sub suffix()
{
    my( $self ) = @_;
    my $name = ChemOnomatopist::Name->new;
    if( $self->{ketone_element} ne 'O' ) {
        $name->append_element( $elements{$self->{ketone_element}}->{prefix} );
        $name->[-1]->{value} =~ s/a$/o/;
    }
    return $name . 'urea';
}

1;
