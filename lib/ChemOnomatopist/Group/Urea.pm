package ChemOnomatopist::Group::Urea;

# ABSTRACT: Urea group
# VERSION

use strict;
use warnings;

use ChemOnomatopist::Elements qw( %elements );
use ChemOnomatopist::Group::Amide;
use ChemOnomatopist::Name;
use List::Util qw( any first );
use Scalar::Util qw( blessed );

use parent ChemOnomatopist::Group::, ChemOnomatopist::Chain::;

sub new
{
    my( $class, $graph, @vertices ) = @_;
    my $ketone = first { blessed $_ && $_->isa( ChemOnomatopist::Group::Ketone:: ) } @vertices;
    if( $graph->degree( $vertices[3] ) > $graph->degree( $vertices[2] ) ) {
        @vertices[2..3] = reverse @vertices[2..3];
    }
    return bless { graph => $graph,
                   vertices => \@vertices,
                   ketone_element => $ketone->{element} },
                 $class;
}

sub heteroatom_positions() { my @not_applicable }

sub needs_substituent_locants()
{
    my( $self ) = @_;
    return 1 if $self->number_of_branches > 1 && $self->number_of_branches < 4;
    return 1 if any { $_->has_substituent_locant } map { @$_ } $self->locant_names;
    return '';
}

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
