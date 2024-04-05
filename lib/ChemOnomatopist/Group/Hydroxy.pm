package ChemOnomatopist::Group::Hydroxy;

# ABSTRACT: Hydroxy group
# VERSION

use strict;
use warnings;

use ChemOnomatopist::Name;

use parent ChemOnomatopist::Group::;

sub new
{
    my( $class, $atom ) = @_;
    return bless { atom => $atom }, $class;
}

sub element() { $_[0]->{atom}{symbol} }

# From BBv2 P-63.1.5
my %prefixes = ( O => 'hydroxy', S => 'sulfanyl', Se => 'selanyl', Te => 'tellanyl' );
my %suffixes = ( O => 'ol', S => 'thiol', Se => 'selenol', Te => 'tellurol' );

sub prefix
{
    my( $self ) = @_;
    return ChemOnomatopist::Name->new( $prefixes{$self->element} );
}

sub suffix
{
    my( $self ) = @_;

    my $suffix = '';
    # FIXME: Isotopes have to come inside the same parenthesis
    if( exists $self->{atom}{isotope} ) {
        $suffix = '(' . $self->{atom}{isotope} . $self->element . ')';
    }
    if( @{$self->{atom}{h_isotope}} && defined $self->{atom}{h_isotope}[0] ) {
        $suffix = '(' . $self->{atom}{h_isotope}[0] . 'H)';
    }

    return $suffix . $suffixes{$self->element};
}

sub _cmp_instances
{
    my( $A, $B ) = @_;
    return $A->element cmp $B->element
}

1;
