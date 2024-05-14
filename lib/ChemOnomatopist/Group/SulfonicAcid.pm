package ChemOnomatopist::Group::SulfonicAcid;

# ABSTRACT: Sulfonic acid group
# VERSION

use strict;
use warnings;

use parent ChemOnomatopist::Group::;

use ChemOnomatopist::Name;
use List::Util qw( all first );

sub new
{
    my( $class, @attachments ) = @_;
    return bless { attachments => \@attachments }, $class;
}

sub element() { 'S' }

my %suffixes = ( O => '', S => 'thio', Se => 'seleno', Te => 'telluro' );

# From BBv2 P-65.3.0 and Table 4.3
sub prefix() { ChemOnomatopist::Name->new( 'sulfo' ) }
sub suffix()
{
    my( $self ) = @_;

    my @attachments = @{$self->{attachments}};
    my $hydroxy = first {  $_->isa( ChemOnomatopist::Group::Hydroxy:: ) } @attachments;
    my @ketones = grep  { !$_->isa( ChemOnomatopist::Group::Hydroxy:: ) } @attachments;
    my @elements = sort grep { $_ ne 'O' } map { $_->element } @attachments;
    if( @ketones == 2 && !@elements ) {
        return ChemOnomatopist::Name->new( 'sulfonic acid' );
    }

    my $name = 'sulfono' . join( '', map { $suffixes{$_} } @elements ) . 'ic ' .
               $hydroxy->element . '-acid';
    return ChemOnomatopist::Name->new( $name );
}

1;
