package ChemOnomatopist::Group::SulfonicAcid;

# ABSTRACT: Sulfonic acid group
# VERSION

use strict;
use warnings;

use parent ChemOnomatopist::Group::;

use ChemOnomatopist::Name;
use List::Util qw( first );

sub new
{
    my( $class, @attachments ) = @_;
    return bless { attachments => \@attachments }, $class;
}

sub element() { 'S' }

# From BBv2 P-65.3.0 and Table 4.3
sub prefix() { ChemOnomatopist::Name->new( 'sulfo' ) }
sub suffix()
{
    my( $self ) = @_;

    my $hydroxy = first { $_->isa( ChemOnomatopist::Group::Hydroxy:: ) } @{$self->{attachments}};
    if( $hydroxy->element eq 'S' ) {
        return ChemOnomatopist::Name->new( 'sulfonothioic S-acid' );
    } else {
        return ChemOnomatopist::Name->new( 'sulfonic acid' );
    }
}

1;
