package ChemOnomatopist::Group::Sulfonamide;

# ABSTRACT: Sulfonamide group or its Se/Te equivalent
# VERSION

use strict;
use warnings;

use parent ChemOnomatopist::Group::;

use ChemOnomatopist::Name;
use List::Util qw( all uniq );

my %prefixes = (
    S  => 'sulfonamide',
    Se => 'selenonamide',
    Te => 'telluronamide',
);

my %ketone_infixes = (
    O  => '',
    S  => 'thio',
    Se => 'seleno',
    Te => 'telluro',
);

sub new
{
    my( $class, $element, @ketones ) = @_;
    return bless { element => $element, ketones => \@ketones }, $class;
}

sub prefix { ChemOnomatopist::Name->new( 'sulfamoyl' ) } # FIXME: May be incorrect

sub suffix
{
    my( $self ) = @_;
    if( all { $_->element eq 'O' } @{$self->{ketones}} ) {
        return ChemOnomatopist::Name->new( $prefixes{$self->element} );
    }

    my $name = $prefixes{$self->element};
    $name =~ s/amide$/o/;

    if( uniq( map { $_->element } @{$self->{ketones}} ) == 1 ) {
        $name .= 'di' . $ketone_infixes{$self->{ketones}[0]->element};
    } else {
        $name .= join '', sort map { $ketone_infixes{$_->element} } @{$self->{ketones}};
    }
    $name .= 'amide';

    return ChemOnomatopist::Name->new( $name );
}

1;
