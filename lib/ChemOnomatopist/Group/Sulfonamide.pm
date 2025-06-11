package ChemOnomatopist::Group::Sulfonamide;

# ABSTRACT: Sulfonamide group or its Se/Te equivalent
# VERSION

use strict;
use warnings;

use parent ChemOnomatopist::Group::;

use ChemOnomatopist::Name;

my %prefixes = (
    S  => 'sulfonamide',
    Se => 'selenonamide',
    Te => 'telluronamide',
);

sub new
{
    my( $class, $element, @ketones ) = @_;
    return bless { element => $element, ketones => \@ketones }, $class;
}

sub prefix { ChemOnomatopist::Name->new( 'sulfamoyl' ) } # FIXME: May be incorrect
sub suffix { ChemOnomatopist::Name->new( $prefixes{$_[0]->element} ) }

1;
