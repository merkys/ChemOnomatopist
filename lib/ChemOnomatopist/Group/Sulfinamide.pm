package ChemOnomatopist::Group::Sulfinamide;

# ABSTRACT: Sulfinamide group or its Se/Te equivalent
# VERSION

use strict;
use warnings;

use parent ChemOnomatopist::Group::;

use ChemOnomatopist::Name;

my %prefixes = (
    S  => 'sulfinamide',
    Se => 'seleninamide',
    Te => 'tellurinamide',
);

sub new
{
    my( $class, $element, $ketone ) = @_;
    return bless { element => $element, ketone => $ketone }, $class;
}

sub prefix { ChemOnomatopist::Name->new( 'sulfinamido' ) } # FIXME: May be incorrect
sub suffix
{
    my( $self ) = @_;
    if( $self->{ketone}->element eq 'O' ) {
        return ChemOnomatopist::Name->new( $prefixes{$self->element} );
    }

    my $name = $prefixes{$self->element};
    $name =~ s/amide$/o/;
    $name .= $prefixes{$self->{ketone}->element};
    $name =~ s/inamide$/oamide/;
    return ChemOnomatopist::Name->new( $name );
}

1;
