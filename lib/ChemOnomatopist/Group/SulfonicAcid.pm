package ChemOnomatopist::Group::SulfonicAcid;

# ABSTRACT: Sulfonic acid group
# VERSION

use strict;
use warnings;

use parent ChemOnomatopist::Group::;

use ChemOnomatopist::Name;
use ChemOnomatopist::Util qw( array_frequencies );
use List::Util qw( all first );
use Scalar::Util qw( blessed );

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
    my $hydroxy = first { blessed $_ && ( $_->isa( ChemOnomatopist::Group::Hydroxy:: ) ||
                                          $_->isa( ChemOnomatopist::Group::Hydroperoxide:: ) ) }
                        @attachments;
    my @ketones = grep { $_ != $hydroxy } @attachments;
    my @elements = sort grep { $_ ne 'O' } map { ChemOnomatopist::element( $_ ) } @attachments;
    if( @ketones == 2 && !@elements ) {
        return ChemOnomatopist::Name->new( 'sulfonic acid' ) if $hydroxy->isa( ChemOnomatopist::Group::Hydroxy:: );
        return ChemOnomatopist::Name->new( 'sulfonoperoxoic acid' );
    }

    my %elements = array_frequencies( @elements ); use Data::Dumper; print Dumper \%elements;
    @elements = ();
    for my $element (sort keys %elements) {
        if( $elements{$element} > 1 ) {
            push @elements, ChemOnomatopist::IUPAC_numerical_multiplier( $elements{$element} );
        }
        push @elements, $suffixes{$element};
    }

    my @nonketone_elements = $hydroxy->isa( ChemOnomatopist::Group::Hydroxy:: )
                                ? ( $hydroxy->element )
                                : map { ChemOnomatopist::element( $_ ) } @{$hydroxy->{atoms}};

    local $" = '';
    my $name = 'sulfono';
    $name .= 'peroxo' if $hydroxy->isa( ChemOnomatopist::Group::Hydroperoxide:: );
    $name .= "@{elements}ic @nonketone_elements-acid";
    return ChemOnomatopist::Name->new( $name );
}

1;
