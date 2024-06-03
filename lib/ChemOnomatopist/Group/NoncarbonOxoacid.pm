package ChemOnomatopist::Group::NoncarbonOxoacid;

# ABSTRACT: Mononuclear noncarbon oxoacid as per BBv3 P-67.1.1.1
# VERSION

use strict;
use warnings;

use parent ChemOnomatopist::Group::;

use ChemOnomatopist;
use ChemOnomatopist::Group::Hydroxy;
use ChemOnomatopist::Group::Ketone;
use ChemOnomatopist::Group::XO3;
use Scalar::Util qw( blessed );

sub new
{
    my( $class, $atom, @attachments ) = @_;
    return bless { attachments => \@attachments, atom => $atom }, $class;
}

sub element() { ChemOnomatopist::element( $_[0]->{atom} ) }

sub is_part_of_chain() { 1 }

my %elements = (
    As => 'ars',
    N  => 'az',
    B  => 'bor',
    Br => 'brom',
    Cl => 'chlor',
    F  => 'fluor',
    I  => 'iod',
    P  => 'phosph',
    Se => 'sele',
    Si => 'silic',
    Sb => 'stib',
    S  => 'sulf',
    Te => 'tellur',
);

my %suffixes = (
    0 => { 1 => 'inous', 2 => 'onous', 3 => 'orous', 4 => 'ic' },
    1 => { 1 => 'inic',  2 => 'onic',  3 => 'oric'  },
    2 => { 1 => 'ic' },
);

sub suffix
{
    my( $self ) = @_;
    my $name = $elements{$self->element};
    if( blessed $self->{atom} && $self->{atom}->isa( ChemOnomatopist::Group::XO3:: ) ) {
        $name = 'per' . $name . 'ic';
    } else {
        my @attachments = @{$self->{attachments}};
        my $hydroxy = grep { blessed $_ && $_->isa( ChemOnomatopist::Group::Hydroxy:: ) } @attachments;
        my $ketones = grep { blessed $_ && $_->isa( ChemOnomatopist::Group::Ketone:: ) }  @attachments;

        $name .= $suffixes{$ketones}->{$hydroxy};
    }
    return $name . ' acid';
}

1;
