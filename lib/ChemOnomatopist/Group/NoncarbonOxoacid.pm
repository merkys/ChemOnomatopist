package ChemOnomatopist::Group::NoncarbonOxoacid;

# ABSTRACT: Mononuclear noncarbon oxoacid as per BBv3 P-67.1.1.1
# VERSION

use strict;
use warnings;

use parent ChemOnomatopist::Group::;

use ChemOnomatopist::Group::Hydroxy;
use ChemOnomatopist::Group::Ketone;
use ChemOnomatopist::Group::Nitro;
use ChemOnomatopist::Group::Sulfinyl;
use ChemOnomatopist::Group::Sulfonyl;
use ChemOnomatopist::Group::XO3;
use ChemOnomatopist::Util;
use Scalar::Util qw( blessed );

sub new
{
    my( $class, $atom, @attachments ) = @_;
    return bless { attachments => \@attachments, atom => $atom }, $class;
}

sub element() { ChemOnomatopist::Util::element( $_[0]->{atom} ) }

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
    Se => 'selen',
    Si => 'silic',
    Sb => 'stib',
    S  => 'sulfur',
    Te => 'tellur',
);

my %suffixes = (
    0 => { 1 => 'inous', 2 => 'onous', 3 => 'orous', 4 => 'ic' },
    1 => { 1 => 'inic',  2 => 'onic',  3 => 'oric'  },
    2 => { 1 => 'ic',    2 => 'ic' },
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
        $hydroxy += 1 if blessed $self->{atom} && $self->{atom}->isa( ChemOnomatopist::Group::Nitro:: );
        $ketones += 1 if blessed $self->{atom} && $self->{atom}->isa( ChemOnomatopist::Group::Sulfinyl:: );
        $ketones += 2 if blessed $self->{atom} && $self->{atom}->isa( ChemOnomatopist::Group::Sulfonyl:: );

        if(      $self->element eq 'B' ) {
            $name .= 'inic' if $hydroxy == 1;
            $name .= 'onic' if $hydroxy == 2;
            $name .= 'ic'   if $hydroxy == 3;
        } elsif( $self->element eq 'N' ) {
            $name .= 'inic'     if $hydroxy == 2 && !$ketones;
            $name  = 'nitric'   if $hydroxy == 2 &&  $ketones == 1;
            $name .= 'onic'     if $hydroxy == 3;
            $name  = 'nitroric' if $hydroxy == 4;
        } elsif( $self->element =~ /^(F|Cl|Br|I)$/ ) {
            $name .= 'ous' if $hydroxy == 2 && $ketones == 1;
            $name .= 'ic'  if $hydroxy == 1 && $ketones == 2;
            $name = 'hypo' . $name . 'ous' if $hydroxy == 1 && !$ketones;
        } elsif( $self->element =~ /^(S|Se|Te)$/ ) {
            $name .= 'ous' if $hydroxy == 2 && $ketones == 1;
            $name .= 'ic'  if $hydroxy == 2 && $ketones == 2;
        } else {
            $name .= $suffixes{$ketones}->{$hydroxy};
        }
    }
    return $name . ' acid';
}

1;
