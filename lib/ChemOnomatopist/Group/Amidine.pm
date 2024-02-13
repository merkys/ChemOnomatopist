package ChemOnomatopist::Group::Amidine;

# ABSTRACT: Amidine group
# VERSION

use strict;
use warnings;

use parent ChemOnomatopist::Group::, ChemOnomatopist::Chain::;

use ChemOnomatopist::Name::Part::Locants;
use Scalar::Util qw( blessed );

sub new
{
    my( $class, $graph, @vertices ) = @_;
    return bless { graph => $graph, vertices => \@vertices }, $class;
}

sub nonstandard_valence_positions()
{
    my( $self ) = @_;

    return @{$self->{nonstandard_valence_positions}} if $self->{nonstandard_valence_positions};

    my @vertices = $self->vertices;
    my @nonstandard_valence_positions;
    for (1..$#vertices) { # Nonstandard valence of the central atom is not important, hence skipped
        next if blessed $vertices[$_];
        next if ChemOnomatopist::is_element( $vertices[$_], 'C' );
        next unless exists $vertices[$_]->{valence};
        push @nonstandard_valence_positions, $_;
    }

    $self->{nonstandard_valence_positions} = \@nonstandard_valence_positions;
    return @nonstandard_valence_positions;
}

sub needs_heteroatom_locants() { '' }
sub needs_heteroatom_names() { '' }
sub needs_substituent_locants { '' }

my %prefixes = ( C => '', S => 'sulf', Se => 'selen', Te => 'tellur' ); # CHECKME: Is Te correct?

sub prefix()
{
    my( $self ) = @_;
    my( $central_atom, @others ) = $self->vertices;
    return 'carbamimidoyl' unless $central_atom->{symbol} eq 'S';

    my $N = grep { ChemOnomatopist::element( $_ ) eq 'N' } @others;
    my $O = grep { ChemOnomatopist::element( $_ ) eq 'O' } @others;

    my $name = ChemOnomatopist::Name::Part::Locants->new( 'S-' )->to_name;
    $name .= 'amino';

    $name .= 'sulfon'    if $N == 2 &&  $O == 1;
    $name .= 'sulfonodi' if $N == 3 && !$O;
    $name .= 'sulfin';

    $name .= 'imidoyl';
    return $name;
}

sub suffix()
{
    my( $self ) = @_;
    my( $central_atom, @others ) = $self->vertices;

    my $name = $prefixes{$central_atom->{symbol}};
    if( $central_atom->{symbol} ne 'C' && @others == 3 ) {
        my $N = grep { ChemOnomatopist::element( $_ ) eq 'N' } @others;
        my $O = grep { ChemOnomatopist::element( $_ ) eq 'O' } @others;
        $name .= 'on'    if $N == 2 && $O == 1;
        $name .= 'onodi' if $N == 3;
    }

    return $name . 'imidamide';
}

1;
