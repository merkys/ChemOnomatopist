package ChemOnomatopist::Chain::Sulfimide;

# ABSTRACT: Sulfimide chain
# VERSION

use strict;
use warnings;

use ChemOnomatopist::Name;

use parent ChemOnomatopist::Group::, ChemOnomatopist::Chain::;

sub new
{
    my( $class, $graph, @vertices ) = @_;
    return bless { graph => $graph, vertices => \@vertices }, $class;
}

sub needs_heteroatom_locants() { '' }
sub needs_heteroatom_names() { '' }

sub locants($@)
{
    my $self = shift;
    map { $_ ? $self->{vertices}[1]->{symbol} : 'N' } @_;
}

# CHECKME: Se and Te suffixes may be incorrect, BBv3 has no examples
my %suffix = ( S => 'sulfan', Se => 'selan', Te => 'telan' );
sub suffix() { ChemOnomatopist::Name->new( $suffix{$_[0]->{vertices}[1]->{symbol}} . 'imine' ) }

1;
