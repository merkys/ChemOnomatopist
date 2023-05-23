package ChemOnomatopist::Chain::Circular;

use strict;
use warnings;

# ABSTRACT: Circular chain
# VERSION

use ChemOnomatopist;
use ChemOnomatopist::ChainHalf; # FIXME: Not sure why it is needed
use ChemOnomatopist::Util::SMILES qw( cycle_SMILES );
use List::Util qw( all any );

use parent ChemOnomatopist::ChainHalf::;

# From BBv2 P-22.2.1
our %names = (
    OCNCC => '1,3-oxazolidine',
    ONCCC => '1,2-oxazolidine',
    NCCCC => 'pyrrolidine',
    NNCCC => 'pyrazolidine',
    NCNCC => 'imidazolidine',

    NCCCCC => 'piperidine',
    NCCNCC => 'piperazine',
    OCCNCC => 'morpholine',

    # 5-membered aromatic
    'OC=CC=C'  => 'furan',
    'NC=NC=C'  => '1H-imidazole', # FIXME: Adjust for isomerism
    'OC=NC=C'  => '1,3-oxazole',
    'ON=CC=C'  => '1,2-oxazole',
    'NN=CC=C'  => '1H-pyrazole', # FIXME: Adjust for isomerism
    'NC=CC=C'  => '1H-pyrole', # FIXME: Adjust for isomerism
    'SeC=CC=C' => 'selenophene',
    'TeC=CC=C' => 'tellurophene',
    'SC=CC=C'  => 'thiophene',

    # 6-membered aromatic
    'C=CC=CC=C' => 'benzene',
    'OCC=CC=C'  => '2H-pyran', # FIXME: Adjust for isomerism
    'N=CC=NC=C' => 'pyrazine',
    'N=NC=CC=C' => 'pyridazine',
    'N=CC=CC=C' => 'pyridine',
    'N=CN=CC=C' => 'pyrimidine',
);

sub new
{
    my( $class, $graph, @vertices ) = @_;
    return bless { graph => $graph, vertices => \@vertices }, $class;
}

# Need to override as the terminal bond is also important
sub backbone_SMILES()
{
    my( $self ) = @_;
    return cycle_SMILES( $self->{graph}, $self->vertices );
}

sub is_benzene()
{
    my( $self ) = @_;
    return '' unless $self->length == 6;
    return $self->backbone_SMILES =~ /^(C=CC=CC=C|CC=CC=CC=)$/;
}

sub name()
{
    my( $self ) = @_;

    my $graph = $self->{graph};
    my $SMILES = $self->backbone_SMILES;

    # Check the preserved names
    return $names{$SMILES} if exists $names{$SMILES};

    # FIXME: This is an exception, but a proper fix should be implemented instead
    return 'benzene' if $self->is_benzene;

    # Check for cycloalkanes
    if( all { $_->{symbol} eq 'C' } $self->vertices ) {
        return 'cyclo' . ChemOnomatopist::alkane_chain_name( $self->length ) . 'ane';
    }

    # Check for annulenes
    # FIXME: Check for kekulized compounds
    if( ( all { $_->{symbol} eq 'c' } $self->vertices ) &&
        $self->length =~ /^(4|6|8|10|12|14|16)$/ ) {
        # Annulene detected
        return 'cyclo' .
               ChemOnomatopist::IUPAC_numerical_multiplier( $self->length, 1 ) .
               ChemOnomatopist::IUPAC_numerical_multiplier( $self->length / 2, 1 ) . 'ene';
    }

    # No other types of graphs with cycles can be processed for now
    die "cannot handle complicated monocycles for now\n";
}

# sub branch_positions() # TODO: Maybe need to add 1 to all returned positions?

sub _disconnected_chain_graph()
{
    my( $self ) = @_;

    my $graph = $self->{graph}->copy;
    $graph->delete_cycle( $self->vertices );

    return $graph;
}

1;
