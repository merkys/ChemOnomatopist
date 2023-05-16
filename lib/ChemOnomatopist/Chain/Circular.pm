package ChemOnomatopist::Chain::Circular;

use strict;
use warnings;

# ABSTRACT: Circular chain
# VERSION

use ChemOnomatopist::ChainHalf; # FIXME: Not sure why it is needed

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

    my @vertices = $self->vertices;
    my $SMILES = $self->SUPER::backbone_SMILES;

    if( $self->{graph}->has_edge_attribute( $vertices[0], $vertices[$#vertices], 'bond' ) ) {
        $SMILES .= $self->{graph}->get_edge_attribute( $vertices[0], $vertices[$#vertices], 'bond' );
    }

    return $SMILES;
}

sub is_benzene()
{
    my( $self ) = @_;
    return '' unless $self->length == 6;
    return $self->backbone_SMILES =~ /^(C=CC=CC=C|CC=CC=CC=)$/;
}

# sub branch_positions() # TODO: Maybe need to add 1 to all returned positions?

sub _disconnected_chain_graph()
{
    my( $self ) = @_;

    my $graph = $self->{graph}->copy;
    $graph->delete_path( $self->vertices );

    return $graph;
}

1;
