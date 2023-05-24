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
    CCNCO => '1,3-oxazolidine',
    CCCNO => '1,2-oxazolidine',
    CCCCN => 'pyrrolidine',
    CCCNN => 'pyrazolidine',
    CCNCN => 'imidazolidine',

    CCCCCN => 'piperidine',
    CCNCCN => 'piperazine',
    CCNCCO => 'morpholine',

    # 5-membered aromatic
    'C=CC=CO'  => 'furan',
    'NC=NC=C'  => '1H-imidazole', # FIXME: Adjust for isomerism
    'OC=NC=C'  => '1,3-oxazole',
    'C=CC=NO'  => '1,2-oxazole',
    'C=CC=NN'  => '1H-pyrazole', # FIXME: Adjust for isomerism
    'C=CC=CN'  => '1H-pyrole', # FIXME: Adjust for isomerism
    'SeC=CC=C' => 'selenophene',
    'TeC=CC=C' => 'tellurophene',
    'C=CC=CS'  => 'thiophene',

    # 6-membered aromatic
    'C=CC=CC=C' => 'benzene',
    'C=CC=CCO'  => '2H-pyran', # FIXME: Adjust for isomerism
    'N=CC=NC=C' => 'pyrazine',
    'C=CC=CN=N' => 'pyridazine',
    'C=CC=CC=N' => 'pyridine',
    'N=CN=CC=C' => 'pyrimidine',
);

sub new
{
    my( $class, $graph, @vertices ) = @_;
    return bless { graph => $graph, vertices => \@vertices }, $class;
}

# Selecting the candidate with the lowest alphabetical order
sub backbone_SMILES()
{
    my( $self ) = @_;

    my @vertices = $self->vertices;
    my @candidates;
    for (0..$#vertices) {
        push @candidates, cycle_SMILES( $self->{graph}, @vertices );
        push @vertices, shift @vertices;
    }
    @vertices = reverse @vertices;
    for (0..$#vertices) {
        push @candidates, cycle_SMILES( $self->{graph}, @vertices );
        push @vertices, shift @vertices;
    }

    my( $SMILES ) = sort @candidates;
    return $SMILES;
}

sub name()
{
    my( $self ) = @_;

    my $graph = $self->{graph};
    my $SMILES = $self->backbone_SMILES;

    # Check the preserved names
    return $names{$SMILES} if exists $names{$SMILES};

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
