package ChemOnomatopist::Chain::Circular;

use strict;
use warnings;

# ABSTRACT: Circular chain
# VERSION

use ChemOnomatopist;
use ChemOnomatopist::ChainHalf; # FIXME: Not sure why it is needed
use ChemOnomatopist::Elements qw( %elements );
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
    'C=CC=NC'  => '1H-imidazole', # FIXME: Adjust for isomerism
    'C=CN=CO'  => '1,3-oxazole',
    'C=CC=NO'  => '1,2-oxazole',
    'C=CC=NN'  => '1H-pyrazole', # FIXME: Adjust for isomerism
    'C=CC=CN'  => '1H-pyrole', # FIXME: Adjust for isomerism
    'C=CC=C[Se]' => 'selenophene',
    'C=CC=C[Te]' => 'tellurophene',
    'C=CC=CS'  => 'thiophene',

    # 6-membered aromatic
    'C=CC=CC=C' => 'benzene',
    'C=CC=CCO'  => '2H-pyran', # FIXME: Adjust for isomerism
    'C=CN=CC=N' => 'pyrazine',
    'C=CC=CN=N' => 'pyridazine',
    'C=CC=CC=N' => 'pyridine',
    'C=CC=NC=N' => 'pyrimidine',
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

    # Check for aromatic notation
    if( $SMILES =~ /:/ ) {
        $SMILES =~ s/([a-z]):/'' . uc( $1 )/ge;
        $SMILES =~ s/\[([a-z]{1,2})\]:/'[' . uc( $1 ) . ']'/ge;
        for my $SMILES_for_name (keys %names) {
            next unless $SMILES_for_name =~ /=/;
            my $name = $names{$SMILES_for_name};
            $SMILES_for_name =~ s/=//g;
            return $name if $SMILES eq $SMILES_for_name;
        }
    }

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

    if( $self->length >= 3 && $self->length <= 10 &&
        any { $_->{symbol} =~ /[cC]/ } $self->vertices ) {
        # Hantzsch-Widman names
        my @hetero = grep { $_->{symbol} !~ /[cC]/ } $self->vertices;
        die "cannot handle complicated monocycles for now\n" unless @hetero == 1; # TODO
        if( $self->length >= 7 ) {
            my @stems = ( 'epane', 'ocane', 'onane', 'ecane' );
            my $name = $elements{ucfirst $hetero[0]->{symbol}}->{prefix};
            $name =~ s/a$// if $stems[$self->length - 7] =~ /^[aeiou]/;
            return $name . $stems[$self->length - 7];
        }
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
