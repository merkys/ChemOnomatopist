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
        # Hantzsch-Widman names (BBv2 P-22.2.2.1)

        # Select the best numbering for heteroatoms
        my @chains;
        my @vertices = $self->vertices;
        my $cycle = $self->{graph}->subgraph( \@vertices ); # TODO: Add attributes
        for (0..$#vertices) {
            push @chains,
                 ChemOnomatopist::Chain::Circular->new( $cycle, @vertices );
            push @vertices, shift @vertices;
        }
        @vertices = reverse @vertices;
        for (0..$#vertices) {
            push @chains,
                 ChemOnomatopist::Chain::Circular->new( $cycle, @vertices );
            push @vertices, shift @vertices;
        }
        my $chain = ChemOnomatopist::filter_chains( @chains );
        @vertices = $chain->vertices;

        # Collect the types of heteroatoms and their attachment positions
        my %heteroatoms;
        for my $i (0..$#vertices) {
            next if ChemOnomatopist::is_element( $vertices[$i], 'C' );
            my $symbol = ucfirst $vertices[$i]->{symbol};
            $heteroatoms{$symbol} = [] unless $heteroatoms{$symbol};
            push @{$heteroatoms{$symbol}}, $i;
        }

        my $name = ChemOnomatopist::Name->new;
        my $least_senior_element;
        for my $element (sort { $elements{$a}->{seniority} <=> $elements{$b}->{seniority} }
                              keys %heteroatoms) {
            unless( scalar keys %heteroatoms == 1 &&
                    (@{$heteroatoms{$element}} == 1 ||
                     @{$heteroatoms{$element}} == $self->length - 1) ) {
                $name->append_locants( map { $_ + 1 } @{$heteroatoms{$element}} );
            }
            if( @{$heteroatoms{$element}} > 1 ) {
                $name->append_multiplier( ChemOnomatopist::IUPAC_numerical_multiplier( scalar @{$heteroatoms{$element}} ) );
            }
            $name->append_element( exists $elements{$element}->{HantzschWidman} ? $elements{$element}->{HantzschWidman} : $elements{$element}->{prefix} );
            $least_senior_element = $element;
        }
        $name->{name} =~ s/a$//;

        if(      $self->length <= 5 ) {
            my @stems = ( 'ir', 'et', 'ol' );
            $name .= $stems[$self->length - 3];
            $name .= $heteroatoms{N} ? 'idine' : 'ane';
            return $name;
        } elsif( $self->length == 6 ) {
            if( ($elements{$least_senior_element}->{seniority} >= 5 &&
                 $elements{$least_senior_element}->{seniority} <= 8) || $least_senior_element eq 'Bi' ) {
                return $name . 'ane';
            } else {
                return $name . 'inane';
            }
        } elsif( $self->length >= 7 ) {
            my @stems = ( 'epane', 'ocane', 'onane', 'ecane' );
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
