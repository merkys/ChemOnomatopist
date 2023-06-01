package ChemOnomatopist::Group::Monocycle;

use strict;
use warnings;

# ABSTRACT: Monocyclic group
# VERSION

use parent ChemOnomatopist::Group::, ChemOnomatopist::Chain::Circular::;

use ChemOnomatopist;
use List::Util qw( all );

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

    # FIXME: For now we generate all possible traversals of the same cycle.
    #        This is not optimal, some caching could be introduced.
    my @chains;
    for (0..$#vertices) {
        push @chains, ChemOnomatopist::Chain::Circular->new( $graph, @vertices );
        push @vertices, shift @vertices;
    }
    @vertices = reverse @vertices;
    for (0..$#vertices) {
        push @chains, ChemOnomatopist::Chain::Circular->new( $graph, @vertices );
        push @vertices, shift @vertices;
    }

    my( $chain ) = ChemOnomatopist::filter_chains( @chains );

    return bless { graph => $graph, vertices => [ $chain->vertices ] }, $class;
}

sub needs_heteroatom_names()
{
    my( $self ) = @_;
    return $self->length < 3 || $self->length > 10 || all { $_->{symbol} !~ /^[cC]$/ } $self->vertices;
}

sub prefix()
{
    my( $self ) = @_;
    my $name = $self->name;
    $name =~ s/ane$/yl/;
    return $name;
}

# FIXME: This is a bit strange: class and object method with the same name
sub suffix()
{
    my( $self ) = @_;
    return '' unless ref $self;
    return $self->name;
}

1;
