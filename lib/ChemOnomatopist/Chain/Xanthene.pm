package ChemOnomatopist::Chain::Xanthene;

use strict;
use warnings;

# ABSTRACT: Xanthene or its close derivative
# VERSION

use parent ChemOnomatopist::Chain::Circular::;

use ChemOnomatopist::Elements qw( %elements );
use ChemOnomatopist::Util::Graph qw( subgraph );
use List::Util qw( any uniq );

sub new
{
    my( $class, $graph, @cycles ) = @_;

    my @benzenes = grep {  $_->is_benzene } @cycles;
    my( $other ) = grep { !$_->is_benzene } @cycles;

    # Find the correct vertex order
    my $subgraph = subgraph( $graph, map { $_->vertices } @cycles );
    my @bridges = grep { $subgraph->degree( $_->[0] ) == 3 &&
                         $subgraph->degree( $_->[1] ) == 3 }
                         $subgraph->edges;
    $subgraph->delete_edges( map { @$_ } @bridges );

    my @heteroatom_positions = $other->heteroatom_positions;
    my @heteroatoms = $other->heteroatoms;

    if( uniq( @heteroatoms ) == 2 ) {
        @heteroatom_positions = reverse @heteroatom_positions if $heteroatoms[0] eq 'N';
        @heteroatom_positions = reverse @heteroatom_positions if $heteroatoms[1] eq 'O';
        @heteroatom_positions = reverse @heteroatom_positions if join( ',', @heteroatoms ) eq 'As,S';
    }

    my @other_vertices = $other->vertices;
    $subgraph->delete_vertex( $other_vertices[($heteroatom_positions[0] + 2) % 6] );
    my( $start ) = grep { $subgraph->has_vertex( $_ ) && $subgraph->degree( $_ ) == 1 }
                   map  { $_->vertices } @benzenes;
    my @vertices = ( reverse( Graph::Traversal::DFS->new( $subgraph, start => $start )->dfs ),
                     $other_vertices[($heteroatom_positions[0] + 2) % 6] );

    return bless { graph => $graph, vertices => \@vertices }, $class;
}

sub candidates()
{
    my( $self ) = @_;
    my @candidates = ( $self );
    my @vertices = $self->vertices;

    push @candidates,
         bless { graph => $self->graph,
                 vertices => [ reverse @vertices[11..13], @vertices[0..10] ],
                 candidate_for => $self };

    if( $self->number_of_heteroatoms == 2 &&
        uniq( $self->heteroatoms ) == 1 ) { # TODO: Add two more candidates
    }

    return @candidates;
}

sub locants(@)
{
    my $self = shift;
    my @locant_map;
    if( $self->number_of_heteroatoms == 1 ) {
        @locant_map = ( 1..4, '4a', 10, '10a', 5..8, '8a',  9,  '9a' );
    } else {
        @locant_map = ( 1..4, '4a',  5,  '5a', 6..9, '9a', 10, '10a' );
    }
    return map { $locant_map[$_] } @_;
}

sub ideal_graph()
{
    return ChemOnomatopist::Chain::Polyacene->ideal_graph( 14 );
}

sub needs_heteroatom_locants() { return '' }
sub needs_heteroatom_names() { return '' }

sub prefix()
{
    my( $self ) = @_;

    my $name = $self->suffix;
    $name->[-1]{value} =~ s/e$//;

    if( $self->parent ) {
        my @vertices = $self->vertices;
        my( $position ) = grep { $self->graph->has_edge( $self->parent, $vertices[$_] ) } 0..$#vertices;
        die "unknown locant in multicyclic compound\n" unless defined $position;
        $name->append_substituent_locant( $self->locants( $position ) );
    }

    $name .= 'yl';
    return $name;
}

sub suffix()
{
    my( $self ) = @_;

    my  @heteroatoms = $self->heteroatoms;
    if( @heteroatoms == 1 ) {
        my $name = ChemOnomatopist::Name::Part::Locants->new( '9H-' )->to_name;
        my $stem = '';
        if( $heteroatoms[0] ne 'O' ) {
            $stem .= $elements{$heteroatoms[0]}->{prefix};
            $stem =~ s/a$/o/;
        }
        $stem .= 'xanthene';
        return $name->append_stem( $stem );
    } elsif( @heteroatoms == 2 && $heteroatoms[0] eq $heteroatoms[1] &&  $heteroatoms[0] eq 'N' ) {
        return ChemOnomatopist::Name::Part::Stem->new( 'phenazine' )->to_name;
    } elsif( @heteroatoms == 2 && $heteroatoms[0] eq $heteroatoms[1] ) {
        my $name = $elements{$heteroatoms[0]}->{prefix};
        $name =~ s/a$//;
        return ChemOnomatopist::Name::Part::Stem->new( $name . 'anthrene' )->to_name;
    } elsif( @heteroatoms == 2 && $heteroatoms[1] eq 'N' ) {
        # BBv2 P-25.2.2.3
        my $name = ChemOnomatopist::Name::Part::Locants->new( '10H-' )->to_name;
        my $stem = 'pheno';
        $stem =~ s/o$// if $elements{$heteroatoms[0]}->{prefix} =~ /^o/;
        $name->append_stem( $stem . $elements{$heteroatoms[0]}->{prefix} . 'zine' );
        return $name;
    } elsif( @heteroatoms == 2 && $heteroatoms[0] eq 'O' ) {
        # BBv2 P-25.2.2.3
        return 'phenoxathiine' if $heteroatoms[1] eq 'S';
        if(      any { $heteroatoms[1] eq $_ } qw( Se Te ) ) {
            my $stem = 'phenoxa' . $elements{$heteroatoms[1]}->{prefix};
            $stem =~ s/a$/ine/;
            return ChemOnomatopist::Name::Part::Stem->new( $stem )->to_name;
        } elsif( any { $heteroatoms[1] eq $_ } qw( P As Sb ) ) {
            my $stem = 'phenoxa' . $elements{$heteroatoms[1]}->{prefix};
            $stem =~ s/a$/inine/;
            return ChemOnomatopist::Name::Part::Stem->new( $stem )->to_name;
        }
    } elsif( @heteroatoms == 2 && $heteroatoms[0] eq 'S' && $heteroatoms[1] eq 'As' ) {
        # BBv2 P-25.2.2.3
        return ChemOnomatopist::Name::Part::Stem->new( 'phenothiarsinine' )->to_name;
    }

    die "cannot name xanthene derivative\n";
}

1;
