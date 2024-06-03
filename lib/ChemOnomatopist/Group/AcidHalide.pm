package ChemOnomatopist::Group::AcidHalide;

# ABSTRACT: Acid halide group
# VERSION

use parent ChemOnomatopist::Group::;

use ChemOnomatopist::Elements qw( %elements );

sub new
{
    my( $class, $group, $element ) = @_;
    return bless { group => $group, element => $element }, $class;
}

sub element() { $_[0]->{group}->element }

sub prefix() { $_[0]->suffix }
sub suffix()
{
    my( $self ) = @_;
    my $name = $self->{group}->prefix . ' ';
    my $halide = $elements{$self->{element}}->{prefix};
    $halide =~ s/a$/ide/;
    return $name . $halide;
}

1;
