package ChemOnomatopist::Elements;

use strict;
use warnings;

# ABSTRACT: Element properties from IUPAC Blue Book
# VERSION

use parent Exporter::;
our @EXPORT_OK = qw( %elements );

# Replacement prefixes and seniorities are taken from BBv2 Appendix 1
our %elements = (
    F => {
        prefix => 'fluora',
        seniority => 0,
    },
    Cl => {
        prefix => 'chlora',
        seniority => 1,
    },
    Br => {
        prefix => 'broma',
        seniority => 2,
    },
    I => {
        prefix => 'ioda',
        seniority => 3,
    },
    At => {
        prefix => 'astata',
        seniority => 4,
    },
    O => {
        prefix => 'oxa',
        seniority => 5,
    },
    S => {
        prefix => 'thia',
        seniority => 6,
    },
    Se => {
        prefix => 'selena',
        seniority => 7,
    },
    Te => {
        prefix => 'tellura',
        seniority => 8,
    },
    Po => {
        prefix => 'polona',
        seniority => 9,
    },
    Lv => {
        prefix => 'livermora',
        seniority => 10,
    },
    N => {
        prefix => 'aza',
        seniority => 11,
    },
    P => {
        prefix => 'phospha',
        seniority => 12,
    },
    As => {
        prefix => 'arsa',
        seniority => 13,
    },
    Sb => {
        prefix => 'stiba',
        seniority => 14,
    },
    Bi => {
        prefix => 'bisma',
        seniority => 15,
    },
    C => {
        prefix => 'carba',
        seniority => 16,
    },
    Si => {
        prefix => 'sila',
        seniority => 17,
    },
    Ge => {
        prefix => 'germa',
        seniority => 18,
    },
    Sn => {
        prefix => 'stanna',
        seniority => 19,
    },
    Pb => {
        prefix => 'plumba',
        seniority => 20,
    },
    Fl => {
        prefix => 'flerova',
        seniority => 21,
    },
    B => {
        prefix => 'bora',
        seniority => 22,
    },
    Al => {
        prefix => 'alumina',
        seniority => 23,
    },
    Ga => {
        prefix => 'galla',
        seniority => 24,
    },
    In => {
        prefix => 'inda',
        seniority => 25,
    },
    Tl => {
        prefix => 'thalla',
        seniority => 26,
    },
    Zn => {
        prefix => 'zinca',
        seniority => 27,
    },
    Cd => {
        prefix => 'cadma',
        seniority => 28,
    },
    Hg => {
        prefix => 'mercura',
        seniority => 29,
    },
    Cn => {
        prefix => 'copernica',
        seniority => 30,
    },
    Cu => {
        prefix => 'cupra',
        seniority => 31,
    },
    Ag => {
        prefix => 'argenta',
        seniority => 32,
    },
    Au => {
        prefix => 'aura',
        seniority => 33,
    },
    Rg => {
        prefix => 'roentgena',
        seniority => 34,
    },
    Ni => {
        prefix => 'nickela',
        seniority => 35,
    },
    Pd => {
        prefix => 'pallada',
        seniority => 36,
    },
    Pt => {
        prefix => 'platina',
        seniority => 37,
    },
    Ds => {
        prefix => 'darmstadta',
        seniority => 38,
    },
    Co => {
        prefix => 'cobalta',
        seniority => 39,
    },
    Rh => {
        prefix => 'rhoda',
        seniority => 40,
    },
    Ir => {
        prefix => 'irida',
        seniority => 41,
    },
    Mt => {
        prefix => 'meitnera',
        seniority => 42,
    },
    Fe => {
        prefix => 'ferra',
        seniority => 43,
    },
    Ru => {
        prefix => 'ruthena',
        seniority => 44,
    },
    Os => {
        prefix => 'osma',
        seniority => 45,
    },
    Hs => {
        prefix => 'hassa',
        seniority => 46,
    },
    Mn => {
        prefix => 'mangana',
        seniority => 47,
    },
    Te => {
        prefix => 'techneta',
        seniority => 48,
    },
    Re => {
        prefix => 'rhena',
        seniority => 49,
    },
    Bh => {
        prefix => 'bohra',
        seniority => 50,
    },
    Cr => {
        prefix => 'chroma',
        seniority => 51,
    },
    Mo => {
        prefix => 'molybda',
        seniority => 52,
    },
    W => {
        prefix => 'tungsta',
        seniority => 53,
    },
    Sg => {
        prefix => 'seaborga',
        seniority => 54,
    },
    V => {
        prefix => 'vanada',
        seniority => 55,
    },
    Nb => {
        prefix => 'nioba',
        seniority => 56,
    },
    Ta => {
        prefix => 'tantala',
        seniority => 57,
    },
    Db => {
        prefix => 'dubna',
        seniority => 58,
    },
    Ti => {
        prefix => 'titana',
        seniority => 59,
    },
    Zr => {
        prefix => 'zircona',
        seniority => 60,
    },
    Hf => {
        prefix => 'hafna',
        seniority => 61,
    },
    Rf => {
        prefix => 'rutherforda',
        seniority => 62,
    },
    Sc => {
        prefix => 'scanda',
        seniority => 63,
    },
    Y => {
        prefix => 'yttra',
        seniority => 64,
    },
    La => {
        prefix => 'lanthana',
        seniority => 65,
    },
    Ce => {
        prefix => 'cera',
        seniority => 66,
    },
    Pr => {
        prefix => 'praseodyma',
        seniority => 67,
    },
    Nd => {
        prefix => 'neodyma',
        seniority => 68,
    },
    Pm => {
        prefix => 'prometha',
        seniority => 69,
    },
    Sm => {
        prefix => 'samara',
        seniority => 70,
    },
    Eu => {
        prefix => 'europa',
        seniority => 71,
    },
    Gd => {
        prefix => 'gadolina',
        seniority => 72,
    },
    Tb => {
        prefix => 'terba',
        seniority => 73,
    },
    Dy => {
        prefix => 'dysprosa',
        seniority => 74,
    },
    Ho => {
        prefix => 'holma',
        seniority => 75,
    },
    Er => {
        prefix => 'erba',
        seniority => 76,
    },
    Tm => {
        prefix => 'thula',
        seniority => 77,
    },
    Yb => {
        prefix => 'ytterba',
        seniority => 78,
    },
    Lu => {
        prefix => 'luteta',
        seniority => 79,
    },
    Ac => {
        prefix => 'actina',
        seniority => 80,
    },
    Th => {
        prefix => 'thora',
        seniority => 81,
    },
    Pa => {
        prefix => 'protactina',
        seniority => 82,
    },
    U => {
        prefix => 'urana',
        seniority => 83,
    },
    Np => {
        prefix => 'neptuna',
        seniority => 84,
    },
    Pu => {
        prefix => 'plutona',
        seniority => 85,
    },
    Am => {
        prefix => 'america',
        seniority => 86,
    },
    Cm => {
        prefix => 'cura',
        seniority => 87,
    },
    Bk => {
        prefix => 'berkela',
        seniority => 88,
    },
    Cf => {
        prefix => 'californa',
        seniority => 89,
    },
    Es => {
        prefix => 'einsteina',
        seniority => 90,
    },
    Fm => {
        prefix => 'ferma',
        seniority => 91,
    },
    Md => {
        prefix => 'mendeleva',
        seniority => 92,
    },
    No => {
        prefix => 'nobela',
        seniority => 93,
    },
    Lr => {
        prefix => 'lawrenca',
        seniority => 94,
    },
    Be => {
        prefix => 'berylla',
        seniority => 95,
    },
    Mg => {
        prefix => 'magnesa',
        seniority => 96,
    },
    Ca => {
        prefix => 'calca',
        seniority => 97,
    },
    Sr => {
        prefix => 'stronta',
        seniority => 98,
    },
    Ba => {
        prefix => 'bara',
        seniority => 99,
    },
    Ra => {
        prefix => 'rada',
        seniority => 100,
    },
    Li => {
        prefix => 'litha',
        seniority => 101,
    },
    Na => {
        prefix => 'soda',
        seniority => 102,
    },
    K => {
        prefix => 'potassa',
        seniority => 103,
    },
    Rb => {
        prefix => 'rubida',
        seniority => 104,
    },
    Cs => {
        prefix => 'caesa',
        seniority => 105,
    },
    Fr => {
        prefix => 'franca',
        seniority => 106,
    },
    He => {
        prefix => 'hela',
        seniority => 107,
    },
    Ne => {
        prefix => 'neona',
        seniority => 108,
    },
    Ar => {
        prefix => 'argona',
        seniority => 109,
    },
    Kr => {
        prefix => 'kryptona',
        seniority => 110,
    },
    Xe => {
        prefix => 'xenona',
        seniority => 111,
    },
    Rn => {
        prefix => 'radona',
        seniority => 112,
    },
);

# Hantzsch-Widman system prefixes, where different, taken from BBv2 P-22.2.2.1.1, Table 2.4
my %old_elements = (
    B => {
        standard_bonding_number => 3,
    },
    C => {
        standard_bonding_number => 4,
    },
    N => {
        standard_bonding_number => 3,
    },
    O => {
        standard_bonding_number => 2,
    },
    F => {
        standard_bonding_number => 1,
    },


    Al => {
        HantzschWidman => 'aluma',
        standard_bonding_number => 3,
    },
    Si => {
        standard_bonding_number => 4,
    },
    P => {
        standard_bonding_number => 3,
    },
    S => {
        standard_bonding_number => 2,
    },
    Cl => {
        standard_bonding_number => 1,
    },

    Ga => {
        standard_bonding_number => 3,
    },
    Ge => {
        standard_bonding_number => 4,
    },
    As => {
        standard_bonding_number => 3,
    },
    Se => {
        standard_bonding_number => 2,
    },
    Br => {
        standard_bonding_number => 1,
    },

    In => {
        HantzschWidman => 'indiga',
        standard_bonding_number => 3,
    },
    Sn => {
        standard_bonding_number => 4,
    },
    Sb => {
        standard_bonding_number => 3,
    },
    Te => {
        standard_bonding_number => 2,
    },
    I => {
        standard_bonding_number => 1,
    },

    Tl => {
        standard_bonding_number => 3,
    },
    Pb => {
        standard_bonding_number => 4,
    },
    Bi => {
        standard_bonding_number => 3,
    },
    Po => {
        standard_bonding_number => 2,
    },
    At => {
        standard_bonding_number => 1,
    },
);

for my $element (keys %old_elements) {
    for my $key (qw( HantzschWidman standard_bonding_number )) {
        next unless exists $old_elements{$element}->{$key};
        $elements{$element}->{$key} = $old_elements{$element}->{$key};
    }
}

1;
