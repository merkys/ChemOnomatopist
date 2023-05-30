# ChemOnomatopist

*ChemOnomatopist* is a tool to derive IUPAC systematic names for chemical structures:

    $ echo C=C1C=CC=C1 | chemonomatopist
    5-methylidenecyclopenta-1,3-diene

*ChemOnomatopist* analyses chemical graphs to determine IUPAC names according to the [Nomenclature of Organic Chemistry. IUPAC Recommendations
and Preferred Names 2013](https://iupac.qmul.ac.uk/BlueBook/PDF/BlueBookV2.pdf).

## Installation

*ChemOnomatopist* is written in Perl.
The easiest way to install it is using [Dist::Zilla](https://metacpan.org/release/Dist-Zilla):

    $ git clone https://github.com/merkys/ChemOnomatopist
    $ dzil install

## License

*ChemOnomatopist* is free software licensed under BSD-3-Clause license.
