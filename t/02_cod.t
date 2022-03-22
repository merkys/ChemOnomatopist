#!/usr/bin/perl

use strict;
use warnings;

use ChemOnomatopist;
use DBI;
use IPC::Run3;
use Test::More;

my $dbh = db_connect('mysql', 'www.crystallography.net', 'cod', 3306, 'cod_reader', '');

# FIXME: Skip tests if connection is unsuccessful. (A.M.)
# FIXME: Skip tests if OPSIN is not installed. (A.M.)

my $sth = $dbh->prepare( 'SELECT chemname, value AS smiles FROM data JOIN smiles ON file = cod_id WHERE chemname IS NOT NULL' );
$sth->execute;

my %tests;
while (my $item = $sth->fetchrow_hashref) {
    if ($item->{'smiles'} =~ /\A[CchH\[\]\(\)]*\z/) {
        $tests{$item->{'chemname'}} = $item->{'smiles'};
    }
}

my %opsin_approved;
for my $compound (keys %tests) {
    my( $smiles, $stderr );
    run3 'java -jar /usr/share/java/opsin.jar', \"$compound\n", \$smiles, \$stderr;
    chomp $smiles;
    next unless $smiles; # Skipping names that were not understood by OPSIN
    if( $smiles eq $tests{$compound} ) { # Ensuring OPSIN approval
        $opsin_approved{$smiles} = $compound;
    }
}

plan tests => scalar keys %opsin_approved;

for my $case (keys %opsin_approved) {
    # FIXME: Chemical name may have initial letter uppercased, but it may
    #        not be the right choice to lowercase it before comparison.
    #        Need to think a bit more on how to deal with it. (A.M.)
    is( ChemOnomatopist::get_name( $case ), $opsin_approved{$case} );
}

sub db_connect
{
    my ($db_platform, $db_host, $db_name, $db_port, $db_user, $db_pass) = @_;
    my $dsn = "dbi:$db_platform:hostname=$db_host;dbname=$db_name;" . 
              "user=$db_user;password=$db_pass";
    my $options = { PrintError => 0, mysql_enable_utf8 => 1 };
    my $dbh = DBI->connect( $dsn, $db_user, $db_pass, $options );
    die 'could not connect to the database - ' . lcfirst( $DBI::errstr ) unless $dbh;
    return $dbh
}
