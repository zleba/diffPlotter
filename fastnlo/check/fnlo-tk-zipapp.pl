#!/usr/bin/env perl
#
# Unit test to read a fastNLO table
# Version:
#
# created by K. Rabbertz: 20.10.2016
#
#-----------------------------------------------------------------------
#
use Cwd;
use English;
use Getopt::Std;
use strict;
use warnings;

# Define source dir as required for VPATH builds like in make distcheck
my $src = ".";

# Remove potentially left-over temporary files
my $tab   = "fnl1014_I902309_2";
my $tabfl = ${tab}.".tab";
my $tabgz = ${tabfl}.".gz";
foreach my $file ( "zipapp.tab.gz", "zipapp.log", "zipcpp0.log", "zipcpp1.log", "zipappdiff.log", "zipcppdiff.log", "${tabfl}", "${tabgz}" ) {
    if ( -e $file ) {
        system("rm -f $file");
    }
}

# Prepare test setup
# Copy gzipped default tables from data/check storage
my $cmd = "cp -f ${src}/../data/check/${tabgz} .";
print "Executing command: $cmd\n";
my $ret = system("$cmd");
if ( $ret ) {die "fnlo-tk-zipapp: Copying test table ${tabgz} failed: $ret, aborted!\n";}

# Use append to normalise gzipped table directly and write out gzipped table
$cmd = "../src/fnlo-tk-append ${tabgz} zipapp.tab.gz > zipapp.log";
print "Executing command: $cmd\n";
$ret = system("$cmd");
if ( $ret ) {die "fnlo-tk-zipapp: Appending with gzipped test table ${tabgz} failed: $ret, aborted!\n";}

# Evaluate both tables, the only difference should be the NEVENT normalisation
$cmd = "../src/fnlo-tk-cppread ${tabgz} > zipcpp0.log";
print "Executing command: $cmd\n";
$ret = system("$cmd");
if ( $ret ) {die "fnlo-tk-zipapp: Evaluation gzipped test table ${tabgz} failed: $ret, aborted!\n";}
$cmd = "../src/fnlo-tk-cppread zipapp.tab.gz > zipcpp1.log";
print "Executing command: $cmd\n";
$ret = system("$cmd");
if ( $ret ) {die "fnlo-tk-zipapp: Evaluation gzipped test table zipapp.tab.gz failed: $ret, aborted!\n";}

# Determine difference between the two outputs and compare to expectation
$cmd = "diff zipcpp0.log zipcpp1.log | grep -v '^[0-9]' > zipappdiff.log";
print "Executing command: $cmd\n";
$ret = system("$cmd");
# Do not check this $ret here, since a difference is EXPECTED!
$cmd = "diff ${src}/../data/check/zipappdiff.log zipappdiff.log > zipcppdiff.log";
print "Executing command: $cmd\n";
$ret = system("$cmd");
# Now check, since here the difference should be gone
if ( $ret ) {die "fnlo-tk-zipapp: Result comparison with 'diff' failed: $ret, aborted!\n";}

# The diff.log must be empty
if ( ! -z "zipcppdiff.log" ) {
    print "fnlo-tk-zipapp: Evaluation of gzipped test and of normalised gzipped test table differs from default:\n";
    $cmd = "cat zipcppdiff.log";
    print "Executing command: $cmd\n";
    $ret = system("$cmd");
    print "fnlo-tk-zipapp: Do you use LHAPDF version 6? Do you have the CT10nlo PDF set installed? Is it found?\n";
    die "fnlo-tk-zipapp: Table evaluation unit test failed, please fix!\n";
}

print "fnlo-tk-zipapp: Table evaluation unit test passed.\n";

exit 0;
