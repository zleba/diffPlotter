#!/usr/bin/env perl
#
# Unit test to catenate two fastNLO tables
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
my @tabs = ("InclusiveNJetEvents_fnl5662i_v23_fix_2-2", "InclusiveNJetEvents_fnl5662i_v23_fix_0-0-2", "InclusiveNJetEvents_fnl5662i_v23_fix_2-2-2");
my @tabd = ("NJetEvents_2-2.tab.gz", "NJetEvents_0-0-2.tab.gz", "NJetEvents_2-2-2.tab.gz");
my @tabfls;
my @tabgzs;
foreach my $tab ( @tabs ) {
    push @tabfls, ${tab}.".tab";
    push @tabgzs, ${tab}.".tab".".gz";
}
foreach my $file ( "cattest.tab", "catdiff.log", @tabfls, @tabgzs ) {
    if ( -e $file ) {
        system("rm -f $file");
    }
}

# Prepare test setup
# Copy gzipped default tables from data/check storage
for ( my $i=0; $i <= $#tabd ; $i++ ) {
    my $cmd = "cp -f ${src}/../data/check/$tabd[$i] $tabgzs[$i]";
    print "Executing command: $cmd\n";
    my $ret = system("$cmd");
    if ( $ret ) {die "fnlo-tk-cattest: Copying test table $tabd[$i] failed: $ret, aborted!\n";}
    $cmd = "gunzip $tabgzs[$i]";
    print "Executing command: $cmd\n";
    $ret = system("$cmd");
    if ( $ret ) {die "fnlo-tk-cattest: Ungzipping test table $tabgzs[$i] failed: $ret, aborted!\n";}
}

# Catenate
my $cmd = "../src/fnlo-tk-cat $tabfls[0] $tabfls[1] cattest.tab";
print "Executing command: $cmd\n";
my $ret = system("$cmd");
if ( $ret ) {die "fnlo-tk-cattest: Catenating test tables $tabfls[0] and $tabfls[1] failed: $ret, aborted!\n";}

# Determine difference to default catenated table
$cmd = "diff $tabfls[2] cattest.tab > catdiff.log";
print "Executing command: $cmd\n";
$ret = system("$cmd");
if ( $ret ) {die "fnlo-tk-cattest: Result comparison with 'diff' failed: $ret, aborted!\n";}

# The diff.log must be empty
if ( ! -z "catdiff.log" ) {
    print "fnlo-tk-cattest: Catenated test table differs from default:\n";
    $cmd = "cat catdiff.log";
    print "Executing command: $cmd\n";
    $ret = system("$cmd");
    die "fnlo-tk-cattest: Catenation unit test failed, please fix!\n";
}

print "fnlo-tk-cattest: Catenation unit test passed.\n";

exit 0;
