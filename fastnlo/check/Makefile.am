#                                               -*- Automake -*-
# -----------------------------------------------------------------------------
#
# K. Rabbertz, 20.10.2016
#
# fastNLO_toolkit
#
# -----------------------------------------------------------------------------
# Use GNU standards
AUTOMAKE_OPTIONS = gnu

#
# Targets
#
# perl script targets to check the package (distributed, but not installed)
#
dist_check_SCRIPTS = fnlo-tk-cpptest.pl fnlo-tk-cppnormtest.pl fnlo-tk-cppscalestest.pl fnlo-tk-cattest.pl \
	fnlo-tk-erasetest.pl fnlo-tk-erasetest2.pl fnlo-tk-mergetest.pl fnlo-tk-stattest.pl
if HAVE_LIBZ
   dist_check_SCRIPTS += fnlo-tk-ziptest.pl fnlo-tk-zipapp.pl
endif
TESTS = $(dist_check_SCRIPTS)

# Additional clean up for check targets
CLEANFILES = *.log *.tab *.tab.gz *.str *.txt *.yoda

# Files that are required to build the project but that can be
# recreated IF the necessary autotools etc. are available by:
# "autoreconf -i; ./configure"
# Use "make maintainer-clean" (on a build machine) before svn commits,
# but not for distribution
MAINTAINERCLEANFILES = Makefile.in
