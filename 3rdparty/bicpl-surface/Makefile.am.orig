SUBDIRS = Include Data_structures Transforms Geometry Numerical \
	Objects Prog_utils bicpl_clapack .

AUTOMAKE_OPTIONS = check-news
ACLOCAL_AMFLAGS = -I m4

EXTRA_DIST = README.release \
	autogen.sh \
	epm-header.in \
	Makefile.msvc-win32 \
	config.h.msvc-win32 \
	$(m4_files)

m4_files = m4/mni_REQUIRE_LIB.m4                \
	m4/mni_REQUIRE_OPENINVENTOR.m4          \
	m4/mni_cxx_have_koenig_lookup.m4        \
	m4/smr_CGAL_MAKEFILE.m4                 \
	m4/smr_OPTIONAL_LIB.m4                  \
	m4/smr_REQUIRED_LIB.m4                  \
	m4/smr_WITH_BUILD_PATH.m4

# version info argument is CURRENT[:REVISION[:AGE]]
# see README.release for update instructions
# DO NOT SET THIS TO BE THE SAME AS THE PACKAGE VERSION!!!
libbicpl_la_LDFLAGS = -version-info 3:0:2

lib_LTLIBRARIES = libbicpl.la

libbicpl_la_SOURCES = version.c
libbicpl_la_LIBADD = \
	$(top_builddir)/Data_structures/libbicpl_ds.la \
	$(top_builddir)/Numerical/libbicpl_n.la \
	$(top_builddir)/Prog_utils/libbicpl_pu.la \
	$(top_builddir)/Geometry/libbicpl_g.la \
	$(top_builddir)/Objects/libbicpl_o.la \
	$(top_builddir)/Transforms/libbicpl_t.la \
	$(top_builddir)/bicpl_clapack/libbicpl_cla.la


