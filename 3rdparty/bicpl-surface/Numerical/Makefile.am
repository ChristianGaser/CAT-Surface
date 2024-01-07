AM_CPPFLAGS = -I$(top_srcdir)/Include

noinst_LTLIBRARIES = libbicpl_n.la
libbicpl_n_la_SOURCES = \
                 amoeba.c \
                 gaussian.c \
                 gradient_minimize.c \
                 histogram.c \
                 least_squares.c \
                 matrix_svd.c \
                 minimize_lsq.c \
                 minimize_lsq_float.c \
                 numerical.c \
                 quadratic.c \
                 real_quadratic.c \
                 statistics.c \
                 t_stat.c


# Despite the name ending in '.c', this is an #included file.
noinst_HEADERS = minimize_lsq_include.c
