# This is organized as a convenience library
AM_CPPFLAGS = -I$(top_srcdir)/libcore

noinst_LTLIBRARIES = libpso.la

libpso_la_SOURCES = \
	gbestpso.c \
	inspiral_pso_fitness.c \
	inspiral_pso_fitness.h \
	lbestpso.c \
	parallel.c \
	parallel.h \
	ptapso_maxphase.c \
	ptapso_maxphase.h \
	pso.c \
	pso.h \
	spso.c

libpso_la_LDFLAGS = 

if HAVE_OPENMP
libpso_la_LDFLAGS += $OPENMP_CFLAGS
endif

libpso_la_LIBADD = ../libcore/libcore.la -lgsl -lgslcblas -lhdf5 -lhdf5_hl -lm
