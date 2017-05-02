# This is organized as a convenience library
AM_CPPFLAGS = -I$(top_srcdir)/libcore

noinst_LTLIBRARIES = libpso.la

libpso_la_SOURCES = \
	gbestpso.c \
	inspiral_pso_fitness.c \
	inspiral_pso_fitness.h \
	lbestpso.c \
	ptapso_maxphase.c \
	ptapso_maxphase.h \
	pso.c \
	pso.h

libpso_la_LDFLAGS = -fopenmp
libpso_la_LIBADD = ../libcore/libcore.la -lgsl -lgslcblas -lhdf5 -lhdf5_hl -lm
