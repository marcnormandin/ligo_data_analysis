# This is organized as a convenience library
AM_CPPFLAGS = -I$(top_srcdir)/libcore

noinst_LTLIBRARIES = libpso.la

libpso_la_SOURCES = \
	inspiral_pso_fitness.c \
	inspiral_pso_fitness.h \
	ptapso_maxphase.c \
	ptapso_maxphase.h \
	ptapso.c \
	ptapso.h

libpso_la_LIBADD = -lgsl -lgslcblas -lhdf5 -lhdf5_hl -lm
