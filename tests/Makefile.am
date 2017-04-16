AM_CPPFLAGS = -I$(top_srcdir)/libcore -I$(top_srcdir)/libpso

bin_PROGRAMS = 

if HAVE_GTEST
bin_PROGRAMS += test_libcore
test_libcore_SOURCES = libcore_test.cpp
test_libcore_LDADD = -lgtest_main -lgtest -lhdf5 -lhdf5_hl
endif
