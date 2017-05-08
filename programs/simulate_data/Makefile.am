AM_CPPFLAGS = -I$(top_srcdir)/libcore -I$(top_srcdir)/libpso

bin_PROGRAMS = lda_simulate_data

lda_simulate_data_LDADD = ../../libcore/libcore.la ../../libpso/libpso.la

lda_simulate_data_SOURCES = \
	inspiral_signal.c \
	inspiral_signal.h \
	lda_simulate_data_main.c \
	simulation_options.c \
	simulation_options.h \
	simulate_noise.c \
	simulate_noise.h \
	simulation_file.c \
	simulation_file.h \
	simulation_settings.c \
	simulation_settings.h 
