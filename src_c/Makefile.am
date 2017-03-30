dist_bin_SCRIPTS = \
	strain.txt 

common_sources = \
	common.h \
	detector_antenna_patterns.c detector_antenna_patterns.h \
	detector_network.c detector_network.h \
	detector_time_delay.c detector_time_delay.h \
	detector.c detector.h \
	inspiral_chirp.c inspiral_chirp.h \
	inspiral_network_statistic.c inspiral_network_statistic.h \
	inspiral_pso_fitness.c inspiral_pso_fitness.h \
	inspiral_signal.c inspiral_signal.h \
	inspiral_source.c inspiral_source.h \
	inspiral_stationary_phase.c inspiral_stationary_phase.h \
	options.c options.h \
	ptapso_maxphase.c ptapso_maxphase.h \
	ptapso.c ptapso.h \
	random.c random.h \
	sampling_system.c sampling_system.h \
	settings.c settings.h \
	signal.c signal.h \
	simulate_noise.c simulate_noise.h \
	sky.h \
	strain_interpolate.c strain_interpolate.h \
	strain.c strain.h	
	
bin_PROGRAMS = lda_test_pso lda_mpi_pso lda_test_sky_map lda_test_snr_hist

lda_test_pso_LDADD = -lgsl -lgslcblas
lda_test_pso_SOURCES = lda_test_pso_main.c $(common_sources)

lda_mpi_pso_mpi_CFLAGS=-DHAVE_MPI
lda_mpi_pso_CXXFLAGS=-DHAVE_MPI
lda_mpi_pso_CPPFLAGS=-DHAVE_MPI

lda_mpi_pso_LDADD = -lgsl -lgslcblas -lmpi
lda_mpi_pso_SOURCES = lda_test_pso_main.c $(common_sources)


lda_test_sky_map_LDADD = -lgsl -lgslcblas
lda_test_sky_map_SOURCES = lda_test_sky_map_main.c $(common_sources)

lda_test_snr_hist_LDADD = -lgsl -lgslcblas
lda_test_snr_hist_SOURCES = lda_test_snr_hist_main.c $(common_sources)
