dist_bin_SCRIPTS = \
	strain.txt 

common_sources = \
	antenna_patterns.c antenna_patterns.h \
	array.c array.h \
	chirp.c chirp.h \
	common.h \
	detector_network.c detector_network.h \
	detector.c detector.h \
	network_analysis.c network_analysis.h \
	options.c options.h \
	ptapso_estimate.c ptapso_estimate.h \
	ptapso_func.c ptapso_func.h \
	ptapso_maxphase.c ptapso_maxphase.h \
	ptapso.c ptapso.h \
	random.c random.h \
	sampling_system.c sampling_system.h \
	signal.c signal.h \
	simulate_inspiral.c simulate_inspiral.h \
	simulate_noise.c simulate_noise.h \
	sky.h \
	snr_sky_map.c snr_sky_map.h \
	source.c source.h \
	stationary_phase.c stationary_phase.h \
	strain_interpolate.c strain_interpolate.h \
	strain.c strain.h \
	time_delay.c time_delay.h
	
	
bin_PROGRAMS = lda_test_pso lda_test_sky_map lda_test_snr_hist

lda_test_pso_LDADD = -lgsl -lgslcblas
lda_test_pso_SOURCES = lda_test_pso_main.c $(common_sources)

lda_test_sky_map_LDADD = -lgsl -lgslcblas
lda_test_sky_map_SOURCES = lda_test_sky_map_main.c $(common_sources)

lda_test_snr_hist_LDADD = -lgsl -lgslcblas
lda_test_snr_hist_SOURCES = lda_test_snr_hist_main.c $(common_sources)
