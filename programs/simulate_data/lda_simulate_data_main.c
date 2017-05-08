#ifdef HAVE_CONFIG_H
	#include "config.h"
#endif

#include <assert.h>

#include "detector_network.h"
#include "simulation_file.h"
#include "simulation_settings.h"

int main(int argc, char* argv[])
{
	size_t i, j;
	simulation_settings_t ps;
	simulation_settings_init( argc, argv, &ps );

	simulated_strain_file_create( ps.output_filename );
	simulated_strain_file_save_settings( ps.output_filename, &ps );

	detector_network_t *net = Detector_Network_load( ps.detector_mapping_filename,
			ps.num_time_samples, ps.sampling_frequency, ps.f_low, ps.f_high );
	simulated_strain_file_save_detector_network( ps.output_filename, net );

	simulate( &ps, net );

	Detector_Network_free( net );

	return 0;
}
