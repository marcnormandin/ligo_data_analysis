/*
 * options.c
 *
 *  Created on: Mar 20, 2017
 *      Author: marcnormandin
 */

#include <stddef.h>
#include <stdio.h>
#include <stdlib.h>
#include <getopt.h>
#include <gsl/gsl_const_mksa.h>

#include "random.h"
#include "inspiral_signal.h"

/* Flag set by ‘--verbose’. */
static int verbose_flag = 0;
static int usage_flag = 0;

void print_usage() {
    printf("OVERVIEW: inspiral network analysis\n\n");
    printf("USAGE: lda [options] <inputs>\n\n");
    printf("OPTIONS:\n");
    printf("--help, --usage\t Print this help message\n");
    printf("--inspiral_ra, -r\t Set the inspiral sky right-ascension (radians -pi to pi)\n");
    printf("--inspiral_dec, -d\t Set the inspiral sky declination\n (radians -0.5*pi to 0.5*pi)\n");
    printf("--inspiral_polarization, -p\t Set the inspiral polarization\n");
    printf("--inspiral_coalesce_phase, -c\t Set the inspiral coalesce phase\n");
    printf("--inspiral_inclination, -i\t Set the inspiral inclination angle\n");
    printf("--inspiral_mass_one, -o\t Set the inspiral mass one (solar mass units)\n");
    printf("--inspiral_mass_two, -t\t Set the inspiral mass two (solar mass units)\n");
    printf("--inspiral_snr, -s\t Set the inspiral's network SNR\n");
    printf("--rng_seed, -z\t Set the random number generator seed. Use '0' for a random seed.\n");
}

int
process_command_options (int argc, char **argv, source_t *s, gslseed_t *seed, int *last_index)
{
  int c;

    static struct option long_options[] =
    {
        /* These options set a flag. */
        {"verbose", no_argument,       &verbose_flag, 1},
        {"brief",   no_argument,       &verbose_flag, 0},
        /* These options don’t set a flag.
         We distinguish them by their indices. */
        {"inspiral_ra",             required_argument, 0, 'r'},
        {"inspiral_dec",            required_argument, 0, 'd'},
        {"inspiral_polarization",   required_argument, 0, 'p'},
        {"inspiral_coalesce_phase", required_argument, 0, 'c'},
        {"inspiral_inclination",    required_argument, 0, 'i'},
        {"inspiral_mass_one",       required_argument, 0, 'o'},
        {"inspiral_mass_two",       required_argument, 0, 't'},
        {"inspiral_snr",            required_argument, 0, 's'},
		{"rng_seed",			    required_argument, 0, 'z'},
        {"help", no_argument, &usage_flag, 1},
        {"usage", no_argument, &usage_flag, 1},

        {0, 0, 0, 0}
    };

    int option_index = 0;

  while (1)
    {
      /* getopt_long stores the option index here. */


      /*c = getopt_long (argc, argv, "r:d:p:c:i:o:t:s:z:",
                       long_options, &option_index);*/

        c = getopt_long (argc, argv, "r:d:p:c:i:o:t:s:z:",
                         long_options, &option_index);

      /* Detect the end of the options. */
      if (c == -1)
        break;

      switch (c)
        {
        case 0:
          /* If this option set a flag, do nothing else now. */
                if (long_options[option_index].flag != 0)
            break;
          printf ("option %s", long_options[option_index].name);
          if (optarg)
            printf (" with arg %s", optarg);
          printf ("\n");
          break;

        case 'r':
          printf ("option --inspiral_ra, -r with value `%s'\n", optarg);
                s->sky.ra = atof(optarg);

          break;

        case 'd':
          printf ("option --inspiral_dec, -d with value `%s'\n", optarg);
                s->sky.dec = atof(optarg);
          break;

        case 'p':
          printf ("option --inspiral_polarization, -p with value `%s'\n", optarg);
                s->polarization_angle = atof(optarg);
          break;

            case 'c':
                printf ("option --inspiral_coalesce_phase, -c with value `%s'\n", optarg);
                s->coalesce_phase = atof(optarg);
                break;
            case 'i':
                printf ("option --inspiral_inclination, -i with value `%s'\n", optarg);
                s->inclination_angle = atof(optarg);
                break;
            case 'o':
                printf ("option --inspiral_mass_one, -o with value `%s'\n", optarg);
                s->m1 = atof(optarg) * GSL_CONST_MKSA_SOLAR_MASS;
                break;
            case 't':
                printf ("option --inspiral_mass_two, -t with value `%s'\n", optarg);
                s->m2 = atof(optarg) * GSL_CONST_MKSA_SOLAR_MASS;
                break;
            case 's':
                s->snr = atof(optarg);

                printf ("option --inspiral_snr, -s with value `%s'\n", optarg);
                break;
            case 'z':
            	*seed = atoi(optarg);
            	printf ("option --rng_seed, -z with value `%s'\n", optarg);
            	break;

        case '?':
          /* getopt_long already printed an error message. */
                print_usage();
          break;

        default:
          abort ();
        }
    }

  /* Instead of reporting ‘--verbose’
     and ‘--brief’ as they are encountered,
     we report the final status resulting from them. */
  /*if (verbose_flag)
    puts ("verbose flag is set");*/


  /* Print any remaining command line arguments (not options). */
  *last_index = optind;
  /*if (optind < argc)
    {
      printf ("non-option ARGV-elements: ");
      while (optind < argc)
        printf ("%s ", argv[optind++]);
      putchar ('\n');
    }
*/
    if (usage_flag) {
        print_usage();
        return -1;
    }

  return 0;
}




