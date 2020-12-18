
/*********************************************************************************************

    This is public domain software that was developed by or for the U.S. Naval Oceanographic
    Office and/or the U.S. Army Corps of Engineers.

    This is a work of the U.S. Government. In accordance with 17 USC 105, copyright protection
    is not available for any work of the U.S. Government.

    Neither the United States Government, nor any employees of the United States Government,
    nor the author, makes any warranty, express or implied, without even the implied warranty
    of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE, or assumes any liability or
    responsibility for the accuracy, completeness, or usefulness of any information,
    apparatus, product, or process disclosed, or represents that its use would not infringe
    privately-owned rights. Reference herein to any specific commercial products, process,
    or service by trade name, trademark, manufacturer, or otherwise, does not necessarily
    constitute or imply its endorsement, recommendation, or favoring by the United States
    Government. The views and opinions of authors expressed herein do not necessarily state
    or reflect those of the United States Government, and shall not be used for advertising
    or product endorsement purposes.

*********************************************************************************************/

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <memory.h>
#include <errno.h>
#include <string.h>
#include <getopt.h>

#include "nvutility.h"

#include "pfm.h"
#include "chrtr2.h"


/*  Commented out until Dave Fabre fixes the tquantile code.  JCD 05/07/16

#include "tquantile.h"

*/

#include "version.h"


#define  MAX_BEAMS  1024


/*

  pfm_residual

  Jan C. Depner

  February 20th 2000

  This program compares two PFM structures.  It can compare the actual values of one structure against any
  of the surfaces of the other or compare the two surfaces.

*/



void usage ()
{
  fprintf (stderr, "\nUsage: pfm_residual [-M | -X | -A] [-b] PFM1 PFM2\n");
  fprintf (stderr, "Where\n\n");
  fprintf (stderr, "\t-M = compare minimum surface of PFM1 to data points in PFM2\n");
  fprintf (stderr, "\t-X = compare maximum surface of PFM1 to data points in PFM2\n");
  fprintf (stderr, "\t-A = compare average or MISP surface of PFM1 to data points\n");
  fprintf (stderr, "\t     in PFM2 (default)\n");
  fprintf (stderr, "\t-b = compare both surfaces to each other (this is only valid\n");
  fprintf (stderr, "\t     if the PFM structures are exactly the same size in terms\n");
  fprintf (stderr, "\t     of bin size and area).  A CHRTR2 file of the difference\n");
  fprintf (stderr, "\t     surface will be created if you use the -b option.  The file\n");
  fprintf (stderr, "\t     will be named PFM1.ch2.  Specify -M, -X, or -A option.\n\n");
  exit (-1);
}


/*  Commented out until Dave Fabre fixes the tquantile code.  JCD 05/07/16

    This is the sort function for qsort.

static int32_t compare_z (const void *a, const void *b)
{
    float *sa = (float *) (a);
    float *sb = (float *) (b);

    return (sa < sb ? 0 : 1);
}

*/


int32_t main (int32_t argc, char **argv)
{
  PFM_OPEN_ARGS   open_args[2];
  NV_I32_COORD2   coord[2], prev_coord = {-1, -1};
  BIN_RECORD      bin_record[2];
  DEPTH_RECORD    *depth;
  int32_t         i, j, k, m, percent, old_percent, neg_count[MAX_BEAMS], pos_count[MAX_BEAMS], option_index, chrtr2_handle = -1, surface_type = 0,
                  min_beams = MAX_BEAMS + 1, max_beams = -1, pfm_handle[2], recnum, total_bins;
  float           rms[MAX_BEAMS], max_val[MAX_BEAMS], min_val[MAX_BEAMS], neg_percent[MAX_BEAMS], pos_percent[MAX_BEAMS], depthtot[MAX_BEAMS],
                  min_depth[MAX_BEAMS], max_depth[MAX_BEAMS], diff = 0.0;
  double          sum[MAX_BEAMS], sum2[MAX_BEAMS], meandiff[MAX_BEAMS], meandepth[MAX_BEAMS], ss[MAX_BEAMS], var[MAX_BEAMS], stddev[MAX_BEAMS],
                  sddepth[MAX_BEAMS], dep = 0.0;
  NV_F64_COORD2   xy;
  uint8_t         both = NVFalse;
  CHRTR2_HEADER   chrtr2_header;
  CHRTR2_RECORD   chrtr2_record;
  char            chrtr2_file[512], diff_file[512], string[30];
  char            c;
  FILE            *fp2 = NULL, *fp24 = NULL, *fp3 = NULL, *fp4 = NULL;
  extern char     *optarg;
  extern int      optind;


  /*  Commented out until Dave Fabre fixes the tquantile code.  JCD 05/07/16

  int32_t         count = 0, two_count = 0, two4_count = 0, three_count = 0, four_count = 0;
  double          z[2000], median, uncertainty;
  TQ              tq[100];
  int32_t         maxdof = 30;
  double          ppp = 0.975;

  int janinit (TQ *, int);
  int janapply (int, TQ *, int, double *, double , double *, double *);

  */


  fprintf (stderr, "\n\n %s \n\n", VERSION);
  fflush (stderr);


  /*  Commented out until Dave Fabre fixes the tquantile code.  JCD 05/07/16

      janinit (tq, maxdof);

  */
 

  option_index = 0;
  while (NVTrue) 
    {
      static struct option long_options[] = {{0, no_argument, 0, 0}};

      c = (char) getopt_long (argc, argv, "MXAbSst:", long_options, &option_index);
      if (c == -1) break;

      switch (c) 
        {
        case 0:

          switch (option_index)
            {
            case 0:
              break;
            }
          break;

        case 'M':
          surface_type = 0;
          break;

        case 'X':
          surface_type = 1;
          break;

        case 'A':
          surface_type = 2;
          break;


          /*  Commented out until Dave Fabre fixes the tquantile code.  JCD 05/07/16

              This is an unadvertised statistics option.  It is used to report the difference between the minimum and average
              surface of two PFMs.  It outputs a CHRTR2 difference surface and 4 difference files in .pts format for
              display as overlays in pfmView (only differences that exceed the STD times the multiplier by more than one
              meter are output).  In practice you will usually specify the same PFM file for PFM1 and PFM2.

        case 's':
          surface_type = 3;
          both = NVTrue;
          break;


              This is an unadvertised statistics option.  It does the same thing as the -s option described above but it uses the
              median value instead of the average value.

        case 'S':
          surface_type = 4;
          both = NVTrue;
          break;

        case 't':
          surface_type = 5;
          sscanf (optarg, "%d", &maxdof);
          both = NVTrue;
          break;

          */


        case 'b':
          both = NVTrue;
          break;

        default:
          usage ();
          break;
        }
    }


  if (optind >= argc) usage ();


  strcpy (open_args[0].list_path, argv[optind]);
  strcpy (open_args[1].list_path, argv[optind + 1]);


  for (i = 0 ; i < MAX_BEAMS ; i++)
    {
      neg_count[i] = 0;
      pos_count[i] = 0;
      min_val[i] = 99999.0;
      max_val[i] = -99999.0;
      min_depth[i] = 99999.0;
      max_depth[i] = -99999.0;
      sum[i] = 0.0;
      sum2[i] = 0.0;
      depthtot[i] = 0.0;
    }


  /*  Open the files.  */

  open_args[0].checkpoint = open_args[1].checkpoint = 0;
  if ((pfm_handle[0] = open_existing_pfm_file (&open_args[0])) < 0) pfm_error_exit (pfm_error);
  if ((pfm_handle[1] = open_existing_pfm_file (&open_args[1])) < 0) pfm_error_exit (pfm_error);


  /*  Check for matching PFMs if -b or -s option was used.  */

  if (both)
    {
      if (open_args[0].head.bin_size_xy != open_args[1].head.bin_size_xy || 
          open_args[0].head.mbr.min_x != open_args[1].head.mbr.min_x ||
          open_args[0].head.mbr.min_y != open_args[1].head.mbr.min_y ||
          open_args[0].head.mbr.max_x != open_args[1].head.mbr.max_x ||
          open_args[0].head.mbr.max_y != open_args[1].head.mbr.max_y)
        {
          fprintf (stderr, "\n\nCannot use -b option.\n");
          fprintf (stderr, "PFM bin sizes or areas do not match\n\n");
          exit (-1);
        }


      memset (&chrtr2_header, 0, sizeof (CHRTR2_HEADER));


      /*  Generate the chrtr2 file name.  */

      strcpy (chrtr2_file, open_args[0].list_path);
      strcpy (&chrtr2_file[strlen (chrtr2_file) - 4], ".ch2");


      /*  Populate the chrtr2 header prior to creating the file.  */

      strcpy (chrtr2_header.creation_software, VERSION);
      chrtr2_header.z_units = CHRTR2_METERS;
      chrtr2_header.mbr.wlon = open_args[0].head.mbr.min_x;
      chrtr2_header.mbr.slat = open_args[0].head.mbr.min_y;
      chrtr2_header.width = open_args[0].head.bin_width;
      chrtr2_header.height = open_args[0].head.bin_height;
      chrtr2_header.lat_grid_size_degrees = open_args[0].head.y_bin_size_degrees;
      chrtr2_header.lon_grid_size_degrees = open_args[0].head.x_bin_size_degrees;
      chrtr2_header.min_z = -326.00;
      chrtr2_header.max_z = 326.00;
      chrtr2_header.z_scale = 100.0;
      chrtr2_header.max_number_of_points = 0;
      chrtr2_header.min_uncertainty = open_args[0].head.min_standard_dev;
      chrtr2_header.max_uncertainty = open_args[0].head.max_standard_dev;
      chrtr2_header.uncertainty_scale = open_args[0].scale;
      chrtr2_header.horizontal_uncertainty_scale = 0.0;
      chrtr2_header.vertical_uncertainty_scale = 0.0;


      /*  Try to create and open the chrtr2 file.  */

      chrtr2_handle = chrtr2_create_file (chrtr2_file, &chrtr2_header);
      if (chrtr2_handle < 0)
        {
          chrtr2_perror ();
          exit (-1);
        }


      /*  Generate the diff file names and open the diff files if we're using the -s option.  */

      if (surface_type >= 3)
        {
          if (surface_type == 3)
            {
              strcpy (string, "_avg_diff_");
            }
          else
            {
              strcpy (string, "_med_diff_");
            }

          strcpy (diff_file, open_args[0].list_path);
          strcpy (&diff_file[strlen (diff_file) - 4], string);
          strcat (diff_file, "1.96.pts");
          if ((fp2 = fopen (diff_file, "w")) == NULL)
            {
              perror (diff_file);
              exit (-1);
            }

          strcpy (diff_file, open_args[0].list_path);
          strcpy (&diff_file[strlen (diff_file) - 4], string);
          strcat (diff_file, "2.40.pts");
          if ((fp24 = fopen (diff_file, "w")) == NULL)
            {
              perror (diff_file);
              exit (-1);
            }

          strcpy (diff_file, open_args[0].list_path);
          strcpy (&diff_file[strlen (diff_file) - 4], string);
          strcat (diff_file, "3.00.pts");
          if ((fp3 = fopen (diff_file, "w")) == NULL)
            {
              perror (diff_file);
              exit (-1);
            }

          strcpy (diff_file, open_args[0].list_path);
          strcpy (&diff_file[strlen (diff_file) - 4], string);
          strcat (diff_file, "4.00.pts");
          if ((fp4 = fopen (diff_file, "w")) == NULL)
            {
              perror (diff_file);
              exit (-1);
            }
        }
    }


  fprintf(stderr, "\n\n");
  fflush (stderr);


  percent = 0;
  old_percent = -1;


  total_bins = open_args[1].head.bin_height * open_args[1].head.bin_width;


  /* Process all records in the PFM index file */

  for (i = 0 ; i < open_args[1].head.bin_height ; i++)
    {
      coord[1].y = i;

      for (j = 0 ; j < open_args[1].head.bin_width ; j++)
        {
          coord[1].x = j;


          if (both)
            {
              read_bin_record_index (pfm_handle[0], coord[1], &bin_record[0]);
              read_bin_record_index (pfm_handle[1], coord[1], &bin_record[1]);


              if ((bin_record[0].validity & PFM_DATA) && (bin_record[1].validity & PFM_DATA))
                {
                  switch (surface_type)
                    {
                    case 0:
                      dep = bin_record[1].min_filtered_depth;
                      diff = bin_record[0].min_filtered_depth - bin_record[1].min_filtered_depth;
                      break;

                    case 1:
                      dep = bin_record[1].max_filtered_depth;
                      diff = bin_record[0].max_filtered_depth - bin_record[1].max_filtered_depth;
                      break;

                    case 2:
                      dep = bin_record[1].avg_filtered_depth;
                      diff = bin_record[0].avg_filtered_depth - bin_record[1].avg_filtered_depth;
                      break;


                      /*  Commented out until Dave Fabre fixes the tquantile code.  JCD 05/07/16

                    case 3:
                      dep = bin_record[0].avg_filtered_depth;
                      diff = bin_record[0].avg_filtered_depth - bin_record[1].min_filtered_depth;
                      if (diff > bin_record[1].standard_dev * 4.0)
                        {
                          four_count++;
                          if (diff - (bin_record[1].standard_dev * 4.0) > 1.0) 
                            fprintf (fp4, "%.9f,%.9f,4.00/%.2f\n", bin_record[0].xy.y, bin_record[0].xy.x, diff - bin_record[1].standard_dev * 4.00);
                        }
                      if (diff > bin_record[1].standard_dev * 3.0)
                        {
                          three_count++;
                          if (diff - (bin_record[1].standard_dev * 3.0) > 1.0)
                            fprintf (fp3, "%.9f,%.9f,3.00/%.2f\n", bin_record[0].xy.y, bin_record[0].xy.x, diff - bin_record[1].standard_dev * 3.00);
                        }
                      if (diff > bin_record[1].standard_dev * 2.4)
                        {
                          two4_count++;
                          if (diff - (bin_record[1].standard_dev * 2.4) > 1.0)
                            fprintf (fp24, "%.9f,%.9f,2.40/%.2f\n", bin_record[0].xy.y, bin_record[0].xy.x, diff - bin_record[1].standard_dev * 2.40);
                        }
                      if (diff > bin_record[1].standard_dev * 1.96)
                        {
                          two_count++;
                          if (diff - (bin_record[1].standard_dev * 1.96) > 1.0)
                            fprintf (fp2, "%.9f,%.9f,1.96/%.2f\n", bin_record[0].xy.y, bin_record[0].xy.x, diff - bin_record[1].standard_dev * 1.96);
                        }
                      count++;
                      break;

                    case 4:
                      if (!read_depth_array_index (pfm_handle[1], coord[1], &depth, &recnum))
                        {
                          k = 0;
                          for (m = 0 ; m < recnum ; m++)
                            {
                              if ((!(depth[m].validity & (PFM_DELETED | PFM_INVAL | PFM_REFERENCE)) && dep < open_args[1].head.null_depth))
                                {
                                  z[k] = depth[m].xyz.z;
                                  k++;
                                }
                            }
                          free (depth);


                          if (k)
                            {
                              qsort (z, k, sizeof (float), compare_z);

                              if (k % 2)
                                {
                                  dep = z[k / 2];
                                }
                              else
                                {
                                  dep = (z[k / 2 - 1] + z[k / 2]) / 2.0;
                                }

                              diff = dep - bin_record[1].min_filtered_depth;
                              if (diff > bin_record[1].standard_dev * 4.0)
                                {
                                  four_count++;
                                  if (diff - (bin_record[1].standard_dev * 4.0) > 1.0) 
                                    fprintf (fp4, "%.9f,%.9f,4.00/%.2f\n", bin_record[0].xy.y, bin_record[0].xy.x, diff - bin_record[1].standard_dev * 4.00);
                                }
                              if (diff > bin_record[1].standard_dev * 3.0)
                                {
                                  three_count++;
                                  if (diff - (bin_record[1].standard_dev * 3.0) > 1.0)
                                    fprintf (fp3, "%.9f,%.9f,3.00/%.2f\n", bin_record[0].xy.y, bin_record[0].xy.x, diff - bin_record[1].standard_dev * 3.00);
                                }
                              if (diff > bin_record[1].standard_dev * 2.4)
                                {
                                  two4_count++;
                                  if (diff - (bin_record[1].standard_dev * 2.4) > 1.0)
                                    fprintf (fp24, "%.9f,%.9f,2.40/%.2f\n", bin_record[0].xy.y, bin_record[0].xy.x, diff - bin_record[1].standard_dev * 2.40);
                                }
                              if (diff > bin_record[1].standard_dev * 1.96)
                                {
                                  two_count++;
                                  if (diff - (bin_record[1].standard_dev * 1.96) > 1.0)
                                    fprintf (fp2, "%.9f,%.9f,1.96/%.2f\n", bin_record[0].xy.y, bin_record[0].xy.x, diff - bin_record[1].standard_dev * 1.96);
                                }
                              count++;
                            }
                        }
                      break;

                    case 5:
                      if (!read_depth_array_index (pfm_handle[1], coord[1], &depth, &recnum))
                        {
                          k = 0;
                          for (m = 0 ; m < recnum ; m++)
                            {
                              if ((!(depth[m].validity & (PFM_DELETED | PFM_INVAL | PFM_REFERENCE)) && dep < open_args[1].head.null_depth))
                                {
                                  z[k] = depth[m].xyz.z;
                                  k++;
                                }
                            }
                          free (depth);


                          if (k > 2)
                            {
                              janapply( maxdof, tq, k, z, ppp, &median, &uncertainty );

                              diff = bin_record[1].avg_filtered_depth - bin_record[1].min_filtered_depth;
                              if (diff - uncertainty > 1.0)
                                {
                                  two_count++;
                                  fprintf (fp2, "%.9f,%.9f,TQ %.2f\n", bin_record[0].xy.y, bin_record[0].xy.x, diff);
                                }

                              fprintf (fp24, "%.9f %.9f %.2f %.2f %.2f %.4f %.4f %d\n", bin_record[0].xy.y, bin_record[0].xy.x,
                                       bin_record[1].min_filtered_depth, median, bin_record[1].avg_filtered_depth, bin_record[1].standard_dev * 1.96,
                                       uncertainty, k);

                              count++;
                            }
                        }
                      break;

                      */
                    }


                  memset (&chrtr2_record, 0, sizeof (CHRTR2_RECORD));

                  chrtr2_record.z = diff;
                  chrtr2_record.uncertainty = bin_record[1].standard_dev;
                  chrtr2_record.status = CHRTR2_REAL;

                  chrtr2_write_record_row_col (chrtr2_handle, i, j, chrtr2_record);


                  if (dep < min_depth[0]) min_depth[0] = dep;

                  if (dep > max_depth[0]) max_depth[0] = dep;

                  depthtot[0] += dep;
                  sum[0] += diff;
                  sum2[0] += diff * diff;

                  if (diff < 0.0)
                    {
                      neg_count[0]++;
                    }
                  else
                    {
                      pos_count[0]++;
                    }

                  if (fabs((double) diff) < min_val[0]) min_val[0] = fabs((double) diff);

                  if (fabs((double) diff) > max_val[0]) max_val[0] = fabs((double) diff);
                }
            }
          else
            {
              if (!read_depth_array_index (pfm_handle[1], coord[1], &depth, &recnum))
                {
                  for (m = 0 ; m < recnum ; m++)
                    {
                      dep = depth[m].xyz.z;

                      if ((!(depth[m].validity & (PFM_DELETED | PFM_INVAL | PFM_REFERENCE)) && dep < open_args[1].head.null_depth))
                        {
                          k = depth[m].beam_number;

                          xy.y = depth[m].xyz.y;
                          xy.x = depth[m].xyz.x;

                          if (bin_inside_ptr (&open_args[0].head, xy))
                            {
                              compute_index_ptr (xy, &coord[0], &open_args[0].head);

                              if (prev_coord.x != coord[0].x || prev_coord.y != coord[0].y)
                                {
                                  read_bin_record_index (pfm_handle[0], coord[0], &bin_record[0]);
                                  prev_coord = coord[0];
                                }

                              if (bin_record[0].validity & PFM_DATA)
                                {
                                  switch (surface_type)
                                    {
                                    case 0:
                                      diff = bin_record[0].min_filtered_depth - dep;
                                      break;

                                    case 1:
                                      diff = bin_record[0].max_filtered_depth - dep;
                                      break;

                                    case 2:
                                      diff = bin_record[0].avg_filtered_depth - dep;
                                      break;
                                    }


                                  if (dep < min_depth[k]) min_depth[k] = dep;

                                  if (dep > max_depth[k]) max_depth[k] = dep;

                                  depthtot[k] += dep;
                                  sum[k] += diff;
                                  sum2[k] += diff * diff;

                                  if (diff < 0.0)
                                    {
                                      neg_count[k]++;
                                    }
                                  else
                                    {
                                      pos_count[k]++;
                                    }

                                  if (fabs((double) diff) < min_val[k]) min_val[k] = fabs((double) diff);

                                  if (fabs((double) diff) > max_val[k]) max_val[k] = fabs((double) diff);

                                  if (k > max_beams) max_beams = k;
                                  if (k < min_beams) min_beams = k;
                                }
                            }
                        }
                    }
                  free (depth);
                }
            }
        }

      percent = ((float) (i * open_args[1].head.bin_width + j) / (float) total_bins) * 100.0;
      if (old_percent != percent)
        {
          fprintf (stderr, "%03d%% processed     \r", percent);
          fflush (stderr);
          old_percent = percent;
        }
    }

  percent = 100;
  fprintf (stderr, "%03d%% processed        \n\n", percent);
  fflush (stderr);


  /* print out results */

  switch (surface_type)
    {
    case 0:
      printf ("#\n#SURFACE TYPE = MINIMUM\n");
      if (both)
        {
          printf ("#\n#This data represents the reference minimum filtered bin value\n");
          printf ("#from the first PFM file minus the reference minimum filtered bin value\n");
          printf ("#from the second PFM file.\n#\n");
        }
      else
        {
          printf ("#\n#This data represents the reference minimum filtered bin value\n");
          printf ("#from the first PFM file minus the depth values from the\n");
          printf ("#second PFM file.\n#\n");
        }
      break;

    case 1:
      printf ("#\n#SURFACE TYPE = MAXIMUM\n");
      if (both)
        {
          printf ("#\n#This data represents the reference maximum filtered bin value\n");
          printf ("#from the first PFM file minus the reference maximum filtered bin value\n");
          printf ("#from the second PFM file.\n#\n");
        }
      else
        {
          printf ("#\n#This data represents the reference maximum filtered bin value\n");
          printf ("#from the first PFM file minus the depth values from the\n");
          printf ("#second PFM file.\n#\n");
        }
      break;

    case 2:
      printf ("#\n#SURFACE TYPE = %s\n", open_args[0].head.average_filt_name);
      if (both)
        {
          printf ("#\n#This data represents the %s bin value\n", open_args[0].head.average_filt_name);
          printf ("#from the first PFM file minus the %s bin value\n", open_args[0].head.average_filt_name);
          printf ("#from the second PFM file.\n#\n");
        }
      else
        {
          printf ("#\n#This data represents the %s bin value\n", open_args[0].head.average_filt_name);
          printf ("#from the first PFM file minus the depth values from the\n");
          printf ("#second PFM file.\n#\n");
        }
      break;


      /*  Commented out until Dave Fabre fixes the tquantile code.  JCD 05/07/16

    case 3:
      printf ("#\n#This data represents the %s bin value\n", open_args[1].head.average_filt_name);
      printf ("#from the second PFM file minus the minimum bin value\n");
      printf ("#from the first PFM file.\n#\n");
      printf ("Points exceeding 2.0 standard deviations : %d\n", two_count);
      printf ("Points exceeding 2.4 standard deviations : %d\n", two4_count);
      printf ("Points exceeding 3.0 standard deviations : %d\n", three_count);
      printf ("Points exceeding 4.0 standard deviations : %d\n", four_count);
      printf ("Points : %d\n", count);
      break;

    case 4:
    case 5:
      printf ("#\n#This data represents the median bin value\n");
      printf ("#from the second PFM file minus the minimum bin value\n");
      printf ("#from the first PFM file.\n#\n");
      printf ("Points exceeding 2.0 standard deviations : %d\n", two_count);
      printf ("Points exceeding 2.4 standard deviations : %d\n", two4_count);
      printf ("Points exceeding 3.0 standard deviations : %d\n", three_count);
      printf ("Points exceeding 4.0 standard deviations : %d\n", four_count);
      printf ("Points : %d\n", count);
      break;

      */
    }

  printf ("#FIRST PFM file  : %s\n", open_args[0].list_path);
  printf ("#SECOND PFM file : %s\n#\n", open_args[1].list_path);


  if (both)
    {
      chrtr2_close_file (chrtr2_handle);


      printf ("#       RMS       MEAN DIFF          STD             STD%%    NEG%%   POS%%      MAX RESID    MEAN DEPTH    # POINTS\n#\n");

      meandiff[0] = sum[0] / (float) (neg_count[0] + pos_count[0]);
      meandepth[0] = depthtot[0] / (float) (neg_count[0] + pos_count[0]);
      ss[0] = sum2[0] - (sum[0] * meandiff[0]);
      var[0] = ss[0] / ((neg_count[0] + pos_count[0]) - 1);
      stddev[0] = sqrt (var[0]);
      sddepth[0] = (stddev[0] / meandepth[0]) * 100;
      rms[0] = sqrt((double) (sum2[0] / (float) (neg_count[0] + pos_count[0])));
      neg_percent[0] = ((float) neg_count[0] / (float) (neg_count[0] + pos_count[0])) * 100.0;
      pos_percent[0] = ((float) pos_count[0] / (float) (neg_count[0] + pos_count[0])) * 100.0;

      if (sum[0] != 0.0)
        {
          printf(" %10.3f   %10.3f      %10.3f      %10.4f    %03d    %03d   %10.3f    %10.3f  %12d\n", 
                 rms[0], meandiff[0], stddev[0], sddepth[0], NINT (neg_percent[0]), NINT (pos_percent[0]), 
                 max_val[0], meandepth[0], (neg_count[0] + pos_count[0]));
        }
    }
  else
    {
      printf ("# BEAM #     RMS       MEAN DIFF          STD             STD%%    NEG%%   POS%%      MAX RESID    MEAN DEPTH    # POINTS\n#\n");

      for (i = min_beams ; i <= max_beams ; i++)
        {
          meandiff[i] = sum[i] / (float) (neg_count[i] + pos_count[i]);
          meandepth[i] = depthtot[i] / (float) (neg_count[i] + pos_count[i]);
          ss[i] = sum2[i] - (sum[i] * meandiff[i]);
          var[i] = ss[i] / ((neg_count[i] + pos_count[i]) - 1);
          stddev[i] = sqrt (var[i]);
          sddepth[i] = (stddev[i] / meandepth[i]) * 100;
          rms[i] = sqrt((double) (sum2[i] / (float) (neg_count[i] + pos_count[i])));
          neg_percent[i] = ((float) neg_count[i] / (float) (neg_count[i] + pos_count[i])) * 100.0;
          pos_percent[i] = ((float) pos_count[i] / (float) (neg_count[i] + pos_count[i])) * 100.0;

          if (sum[i] != 0.0)
            {
              printf(" %3d   %10.3f   %10.3f      %10.3f      %10.4f    %03d    %03d   %10.3f    %10.3f  %12d\n", 
                     i, rms[i], meandiff[i], stddev[i], sddepth[i], NINT (neg_percent[i]), NINT (pos_percent[i]), 
                     max_val[i], meandepth[i], (neg_count[i] + pos_count[i]));
            }
        }
    }

  printf ("\n\n\n");


  if (both && surface_type >= 3)
    {
      fclose (fp2);
      fclose (fp24);
      fclose (fp3);
      fclose (fp4);
    }

  return (0);
}
