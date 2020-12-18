
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


/*********************************************************************************************

    This program is public domain software that was developed by 
    the U.S. Naval Oceanographic Office.

    This is a work of the US Government. In accordance with 17 USC 105,
    copyright protection is not available for any work of the US Government.

    This software is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.

*********************************************************************************************/

#ifndef VERSION

#define     VERSION     "PFM Software - pfm_residual V3.00 - 05/07/14"

#endif

/*

    Version 1.0
    Jan C. Depner
    02/20/00


    Version 1.1
    Jan C. Depner
    03/29/00

    Uses version 2.1 of the pfm library.


    Version 1.2
    Jan C. Depner
    09/04/00

    Replaced call to read_depth_record_index with read_depth_array_index.


    Version 1.3
    Jan C. Depner
    02/21/01

    Passing scale to open_pfm_file as a pointer.


    Version 1.4
    Jan C. Depner
    06/21/01

    Passing structure args to open_pfm_file.


    Version 1.5
    Jan C. Depner
    07/19/01
 
    4.0 PFM library changes.


    Version 1.6
    Jan C. Depner
    10/07/01
 
    Precede all comments with # for gnuplot.


    Version 1.7
    Jan C. Depner
    10/11/01
 
    Increased limits to 256 beams.


    Version 1.8
    Jan C. Depner
    10/25/01
 
    Bug in test for null depth.


    Version 1.81
    Jan C. Depner
    12/16/04
 
    Changed Usage message for PFM 4.5 directory input.


    Version 1.82
    Jan C. Depner
    02/25/05

    Switched to open_existing_pfm_file from open_pfm_file.


    Version 1.83
    Jan C. Depner
    03/04/05

    Fix return from open_existing_pfm_file.


    Version 1.84
    Jan C. Depner
    10/26/05

    Changed usage for PFM 4.6 handle file use.


    Version 1.85
    Jan C. Depner
    08/06/07

    Cleaned things up a bit and found a serious bug in PFM.


    Version 2.0
    Jan C. Depner
    08/21/07

    Added ability to use different surfaces and to compare two surfaces.


    Version 2.01
    Jan C. Depner
    09/17/07

    Replaced compute_index with compute_index_ptr.


    Version 2.02
    Jan C. Depner
    10/22/07

    Added fflush calls after prints to stderr since flush is not automatic in Windows.


    Version 2.03
    Jan C. Depner
    04/07/08

    Replaced single .h files from utility library with include of nvutility.h


    Version 2.04
    Jan C. Depner
    01/29/09

    Set checkpoint to 0 prior to calling open_existing_pfm_file.


    Version 2.05
    Jan C. Depner
    07/12/10

    Added output of CHRTR2 file when using the -b option.


    Version 2.06
    Jan C. Depner
    08/23/10

    Added a hidden option to produce some statistics for analyzing the difference between the minimum and
    average surfaces.


    Version 2.07
    Jan C. Depner
    05/06/11

    Fixed problem with getopt that only happens on Windows.


    Version 2.08
    Jan C. Depner
    09/21/11

    Replaced bin_inside calls with bin_inside_ptr calls.

   Version 2.09
   Stacy Johnson

   Switched %ld and %lf to %d and %f respectively (only in printf type statements, they're valid, and needed, in
   scanf type statements).


   Version 2.10
   Jan C. Depner (PFM Software)
   07/23/14

   - Switched from using the old NV_INT64 and NV_U_INT32 type definitions to the C99 standard stdint.h and
     inttypes.h sized data types (e.g. int64_t and uint32_t).


   Version 2.11
   Jan C. Depner (PFM Software)
   07/29/14

   - Fixed some of the errors discovered by cppcheck.  Please look at tquantile.c - it is FUBAR.


   Version 3.00
   Jan C. Depner (PFM Software)
   05/07/16

   - I have removed all of the tquantile code since the original author has not fixed it (in two years) and it is now breaking
     the program.  If he wants to fix it I'll put it back in at some point in the future.  The options it was used for were not
     documented anyway so I doubt that anyone will care.

*/
