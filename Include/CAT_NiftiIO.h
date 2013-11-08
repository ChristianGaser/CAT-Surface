/* Christian Gaser - christian.gaser@uni-jena.de
 * Department of Psychiatry
 * University of Jena
 *
 * Copyright Christian Gaser, University of Jena.
 * $Id: CAT_SurfaceIO.h 304 2013-11-07 15:16:36Z gaser $
 *
*/

#include  <bicpl.h>
#include <limits.h>
#include <float.h>

#include "niftilib/nifti1_io.h"
#include "niftilib/nifti1_local.h"
#include "niftilib/analyze75.h"

Status input_volume_all(char *, int, char **, nc_type, BOOLEAN,
                 double, double, BOOLEAN, Volume *, minc_input_options *);
