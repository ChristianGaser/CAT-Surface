#!/usr/bin/perl

# converts .obj files to .vtk format
# these are mesh files for polygonal data.
# currently we only process v (vertex) and f (face), all other .obj tags are ignored.
# obj file format taken from: http://www.royriggs.com/obj.html

# check for proper argument
if (scalar(@ARGV) != 1 || ($infile = $ARGV[0]) !~ /\.obj$/) {
  print "Usage: obj2vtk filename.obj\n";
  exit (0);
}
$outfile = $infile;
# determine outfile name
$outfile =~ s/\.obj$/\.vtk/;

@points = ();
@polygons = ();
$numpolyvals = 0;

open (INFILE, $infile);
while ($line = <INFILE>)
{
	# parse each line
	chop($line);
	@parts = split(/ /, $line);
	# determine type of line
	$op = shift(@parts);
	if ($op eq "v") {
		# proper vertices have 3 values
  		if (scalar(@parts) == 3) { push (@points, join(" ", @parts) ); }
	} elsif ($op eq "f") {
		# remove all but the polygon number (normals, textures)
		for ($i = 0; $i < scalar(@parts); $i++) {
			($parts[$i],@junk) = split(/\//, $parts[$i]); 
			# also decrement poly vals by 1 - start with 0 instead of 1
			$parts[$i]--;
		}
		# put the number of points/sides on first
		unshift (@parts, scalar(@parts));
		# count the number of values in the polygon list
		$numpolyvals += scalar(@parts);
		# save the polygon
		push (@polygons, join(" ", @parts));
	}

}

close(INFILE);

open (OUTFILE, ">".$outfile);
print OUTFILE "# vtk DataFile Version 2.0\nvtk output\nASCII\nDATASET POLYDATA\n";
print OUTFILE "POINTS " . scalar(@points) . " float\n";
$" = "\n";
print OUTFILE "@points";
print OUTFILE "\nPOLYGONS " . scalar(@polygons) . " $numpolyvals\n";
$" = "\n";
print OUTFILE "@polygons";

close (OUTFILE);
