#############################################
Sebastian E. Ramos-Onsins
Centre for Research in Agricultural Genomics
February 18, 2014
#############################################

The script "Recombination_spline_.R" calculates the
recombination rates per position given Genetic and Physic
maps:
First, the genetic map was plotted against the physical map
and incompatible values were discarded (that is, those
points that in the cumulative genetic map were not equal or
higher than the previous values). Cumulative recombination
curve for each chromosome were estimated from the Genetic
and physic maps using a monotonic cubic spline interpolation
method (implemented in the standard library of R,
http://www.r-project.org) using Hyman filtering (Hyman,
1983) over windows of sizes 50Kb, 100Kb, 500Kb, 1Mb and 5Mb.
The recombination value per position was obtained
calculating the slope per window (that is, the derivative).

References:
The R project for Statistical Computing.
http://www.r-project.org
Accurate Monotonicity preserving cubic interpolation. Hyman,
James M. SIAM Journal of Science Statistical Computating.
4(4). 1983.

Usage of the script "Recombination_spline_.R"

The script can be used on command-line with the following
arguments:
first argument: nchromosomes
second argument: file with the total length of each
chromosome, on rows.
third argument: file with the name of marker, number of
chromosome, cM and Physical position (in bp)
fourth and rest arguments: size or sizes (in bp) of the
windows.

Example:
R --no-save --args 12 ./lengths_chrs_melon.txt
./marker_melon_LG.txt 5e4 1e5 5e5 1e6 5e6 <
./Recombination_spline_.R


File examples are included in the release.

Sebastian


