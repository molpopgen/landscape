# Landscape simulations using [fwdpp](https://github.com/molpopgen/fwdpp)

## License

GPLv2 or later. 

## Based on:

* [fwdpp](https://github.com/molpopgen/fwdpp)
* R-trees from Boost [geometry](http://www.boost.org/doc/libs/1_57_0/libs/geometry/doc/html/index.html) library
* [GSL](http://gnu.org/software/gsl)

## Programs in src/ are:

__NOTE:__ these programs are compiled in "debug" mode, allowing extensive checking of correctness of the data structures
behind the scenes.  Debug mode imposes a big runtime penalty.  "Production" code should be compiled using -DNDEBUG to
get rid of the debugging blocks.

### rtree_example.cc: 

This program will make some random points, then query for how many points are found in a certain box. This example
   does suggest that we should simply add pair<point,std::size_t> as a data member to a custom diploid for this sort of
   simulation.

Lessons learned:

1. boost < 1.5.9 has no method to iterate over an rtree without doing a query.  That is unfortunate--Ubuntu/Debian
   currently have 1.5.8 by default.

### wflandscape.cc

An implementation of a simple landscape model + Wright-Fisher sampling. This example serves to demonstrate how to
manipulate diploids, add new ones to landscapes, etc., and still be able to fwdpp's existing machinery.

Usage:

1. N = population size
2. theta = 4Nu
3. rho = 4Nr
4. s = selection coefficient for selected sites
6. h = dominance of selected mutations
7. mu_s = mutation rate to selected sites
8. radius = radius around individual to look for mates
9. dispersal = std. dev. of Gaussian disperal of offspring
10. seed = random number seed.

The model in brief:

* An individual is chosen proportional to fitness from the entire population
* Possible mates are discovered within a radius of the first individual.
* If no mates are found, the individual selfs.
* Othwerwise, the mate is chosen according to fitness.  Selfing can occur here, too.
* An offsprings location in x,y space is the midpoint of the parents + a Gaussian noise term added independently to each
  coordinate.
* The "landscape" is a square from [0,0] to [1,1].
* Initially, 1/2 the population is placed in the upper left and lower right quadrants.
* Selection coefficients have the opposite sign in the bottom left quadrant.

Implementation details:

* Custom rules class handles the landscape details.  The rules class design conforms to the specifications of fwdpp's
  "experimenta" API.
* Custom fitness function.
