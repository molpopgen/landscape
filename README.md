# Landscape simulations using [fwdpp](https://github.com/molpopgen/fwdpp)

## License

GPLv2 or later. 

## Based on:

* [fwdpp](https://github.com/molpopgen/fwdpp)
* R-trees from Boost [geometry](http://www.boost.org/doc/libs/1_57_0/libs/geometry/doc/html/index.html) library
* [GSL](http://gnu.org/software/gsl)

## Programs in src/ are:

### rtree_example.cc: 

This program will make some random points, then query for how many points are found in a certain box. This example
   does suggest that we should simply add pair<point,std::size_t> as a data member to a custom diploid for this sort of
   simulation.

Lessons learned:

1. boost < 1.5.9 has no method to iterate over an rtree without doing a query.  That is unfortunate--Ubuntu/Debian
   currently have 1.5.8 by default.
