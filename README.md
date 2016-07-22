# Landscape simulations using [fwdpp](https://github.com/molpopgen/fwdpp)

## License

GPLv2 or later. 

## Based on:

* [fwdpp](https://github.com/molpopgen/fwdpp)
* R-trees from Boost [geometry](http://www.boost.org/doc/libs/1_57_0/libs/geometry/doc/html/index.html) library
* [GSL](http://gnu.org/software/gsl)

## Programs in src/ are:

1. rtree_example.cc: make some random points, then query for how many points are found in a certain box. This example
   does suggest that we should simply add pair<point,std::size_t> as a data member to a custom diploid for this sort of
   simulation
