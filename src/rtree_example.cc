/*
 * This is an example of using boost's rtree
 * data types from the boost::geometry library
 * */

#include <boost/geometry.hpp>
#include <boost/geometry/geometries/point.hpp>
#include <boost/geometry/geometries/box.hpp>

#include <boost/geometry/index/rtree.hpp>

// to store queries results
#include <vector>
#include <utility>

// just for output
#include <iostream>
#include <boost/foreach.hpp>

//use fwdpp's smart pointer around gsl_rng
#include <fwdpp/sugar/GSLrng_t.hpp>

#include <gsl/gsl_randist.h>

namespace bg = boost::geometry;
namespace bgi = boost::geometry::index;

using namespace std;

//typedefs for simplity
using point = bg::model::point<double, 2, bg::cs::cartesian>;
using box = bg::model::box<point>;
using value = std::pair<box, unsigned>;

int main(int argc, char ** argv)
{
    //init random number generator
    KTfwd::GSLrng_t<KTfwd::GSL_RNG_MT19937> rng(101);

    // create the rtree using default constructor
    bgi::rtree< point, bgi::quadratic<16> > rtree;

    for(unsigned i=0; i<1000; ++i)
    {
        point p(2+gsl_ran_gaussian(rng.get(),0.025),
                2+gsl_ran_gaussian(rng.get(),0.025));
        rtree.insert(p);
    }
    vector<point> values;
    box region(point(1.8,1.8),point(2.25,2.25));
    rtree.query(bgi::covered_by(region),std::back_inserter(values));
    for(auto v : values)
    {
        cout << v.get<0>() << ' ' << v.get<1>() << '\n';
    }
}
