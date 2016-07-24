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

/*
 * Concepts intended:
 * A point is a point in xy space.  Kinda obvious...
 * A box is a region in which we will search
 * A value associates a point with some other label.
 *
 * The notion of a value bring up two possibilities:
 * 1. We can add values like this to diploids in fwdpp,
 *    where the value represents the index of the diploid
 *    at that point in the diploids container.
 * 2. Try to use an rtree as a container of diploids.  
 *    This doesn't seem practical--there's no random-access
 *    to the objects. See below where we try to iterate.
 *
 * So it seems like a diploid should just have a "value"
 * as a data member...
 */
using point = bg::model::point<double, 2, bg::cs::cartesian>;
using box = bg::model::box<point>;
using value = std::pair<point, std::size_t>;

int main(int argc, char ** argv)
{
    std::cout << "sizeof point = " << sizeof(point) 
        << ", sizeof box = " << sizeof(box) 
        << ", sizeof value = " << sizeof(value) << '\n';

    //init random number generator
    KTfwd::GSLrng_t<KTfwd::GSL_RNG_MT19937> rng(101);

    // create the rtree using default constructor
    bgi::rtree< value, bgi::quadratic<16> > rtree;

    for(unsigned i=0; i<1000; ++i)
    {
        point p(2+gsl_ran_gaussian(rng.get(),0.025),
                2+gsl_ran_gaussian(rng.get(),0.025));
        rtree.insert(make_pair(p,i));
    }
    //So we can iterate over an rtree, but we have no operator[],
    //which makes some fwdpp-esque things quite inconvenient
    cout << "iterate over our tree:\n";
    for( const auto & v : rtree ) cout << v.first.get<0>() << ' ' << v.first.get<1>() << ' ' << v.second << '\n';
    //Do a search in our tree:
    vector<value> values;
    box region(point(1.8,1.8),point(2.25,2.25));
    rtree.query(bgi::covered_by(region),std::back_inserter(values));
    cout << "Search for elements within the box found " << values.size() << " items:\n";
    for(auto v : values)
    {
        cout << v.first.get<0>() << ' ' << v.first.get<1>() << ' ' << v.second << '\n';
    }
    //Try to remove elements from a tree:
    vector<value> values2;
    box region2remove(point(1.9,1.9),point(1.95,1.95));
    rtree.query(bgi::covered_by(region2remove),std::back_inserter(values2));
    for(auto v : values2) rtree.remove(v);

    //redo our initial search:
    values.clear();
    rtree.query(bgi::covered_by(region),std::back_inserter(values));
    cout << "found " << values.size() << " items\n";

}
