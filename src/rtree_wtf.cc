/*
 * Why does the output of my simulation change
 * depending on how the rtree is constructed/filled?
 *
 * Are the search results affected?
 */

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
    bgi::rtree<value,bgi::quadratic<32>> rtree2;
    std::vector<value> temp;
    for(unsigned i=0; i<1000; ++i)
    {
        point p(2+gsl_ran_gaussian(rng.get(),0.025),
                2+gsl_ran_gaussian(rng.get(),0.025));
        rtree.insert(make_pair(p,i));
        rtree2.insert(make_pair(p,i));
        temp.push_back(make_pair(p,i));
    }
    decltype(rtree) rtree3(temp);
    decltype(rtree2) rtree4(temp);

    //OK, we now have 4 trees that have the same data, right?
    std::cout << rtree.size() << ' ' << rtree2.size() << ' '
              << rtree3.size() << ' ' << rtree4.size() << '\n';

    //Let's try a query based on finding all points in a box.
    box region(point(1.8,1.8),point(2.0,2.0));
	vector<value> v1,v2,v3,v4;	
	rtree.query(bgi::covered_by(region),std::back_inserter(v1));
	rtree2.query(bgi::covered_by(region),std::back_inserter(v2));
	rtree3.query(bgi::covered_by(region),std::back_inserter(v3));
	rtree4.query(bgi::covered_by(region),std::back_inserter(v4));

	std::cout << v1.size() << ' ' << v2.size() << ' ' << v3.size() << ' ' << v4.size() << '\n';
	for(std::size_t i=0;i<v1.size();++i)
	{
		std::cout << v1[i].first.get<0>() << ' ' << v1[i].first.get<1>() << ' ' << v1[i].second << '|';
		std::cout << v2[i].first.get<0>() << ' ' << v2[i].first.get<1>() << ' ' << v2[i].second << '|';
		std::cout << v3[i].first.get<0>() << ' ' << v3[i].first.get<1>() << ' ' << v3[i].second << '|';
		std::cout << v4[i].first.get<0>() << ' ' << v4[i].first.get<1>() << ' ' << v4[i].second << '\n';

	}
}
