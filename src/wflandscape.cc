/*
 * Simple landscape model under Wright-Fisher life cycle.
 * This example is handy because it shows us how to use
 * boost::geometry to implement landscape concepts while
 * being able to use fwdpp "as is".  Overlapping generations
 * will require some custom code b/c fwdpp doesn't currently
 * provide support for such life cycles
 */
#include "simtypes.hpp"
#include "wfrules.hpp"
#include <cassert> //fwdpp has this missing in one of its headers...
#include <cstdlib>
#include <functional>
#include <iostream>
#include <fwdpp/diploid.hh>  //Main fwdpp library header
#include <fwdpp/experimental/sample_diploid.hpp> //"Experimental" version of WF sampling function
#include <fwdpp/sugar/infsites.hpp> //Infinitely-many sites mutation scheme.  This works with popgenmut out of the box.
#include <fwdpp/sugar/GSLrng_t.hpp> //Smart-pointer wrapper around the gsl_rng *
#include <gsl/gsl_randist.h>
#include <boost/geometry/index/rtree.hpp>

namespace bg = boost::geometry;
namespace bgi = boost::geometry::index;

using rtree_type = bgi::rtree< landscape::csdiploid::value, bgi::quadratic<16> >;
using rules_type = landscape::WFLandscapeRules<rtree_type>;

int main(int argc, char ** argv)
{
    int argn = 1;
    const unsigned N = atoi(argv[argn++]);
    const double theta = atof(argv[argn++]);
    const double rho = atof(argv[argn++]);
    const double s = atof(argv[argn++]);       //selection coefficient
    const double h = atof(argv[argn++]);       //dominance.  Fitnesses will be 1,1+sh,1+2s, so h=1=additive.
    const double mu = atof(argv[argn++]);      //mutation rate to selected variants
    const double radius = atof(argv[argn++]);  //Radius in which to search for mates.
	const double dispersal = atof(argv[argn++]); //std. deviation in offspring dispersal
    const unsigned seed = atoi(argv[argn++]);  //RNG seed.

    //per-generation rates
    const double mu_n = theta/double(4*N);
    const double littler = rho/double(4*N);

    //This is our random number generator.
    //It is a simple wrapper around the gsl_rng *.
    //We do not need to free it.  Deletion of the object
    //takes care of that for us.
    KTfwd::GSLrng_t<KTfwd::GSL_RNG_MT19937> rng(seed);

    //Our genetic map is uniform on the interval (0,1]
    std::function<double(void)> recmap = std::bind(gsl_rng_uniform,rng.get());

    //This is our population:
    landscape::poptype pop(N);

    //Assign random points in space to our diploids.
    //The geometry is a square (0,0) to (1,1)
	rtree_type rtree;
	for(std::size_t i=0; i<N; ++i)
    {
        //Yes, I am aware that this does not truly
        //sample uniformly in 2d...
        double x = gsl_ran_flat(rng.get(),0.,1.);
        double y = gsl_ran_flat(rng.get(),0.,1.);
        pop.diploids[i].v = landscape::csdiploid::value(std::make_pair(landscape::csdiploid::point(x,y),i));
		rtree.insert(pop.diploids[i].v);
    }

    //pre-allocate space for a good guess as to the total # mutations
    //expected at equilibrium.
    pop.mutations.reserve(size_t(std::ceil(std::log(2*N)*theta+0.667*theta)));

	rules_type rules(std::move(rtree),radius,dispersal);
    for(unsigned generation = 0 ; generation < 10*N ; ++generation )
    {
        double wbar = KTfwd::experimental::sample_diploid(rng.get(),
                      pop.gametes,pop.diploids,pop.mutations,pop.mcounts,
                      N,mu_n+mu,
                      std::bind(KTfwd::infsites(),std::placeholders::_1,std::placeholders::_2,rng.get(),std::ref(pop.mut_lookup),generation,
        mu_n,mu,[&rng]() {
            return gsl_rng_uniform(rng.get());
        },[&s]() {
            return s;
        },[&h]() {
            return h;
        }),
        std::bind(KTfwd::poisson_xover(),rng.get(),littler,0.,1.,
                  std::placeholders::_1,std::placeholders::_2,std::placeholders::_3),
        std::bind(KTfwd::multiplicative_diploid(),std::placeholders::_1,std::placeholders::_2,
                  std::placeholders::_3,2.),
        pop.neutral,pop.selected,0,rules);
    }
}