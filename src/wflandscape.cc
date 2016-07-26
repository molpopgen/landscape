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
#include <fwdpp/sugar/sampling.hpp>
#include <Sequence/SimData.hpp>
#include <gsl/gsl_randist.h>
#include <boost/geometry/index/rtree.hpp>

namespace bg = boost::geometry;
namespace bgi = boost::geometry::index;

//typedefs to simplify life
using rtree_type = bgi::rtree< landscape::csdiploid::value, bgi::quadratic<16> >;
using rules_type = landscape::WFLandscapeRules<rtree_type>;

//Arbitrary model for fitness.
//Treat s as -s in the
//lower left quandrant of the landscape,
//otherwise as s.
//Fitness is multiplicative across sites.
struct spatial_fitness
{
    /* This function makes a spatial fitness object
     * behave as a function.
     * Note: if we had an additional landscape that reflected
     * how the landscape modified genetic values of fitness,
     * we could bind it here.
     */
    inline double operator()(const landscape::csdiploid & dip,
                             const std::vector<KTfwd::gamete> & gametes,
                             const std::vector<KTfwd::popgenmut> & mutations) const
    {
        double x = boost::geometry::get<0>(dip.v.first);
        double y = boost::geometry::get<1>(dip.v.first);
        KTfwd::site_dependent_fitness s;
        double geographic_factor = 1.0;
        if(x<=0.5 && y <= 0.5) geographic_factor = -1.0;

        return std::max(0.0,
                        s(dip,gametes,mutations,
        [&geographic_factor](double & w,const KTfwd::popgenmut & m) {
            w *= (1.0 + geographic_factor*2.0*m.s);
        },
        [&geographic_factor](double & w,const KTfwd::popgenmut & m) {
            w *= (1.0 + geographic_factor*m.h*m.s);
        },
        1.0));
    }
};

int main(int argc, char ** argv)
{
    if(argc!=11)
    {
        std::cerr << "Incorrect number of arguments.\n"
                  << "Usage:\n"
                  << argv[0] << ' '
                  << "N "
                  << "theta "
                  << "rho "
                  << "s "
                  << "h "
                  << "mutrate_to_selected "
                  << "radius "
                  << "dispersal "
                  << "seed "
                  << "format\n"
                  << "\n"
                  << "Note: format = 0 means list of diploids + selected mutations\n"
                  << "format = nsam > 0  = ms-style output of nsam diploids + their geographic locations\n"
				  << "format = N = output info for whole population\n"
				  << "format > N = bad bad bad\n";
        exit(0);
    }
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
    const unsigned format = atoi(argv[argn++]);

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
    //The geometry is a square (0,0) to (1,1).
    //We assign 1/2 of diploids to upper left,
    //and 1/2 to lower right of the landscape initially.
    rtree_type rtree;
    for(std::size_t i=0; i<N; ++i)
    {
        double x=0.0,y=0.0;
        if(i<N/2)
        {
            //Yes, I am aware that this does not truly
            //sample uniformly in 2d...
            x = gsl_ran_flat(rng.get(),0.,0.5);
            y = gsl_ran_flat(rng.get(),0.5,1.);
        }
        else
        {
            x = gsl_ran_flat(rng.get(),0.5,1);
            y = gsl_ran_flat(rng.get(),0.,0.5);
        }
        pop.diploids[i].v = landscape::csdiploid::value(std::make_pair(landscape::csdiploid::point(x,y),i));
        rtree.insert(pop.diploids[i].v);
    }

    //pre-allocate space for a good guess as to the total # mutations
    //expected at equilibrium.
    pop.mutations.reserve(size_t(std::ceil(std::log(2*N)*theta+0.667*theta)));

    /* Our rules type (defined in wfrules.hpp),
     * gets initialized with the rtree from above,
     * the "mating radius" and the "dispersal radius"
     */
    rules_type rules(std::move(rtree),radius,dispersal);

    /* Now, we define our recombination,
     * fitness, and mutation models.
     *
     * Recombination is uniform on the 1/2-open interval [0,1),
     * which is what the 0 and 1 are in the call below:
     */
    auto recombination_model=std::bind(KTfwd::poisson_xover(),rng.get(),littler,0.,1.,
                                       std::placeholders::_1,std::placeholders::_2,std::placeholders::_3);
    /* Fitness is multiplicative, but with s treated as -s in a square
     * bounded by (0,0) to (0.5,0.5)
     */
    auto fitness_model = std::bind(spatial_fitness(),std::placeholders::_1,std::placeholders::_2,
                                   std::placeholders::_3);
    //We're going to initialized our generation here...
    unsigned generation=0;

    /* The mutation model is infinitely-many sites.
     * This code uses a function from the "sugar"
     * part of fwdpp.  The mutation type is auto-detected
     * by the compiler to be KTfwd::popgenmut.
     *
     * We'll cover the details as we move through creating
     * this object.
     *
     * Note that the expressions below are easily modified
     * to allow distributions on s,h, etc., and it is
     * not hard to have variation in mutation and recombination rates,
     * but that is just beyond the scope here.
     */

    //Mutation positions are uniform on (0,1]
    auto mutation_positions = [&rng] { return gsl_rng_uniform(rng.get()); };
    //Every mutation has the same selection coefficient...
    auto selection_coefficients = [&s] { return s; };
    //...and the same dominance
    auto dominance = [&h] {return h;};
    auto mutation_model = std::bind(KTfwd::infsites(),
                                    std::placeholders::_1,
                                    std::placeholders::_2,
                                    rng.get(),
                                    std::ref(pop.mut_lookup),
                                    &generation,//this pointer to generation ensures that origin time of each mutation is recored
                                    mu_n,
                                    mu,
                                    mutation_positions,
                                    selection_coefficients,
                                    dominance);
    /* This is the main workhorse.
     * We'll evolve for 10N generations.
     * Each generation, a call to fwdpp's experimental
     * version of sample_diploid applies the rules class
     * for calculating fitnesses, picking parents, and updating
     * offspring information.
     *
     * After each iteration, we remove fixations for efficiency.
     * Here, a fixation would mean a variant fixed in the entire
     * "pop" object.  Here, all neutral & selected variants get
     * removed, but that is optional--fwdpp has a few variants of that
     * function
     */
    for( ; generation < 10*N ; ++generation )
    {
        double wbar = KTfwd::experimental::sample_diploid(rng.get(),
                      pop.gametes,
                      pop.diploids,
                      pop.mutations,
                      pop.mcounts,
                      N, //Population size will be constant.  If changing, we also pass in the next size
                      mu_n+mu, //TOTAL mutation rate = neutral + selected mutation rates
                      mutation_model,
                      recombination_model,
                      fitness_model,
                      pop.neutral,pop.selected,
                      0,//This is probability(selfing). It defaults to 0, but we need to pass it
                      //so that we can pass our "rules" along on the next line
                      rules);
        //Take any fixed variants, transfer them out of population and into fixation time containers
        KTfwd::update_mutations(pop.mutations,pop.fixations,pop.fixation_times,pop.mut_lookup,pop.mcounts,generation,2*N);
    }
    if(!format)
    {
        //At this point, we would do some analysis...
        //Here, we'll print out each diploid, and
        //the position + s for each mutation on each chromosome,
        //plus its coordinate.  Output will be "tidy",
        //e.g. ready for dplyr.
        std::cout << "dip x y chrom pos s\n";
        for(std::size_t i=0; i<pop.diploids.size(); ++i)
        {
            auto x = pop.diploids[i].v.first.get<0>();
            auto y = pop.diploids[i].v.first.get<1>();
            if(pop.gametes[pop.diploids[i].first].smutations.empty())
            {
                std::cout << i << ' ' << x << ' ' << y << " 0 " <<"NA NA" << '\n';
            }
            else
            {
                for(const auto & m : pop.gametes[pop.diploids[i].first].smutations)
                {
                    std::cout << i << ' ' << x << ' ' << y << " 0 "
                              << pop.mutations[m].pos << ' '
                              << pop.mutations[m].s << '\n';
                }
            }
            if(pop.gametes[pop.diploids[i].second].smutations.empty())
            {
                std::cout << i << ' ' << x << ' ' << y << " 1 " << "NA NA" << '\n';
            }
            else
            {
                for(const auto & m : pop.gametes[pop.diploids[i].second].smutations)
                {
                    std::cout << i << ' ' << x << ' ' << y << " 1 "
                              << pop.mutations[m].pos << ' '
                              << pop.mutations[m].s << '\n';
                }
            }
        }
    }
    else
    {
        /* Sample "format" random diploids.  We want to get their
         * geographic info, so we'll randomly choose individuals
         * w/o replacement, and use fwdpp to get an "ms" block
         * from the population
         */
        std::vector<unsigned> diploids2sample;
        if( format < N )
        {
            for(unsigned i=0; i<format; ++i)
            {
                auto ind = unsigned(gsl_ran_flat(rng.get(),0.0,double(N)));
                while(find(diploids2sample.begin(),diploids2sample.end(),ind)!=diploids2sample.end())
                {
                    ind = unsigned(gsl_ran_flat(rng.get(),0.0,double(N)));
                }
                diploids2sample.push_back(ind);
            }
        } else
        {
            diploids2sample.resize(N);
            unsigned i=0;
            std::generate(diploids2sample.begin(),diploids2sample.end(),[&i] {return i++;});
        }

        auto popsample = KTfwd::sample_separate(pop,diploids2sample,true);//true = do not include variants fixed in the sample
        /*
         * Print out geographic info for each individual,
         * then neutral genotypes, then selected genotypes
         */
        for(auto & d : diploids2sample)
        {
            std::cout << pop.diploids[d].v.first.get<0>() << ' ' << pop.diploids[d].v.first.get<1>() << '\n';
        }
        //Now we need libsequence
        Sequence::SimData neutral(popsample.first.begin(),popsample.first.end()),
                 selected(popsample.second.begin(),popsample.second.end());
        std::cout << neutral << '\n' << selected << '\n';
    }
}
