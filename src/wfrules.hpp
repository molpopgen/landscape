#ifndef WFRULES_HPP
#define WFRULES_HPP
#include <vector>
#include <cmath>
#include <fwdpp/internal/gsl_discrete.hpp>
#include <fwdpp/type_traits.hpp>
#include <boost/geometry/index/rtree.hpp>
// Boost.Range
#include <boost/range.hpp>
// adaptors
#include <boost/range/adaptor/indexed.hpp>
#include <boost/range/adaptor/transformed.hpp>

namespace landscape
{
/* A "rules" class must define the following functions:
 * "w", which is the first function to be called each generation.
 * Typically, "w" calculated mean fitness and generates a lookup table
 * of fitnesses.
 *
 * "pick1" is a function to pick the first parent.
 *
 * "pick2" is a function to pick the second parent.  The first parent
 * is passed in.
 *
 * "update" is the last function to be called, and it allows any extra
 * data in the offspring to be modified.
 *
 * Thus, a W-F generation is:
 * 1. call rules.w()
 * 2. call rules.pick1() and rules.pick2()
 * 3. mutate and recombine gametes from results of pick1 and pick2
 *    to generate offspring.  These details are handled by fwdpp
 *    internally.
 * 4. call rules.update().
 *
 * The rules class is a template.  The template type
 * must be something with the API of a boost::geometry::rtree.
 */
template<typename rtree_type>
struct WFLandscapeRules
{
    //These are data that our rules class
    //will have access to
    double wbar,radius,dispersal;
    std::size_t dipindex;
    std::vector<double> fitnesses,fitnesses_temp;
    //These are smart pointer wrappers around
    //gsl_ran_discrete_t
    KTfwd::fwdpp_internal::gsl_ran_discrete_t_ptr lookup,lookup2;
    rtree_type parental_rtree;
    using rtree_type_value_t = typename rtree_type::value_type;
    using offspring_vec = std::vector<rtree_type_value_t>;
    offspring_vec offspring_locations;

    //"Constructor" function initialized the object.
    //We need an initial rtree, the "mating radius",
    //and the dispersal radius.  The initial rtree
    //gets moved in instead of copied--it will be left
    //in an invalid state in the calling environment.
    WFLandscapeRules(rtree_type && r,double radius_,double dispersal_) :
        wbar(0.),radius(radius_),dispersal(dispersal_),dipindex(0),
        fitnesses(std::vector<double>()),
        lookup(KTfwd::fwdpp_internal::gsl_ran_discrete_t_ptr(nullptr)),
        lookup2(KTfwd::fwdpp_internal::gsl_ran_discrete_t_ptr(nullptr)),
        parental_rtree(std::move(r)), //move construct from the input data
        offspring_locations(std::vector<rtree_type_value_t>())
    {
    }
    //Taken from http://www.boost.org/doc/libs/1_61_0/libs/geometry/doc/html/geometry/spatial_indexes/rtree_examples/range_adaptors.html
    template <typename First, typename Second>
    struct pair_maker
    {
        typedef std::pair<First, Second> result_type;
        template<typename T>
        inline result_type operator()(T const& v) const
        {
			return result_type(v.value().first,v.value().second);
		}
    };

    //Get fitnesses for each diploid, tally current mean fitness.
    //Create fast lookup table for individuals based on fitness
    //Because this fxn is called first, it plays a role as a "setup"
    //function for the above data each generation.
    template<typename dipcont_t,
             typename gcont_t,
             typename mcont_t,
             typename fitness_func>
    void w(const dipcont_t & diploids,
           gcont_t & gametes,
           const mcont_t & mutations,
           const fitness_func & ff)
    {
        if(!offspring_locations.empty())//then we've been through at least 1 generation...
        {
			using point_t = typename dipcont_t::value_type::value::first_type;
            parental_rtree = rtree_type(offspring_locations | boost::adaptors::indexed()
                                        | boost::adaptors::transformed(pair_maker<point_t,std::size_t>()));
			offspring_locations.clear();
        }
        offspring_locations.reserve(diploids.size());
        //set "dipindex to 0.
        dipindex=0;
        //Debug loop.  Will not be executed if compiled
        //with -DNDEBUG.
        //Makes sure KRT hasn't messed things up.
#ifndef NDEBUG
        for(std::size_t i = 0 ; i < diploids.size() ; ++i)
        {
            std::vector<typename dipcont_t::value_type::value> v;
            parental_rtree.query(boost::geometry::index::satisfies([&diploids,i](const typename dipcont_t::value_type::value & vi)
            {
                auto x = boost::geometry::get<0>(diploids[i].v.first);
                auto y = boost::geometry::get<1>(diploids[i].v.first);
                auto x2 = boost::geometry::get<0>(vi.first);
                auto y2 = boost::geometry::get<1>(vi.first);
                return x==x2&&y==y2;
            }),
            std::back_inserter(v));
            assert(!v.empty());
            bool found=false;
            for(auto & vi : v)
            {
                if(vi.second==diploids[i].v.second) found = true;
            }
            assert(found);
        }
#endif
        unsigned N_curr = diploids.size();
        //Rules classes are handy, as we can re-use
        //allocated RAM each generation:
        if(fitnesses.size() < N_curr) fitnesses.resize(N_curr);
        for(std::size_t i = 0 ; i < diploids.size() ; ++i)
        {
            gametes[diploids[i].first].n=gametes[diploids[i].second].n=0; //set gamete counts to zero!!!!!
            fitnesses[i]=ff(diploids[i],gametes,mutations); //calc fitness of i-th diploid
            wbar+=fitnesses[i]; //keep track of mean fitness
        }
        wbar /= double(diploids.size());

        //this lookup table now allows picking a diploid in O(1) time!  Yay.
        lookup = KTfwd::fwdpp_internal::gsl_ran_discrete_t_ptr(gsl_ran_discrete_preproc(N_curr,&fitnesses[0]));
    }

    //Pick parent 1 according to fitness
    //from the ENTIRE landscape
    inline size_t pick1(const gsl_rng * r) const
    {
        return gsl_ran_discrete(r,lookup.get());
    }

    //Pick parent two in a radius centered on parent 1.
    //If the only diploid in the radius is parent 1, then
    //parent 1 is chosen as the second parent, and hence
    //selfing occurs.  Othersise, we choose parent 2
    //based on fitnesses within the radius.
    template<typename diploid_t,typename gcont_t,typename mcont_t>
    inline size_t pick2(const gsl_rng * r, const size_t & p1, const double & ,
                        diploid_t & parent1, const gcont_t &, const mcont_t &)
    {
        using value_t = typename diploid_t::value;
        std::vector<value_t> possible_mates;
		double p1x=boost::geometry::get<0>(parent1.v.first);
		double p1y=boost::geometry::get<1>(parent1.v.first);
		using point = typename diploid_t::point;
		using box = boost::geometry::model::box<point>;
		box region(point(p1x-radius,p1y-radius),point(p1x+radius,p1y+radius));
		parental_rtree.query(boost::geometry::index::covered_by(region) &&
				boost::geometry::index::satisfies([&parent1,this](const value_t & v) {
					return boost::geometry::distance(v.first,parent1.v.first)<=radius;
					}),std::back_inserter(possible_mates));
        //find all individuals in population whose Euclidiean distance
        //from parent1 is <= radius.  The "point" info fill up
        //the possible_mates vector.
		//This is a more "idiomatic" version based on
		//http://stackoverflow.com/questions/22909171/boostgeometry-nearest-neighbors-using-a-circle,
		//which shows that lambda expressions as query objects are just a lot slower...
		//cout-debugging confirms that this is equivalent to the commented-out bit below:
		/*
		parental_rtree.query(boost::geometry::index::satisfies([&parent1,this](const value_t & v) { 
				return boost::geometry::distance(v.first,parent1.v.first) <= radius;
				}),std::back_inserter(possible_mates));
		*/
		/*	
        parental_rtree.query(boost::geometry::index::satisfies([&parent1,this](const value_t & v) {
            double p1x=boost::geometry::get<0>(parent1.v.first);
            double p1y=boost::geometry::get<1>(parent1.v.first);
            double p2x=boost::geometry::get<0>(v.first);
            double p2y=boost::geometry::get<1>(v.first);
            double euclid = std::sqrt(std::pow(p1x-p2x,2.0)+std::pow(p1y-p2y,2.0));
            return euclid <= radius;
        }),
        std::back_inserter(possible_mates));
        */
		if(possible_mates.size()==1) return p1; //only possible mate was itself, so we self-fertilize

        //build lookup table of possible mates.
        //selfing still allowed...
        if(fitnesses_temp.size() < possible_mates.size()) fitnesses_temp.resize(possible_mates.size());
        double sumw=0.0;
        for(std::size_t i = 0 ; i < possible_mates.size() ; ++i)
        {
            fitnesses_temp[i]=fitnesses[possible_mates[i].second];
            sumw += fitnesses_temp[i];
        }
        double uni = gsl_ran_flat(r,0.0,sumw);
        double sum=0.0;
        for(std::size_t i=0; i<possible_mates.size(); ++i)
        {
            sum+=fitnesses_temp[i];
            if(uni < sum)
            {
                return possible_mates[i].second;
            }
        }
        //should never (?) get here...
        return possible_mates.back().second;

        /* This next code uses the GSL lookup idea
         * to pick mates according to fitness.
         * This is slower than the above, taking another
         * O(k) step to preprocess...
         */
        //build another one of these fast fitness lookups
        //and return a value from it.
        //for(std::size_t i=0; i<possible_mates.size(); ++i)
        //{
        //    //Each diploid's point data contains both its
        //    //xy coords AND there it is in diploids,
        //    //and thus in fitnesses.
        //    //(see simtypes.hpp)
        //    fitnesses_temp[i]=fitnesses[possible_mates[i].second];
        //}
        //lookup2 = KTfwd::fwdpp_internal::gsl_ran_discrete_t_ptr(gsl_ran_discrete_preproc(possible_mates.size(),fitnesses_temp.data()));
        //return possible_mates[gsl_ran_discrete(r,lookup2.get())].second;
    }

    //! \brief Update some property of the offspring based on properties of the parents
    template<typename diploid_t,typename gcont_t,typename mcont_t>
    void update(const gsl_rng * r, diploid_t & offspring,const diploid_t & parent1,
                const diploid_t & parent2,
                const gcont_t &,
                const mcont_t &)
    {
        //Get coordinates for offspring, based on midpoint of parents +
        //Gaussian dispersal independently along each axis.
        //Another option for linear dispersal is
        //https://www.gnu.org/software/gsl/manual/html_node/Spherical-Vector-Distributions.html
        double x = (boost::geometry::get<0>(parent1.v.first)+boost::geometry::get<0>(parent2.v.first))/2.0 + gsl_ran_gaussian(r,dispersal);
        if (x<0.)x=0.;
        if (x>1.)x=1.;
        double y = (boost::geometry::get<1>(parent1.v.first)+boost::geometry::get<1>(parent2.v.first))/2.0 + gsl_ran_gaussian(r,dispersal);
        if (y<0.)y=0.;
        if (y>1.)y=1.;
        //"Label" the offspring with its coordinates.
        //b/c fwdpp guarantees filling diploids from 0 to N-1,
        //we use dipindex here to record where this offspring is
        //in the diploids container.
        offspring.v = typename diploid_t::value(std::make_pair(typename diploid_t::point(x,y),dipindex++));
        offspring_locations.push_back(offspring.v);
    }
};
}
#endif
