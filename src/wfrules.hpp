#ifndef WFRULES_HPP
#define WFRULES_HPP
#include <vector>
#include <cmath>
#include <fwdpp/internal/gsl_discrete.hpp>
#include <fwdpp/type_traits.hpp>
#include <boost/geometry/index/rtree.hpp>

namespace landscape
{
template<typename rtree_type>
struct WFLandscapeRules
{
    double wbar,radius,dispersal;
    std::size_t dipindex;
    std::vector<double> fitnesses,fitnesses_temp;
    KTfwd::fwdpp_internal::gsl_ran_discrete_t_ptr lookup,lookup2;
    rtree_type parental_rtree,offspring_rtree;
    //! \brief Constructor
    WFLandscapeRules(rtree_type && r,double radius_,double dispersal_) :
        wbar(0.),radius(radius_),dispersal(dispersal_),dipindex(0),
        fitnesses(std::vector<double>()),
        lookup(KTfwd::fwdpp_internal::gsl_ran_discrete_t_ptr(nullptr)),
        lookup2(KTfwd::fwdpp_internal::gsl_ran_discrete_t_ptr(nullptr)),
        parental_rtree(rtree_type()),
        offspring_rtree(std::move(r))
    {
    }

    //Get fitnesses for each diploid, tally current mean fitness.
    //Create fast lookup table for individuals based on fitness
    template<typename dipcont_t,
             typename gcont_t,
             typename mcont_t,
             typename fitness_func>
    void w(const dipcont_t & diploids,
           gcont_t & gametes,
           const mcont_t & mutations,
           const fitness_func & ff)
    {
        parental_rtree = std::move(offspring_rtree);
        offspring_rtree=rtree_type();
        dipindex=0;
#ifndef NEBUG
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
        if(fitnesses.size() < N_curr) fitnesses.resize(N_curr);
        for(std::size_t i = 0 ; i < diploids.size() ; ++i)
        {
            gametes[diploids[i].first].n=gametes[diploids[i].second].n=0; //set gamete counts to zero!!!!!
            fitnesses[i]=ff(diploids[i],gametes,mutations);
            wbar+=fitnesses[i];
        }
        lookup = KTfwd::fwdpp_internal::gsl_ran_discrete_t_ptr(gsl_ran_discrete_preproc(N_curr,&fitnesses[0]));
    }

    //Pick parent 1 according to fitness
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
    inline size_t pick2(const gsl_rng * r, const size_t & p1, const double & f,
                        diploid_t & parent1, const gcont_t &, const mcont_t &)
    {
        using value_t = typename diploid_t::value;
        std::vector<value_t> possible_mates;
        parental_rtree.query(boost::geometry::index::satisfies([&parent1,this](const value_t & v) {
            double p1x=boost::geometry::get<0>(parent1.v.first);
            double p1y=boost::geometry::get<1>(parent1.v.first);
            double p2x=boost::geometry::get<0>(v.first);
            double p2y=boost::geometry::get<1>(v.first);
            return std::fabs(p1x-p2x)<=this->radius && std::fabs(p1y-p2y)<=this->radius;
        }),
        std::back_inserter(possible_mates));
        if(possible_mates.size()==1) return p1; //only possible mate was itself

		//build lookup table of possible mates.
		//selfing still allowed...
        if(fitnesses_temp.size() < possible_mates.size()) fitnesses_temp.resize(possible_mates.size());

        for(std::size_t i=0; i<possible_mates.size(); ++i)
        {
            fitnesses_temp[i]=fitnesses[possible_mates[i].second];
        }
        lookup2 = KTfwd::fwdpp_internal::gsl_ran_discrete_t_ptr(gsl_ran_discrete_preproc(possible_mates.size(),fitnesses_temp.data()));
        return gsl_ran_discrete(r,lookup2.get());
    }

    //! \brief Update some property of the offspring based on properties of the parents
    template<typename diploid_t,typename gcont_t,typename mcont_t>
    void update(const gsl_rng * r, diploid_t & offspring,const diploid_t & parent1,
                const diploid_t & parent2,
                const gcont_t &,
                const mcont_t &)
    {
        double x = (boost::geometry::get<0>(parent1.v.first)+boost::geometry::get<0>(parent2.v.first))/2.0 + gsl_ran_gaussian(r,dispersal);
        //double x = (parent1.v.get<0>()+parent2.v.get<0>())/2.0 + gsl_ran_gaussian(r,dispersal);
        if (x<0.)x=0.;
        if (x<1.)x=1.;
        double y = (boost::geometry::get<1>(parent1.v.first)+boost::geometry::get<1>(parent2.v.first))/2.0 + gsl_ran_gaussian(r,dispersal);
        if (y<0.)y=0.;
        if (y>1.)y=1.;
        offspring.v = typename diploid_t::value(std::make_pair(typename diploid_t::point(x,y),dipindex++));
        offspring_rtree.insert(offspring.v);
    }

};
}
#endif
