#ifndef LANDSCAPE_SIMTYPES_HPP
#define LANDSCAPE_SIMTYPES_HPP

#include <limits>
#include <utility>
#include <boost/geometry/geometries/point.hpp>
#include <boost/geometry/core/cs.hpp>
#include <fwdpp/tags/diploid_tags.hpp>
#include <fwdpp/sugar/singlepop.hpp>
#include <fwdpp/sugar/popgenmut.hpp>



namespace landscape
{
struct csdiploid : public KTfwd::tags::custom_diploid_t
/*
 * Minimal custom diploid in Cartesian space. Inherits fwdpp tag so
 * that it gets "dispatched" properly.
 */
{
    using point = boost::geometry::model::point<double, 2, boost::geometry::cs::cartesian>;
    using value = std::pair<point, std::size_t>;
    using first_type = std::size_t;
    using second_type = std::size_t;
    first_type first; //first gamete
    second_type second;//second gamete
    value v;           //location in space & index in population
    csdiploid() noexcept : v(std::make_pair(point(std::numeric_limits<double>::quiet_NaN(),
                                                std::numeric_limits<double>::quiet_NaN()),
                                                std::numeric_limits<std::size_t>::max()))
    {
    }
    csdiploid(std::size_t i,std::size_t j) : first(i),second(j),
        v(std::make_pair(point(std::numeric_limits<double>::quiet_NaN(),
                               std::numeric_limits<double>::quiet_NaN()),
                         std::numeric_limits<std::size_t>::max()))
    {
    }
};
/*
 * Our population is a single deme.  The mutation type
 * is KTfwd::popgenmut, from fwdpp/sugar/popgenmut.hpp.
 * This mutation is the "standard" mutation type with
 * position, s, h as its main data members.  The origin
 * time of mutations is also recorded.
 *
 * The diploid type is our csdiploid defined above.
 *
 * For more details, see the definition of a "singlepop"
 * in fwdpp/sugar/singlepop and the files included from there.
 *
 * The fwdpp manual online has detailed info, too.
 */
using poptype = KTfwd::singlepop<KTfwd::popgenmut,csdiploid>;
}
#endif
