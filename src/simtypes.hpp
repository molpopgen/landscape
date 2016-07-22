#ifndef LANDSCAPE_SIMTYPES_HPP
#define LANDSCAPE_SIMTYPES_HPP

#include <limits>
#include <utility>
#include <boost/geometry/geometries/point.hpp>
#include <fwdpp/tags/diploid_tags.hpp>

namespace bg = boost::geometry;

using point = bg::model::point<double, 2, bg::cs::cartesian>;
using value = std::pair<point, std::size_t>;

struct diploid : public KTfwd::tags::custom_diploid_t
/*
 * Minimal custom diploid. Inherits fwdpp tag so
 * that it gets "dispatched" properly.
 */
{
    using first_type = std::size_t;
    using second_type = std::size_t;
    first_type first; //first gamete
    second_type second;//second gamete
    value v;           //location in space & index in population
    diploid() noexcept : v(std::make_pair(point(std::numeric_limits<double>::quiet_Nan(),
                                           std::numeric_limits<double>::quiet_Nan()),
                                     std::numeric_limits<std::size_t>::max()))
    {
    }
};

#endif
