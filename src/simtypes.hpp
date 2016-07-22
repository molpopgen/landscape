#ifndef LANDSCAPE_SIMTYPES_HPP
#define LANDSCAPE_SIMTYPES_HPP

#include <boost/geometry/geometries/point.hpp>
#include <fwdpp/tags/diploid_tags.hpp>

namespace bg = boost::geometry;

using point = bg::model::point<double, 2, bg::cs::cartesian>;
using value = std::pair<point, unsigned>;

struct diploid : public KTfwd::tags::custom_diploid_t
/*
 * Minimal custom diploid. Inherits fwdpp tag so
 * that it gets "dispatched" properly.
 */
{
    std::size_t first;
    std::size_t second;
    value v;
};

#endif
