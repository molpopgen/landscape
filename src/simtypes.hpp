#ifndef LANDSCAPE_SIMTYPES_HPP
#define LANDSCAPE_SIMTYPES_HPP

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
    std::size_t first; //first gamete
    std::size_t second;//second gamete
    value v;           //location in space & index in population
};

#endif
