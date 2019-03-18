#include <config.h>

#include <iostream>
#include <algorithm>
#include <initializer_list>

#include <dune/common/exceptions.hh>
#include <dune/common/parallel/mpihelper.hh>
#include <dune/common/fvector.hh>
#include <dune/geometry/multilineargeometry.hh>

#include <dumux/common/geometry/geometryintersection.hh>


#ifndef DOXYGEN
template<int dimworld = 2>
Dune::MultiLinearGeometry<double, 1, dimworld>
makeLine(std::initializer_list<Dune::FieldVector<double, dimworld>>&& c)
{
    return {Dune::GeometryTypes::line, c};
}

template<int dimworld = 2>
bool testIntersection(const Dune::MultiLinearGeometry<double, dimworld, dimworld>& polygon,
                      const Dune::MultiLinearGeometry<double, 1, dimworld>& line,
                      bool foundExpected = true)
{
    using Test = Dumux::GeometryIntersection<Dune::MultiLinearGeometry<double,dimworld,dimworld>,
                                             Dune::MultiLinearGeometry<double,1,dimworld> >;
    typename Test::IntersectionType intersection;
    bool found = Test::intersection(polygon, line, intersection);
    if (!found && foundExpected)
        std::cerr << "Failed detecting intersection with " << line.corner(0) << " " << line.corner(1) << std::endl;
    else if (found && foundExpected)
        std::cout << "Found intersection with " << line.corner(0) << " " << line.corner(1) << std::endl;
    else if (found && !foundExpected)
        std::cerr << "Found false positive: intersection with " << line.corner(0) << " " << line.corner(1) << std::endl;
    else if (!found && !foundExpected)
        std::cout << "No intersection with " << line.corner(0) << " " << line.corner(1) << std::endl;
    return (found == foundExpected);
}
#endif

int main (int argc, char *argv[]) try
{
    // maybe initialize mpi
    Dune::MPIHelper::instance(argc, argv);

    constexpr int dimworld = 2;
    constexpr int dim = 2;

    // we test quadrilateral-line & triangle-line intersections
    using CornerStorage = std::vector<Dune::FieldVector<double, dimworld>>;
    CornerStorage quadCorners({ {0.0, 0.0}, {1.0, 0.0}, {0.0, 1.0}, {1.0, 1.0} });
    CornerStorage triaCorners({ {0.0, 0.0}, {1.0, 0.0}, {0.0, 1.0} });

    using Geometry = Dune::MultiLinearGeometry<double, dim, dimworld>;
    Geometry quad(Dune::GeometryTypes::cube(dimworld), quadCorners);
    Geometry triangle(Dune::GeometryTypes::simplex(dimworld), triaCorners);

    // collect returns to determine exit code
    std::vector<bool> returns;

    // the tests
    returns.push_back(testIntersection(quad, makeLine({ {0.0, 0.0}, {1.0, 0.0} })));
    returns.push_back(testIntersection(quad, makeLine({ {0.0, 0.0}, {0.0, 1.0} })));
    returns.push_back(testIntersection(quad, makeLine({ {0.0, 0.0}, {1.0, 1.0} })));

    returns.push_back(testIntersection(quad, makeLine({ {1.0, 0.0}, {1.0, 1.0} })));
    returns.push_back(testIntersection(quad, makeLine({ {1.0, 1.0}, {0.0, 1.0} })));

    returns.push_back(testIntersection(quad, makeLine({ {0.5, 0.0}, {0.5, 1.0} })));
    returns.push_back(testIntersection(quad, makeLine({ {0.0, 0.5}, {1.0, 0.5} })));
    returns.push_back(testIntersection(quad, makeLine({ {0.5, 0.5}, {0.5, 2.0} })));
    returns.push_back(testIntersection(quad, makeLine({ {0.5, 0.5}, {0.5, -2.0} })));
    returns.push_back(testIntersection(quad, makeLine({ {0.5, 0.5}, {-2.0, 0.5} })));
    returns.push_back(testIntersection(quad, makeLine({ {0.5, 0.5}, {2.0, 0.5} })));

    returns.push_back(testIntersection(quad, makeLine({ {0.5, 0.0}, {0.5, -2.0} }), false));
    returns.push_back(testIntersection(quad, makeLine({ {0.5, 1.0}, {0.5, 2.0} }), false));

    returns.push_back(testIntersection(triangle, makeLine({ {0.0, 0.0}, {1.0, 0.0} })));
    returns.push_back(testIntersection(triangle, makeLine({ {0.0, 0.0}, {0.0, 1.0} })));
    returns.push_back(testIntersection(triangle, makeLine({ {0.0, 0.0}, {1.0, 1.0} })));

    returns.push_back(testIntersection(triangle, makeLine({ {0.5, 0.0}, {0.5, 1.0} })));
    returns.push_back(testIntersection(triangle, makeLine({ {0.0, 0.5}, {1.0, 0.5} })));
    returns.push_back(testIntersection(triangle, makeLine({ {0.5, 0.5}, {0.0, 0.0} })));
    returns.push_back(testIntersection(triangle, makeLine({ {0.0, 0.0}, {0.5, 0.5} })));
    returns.push_back(testIntersection(triangle, makeLine({ {0.5, 0.5}, {-2.0, 0.5} })));
    returns.push_back(testIntersection(triangle, makeLine({ {0.5, 0.5}, {0.5, -2.0} })));

    returns.push_back(testIntersection(triangle, makeLine({ {1.0, 1.0}, {0.0, 1.0} }), false));
    returns.push_back(testIntersection(triangle, makeLine({ {1.0, 0.0}, {1.0, 1.0} }), false));

    // determine the exit code
    if (std::any_of(returns.begin(), returns.end(), [](bool i){ return !i; }))
        return 1;

    std::cout << "All tests passed!" << std::endl;

    return 0;
}
// //////////////////////////////////
//   Error handler
// /////////////////////////////////
catch (const Dune::Exception& e) {
    std::cout << e << std::endl;
    return 1;
}
