dumux_add_test(SOURCES test_0d1d_intersection.cc LABELS unit)
dumux_add_test(SOURCES test_0d3d_intersection.cc LABELS unit)
dumux_add_test(SOURCES test_1d1d_intersection.cc LABELS unit)
dumux_add_test(SOURCES test_1d3d_intersection.cc LABELS unit)
dumux_add_test(SOURCES test_1d2d_intersection.cc LABELS unit)
dumux_add_test(SOURCES test_2d2d_intersection.cc LABELS unit)
dumux_add_test(SOURCES test_2d3d_intersection.cc LABELS unit)
dumux_add_test(SOURCES test_graham_convex_hull.cc LABELS unit)

dune_symlink_to_source_files(FILES ball.msh)
dumux_add_test(SOURCES test_intersectionentityset.cc
               CMAKE_GUARD dune-alugrid_FOUND
               LABELS unit)

dumux_add_test(SOURCES test_makegeometry.cc LABELS unit)
dumux_add_test(SOURCES test_refinementquadraturerule.cc LABELS unit)
