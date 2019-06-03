#include "catch.hpp"

#include "mesh.h"
#include "junction.h"
#include "reservoir.h"
#include "pipe.h"

// Check mesh class
TEST_CASE("Mesh is checked", "[Mesh]") {
    // Tolerance
    const double tolerance = 1.e-12;

    // Mesh index
    const unsigned meshid = 101;

    // Creat a mesh
    auto mesh = std::make_unique<pipenetwork::Mesh>(meshid);

    std::vector<Index> junction_ids{1,2,3,4,5};
    std::vector<double> elevations{5,4,3,2,1};
    std::vector<double> demands{1,2,3,4,5};
    std::vector<double> leak_diameters{0.1,0.4,0.3,0.2,0.1};

    mesh->create_junctions (junction_ids,elevations,demands,leak_diameters);

    std::vector<Index> res_ids{6,7};
    std::vector<double> heads{99,100};

    mesh->create_reservoirs (res_ids,heads);

    std::vector<Index> pipe_ids{1,2,3};
    std::vector<std::pair<Index, Index>> nodeids{std::make_pair (1,2),
                                                 std::make_pair (2,6),
                                                 std::make_pair (4,7)};
    const std::vector<double> length{100,200,300};
    const std::vector<double> diameter{3,4,5};
    const std::vector<double> roughness{0.2,.6,.9};
    const std::vector<Pipe_status> status{open,open,open};

    mesh->create_pipes (pipe_ids,nodeids,length,diameter,roughness,status);

    mesh->print_summary ();

//
//    //! Create Pipe pointers
//    //! \param[in] nodeids pair of end node ids for the pipe
//    //! \param[in] length, length of the pipe
//    //! \param[in] diameter, diameter of the pipe
//    //! \param[in] roughness, roughness of the pipe
//    //! \param[in] status, status of the pipe (open or close)
//    bool create_pipes(const std::vector<Index>& ids,
//                      const std::vector<std::pair<Index, Index>>& nodeids,
//                      const std::vector<double>& length,
//                      const std::vector<double>& diameter,
//                      const std::vector<double>& roughness,
//                      const std::vector<Pipe_status>& status);




}
