#include "config.h"
#include <iostream>
#ifdef HAVE_UG
#include <iomanip>
#include <dune/grid/utility/gridtype.hh>
#include <dune/grid/common/gridinfo.hh>
#include <dune/grid/io/file/dgfparser/dgfparser.hh>
#include <dune/grid/io/file/dgfparser/dgfug.hh>
#include <dune/grid/io/file/dgfparser/dgfs.hh>
#include <dune/grid/io/file/dgfparser/dgfalu.hh>
#include <dune/grid/io/file/dgfparser/dgfalberta.hh>
#include <dune/grid/io/file/dgfparser/dgfyasp.hh>
#include <dune/grid/io/file/vtk/vtkwriter.hh>
#include <dune/istl/io.hh>
#include <dune/common/timer.hh>
#include "lenswithelementid.hh"
#include "dumux/twophase/fv/boxpwsn.hh"
#include "dumux/timedisc/timeloop.hh"
#include "dumux/material/vangenuchtenlaw_deprecated.hh"

int main(int argc, char** argv) 
{
  try{
    // define the problem dimensions  
    const int dim=2;
    typedef double NumberType; 
    Dune::FieldVector<NumberType, dim> outerLowerLeft(0);
    Dune::FieldVector<NumberType, dim> outerUpperRight(500);
    outerUpperRight[1] = 1000;
    Dune::FieldVector<NumberType, dim> innerLowerLeft(1);
    innerLowerLeft[1] = 2;
    Dune::FieldVector<NumberType, dim> innerUpperRight(4);
    innerUpperRight[1] = 3;
    if (argc != 4) {
      std::cout << "usage: test_elementid basefilename tEnd dt" << std::endl;
      return 0;
    }
    	std::string arg1(argv[2]);
	std::istringstream is1(arg1);
	double tEnd;
	is1 >> tEnd;
	std::string arg2(argv[3]);
	std::istringstream is2(arg2);
	double dt;
	is2 >> dt;
     


    // create a grid object
    //typedef Dune::SGrid<dim,dim> GridType; 
    //typedef Dune::ALUSimplexGrid<dim,dim> GridType; 
    //typedef Dune::AlbertaGrid<dim,dim> GridType; 
    //typedef Dune::YaspGrid<dim,dim> GridType; 
//    typedef Dune::UGGrid<dim> GridType; 
//      typedef Dune::ALUCubeGrid<dim,dim> GridType;
	  typedef Dune::ALUSimplexGrid<dim,dim> GridType; 
    // create grid pointer, GridType is defined by gridtype.hh
    Dune::GridPtr<GridType> gridPtr( argv[1] );

    // grid reference 
    GridType& grid = *gridPtr;

// hussam   grid.globalRefine(1); //3 

    Dune::gridinfo(grid);
    
    Air air;
    Water water;
    Dune::VanGenuchtenLaw law(water, air);
    Dune::LensWithElementID<GridType, NumberType> problem(gridPtr, law, outerLowerLeft, outerUpperRight, innerLowerLeft, innerUpperRight);

    typedef Dune::BoxPwSn<GridType, NumberType> TwoPhase;
    TwoPhase twoPhase(grid, problem);
    
    Dune::TimeLoop<GridType, TwoPhase> timeloop(0, tEnd, dt, "layers", 1);
    
    Dune::Timer timer;
    timer.reset();
    timeloop.execute(twoPhase);
    std::cout << "timeloop.execute took " << timer.elapsed() << " seconds" << std::endl;
     
    //printvector(std::cout, *twoPhase.u, "u", "row", 2, 1, 3);

    return 0;
  }
  catch (Dune::Exception &e){
    std::cerr << "Dune reported error: " << e << std::endl;
  }
  catch (...){
    std::cerr << "Unknown exception thrown!" << std::endl;
  }
}
#else 

int main (int argc , char **argv) try
{
  std::cout << "Please install the UG library." << std::endl;

  return 1;
}
catch (...) 
{
    std::cerr << "Generic exception!" << std::endl;
    return 2;
}
#endif 
