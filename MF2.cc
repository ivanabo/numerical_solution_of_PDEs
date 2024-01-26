#ifdef HAVE_CONFIG_H
#include "config.h"
#endif

#include <cmath>
#include <iostream>
#include <vector>
#include <map>
#include <string>
#include <array>
#include <bitset>

#include <dune/common/fvector.hh>
#include <dune/common/exceptions.hh>
#include <dune/common/parallel/mpihelper.hh>

#include <dune/grid/io/file/vtk/subsamplingvtkwriter.hh>
#include <dune/grid/io/file/gmshreader.hh>
#include <dune/grid/yaspgrid.hh>

#include <dune/common/parametertree.hh>
#include <dune/common/parametertreeparser.hh>

#if HAVE_ALBERTA
#include <dune/grid/albertagrid.hh>
#include <dune/grid/albertagrid/dgfparser.hh>
#endif
#if HAVE_ALUGRID
#include <dune/grid/alugrid.hh>
#include <dune/grid/io/file/dgfparser/dgfalu.hh>
#include <dune/grid/io/file/dgfparser/dgfparser.hh>
#endif
#if HAVE_UG
#include <dune/grid/uggrid.hh>
#endif

#include "slvel.hh"
#include "combo_driver.hh"

//===============================================================
// Main program with grid setup
//===============================================================

int main(int argc, char** argv)
{
	//Maybe initialize Mpi
	Dune::MPIHelper& helper = Dune::MPIHelper::instance(argc, argv);

	if (argc!=4)
	{
		if(helper.rank()==0)
			std::cout << "usage: ./MF2 <dt> <tend> <variance>" << std::endl << "Ostali parametri se mijenjaju u datoteci \" MF2.input \"" << std::endl;
		return 1;
	}

	// Pročitaj ulaznu datoteku
	Dune::ParameterTree input_data;

	// Ime ulazne datoteke je ime_programa.input
	std::string filename("MF2.input");
    
	std::cout << "Reading input file \"" << filename << "\"" << std::endl;

	try
	{   // čitamo ulaznu datoteku
		Dune::ParameterTreeParser::readINITree (filename, input_data);
	}
	catch (...)
	{
		std::cerr << "The configuration file \"" << filename << "\" "
						 "could not be read. Exiting..." << std::endl;
		std::exit(1);
	}
	// čitamo podatke koji nam trebaju 
	int level = input_data.get<int>("level"); // razina profinjenja 
	double xLength = input_data.get<double>("xLength"); 
	double yLength = input_data.get<double>("yLength"); 
	int nox = input_data.get<int>("nox"); 
	int noy = input_data.get<int>("noy");
	double mean = input_data.get<double>("mean"); 
	double correlation_length = input_data.get<double>("correlation_length");

	double dt;
	sscanf(argv[1],"%lg",&dt);

	double tend;
	sscanf(argv[2],"%lg",&tend);

	double variance;
	sscanf(argv[3],"%lg",&variance);

	constexpr int dim = 2;

	// sekvencijalna verzija -- kreiraj Grid
   Dune::FieldVector<double,dim> L(1.0);
	L[0] = xLength;
	L[1] = yLength;
   std::array<int,dim>           N{nox,noy};
   std::bitset<dim>              periodic(false);
	int overlap = 0;
   Dune::YaspGrid<dim> grid(L,N,periodic,overlap);

	grid.globalRefine(level);

	// referenca na GridView
	const auto& gv=grid.leafGridView();
	typedef Dune::YaspGrid<dim>::LeafGridView GV;

	Dune::FieldVector<double,GV::Grid::dimension> cl;
	cl = correlation_length;
	koef_sl<GV,double> koef(gv,cl,variance,mean,5000,-1083);
   
	// zovemo driver
	driver(gv,dt,tend,koef);

    return 0;
}
