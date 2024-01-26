#ifndef _DRIVER_MF2_HH_
#define _DRIVER_MF2_HH_

#include <dune/istl/bvector.hh>
#include <dune/istl/io.hh>
#include <dune/istl/operators.hh>
#include <dune/istl/preconditioners.hh>
#include <dune/istl/solvers.hh>
#include <dune/istl/superlu.hh>

#include <dune/pdelab/backend/istl/bcrsmatrixbackend.hh>
#include <dune/pdelab/backend/istl.hh>
#include <dune/pdelab/backend/istlsolverbackend.hh>
#include <dune/pdelab/common/function.hh>
#include <dune/pdelab/common/vtkexport.hh>
#include <dune/pdelab/constraints/conforming.hh>
//#include <dune/pdelab/constraints/constraints.hh>
//#include <dune/pdelab/finiteelementmap/conformingconstraints.hh>
#include <dune/pdelab/finiteelementmap/qkfem.hh>
#include <dune/pdelab/gridfunctionspace/genericdatahandle.hh>
#include <dune/pdelab/gridfunctionspace/gridfunctionspace.hh>
#include <dune/pdelab/gridfunctionspace/gridfunctionspaceutilities.hh>
#include <dune/pdelab/gridfunctionspace/interpolate.hh>
#include <dune/pdelab/gridoperator/gridoperator.hh>
#include <dune/pdelab/gridoperator/onestep.hh>
#include <dune/pdelab/stationary/linearproblem.hh>
#include <dune/pdelab/instationary/onestep.hh>

#include <dune/grid/io/file/gmshreader.hh>
#include <dune/grid/io/file/vtk/subsamplingvtkwriter.hh>

#include "bctype_MF2.hh"
#include "operator_MF2.hh"
#include "space_operator.hh"
#include "time_operator.hh"
#include "slvel.hh"

#include <memory>

template<class GV, typename KOEF>
void driver(const GV& gv, double dt, double tend, KOEF& koef)
{
   using namespace Dune::PDELab;  // Da skratimo imena

	// <<<1>>> Tipovi u domeni i slici
	typedef typename GV::Grid::ctype Coord;
	typedef double Real;
//  const int dim = GV::dimension;

	// Početni trenutak je t=0
	Real time = 0.0;

	// <<<2>>> Zajedničko za tlak i koncentraciju
	const int k = 1;
	typedef ConformingDirichletConstraints                           CONSTRAINTS;
	typedef istl::VectorBackend<>                                    VBE;
	typedef istl::BCRSMatrixBackend<>                                MBE;
	MBE mbe(9);  // traži prosječan broj ne-nul elemenata u redu (=9)

	// <<<3>>> Grid function space za tlak p
	typedef QkLocalFiniteElementMap<GV,Coord,Real,k>                 FEMP;
	FEMP femp(gv);
	typedef GridFunctionSpace<GV,FEMP,CONSTRAINTS,VBE>               GFSP;
	GFSP gfsp(gv,femp);

	typedef typename GFSP::template ConstraintsContainer<Real>::Type CCP;
	CCP ccp;
	DirichletBdry bctypeP;
	constraints(bctypeP, gfsp, ccp);
//	std::cout << "constrained dofs=" << ccp.size() << " of " << gfsp.globalSize() << std::endl;

	// <<<4>>> GridOperator za tlak
	typedef DiffusionLocalOperator<DirichletBdry,koef_sl<GV,double>>  LOPP;
	typedef GridOperator<
		GFSP,GFSP,        /* prostor KE rješenja i test funkcije */
		LOPP,            /* lokalni operator */
		MBE,            /* matrix backend */
		Real,Real,Real, /* tipovi u domeni, slici i jakobijanu */
		CCP,CCP           /* ograničenja za prostor rješenja i test funkcija. */
		> GOP;

	LOPP lopp(bctypeP, koef);
	GOP gop(gfsp,ccp,gfsp,ccp,lopp,mbe);

	// <<<5>>> Konstrukcija rješavača; tlak
	typedef typename GOP::Traits::Domain             P;
	typedef BCExtensionP<GV,Real>                    GP;
	typedef ISTLBackend_SEQ_BCGS_SSOR                LSP; // Linear Solver Pressure
	typedef StationaryLinearProblemSolver<GOP,LSP,P> SLP;

	P p(gfsp,0.0);
	GP gp(gv);
	interpolate(gp,gfsp,p);
	LSP lsp(5000,true);        // max 5000 iteracija, verbosity = true
	SLP slp(gop,lsp,p,1e-10);  // redukcija = 1e-10
	slp.apply();

	// <<<6>>> grafički izlaz (VTK); tlak
	typedef DiscreteGridFunction<GFSP,P> DGFP; // Discrete Grid Function Pressure

	DGFP dgfp(gfsp,p);
	Dune::VTKWriter<GV> vtkwriter(gv,Dune::VTK::conforming);
	vtkwriter.addVertexData(std::make_shared<VTKGridFunctionAdapter<DGFP>>(dgfp,"pressureSolution"));
	vtkwriter.write("pressure",Dune::VTK::ascii); //Dune::VTK::appendedraw);

	// <<<7>>> u = DiscreteGridFunctionGradient P
	typedef DiscreteGridFunctionGradient<GFSP,P> DGFGP; // Discrete Grid Function Gradient Pressure
	DGFGP dgfgp(gfsp,p);

	// <<<8>>> Grid function space za koncentraciju c
	typedef QkLocalFiniteElementMap<GV,Coord,Real,k>                 FEMC;
	FEMC femc(gv);
	typedef Dune::PDELab::GridFunctionSpace<GV,FEMC,CONSTRAINTS,VBE> GFSC;
	GFSC gfsc(gv,femc);

	typedef typename GFSC::template ConstraintsContainer<Real>::Type CCC;
	CCC ccc;
	BCTypeParam bctypeC;
	bctypeC.setTime(time); // b.c. ovisi o vremenu       
	constraints(bctypeC, gfsc, ccc);

	// <<<9>>> Nestacionaran prostorni lokalni operator
	typedef StationaryLocalOperator<BCTypeParam,DGFGP,koef_sl<GV,double> > SLOP; 
	SLOP slop(bctypeC, dgfgp, koef, 4);                                           

	// <<<10>>> Nestacionaran vremenski lokalni operator
	typedef TimeLocalOperator TLOP; 
	TLOP tlop(4);                                                 

	// <<<11>>> Nestacionaran prostorni grid operator
	typedef Dune::PDELab::GridOperator<GFSC,GFSC,SLOP,MBE,Real,Real,Real,CCC,CCC> GO0;
	GO0 go0(gfsc, ccc, gfsc, ccc, slop, mbe);

	// <<<12>>> Nestacionaran vremenski grid operator
	typedef Dune::PDELab::GridOperator<GFSC,GFSC,TLOP,MBE,Real,Real,Real,CCC,CCC> GO1;
	GO1 go1(gfsc, ccc, gfsc, ccc, tlop, mbe);  

	// <<<13>>> Izbor tipa vremenske diskretizacije
	typedef Dune::PDELab::OneStepGridOperator<GO0,GO1> IGO;
	IGO igo(go0,go1);                           // novi grid operator

	// <<<14>>> Interpolacija Dirichletovog rubnog te početnog uvjeta; za koncentraciju c                   
	typedef typename IGO::Traits::Domain C;
	C cold(gfsc,0.0);							// rješenje u t=t^n, početni uvjet je 0
	typedef BCExtensionC<GV,Real> GC;	// rubni uvjet
	GC gc(gv);   								// sada ovisi o vremenu
	gc.setTime(time);                                         
	Dune::PDELab::interpolate(gc,gfsc,cold);	// početni uvjet je dan ovdje

	// <<<15>>> Select a linear solver backend
	typedef Dune::PDELab::ISTLBackend_SEQ_BCGS_SSOR LSC;
	LSC lsc(5000,false);

	// <<<16>>> Solver for linear problem per stage
	typedef Dune::PDELab::StationaryLinearProblemSolver<IGO,LSC,C> PDESOLVER;
	PDESOLVER pdesolver(igo,lsc,1e-10);

	// <<<17>>> Ovdje se bira vremenska diskretizacija (implicitni Euler)
	Dune::PDELab::ImplicitEulerParameter<Real> method;	// koeficijenti metode
	Dune::PDELab::OneStepMethod<Real,IGO,PDESOLVER,C,C> osm(method,igo,pdesolver);
	osm.setVerbosityLevel(0);

	// <<<18>>> ispiši inicijalni uvjet
	Dune::PDELab::FilenameHelper fn("first");			// dodaj broj imenu
	{
		typedef Dune::PDELab::DiscreteGridFunction<GFSC,C> DGFC;
		DGFC dgfc(gfsc,cold);
		Dune::SubsamplingVTKWriter<GV> vtkwriter(gv,3);
		vtkwriter.addVertexData(new Dune::PDELab::VTKGridFunctionAdapter<DGFC>(dgfc,"concentration"));
		vtkwriter.write(fn.getName(), Dune::VTK::ascii);
		fn.increment();					// povećaj broj u imenu datoteke
	}
                                    
	// <<<19>>> vremenska petlja
	C cnew(gfsc,0.0);			// sljedeći vremenski sloj -jednokoračna metoda
	while (time < tend-1e-8) 
	{
		// postavi novo vrijeme u BC klasu
		bctypeC.setTime(time+dt);			// izračunaj Dirichletov r.u
      ccc.clear();							// u ovom vremenskom trenutku
      Dune::PDELab::constraints(bctypeC,gfsc,ccc);
      osm.apply(time,dt,cold,gc,cnew);	// riješi sustav
      
      int noIter = osm.getPDESolver().ls_result().iterations;   
      std::cout << "no of linear iterations = " << noIter << std::endl;
      // kontrola vremenskog koraka
      if(noIter < 7) dt *= 2.0;
      if(noIter > 20) dt /= 2.0;
      // graphics
      typedef Dune::PDELab::DiscreteGridFunction<GFSC,C> DGFC;
      DGFC dgfc(gfsc,cnew);
      Dune::SubsamplingVTKWriter<GV> vtkwriter(gv,3);
      vtkwriter.addVertexData(new Dune::PDELab::VTKGridFunctionAdapter<DGFC>(dgfc,"concentration"));
      vtkwriter.write(fn.getName(),Dune::VTK::ascii);
      fn.increment();

      cold = cnew;			// pripremi sljedeći vremenski korak
      time += dt;
std::cout << "Trenutno vrijeme: " << time << std::endl;
    }

}



#endif
