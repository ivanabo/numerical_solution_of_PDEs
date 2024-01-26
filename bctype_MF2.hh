#ifndef _BCTYPE_MF2_HH_
#define _BCTYPE_MF2_HH_

#include <dune/common/fvector.hh>
#include <dune/pdelab/common/function.hh>
#include <dune/pdelab/constraints/common/constraintsparameters.hh>

/***************************************************************
** Klasa koja određuje je li neka točka na Dirichletovoj granici 
** u zadaći za tlak p
****************************************************************/

class DirichletBdry : public Dune::PDELab::DirichletConstraintsParameters
{
public:
	//  intersection = stranica elementa (u 3D) ili brid elementa (u 2D)
	//  coord        = lokalne koordinate točke na "intersectionu" koja se ispituje
	//  povratna vrijednost: true ako je točka na Dirichletovoj granici
	//                       false ako nije. 
	template<typename I>
	bool isDirichlet(const I & intersection,
	                 const Dune::FieldVector<typename I::ctype, I::dimension-1> & coord
	                ) const
	{
		// Globalne koordinate točke
		Dune::FieldVector<typename I::ctype, I::dimension> xg = intersection.geometry().global( coord );

		if( xg[0]>10.0-1E-6 || xg[0]<1E-6) return true; //  Dirichletov uvjet na lijevoj i desnoj granici

		return false;
	}

};

/***************************************************************
**  Dirichletov rubni uvjet proširen na čitavu domenu; 
**  za zadaću za tlak p.
**  Template parametri:
**     GV = GridView
**     RF = Range Field Type (tip kojim su predstavljeni elementi slike funkcije)
**     
**     Treći parametar u GridFunctionTraits je dimenzija slike funkcije (1 jer su 
**     naše funkcije skalarne). Ta se dimenzija ponavlja u zadnjem parametru 
**     GridFunctionTraits klase.
******************************************************************/

template<typename GV, typename RF>
class BCExtensionP
	: public Dune::PDELab::GridFunctionBase<Dune::PDELab::
	                                        GridFunctionTraits<GV,RF,1,Dune::FieldVector<RF,1> >, 
	                                        BCExtensionP<GV,RF> > {
 // Klasa čuva referencu na GridView objekt.
	const GV& gv;
public:
	typedef Dune::PDELab::GridFunctionTraits<GV,RF,1,Dune::FieldVector<RF,1> > Traits;

	// Konstruktor samo uzima referencu na  GridView objekt. 
	BCExtensionP (const GV& gv_) : gv(gv_) {}

	// Izračunaj Dirichletovu vrijednost na elementu. Ako točka nije na 
	// Dirichletovoj granici, onda funkcija daje proširenje Dirichletovog rubnog
	// uvjeta na čitavu domenu. To je proširenje u osnovi proizvoljno. 
	// e      = element 
	// xlocal = lokalne koordinate točke u kojoj se računa Dirichletova vrijednost
	// y      = izračunata Dirichletova vrijednost
	inline void evaluate (const typename Traits::ElementType& e,
	                      const typename Traits::DomainType& xlocal,
	                      typename Traits::RangeType& y) const
	{
		const int dim = Traits::GridViewType::Grid::dimension;
		typedef typename Traits::GridViewType::Grid::ctype ctype;

		// Pretvori lokalne koordinate u globalne
		Dune::FieldVector<ctype,dim> x = e.geometry().global(xlocal);

		if (x[0]<1E-6)
			y = 1.0; // na lijevoj granici je p=1
		else
			y = 0.0; // na desnoj granici je p=0

		return;
	}

	// Vrati referencu na GridView
	inline const GV& getGridView () {return gv;}
};

/***************************************************************
** Klasa koja određuje je li neka točka na Dirichletovoj granici 
** u zadaći za koncentraciju c
****************************************************************/

class BCTypeParam
	: public Dune::PDELab::DirichletConstraintsParameters 
{
	double time;
public:
	template<typename I>
	bool isDirichlet(const I & intersection,   
	                 const Dune::FieldVector<typename I::ctype, I::dimension-1> & coord
	                 ) const
	{
    
		Dune::FieldVector<typename I::ctype, I::dimension>
		xg = intersection.geometry().global( coord );
    
		if( xg[0]<1E-6 )
			return true; // Dirichlet na lijevoj granici
    
		return false;  // nije Dirichlet na svim ostalim
	}

	//! set time for subsequent evaluation
	void setTime (double t) { time = t; }

};

/***************************************************************
**  Dirichletov rubni uvjet proširen na čitavu domenu; 
**  za zadaću za koncentraciju c.
**  u t=0 daje inicijalni uvjet.
****************************************************************/

template<typename GV, typename RF>
class BCExtensionC
	: public Dune::PDELab::GridFunctionBase<Dune::PDELab::
           GridFunctionTraits<GV,RF,1,Dune::FieldVector<RF,1> >, BCExtensionC<GV,RF> >
{
	const GV& gv;
	RF time;
public:
	typedef Dune::PDELab::GridFunctionTraits<GV,RF,1,Dune::FieldVector<RF,1> > Traits;

	//! construct from grid view
	BCExtensionC (const GV& gv_) : gv(gv_) {}

	//! evaluate extended function on element
	inline void evaluate (const typename Traits::ElementType& e, 
	                      const typename Traits::DomainType& xlocal,
	                      typename Traits::RangeType& y) const
	{  
		const int dim = Traits::GridViewType::Grid::dimension;
		typedef typename Traits::GridViewType::Grid::ctype ctype;
		Dune::FieldVector<ctype,dim> x = e.geometry().global(xlocal);

		if (x[0]<1E-6)
			y = x[1] * (1-x[1]);

		return;
	}
  
	//! get a reference to the grid view
	inline const GV& getGridView () {return gv;}

	//! Postavljanje vremena. Potrebno pozvati prije evaluate() !
	void setTime (double t) {time = t;}
};


#endif
