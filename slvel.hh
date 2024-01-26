#ifndef DUNE_SLVEL_HH
#define DUNE_SLVEL_HH

#include<math.h>
#include <dune/pdelab/common/function.hh>
#include "permeability_generatorMF2.hh"

template<typename GV, typename RF>
class koef_sl 
	: public Dune::PDELab::GridFunctionBase<Dune::PDELab::
				GridFunctionTraits<GV,RF,1,Dune::FieldVector<RF,1> >,
				koef_sl<GV,RF> >
{
public:
//	typedef RF RFType;
	typedef Dune::PDELab::GridFunctionTraits<GV,RF,1,Dune::FieldVector<RF,1> > Traits;
//	typedef Dune::PDELab::GridFunctionBase<Traits,koef_sl<GV,RF> > BaseT;

	koef_sl (const GV& gv_, 
		 Dune::FieldVector<double,GV::dimension> correlation_length,
		 double variance = 1.0, 
		 double mean = 0.0, 
		 long modes = 1000, 
		 long seed = -1083)
	: gv(gv_), is(gv.indexSet()), perm(is.size(0))
	{
		typedef typename GV::Traits::template Codim<0>::Iterator ElementIterator;
		typedef typename Traits::DomainFieldType DF;
		const int dim = GV::dimension;
		double mink=1E100;
		double maxk=-1E100;


		EberhardPermeabilityGenerator<GV::dimension> field(correlation_length,variance,mean,modes,seed);

		for (ElementIterator it = gv.template begin<0>(); it!=gv.template end<0>(); ++it)
		{
			int id = is.index(*it);
			Dune::GeometryType gt = it->geometry().type();
			Dune::FieldVector<DF,dim> localcenter = Dune::ReferenceElements<DF,dim>::general(gt).position(0,0);
			Dune::FieldVector<DF,dim> globalcenter = it->geometry().global(localcenter);
			perm[id]=field.eval(globalcenter);
			mink = std::min(mink,log10(perm[id]));
			maxk = std::max(maxk,log10(perm[id]));
		}
		std::cout << "log10(mink)=" << mink << " log10(maxk)=" << maxk << std::endl;
	}

	koef_sl ( const GV& gv_, const std::vector<RF>& perm_)
	: gv(gv_), is(gv.indexSet()), perm(perm_)
	{}

	inline RF evaluate (const typename Traits::ElementType& e,
                         const typename Traits::DomainType& x) const
	{
		return perm[is.index(e)];
	}

	inline const typename Traits::GridViewType& getGridView () const
	{
		return gv;
	}

private:
	const GV& gv;
	const typename GV::IndexSet& is;
	std::vector<RF> perm;
};

#endif
