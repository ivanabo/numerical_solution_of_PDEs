#ifndef _OPERATOR_HH_
#define _OPERATOR_HH_

#include<dune/geometry/quadraturerules.hh>
#include<dune/geometry/referenceelements.hh>

#include<dune/pdelab/localoperator/defaultimp.hh>
#include<dune/pdelab/localoperator/flags.hh>
#include<dune/pdelab/localoperator/pattern.hh>

/** Lokalni operator za zadaću :
 *
 *   - div(K \nabla p) + a*p = f   in \Omega
 *                  p = g   on \Gamma_D\subseteq\partial\Omega
 *  -K \nabla p \cdot n = j   on \Gamma_N = \partial\Omega\setminus\Gamma_D
 *
 * with conforming finite elements on all types of grids in any dimension
 *
 * \tparam BCType parameter class indicating the type of boundary condition
 */

template<class BCType, class KOEF>
class DiffusionLocalOperator : // derivacijska lista -- jakobijan i pattern računa PDELab
	public Dune::PDELab::NumericalJacobianApplyVolume  <DiffusionLocalOperator<BCType,KOEF> >,
	public Dune::PDELab::NumericalJacobianVolume       <DiffusionLocalOperator<BCType,KOEF> >,
	public Dune::PDELab::NumericalJacobianApplyBoundary<DiffusionLocalOperator<BCType,KOEF> >,
	public Dune::PDELab::NumericalJacobianBoundary     <DiffusionLocalOperator<BCType,KOEF> >,
	public Dune::PDELab::FullVolumePattern,
	public Dune::PDELab::LocalOperatorDefaultFlags
{
public:
	// Zastavice koje signaliziraju da na svakom elementu treba zvati: 
	enum { doPatternVolume = true };  // metodu za računanje patterna (iz volumnih doprinosa)
	enum { doAlphaVolume = true };    // alpha_volume
	enum { doAlphaBoundary = true };  // alpha_boundary         

	DiffusionLocalOperator(const BCType& bctype_, // boundary cond.type
								  const KOEF& koef_,
	                       unsigned int intorder_=2) :
		bctype( bctype_ ), koef( koef_ ), intorder( intorder_ )
	{}

	// Računanje volumnog integrala
	// eg   = element (geometry)
	// lfsu = lokalni prostor funkcija za rješenje
	// lfsv = lokalni prostor funkcija za test funkciju
	// x    = vektor koeficijenata rješenja 
	// r    = lokalni rezidual
	template<typename EG, typename LFSU, typename X, typename LFSV, typename R>
	void alpha_volume (const EG& eg, const LFSU& lfsu, const X& x, const LFSV& lfsv, R& r) const
	{
		// dimenzije
		const int dim = EG::Geometry::dimension;
		const int dimw = EG::Geometry::dimensionworld;

		// tipovi
		typedef typename LFSU::Traits::FiniteElementType::Traits::LocalBasisType::Traits::DomainFieldType DF;
		typedef typename LFSU::Traits::FiniteElementType::Traits::LocalBasisType::Traits::RangeFieldType RF;
		typedef typename LFSU::Traits::FiniteElementType::Traits::LocalBasisType::Traits::JacobianType Jacobian;
		typedef typename LFSU::Traits::FiniteElementType::Traits::LocalBasisType::Traits::RangeType Range;
		typedef Dune::FieldVector<RF,dimw> Gradient;
		typedef typename LFSU::Traits::SizeType size_type;

		// integracijska formula
		auto gt = eg.geometry().type();
		auto& rule = Dune::QuadratureRules<DF,dim>::rule(gt,intorder);

		// petlja po svim integracijskim točkama
		for (auto it=rule.begin(); it!=rule.end(); ++it)
		{
			// računanje baznih funckcija na referentnom elementu
			std::vector<Range> phi(lfsu.size());
			lfsu.finiteElement().localBasis().evaluateFunction(it->position(),phi);

			// rješenje u integracijskoj točki
			RF p=0.0;
			for (size_type i=0; i<lfsu.size(); ++i) p += x(lfsu,i)*phi[i];

			// gradijent baznih funkcija
			std::vector<Jacobian> js(lfsu.size());
			lfsu.finiteElement().localBasis().evaluateJacobian(it->position(),js);

			// transformacija gradijenata s referentnog na fizički element
			const Dune::FieldMatrix<DF,dimw,dim> &jac = eg.geometry().jacobianInverseTransposed(it->position());
			std::vector<Gradient> gradphi(lfsu.size());
			for (size_type i=0; i<lfsu.size(); i++)
				jac.mv(js[i][0],gradphi[i]);

			//gradijent rješenja u integracijskoj točki
			Gradient gradp(0.0);
			for (size_type i=0; i<lfsu.size(); ++i)
				gradp.axpy(x(lfsu,i),gradphi[i]);

			RF f = 0.0;
			RF a = 0.0;
			RF k = koef.evaluate(eg.entity(), it->position());

			// integriramo :  -K * grad p * grad phi_i + a*p*phi_i - f phi_i
			RF factor = it->weight()*eg.geometry().integrationElement(it->position());

			for (size_type i=0; i<lfsu.size(); ++i)
				r.accumulate(lfsu, i, (-k*(gradp*gradphi[i]) + a*p*phi[i] - f*phi[i]) * factor);
			}
		}

	// integral po rubu
	// ig     = intersection (= stranica elementa)
	// lfsu_s = lokalni prostor funkcija na stranici za rješenje
	// lfsu_v = lokalni prostor funkcija na stranici za test funkciju 
	// x_s    = vektor koeficijenata rješenja (na stranici)
	// r_s    = rezidual (na stranici)
	template<typename IG, typename LFSU, typename X, typename LFSV, typename R>
	void alpha_boundary (const IG& ig, const LFSU& lfsu_s, const X& x_s,
	                     const LFSV& lfsv_s, R& r_s) const
	{
		// tipovi
		typedef typename LFSU::Traits::FiniteElementType::Traits::LocalBasisType::Traits::DomainFieldType DF;
		typedef typename LFSU::Traits::FiniteElementType::Traits::LocalBasisType::Traits::RangeFieldType RF;
		typedef typename LFSU::Traits::FiniteElementType::Traits::LocalBasisType::Traits::RangeType Range;
		typedef typename LFSU::Traits::SizeType size_type;

		// dimenzije
		const int dim = IG::dimension;

		// integracijska formula na stranici
		auto gtface = ig.geometryInInside().type();
		const auto& rule = Dune::QuadratureRules<DF,dim-1>::rule(gtface,intorder);

		// petlja po svim integracijskim točkama
		for (auto it=rule.begin(); it!=rule.end(); ++it)
		{
			// Ako smo na Dirichletovoj granici preskačemo petlju
			if ( bctype.isDirichlet( ig, it->position() ) )
				continue;

			// pozicija int. točke u lokalnim koordinatam elementa
			auto local = ig.geometryInInside().global(it->position());

			// izračunaj bazne funkcije u integracijskoj tički
			std::vector<Range> phi(lfsu_s.size());
			lfsu_s.finiteElement().localBasis().evaluateFunction(local,phi);

			// rješenje u integracijskoj točki
			RF p=0.0;
			for (size_type i=0; i<lfsu_s.size(); ++i)
				p += x_s(lfsu_s,i)*phi[i];

			// računanje Neumannovog rubnog uvjeta
			Dune::FieldVector<RF,dim> globalpos = ig.geometry().global(it->position());
			RF j = 0.0; // svugdje gdje je Neumann je nula

			// integracija
			RF factor = it->weight()*ig.geometry().integrationElement(it->position());

			for (size_type i=0; i<lfsu_s.size(); ++i)
				r_s.accumulate(lfsu_s,i, j*phi[i]*factor);
		}
	}

private:
	const BCType& bctype;
	const KOEF& koef;
	unsigned int intorder;
};


#endif
