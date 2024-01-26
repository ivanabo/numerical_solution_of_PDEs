#include<dune/pdelab/localoperator/idefault.hh>

#include<dune/geometry/quadraturerules.hh>
#include<dune/geometry/referenceelements.hh>

#include<dune/pdelab/localoperator/defaultimp.hh>
#include<dune/pdelab/localoperator/flags.hh>
#include<dune/pdelab/localoperator/pattern.hh>

#include <dune/pdelab/gridfunctionspace/gridfunctionspaceutilities.hh>


/** a local operator for solving the equation
 *
 *   - div (D \nabla c) + u \cdot \nabla c = f   in \Omega
 *                                       c = g   on \Gamma_D\subseteq\partial\Omega
 *                       -\nabla c \cdot n = j   on \Gamma_N = \partial\Omega\setminus\Gamma_D
 *
 * \tparam BCType parameter class indicating the type of boundary condition
 * \tparam DGFG discrete grid function gradient
 */
template<class BCType, class DGFG, class KOEF>
class StationaryLocalOperator :
	public Dune::PDELab::NumericalJacobianApplyVolume<StationaryLocalOperator<BCType,DGFG,KOEF> >,
	public Dune::PDELab::NumericalJacobianVolume<StationaryLocalOperator<BCType,DGFG,KOEF> >,
	public Dune::PDELab::NumericalJacobianApplyBoundary<StationaryLocalOperator<BCType,DGFG,KOEF> >,
	public Dune::PDELab::NumericalJacobianBoundary<StationaryLocalOperator<BCType,DGFG,KOEF> >,
	public Dune::PDELab::FullVolumePattern,
	public Dune::PDELab::LocalOperatorDefaultFlags,
	public Dune::PDELab::InstationaryLocalOperatorDefaultMethods<double>         
{
public:
	// pattern assembly flags
	enum { doPatternVolume = true };

	// residual assembly flags
	enum { doAlphaVolume = true };
	enum { doAlphaBoundary = true };                                // assemble boundary

	StationaryLocalOperator(BCType& bctype_, // boundary cond.type
								  const DGFG& dgfg_,
								  const KOEF& koef_,
	                       unsigned int intorder_=2) :
		bctype( bctype_ ), dgfg( dgfg_ ), koef( koef_ ), intorder( intorder_ )
	{}

	// volume integral depending on test and ansatz functions
	template<typename EG, typename LFSU, typename X, typename LFSV, typename R>
	void alpha_volume (const EG& eg, const LFSU& lfsu, const X& x, const LFSV& lfsv, R& r) const
	{
		// assume Galerkin: lfsu == lfsv
		// This yields more efficient code since the local functionspace only
		// needs to be evaluated once, but would be incorrect for a finite volume
		// method

		// dimensions
		const int dim = EG::Geometry::dimension;
		const int dimw = EG::Geometry::dimensionworld;

		// extract some types
		typedef typename LFSU::Traits::FiniteElementType::
			Traits::LocalBasisType::Traits::DomainFieldType DF;
		typedef typename LFSU::Traits::FiniteElementType::
			Traits::LocalBasisType::Traits::RangeFieldType RF;
		typedef typename LFSU::Traits::FiniteElementType::
			Traits::LocalBasisType::Traits::JacobianType Jacobian;
		typedef typename LFSU::Traits::FiniteElementType::
			Traits::LocalBasisType::Traits::RangeType Range;
		typedef Dune::FieldVector<RF,dimw> Gradient;
		typedef typename LFSU::Traits::SizeType size_type;

		// select quadrature rule
		Dune::GeometryType gt = eg.geometry().type();
		const Dune::QuadratureRule<DF,dim>&
		rule = Dune::QuadratureRules<DF,dim>::rule(gt,intorder);

		// loop over quadrature points
		for (typename Dune::QuadratureRule<DF,dim>::const_iterator
				it=rule.begin(); it!=rule.end(); ++it)
		{
			// evaluate basis functions on reference element
			std::vector<Range> phi(lfsu.size());
			lfsu.finiteElement().localBasis().evaluateFunction(it->position(),phi);

			// compute c at integration point
			RF c = 0.0;
			for (size_type i=0; i<lfsu.size(); ++i)
				c += x(lfsu,i)*phi[i];

			// evaluate gradient of basis functions on reference element
			std::vector<Jacobian> js(lfsu.size());
			lfsu.finiteElement().localBasis().evaluateJacobian(it->position(),js);

			// transform gradients from reference element to real element
			const Dune::FieldMatrix<DF,dimw,dim>
				&jac = eg.geometry().jacobianInverseTransposed(it->position());
			std::vector<Gradient> gradphi(lfsu.size());
			for (size_type i=0; i<lfsu.size(); i++)
				jac.mv(js[i][0],gradphi[i]);

			// compute gradient of c
			Gradient gradc(0.0);
			for (size_type i=0; i<lfsu.size(); ++i)
				gradc.axpy(x(lfsu,i),gradphi[i]);

			// evaluate parameters;
			//Dune::FieldVector<RF,dim> globalpos = eg.geometry().global(it->position());
			RF f = 0;
			RF k = koef.evaluate(eg.entity(), it->position());
			// evaluate D(u)
			Dune::FieldVector<RF,dimw> U;
			dgfg.evaluate(eg.entity(), it->position(), U);
			RF D = 1.0 + k*k*U.two_norm2() / ( 1.0 + std::abs(k)*U.two_norm() );

			// integrate D * grad c * grad phi_i + u * grad c * phi_i - f phi_i
			RF factor = it->weight()*eg.geometry().integrationElement(it->position());
			for (size_type i=0; i<lfsu.size(); ++i)
				r.accumulate(lfsu, i, (D*(gradc*gradphi[i]) + k*(U*gradc)*phi[i] - f*phi[i]) * factor);
		}
	}

	// boundary integral
	template<typename IG, typename LFSU, typename X, typename LFSV, typename R>
	void alpha_boundary (const IG& ig, const LFSU& lfsu_s, const X& x_s,
	                     const LFSV& lfsv_s, R& r_s) const
	{
		// assume Galerkin: lfsu_s == lfsv_s
		// This yields more efficient code since the local functionspace only
		// needs to be evaluated once, but would be incorrect for a finite volume
		// method

		// some types
		typedef typename LFSU::Traits::FiniteElementType::
			Traits::LocalBasisType::Traits::DomainFieldType DF;
		typedef typename LFSU::Traits::FiniteElementType::
			Traits::LocalBasisType::Traits::RangeFieldType RF;
		typedef typename LFSU::Traits::FiniteElementType::
			Traits::LocalBasisType::Traits::RangeType Range;
		typedef typename LFSU::Traits::SizeType size_type;

		// dimensions
		const int dim = IG::dimension;

		// select quadrature rule for face
		Dune::GeometryType gtface = ig.geometryInInside().type();
		const Dune::QuadratureRule<DF,dim-1>&
			rule = Dune::QuadratureRules<DF,dim-1>::rule(gtface,intorder);

		// loop over quadrature points and integrate normal flux
		for (typename Dune::QuadratureRule<DF,dim-1>::const_iterator it=rule.begin();
			it!=rule.end(); ++it)
		{
			// skip rest if we are on Dirichlet boundary
			if ( bctype.isDirichlet( ig, it->position() ) )
				continue;

			// position of quadrature point in local coordinates of element
			Dune::FieldVector<DF,dim> local = ig.geometryInInside().global(it->position());

			// evaluate basis functions at integration point
			std::vector<Range> phi(lfsu_s.size());
			lfsu_s.finiteElement().localBasis().evaluateFunction(local,phi);

			// evaluate c (e.g. flux may depend on c)
			RF c=0.0;
			for (size_type i=0; i<lfsu_s.size(); ++i)
				c += x_s(lfsu_s,i)*phi[i];

			// evaluate flux boundary condition
			Dune::FieldVector<RF,dim>
				globalpos = ig.geometry().global(it->position());
			RF j = 0.0;

			// integrate j
			RF factor = it->weight()*ig.geometry().integrationElement(it->position());
			for (size_type i=0; i<lfsu_s.size(); ++i)
				r_s.accumulate(lfsu_s,i,j*phi[i]*factor);
		}
	}

	void preStep (double time, double dt, int stages) {
		bctype.setTime(time+dt); // enable change of boundary condition type
		Dune::PDELab::InstationaryLocalOperatorDefaultMethods<double>::preStep(time,dt,stages);
  }
private:
	BCType& bctype;
	const DGFG& dgfg;
	const KOEF& koef;
	unsigned int intorder;
};
