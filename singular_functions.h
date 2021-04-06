#ifndef OOMPH_SINGULAR_FUNCTIONS_HEADER
#define OOMPH_SINGULAR_FUNCTIONS_HEADER

// the analytic solution for flow around an edge
#include "moffatt_solution.h"
#include "coordinate_conversions.h"
#include "exact_solutions_finite_disk.h"
#include "additional_maths.h"

namespace SingularFunctions
{  
  // how far to shift the radius to compute "infinity" when the actual
  // radius is zero - for output purposes, since integration points are
  // never on the singularity.
  // N.B. not const, since this might be changed from outside
  double dr = 1e-5;

  // tolerance for some floating-point comparisons
  const double tol = 1e-7;

  const double Pi = MathematicalConstants::Pi;

  // store the warped disk object so that we can use it to get
  // surface normals
  CylindricallyWarpedCircularDiskWithAnnularInternalBoundary*
  warped_disk_with_boundary_pt = nullptr;
  
  // typedefs
  // --------------------------------------------------------------------------
  
  // generic function pointer for a singular function
  typedef Vector<double> (*SingularFctPt)(const LagrangianCoordinates&);

  // generic function pointer for a function that computes the derivatives
  // of the cylindrical singular velocity components w.r.t. the oblate spheroidal coords
  typedef Vector<double> (*DSingularFctDoblateSpheroidalPt)(const LagrangianCoordinates&,
							    DenseMatrix<double>&);

  // --------------------------------------------------------------------------

  // ==========================================================================
  // generic function to compute velocity gradient via finite difference
  // ==========================================================================
  void generic_dudx_finite_diff(const LagrangianCoordinates& lagr_coords,
				const SingularFctPt& sing_fct_pt,
				DenseMatrix<double>& du_dx)
  {
    // coordinate increment for finite-difference gradient
    const double fd_dx = 1e-8;
    
    // make sure we've got enough space
    du_dx.resize(3,3, 0.0);
    
    // solution at the current point
    Vector<double> u0 = (*sing_fct_pt)(lagr_coords);

    // get the Cartesian coords of this point
    Vector<double> x(3, 0.0);
    CoordinateConversions::lagrangian_to_eulerian_coordinates(lagr_coords, x);

    /* // compute dzeta/dx at this Lagrangian point */
    /* Vector<double> dzeta_dx = compute_dzeta_dx(edge_coords); */

    // get the sign of the normal coordinate; +ve => upper disk, -ve => lower disk
    int sign_of_xi3 = lagr_coords.sign_of_xi3();
     
    // now compute the solution at a small increment in each Cartesian direction
    for(unsigned j=0; j<3; j++)
    {
      // coordinate vectors with increments applied to their components
      Vector<double> x_plus_dx = x;

      // add the increments in each direction
      x_plus_dx[j] += fd_dx;

      // convert the incremented position into Lagrangian coords
      LagrangianCoordinates lagr_coords_plus_dx;      
      CoordinateConversions::eulerian_to_lagrangian_coordinates(
	x_plus_dx, lagr_coords_plus_dx);

      // get the sign of phi after the dx shift
      int sign_of_xi3_at_dx = lagr_coords_plus_dx.sign_of_xi3();

      // direction of the increment w.r.t. the +ve coordinate axis
      int direction = 1;

      // did the dx increment move the point to the other side of the disk?
      if (sign_of_xi3 != sign_of_xi3_at_dx)
      {
	// if it did, flip the direction sign
	direction = -1;

	// take the FD step in the opposite direction instead
	x_plus_dx[j] -= 2 * fd_dx;

	CoordinateConversions::eulerian_to_lagrangian_coordinates(
	x_plus_dx, lagr_coords_plus_dx);
      }
      
      // get the solution at the incremented position
      Vector<double> u_plus_dx = (*sing_fct_pt)(lagr_coords_plus_dx);
	
      // now do the finite diff to get the velocity gradients
      for(unsigned i=0; i<3; i++)
    	du_dx(i,j) = double(direction) * (u_plus_dx[i] - u0[i]) / fd_dx;
    }
  }
  
  // ==========================================================================
  // ==========================================================================
  // ==========================================================================

  // ==========================================================================
  // Singular functions, really just wrappers for the exact solutions in
  // Lagrangian coordinates
  // ==========================================================================
  
  Vector<double> singular_fct_exact_broadside_translation(const LagrangianCoordinates& lagr_coords)
  {
    // catch the case where we're sat exactly at the singularity
    const double xi1 = (abs(lagr_coords.xi1 - 1.0) < tol) && ((abs(lagr_coords.xi3) < tol)) ?
      1.0 - dr : lagr_coords.xi1;
    
    // convert the Lagrangian curvilinear coordinates into 'flat'
    // oblate spheroidal coordinates
    OblateSpheroidalCoordinates obl_sph_coords(xi1, lagr_coords.xi3);

    // velocity components in 'cylindrical' coordinates (r, theta, z)
    Vector<double> u_lagr(3, 0.0);

    // pressure
    double p = 0.0;
    
    // get the flat disk solution
    FlatDiskExactSolutions::broadside_translation_solution_cylindrical(obl_sph_coords, u_lagr, p);
    
    // now convert to Cartesian
    // ------------------------
    Vector<double> u(3, 0.0);
    CoordinateConversions::lagrangian_to_eulerian_velocity(lagr_coords, u_lagr, u);
    
#ifdef PARANOID    
    // QUEHACERES debug
    bool inf_or_nan = isinf(u[0]) || isinf(u[1]) || isinf(u[2])
      || isnan(u[0]) || isnan(u[1]) || isnan(u[2]);

    if(inf_or_nan)
    {
      oomph_info << "INF or NAN in cartesian coords; singular_fct_exact_asymptotic_broadside()";
      abort();
    }
#endif
    
    // and add the pressure
    u.push_back(p);

    return u;
  }

  // --------------------------------------------------------------------------
  
  Vector<double> singular_fct_exact_broadside_rotation(const LagrangianCoordinates& lagr_coords)
  {
    // catch the case where we're sat exactly at the singularity
    const double xi1 = (abs(lagr_coords.xi1 - 1) < tol) && ((abs(lagr_coords.xi3) < tol)) ?
      1 - dr : lagr_coords.xi1;

    // set the 'azimuthal' angle to where cos(theta) = sin(theta) for slice solution
    const double theta = Pi / 4.0;
    
    // convert the Lagrangian curvilinear coordinates into 'flat'
    // oblate spheroidal coordinates
    OblateSpheroidalCoordinates obl_sph_coords(xi1, theta, lagr_coords.xi3);
    
    // velocity components in 'cylindrical' coordinates (r, theta, z)
    Vector<double> u_lagr(3, 0.0);

    // pressure
    double p = 0.0;

    // get it
    FlatDiskExactSolutions::out_of_plane_rotation_solution_cylindrical(obl_sph_coords, u_lagr, p);
    
    // now convert to Cartesian
    // ------------------------
    Vector<double> u(3, 0.0);
    CoordinateConversions::lagrangian_to_eulerian_velocity(lagr_coords, u_lagr, u);

    // and add the pressure
    u.push_back(p);

    return u;
  }
  
  // --------------------------------------------------------------------------
  
  Vector<double> singular_fct_exact_in_plane_translation(const LagrangianCoordinates& lagr_coords)
  {
    // catch the case where we're sat exactly at the singularity
    const double xi1 = (abs(lagr_coords.xi1 - 1) < tol) && ((abs(lagr_coords.xi3) < tol)) ?
      1 - dr : lagr_coords.xi1;

    // set the 'azimuthal' angle to where cos(theta) = sin(theta) for slice solution
    const double theta = Pi / 4.0;
    
    // convert the Lagrangian curvilinear coordinates into 'flat'
    // oblate spheroidal coordinates
    OblateSpheroidalCoordinates obl_sph_coords(xi1, theta, lagr_coords.xi3);
    
    // velocity components in 'cylindrical' coordinates (r, theta, z)
    Vector<double> u_lagr(3, 0.0);

    // pressure
    double p = 0.0;
    
    // get the flat disk solution
    FlatDiskExactSolutions::in_plane_translation_solution_cylindrical(obl_sph_coords,
								      u_lagr, p);
    
    // now convert to Cartesian
    // ------------------------
    Vector<double> u(3, 0.0);
    CoordinateConversions::lagrangian_to_eulerian_velocity(lagr_coords, u_lagr, u);

#ifdef PARANOID    
    // QUEHACERES debug
    bool inf_or_nan = isinf(u[0]) || isinf(u[1]) || isinf(u[2])
      || isnan(u[0]) || isnan(u[1]) || isnan(u[2]);

    if(inf_or_nan)
    {
      oomph_info << "INF or NAN in cartesian coords; "
		 << "__singular_fct_exact_in_plane_translation()" << std::endl;
      abort();
    }
#endif
    
    // and add the pressure
    u.push_back(p);

    return u;
  }
    
  // --------------------------------------------------------------------------
  
  Vector<double> singular_fct_exact_in_plane_rotation(const LagrangianCoordinates& lagr_coords)
  {
    // catch the case where we're sat exactly at the singularity
    const double xi1 = (abs(lagr_coords.xi1 - 1) < tol) && ((abs(lagr_coords.xi3) < tol)) ?
      1 - dr : lagr_coords.xi1;

    // set the 'azimuthal' angle to where cos(theta) = sin(theta) for slice solution
    const double theta = Pi / 4.0;

    // convert the Lagrangian curvilinear coordinates into 'flat'
    // oblate spheroidal coordinates
    OblateSpheroidalCoordinates obl_sph_coords(xi1, theta, lagr_coords.xi3);
    
    // velocity components in 'cylindrical' coordinates (r, theta, z)
    Vector<double> u_lagr(3, 0.0);

    // pressure
    double p = 0.0;

    // get it
    FlatDiskExactSolutions::in_plane_rotation_solution_cylindrical(obl_sph_coords, u_lagr, p);
      
    // now convert to Cartesian
    // ------------------------
    Vector<double> u(3, 0.0);
    CoordinateConversions::lagrangian_to_eulerian_velocity(lagr_coords, u_lagr, u);

#ifdef PARANOID        
    // QUEHACERES debug
    bool inf_or_nan = isinf(u[0]) || isinf(u[1]) || isinf(u[2])
      || isnan(u[0]) || isnan(u[1]) || isnan(u[2]);

    if(inf_or_nan)
    {
      oomph_info << "INF or NAN in cartesian coords; "
		 << "__singular_fct_exact_in_plane_rotation()" << std::endl;
      abort();
    }
#endif
    
    // and add the pressure
    u.push_back(p);

    return u;
  }

  // ==========================================================================
  // ==========================================================================
  // ==========================================================================

  // ==========================================================================
  /// \short Derivatives of cylindrical velocity components with respect
  /// to oblate spheroidal coordinates (azimuthal derivatives set to zero,
  /// since these are 'sliced' versions of the solution and the azimuthally
  /// varying singular amplitude takes care of the azimuthal derivatives
  // ==========================================================================

  void du_broadside_cylindrical_doblate_sph(const LagrangianCoordinates& lagr_coords,
					    DenseMatrix<double>& du_cyl_doblate_sph)
  {
    // catch the case where we're sat exactly at the singularity
    const double xi1 = (abs(lagr_coords.xi1 - 1.0) < tol) && ((abs(lagr_coords.xi3) < tol)) ?
      1.0 - dr : lagr_coords.xi1;
    
    // convert the Lagrangian curvilinear coordinates into 'flat'
    // oblate spheroidal coordinates
    OblateSpheroidalCoordinates obl_sph_coords(xi1, lagr_coords.xi3);

    // shorthands
    const double lambda = obl_sph_coords.lambda;
    const double zeta   = obl_sph_coords.zeta;
    
    du_cyl_doblate_sph.resize(3, 3, 0.0);

    // du_r/dlambda
    du_cyl_doblate_sph(0,0) = (2*lambda*zeta*Sqrt((1 - Power(zeta,2))/(1 + Power(lambda,2)))
			       *(-Power(lambda,4) + (2 + Power(lambda,2))*Power(zeta,2))) /
      ((1 + Power(lambda,2))*Pi*Power(Power(lambda,2) + Power(zeta,2),2));

    // du_r/dzeta
    du_cyl_doblate_sph(0,1) = (-2*Power(lambda,2)*(Power(zeta,2) + Power(lambda,2)*
						   (-1 + 2*Power(zeta,2)))) /
      ((1 + Power(lambda,2))*Pi*Sqrt((1 - Power(zeta,2))/(1 + Power(lambda,2)))
       *Power(Power(lambda,2) + Power(zeta,2),2));

    // du_r/dtheta = 0
    
    // du_theta/dlambda
    du_cyl_doblate_sph(1,0) = 0.0;

    // du_theta/dzeta
    du_cyl_doblate_sph(1,1) = 0.0;

    // du_theta/dtheta = 0
    
    // du_z/dlambda
    du_cyl_doblate_sph(2,0) = (-2*Power(lambda,2)*(Power(lambda,2) + (3 + Power(lambda,2))
						   *Power(zeta,2) - Power(zeta,4))) /
      ((1 + Power(lambda,2))*Pi*Power(Power(lambda,2) + Power(zeta,2),2));

    // du_z/dzeta
    du_cyl_doblate_sph(2,1) = (4*Power(lambda,3)*zeta)/(Pi*Power(Power(lambda,2) + Power(zeta,2),2));

    // du_z/dtheta = 0
  }

  // --------------------------------------------------------------------------
  
  void du_broadside_rotation_cylindrical_doblate_sph(const LagrangianCoordinates& lagr_coords,
						     DenseMatrix<double>& du_cyl_doblate_sph)
  {
    // catch the case where we're sat exactly at the singularity
    const double xi1 = (abs(lagr_coords.xi1 - 1.0) < tol) && ((abs(lagr_coords.xi3) < tol)) ?
      1.0 - dr : lagr_coords.xi1;
    
    // convert the Lagrangian curvilinear coordinates into 'flat'
    // oblate spheroidal coordinates
    OblateSpheroidalCoordinates obl_sph_coords(xi1, lagr_coords.xi2, lagr_coords.xi3);

    // shorthands
    const double lambda = obl_sph_coords.lambda;
    const double zeta   = obl_sph_coords.zeta;
    const double theta  = obl_sph_coords.theta;
    
    // get the solution (for theta derivs)
    Vector<double> u_cyl(3, 0.0);
    double p_dummy = 0.0;
    FlatDiskExactSolutions::out_of_plane_rotation_solution_cylindrical(obl_sph_coords,
								       u_cyl, p_dummy);
    
    du_cyl_doblate_sph.resize(3, 3, 0.0);

    // du_r/dlambda
    du_cyl_doblate_sph(0,0) = cos(theta) * (2*zeta*((lambda*(-1 + Power(lambda,2)))/Power(1 + Power(lambda,2),2)
				       - (4*Power(lambda,3))/Power(Power(lambda,2) + Power(zeta,2),2)
				       + (4*lambda)/(Power(lambda,2) + Power(zeta,2))
				       - ArcCot(lambda)))/Pi;

    // du_r/dzeta
    du_cyl_doblate_sph(0,1) = cos(theta) * (-2*lambda*(-Power(lambda,5) + lambda*Power(zeta,2)*
					  (2 + Power(zeta,2)) +
					  Power(lambda,3)*(-2 + 4*Power(zeta,2)) +
					  (1 + Power(lambda,2))*Power(Power(lambda,2) +
								      Power(zeta,2),2)*ArcCot(lambda)))/
      ((1 + Power(lambda,2))*Pi*Power(Power(lambda,2) + Power(zeta,2),2));

    // du_r/dtheta
    du_cyl_doblate_sph(0,2) = 0.0; 
    
    // du_theta/dlambda
    du_cyl_doblate_sph(1,0) = sin(theta) * (2*zeta*(-((lambda*(3 + Power(lambda,2))) /
					 Power(1 + Power(lambda,2),2)) + ArcCot(lambda)))/Pi;

    // du_theta/dzeta
    du_cyl_doblate_sph(1,1) = sin(theta) * (2*lambda*(-(lambda/(1 + Power(lambda,2))) +
						      ArcCot(lambda)))/Pi;

    // du_theta/dtheta
    du_cyl_doblate_sph(1,2) = 0.0; 

    // du_z/dlambda
    du_cyl_doblate_sph(2,0) = cos(theta) * (2*(-1 + Power(zeta,2))*(Power(lambda,2)*
						       (Power(lambda,2)*(2 + Power(lambda,2)) +
							6*(1 + Power(lambda,2))*Power(zeta,2) +
							Power(zeta,4)) -
						       (lambda + Power(lambda,3))*
						       Power(Power(lambda,2) +
							     Power(zeta,2),2)*ArcCot(lambda))) /
      (Power(1 + Power(lambda,2),2)*Pi*Sqrt((1 - Power(zeta,2)) /
					    (1 + Power(lambda,2)))*Power(Power(lambda,2)
									 + Power(zeta,2),2));

    // du_z/dzeta
    du_cyl_doblate_sph(2,1) = cos(theta) * (-2*zeta*(-Power(lambda,5) + lambda*Power(zeta,4) +
					4*Power(lambda,3)*(-1 + Power(zeta,2)) +
					(1 + Power(lambda,2))*Power(Power(lambda,2) +
								    Power(zeta,2),2)*ArcCot(lambda))) /
      ((1 + Power(lambda,2))*Pi*Sqrt((1 - Power(zeta,2))/
				     (1 + Power(lambda,2)))*Power(Power(lambda,2) +
								  Power(zeta,2),2));

    // du_z/dtheta
    du_cyl_doblate_sph(2,2) = 0.0; 
  }

  // --------------------------------------------------------------------------

  void du_in_plane_cylindrical_doblate_sph(const LagrangianCoordinates& lagr_coords,
					   DenseMatrix<double>& du_cyl_doblate_sph)
  {
    // catch the case where we're sat exactly at the singularity
    const double xi1 = (abs(lagr_coords.xi1 - 1.0) < tol) && ((abs(lagr_coords.xi3) < tol)) ?
      1.0 - dr : lagr_coords.xi1;
    
    // convert the Lagrangian curvilinear coordinates into 'flat'
    // oblate spheroidal coordinates
    OblateSpheroidalCoordinates obl_sph_coords(xi1, lagr_coords.xi2, lagr_coords.xi3);

    // shorthands
    const double lambda = obl_sph_coords.lambda;
    const double zeta   = obl_sph_coords.zeta;
    const double theta  = obl_sph_coords.theta;

    // get the solution without azimuthal dependence (for theta derivs)
    Vector<double> u_cyl(3, 0.0);
    double p_dummy = 0.0;
    FlatDiskExactSolutions::in_plane_translation_solution_cylindrical(obl_sph_coords,
								      u_cyl, p_dummy);
    
    du_cyl_doblate_sph.resize(3, 3, 0.0);

    // du_r/dlambda
    du_cyl_doblate_sph(0,0) = cos(theta) * (4*(-((2 + 4*Power(lambda,2) + Power(lambda,4)) /
				    Power(1 + Power(lambda,2),2)) - (2*Power(lambda,4)) /
				  Power(Power(lambda,2) + Power(zeta,2),2) +
				  (3*Power(lambda,2))/(Power(lambda,2) +
						       Power(zeta,2))))/(3.*Pi);

    // du_r/dzeta
    du_cyl_doblate_sph(0,1) = cos(theta) * (-8*Power(lambda,3)*zeta)/(3.*Pi*Power(Power(lambda,2) + Power(zeta,2),2));

    // du_r_dtheta
    du_cyl_doblate_sph(0,2) = 0.0; 
    
    // du_theta/dlambda
    du_cyl_doblate_sph(1,0) = sin(theta) * (4*(2 + Power(lambda,2)))/(3.*Power(1 + Power(lambda,2),2)*Pi);

    // du_theta/dzeta
    du_cyl_doblate_sph(1,1) = 0.0;

    // du_theta/dtheta
    du_cyl_doblate_sph(1,2) = 0.0; 
    
    // du_z/dlambda
    du_cyl_doblate_sph(2,0) = cos(theta) * (4*lambda*zeta*Sqrt((1 - Power(zeta,2)) /
						  (1 + Power(lambda,2)))*(-Power(lambda,4) +
									  (2 + Power(lambda,2))*
									  Power(zeta,2)))/
      (3.*(1 + Power(lambda,2))*Pi*Power(Power(lambda,2) + Power(zeta,2),2));

    // du_z/dzeta
    du_cyl_doblate_sph(2,1) = cos(theta) * (-4*Power(lambda,2)*(Power(zeta,2) + Power(lambda,2)*
						   (-1 + 2*Power(zeta,2)))) /
      (3.*(1 + Power(lambda,2))*Pi*Sqrt((1 - Power(zeta,2)) /
					(1 + Power(lambda,2)))*
       Power(Power(lambda,2) + Power(zeta,2),2));

    // du_z/dtheta
    du_cyl_doblate_sph(2,2) = 0.0; 
  }

  // --------------------------------------------------------------------------
  
  void du_in_plane_rotation_cylindrical_doblate_sph(const LagrangianCoordinates& lagr_coords,
						    DenseMatrix<double>& du_cyl_doblate_sph)
  {
    // catch the case where we're sat exactly at the singularity
    const double xi1 = (abs(lagr_coords.xi1 - 1.0) < tol) && ((abs(lagr_coords.xi3) < tol)) ?
      1.0 - dr : lagr_coords.xi1;
    
    // convert the Lagrangian curvilinear coordinates into 'flat'
    // oblate spheroidal coordinates
    OblateSpheroidalCoordinates obl_sph_coords(xi1, lagr_coords.xi3);

    // shorthands
    const double lambda = obl_sph_coords.lambda;
    const double zeta   = obl_sph_coords.zeta;
    
    du_cyl_doblate_sph.resize(3, 3, 0.0);

    // du_r/dlambda
    du_cyl_doblate_sph(0,0) = 0.0;

    // du_r/dzeta
    du_cyl_doblate_sph(0,1) = 0.0;

    // du_theta/dlambda
    du_cyl_doblate_sph(1,0) = (2*Sqrt((1 - Power(zeta,2))/(1 + Power(lambda,2)))*(-2 - Power(lambda,2) + (lambda + Power(lambda,3))*ArcCot(lambda)))/((1 + Power(lambda,2))*Pi);

    // du_theta/dzeta
    du_cyl_doblate_sph(1,1) = (2*zeta*(lambda - (1 + Power(lambda,2))*ArcCot(lambda)))/((1 + Power(lambda,2))*Pi*Sqrt((1 - Power(zeta,2))/(1 + Power(lambda,2))));

    // du_z/dlambda
    du_cyl_doblate_sph(2,0) = 0.0;

    // du_z/dzeta
    du_cyl_doblate_sph(2,1) = 0.0;
  }

  // ==========================================================================
  // ==========================================================================
  // ==========================================================================

  // ==========================================================================
  /// \short Functions to handle the conversion from the natural oblate
  /// spheroidal derivatives of the cylindrical velocity components to
  /// Eulerian derivatives of Eulerian components
  // ==========================================================================
  

  // ==========================================================================
  // Generic function for computing the derivatives of the Cartesian singular
  // velocity components w.r.t. the curvilinear coordinates
  // ==========================================================================
  void du_sing_dx_eulerian(const LagrangianCoordinates& lagr_coords,
			   const Vector<double>& u_cyl,
			   const DenseMatrix<double>& du_cyl_dobl_sph,			       
			   DenseMatrix<double>& du_dx)
  {
#ifdef PARANOID
    
    if(warped_disk_with_boundary_pt == nullptr)
    {
      throw OomphLibError("Error, SingularFunctions::warped_disk_with_boundary_pt "
			  "pointer hasn't been set\n",
			  OOMPH_CURRENT_FUNCTION,
			  OOMPH_EXCEPTION_LOCATION);
    }
    
#endif
    
    // basis vectors
    Vector<double> a1, a2, a3;
    warped_disk_with_boundary_pt->basis_vectors(lagr_coords, a1, a2, a3);

    // curvilinear derivatives of basis vectors
    DenseMatrix<double> da1_dxi, da2_dxi, da3_dxi;
    warped_disk_with_boundary_pt->da_dxi(lagr_coords, da1_dxi, da2_dxi, da3_dxi);
    
    // make vectors for easier indexing
    // --------------------------------

    Vector<Vector<double>> a(3);
    a[0] = a1;
    a[1] = a2;
    a[2] = a3;
    
    Vector<DenseMatrix<double>> da_dxi(3);

    da_dxi[0] = da1_dxi;
    da_dxi[1] = da2_dxi;
    da_dxi[2] = da3_dxi;

    // oblate spheroidal coords for the corresponding flat disk
    OblateSpheroidalCoordinates obl_sph_coords(lagr_coords.xi1, lagr_coords.xi3);

    // derivatives of the oblate spheroidal coordinates w.r.t. the cylindrical
    // coordinates
    DenseMatrix<double> dobl_sph_dcyl(3, 3, 0.0);
    obl_sph_coords.doblate_spheroidal_dcylindrical(dobl_sph_dcyl);

    // compute the derivatives of the Cartesian velocity components w.r.t. the
    // curviliear coordinates
    // --------------------------------------------------
    DenseMatrix<double> du_eulerian_dxi(3, 3, 0.0);
    
    for(unsigned i=0; i<3; i++)
    {
      for(unsigned j=0; j<3; j++)
      {
	for(unsigned k=0; k<3; k++)
	{
	  du_eulerian_dxi(i,j) += da_dxi[k](i,j) * u_cyl[k]; 	    

	  for(unsigned l=0; l<3; l++)
	  {
	    du_eulerian_dxi(i,j) += a[k][i] * du_cyl_dobl_sph(k,l) * dobl_sph_dcyl(l,j);
	  }
	}
      }
    }

    // now compute the Eulerian derivatives
    // ---------------------------------------

    // get the Eulerian derivatives of the curvilinear coordinates
    DenseMatrix<double> dxi_dx = CoordinateConversions::dxi_dx(lagr_coords);

    // zero everything out as we're doing += 
    du_dx.resize(3, 3, 0.0);
    du_dx.initialise(0.0);
    
    for(unsigned i=0; i<3; i++)
    {
      for(unsigned j=0; j<3; j++)
      {
	for(unsigned k=0; k<3; k++)
	{
	  // compute via the chain rule
	  du_dx(i,j) += du_eulerian_dxi(i,k) * dxi_dx(k,j);
	}
      }
    }
  }
  
  // ==========================================================================
  // ==========================================================================

  
  // ==========================================================================
  // Broadside translation Eulerian velocity gradient
  // ==========================================================================
  DenseMatrix<double> gradient_of_singular_fct_exact_broadside_translation(
    const LagrangianCoordinates& lagr_coords)
  {
    // catch the case where we're sat exactly at the singularity
    const double xi1 = (abs(lagr_coords.xi1 - 1.0) < tol) && ((abs(lagr_coords.xi3) < tol)) ?
      1.0 - dr : lagr_coords.xi1;
    
    // convert the Lagrangian curvilinear coordinates into 'flat'
    // oblate spheroidal coordinates
    OblateSpheroidalCoordinates obl_sph_coords(xi1, lagr_coords.xi3);

    // velocity components in 'cylindrical' coordinates (r, theta, z)
    Vector<double> u_lagr(3, 0.0);

    // pressure
    double p = 0.0;
    
    // get the flat disk solution
    FlatDiskExactSolutions::broadside_translation_solution_cylindrical(obl_sph_coords, u_lagr, p);

    // get the derivatives of the 'cylindrical' components of the singular velocity
    // w.r.t. the oblate spheroidal coordinates
    DenseMatrix<double> du_cyl_doblate_sph;
    du_broadside_cylindrical_doblate_sph(lagr_coords, du_cyl_doblate_sph);
	
    // get the singular velocity gradient
    DenseMatrix<double> du_dx_sing(3, 3, 0.0);

    du_sing_dx_eulerian(lagr_coords, u_lagr, du_cyl_doblate_sph, du_dx_sing);

#ifdef PARANOID
    for(unsigned i=0; i<3; i++)
    {
      for(unsigned j=0; j<3; j++)
      {
	if(!isfinite(du_dx_sing(i,j)))
	{
	  ostringstream error_message;

	  error_message << "Error: du_dx_broadside(" << i << "," << j << ") = "
			<< du_dx_sing(i,j) << std::endl;
	    
	  throw OomphLibError(error_message.str().c_str(),
			      OOMPH_CURRENT_FUNCTION,
			      OOMPH_EXCEPTION_LOCATION);
	}
      }
    }
#endif
    
    return du_dx_sing;
  }

  // --------------------------------------------------------------------------
  
  DenseMatrix<double> gradient_of_singular_fct_exact_broadside_rotation(
    const LagrangianCoordinates& lagr_coords)
  {
    // catch the case where we're sat exactly at the singularity
    const double xi1 = (abs(lagr_coords.xi1 - 1.0) < tol) && ((abs(lagr_coords.xi3) < tol)) ?
      1.0 - dr : lagr_coords.xi1;
    
    // set the 'azimuthal' angle to where cos(theta) = sin(theta) for slice solution
    const double theta = Pi / 4.0;
        
    // convert the Lagrangian curvilinear coordinates into 'flat'
    // oblate spheroidal coordinates
    OblateSpheroidalCoordinates obl_sph_coords_slice(xi1, theta, lagr_coords.xi3);

    // velocity components in 'cylindrical' coordinates (r, theta, z)
    Vector<double> u_lagr(3, 0.0);

    // pressure
    double p = 0.0;
    
    // get the flat disk solution
    FlatDiskExactSolutions::out_of_plane_rotation_solution_cylindrical(obl_sph_coords_slice,
								       u_lagr, p);

    // Lagrangian coordinates of the slice where cos(theta) = sin(theta)
    LagrangianCoordinates lagr_coords_slice(lagr_coords.xi1, theta, lagr_coords.xi3);
    
    // get the derivatives of the 'cylindrical' components of the singular velocity
    // w.r.t. the oblate spheroidal coordinates
    DenseMatrix<double> du_cyl_doblate_sph;
    du_broadside_rotation_cylindrical_doblate_sph(lagr_coords_slice, du_cyl_doblate_sph);
    
    // get the singular velocity gradient, passing the actual (non-sliced) Lagrangian coords
    DenseMatrix<double> du_dx_sing(3, 3, 0.0);

    du_sing_dx_eulerian(lagr_coords, u_lagr, du_cyl_doblate_sph, du_dx_sing);

#ifdef PARANOID
    for(unsigned i=0; i<3; i++)
    {
      for(unsigned j=0; j<3; j++)
      {
	if(!isfinite(du_dx_sing(i,j)))
	{
	  ostringstream error_message;

	  error_message << "Error: du_dx_broadside_rot(" << i << "," << j << ") = "
			<< du_dx_sing(i,j) << std::endl;
	    
	  throw OomphLibError(error_message.str().c_str(),
			      OOMPH_CURRENT_FUNCTION,
			      OOMPH_EXCEPTION_LOCATION);
	}
      }
    }
#endif
    
    return du_dx_sing;
  }

  // --------------------------------------------------------------------------
  
  DenseMatrix<double> gradient_of_singular_fct_exact_in_plane_translation(
    const LagrangianCoordinates& lagr_coords)
  {      
    // catch the case where we're sat exactly at the singularity
    const double xi1 = (abs(lagr_coords.xi1 - 1.0) < tol) && ((abs(lagr_coords.xi3) < tol)) ?
      1.0 - dr : lagr_coords.xi1;

    // set the 'azimuthal' angle to where cos(theta) = sin(theta) for slice solution
    const double theta = Pi / 4.0;
     
    // convert the Lagrangian curvilinear coordinates into 'flat'
    // oblate spheroidal coordinates
    OblateSpheroidalCoordinates obl_sph_coords_slice(xi1, theta, lagr_coords.xi3);

    // velocity components in 'cylindrical' coordinates (r, theta, z)
    Vector<double> u_lagr(3, 0.0);

    // pressure
    double p = 0.0;
    
    // get the flat disk solution
    FlatDiskExactSolutions::in_plane_translation_solution_cylindrical(obl_sph_coords_slice,
								      u_lagr, p);

    // Lagrangian coordinates of the slice where cos(theta) = sin(theta)
    LagrangianCoordinates lagr_coords_slice(lagr_coords.xi1, theta, lagr_coords.xi3);
    
    // get the derivatives of the 'cylindrical' components of the singular velocity
    // w.r.t. the oblate spheroidal coordinates
    DenseMatrix<double> du_cyl_doblate_sph;
    du_in_plane_cylindrical_doblate_sph(lagr_coords_slice, du_cyl_doblate_sph);
    
    // get the singular velocity gradient, passing the actual (non-sliced) Lagrangian coords
    DenseMatrix<double> du_dx_sing(3, 3, 0.0);

    du_sing_dx_eulerian(lagr_coords, u_lagr, du_cyl_doblate_sph, du_dx_sing);

#ifdef PARANOID
    for(unsigned i=0; i<3; i++)
    {
      for(unsigned j=0; j<3; j++)
      {
	if(!isfinite(du_dx_sing(i,j)))
	{
	  ostringstream error_message;

	  error_message << "Error: du_dx_in_plane(" << i << "," << j << ") = "
			<< du_dx_sing(i,j) << std::endl;
	    
	  throw OomphLibError(error_message.str().c_str(),
			      OOMPH_CURRENT_FUNCTION,
			      OOMPH_EXCEPTION_LOCATION);
	}
      }
    }
#endif
    return du_dx_sing;
  }

  // --------------------------------------------------------------------------
  
  DenseMatrix<double> gradient_of_singular_fct_exact_in_plane_rotation(
    const LagrangianCoordinates& lagr_coords)
  {    
    // catch the case where we're sat exactly at the singularity
    const double xi1 = (abs(lagr_coords.xi1 - 1.0) < tol) && ((abs(lagr_coords.xi3) < tol)) ?
      1.0 - dr : lagr_coords.xi1;
    
    // convert the Lagrangian curvilinear coordinates into 'flat'
    // oblate spheroidal coordinates
    // N.B. this function is axisymmetric, so can just (implicitly) pass theta=0 here
    OblateSpheroidalCoordinates obl_sph_coords(xi1, lagr_coords.xi3);

    // velocity components in 'cylindrical' coordinates (r, theta, z)
    Vector<double> u_lagr(3, 0.0);

    // pressure
    double p = 0.0;
    
    // get the flat disk solution
    FlatDiskExactSolutions::in_plane_rotation_solution_cylindrical(obl_sph_coords, u_lagr, p);

    // get the derivatives of the 'cylindrical' components of the singular velocity
    // w.r.t. the oblate spheroidal coordinates
    DenseMatrix<double> du_cyl_doblate_sph;
    du_in_plane_rotation_cylindrical_doblate_sph(lagr_coords, du_cyl_doblate_sph);
	
    // get the singular velocity gradient
    DenseMatrix<double> du_dx_sing(3, 3, 0.0);

    du_sing_dx_eulerian(lagr_coords, u_lagr, du_cyl_doblate_sph, du_dx_sing);

#ifdef PARANOID
    for(unsigned i=0; i<3; i++)
    {
      for(unsigned j=0; j<3; j++)
      {
	if(!isfinite(du_dx_sing(i,j)))
	{
	  ostringstream error_message;

	  error_message << "Error: du_dx_in_plane_rot(" << i << "," << j << ") = "
			<< du_dx_sing(i,j) << std::endl;
	    
	  throw OomphLibError(error_message.str().c_str(),
			      OOMPH_CURRENT_FUNCTION,
			      OOMPH_EXCEPTION_LOCATION);
	}
      }
    }
#endif
    
    return du_dx_sing;
  }

  // @@@@@@@@
  // FD versions - keep around for the time being
  
  /* DenseMatrix<double> __gradient_of_singular_fct_exact_broadside_translation( */
  /*   const LagrangianCoordinates& lagr_coords) */
  /* { */
  /*   DenseMatrix<double> du_dx(3, 3, 0.0); */
    
  /*   // compute the velocity gradient via finite difference */
  /*   generic_dudx_finite_diff(lagr_coords, */
  /* 			     &__singular_fct_exact_broadside_translation, */
  /* 			     du_dx); */

  /*   return du_dx; */
  /* } */

  /* DenseMatrix<double> __gradient_of_singular_fct_exact_broadside_rotation( */
  /*   const LagrangianCoordinates& lagr_coords) */
  /* { */
  /*   DenseMatrix<double> du_dx(3, 3, 0.0); */
    
  /*   // compute the velocity gradient via finite difference */
  /*   generic_dudx_finite_diff(lagr_coords, */
  /* 			     &__singular_fct_exact_broadside_rotation, */
  /* 			     du_dx); */

  /*   return du_dx; */
  /* } */

  /* DenseMatrix<double> __gradient_of_singular_fct_exact_in_plane_translation( */
  /*   const LagrangianCoordinates& lagr_coords) */
  /* { */
  /*   DenseMatrix<double> du_dx(3, 3, 0.0); */
    
  /*   // compute the velocity gradient via finite difference */
  /*   generic_dudx_finite_diff(lagr_coords, */
  /* 			     &__singular_fct_exact_in_plane_translation, */
  /* 			     du_dx); */

  /*   return du_dx; */
  /* } */

  /* DenseMatrix<double> __gradient_of_singular_fct_exact_in_plane_rotation( */
  /*   const LagrangianCoordinates& lagr_coords) */
  /* { */
  /*   DenseMatrix<double> du_dx(3, 3, 0.0); */
    
  /*   // compute the velocity gradient via finite difference */
  /*   generic_dudx_finite_diff(lagr_coords, */
  /* 			     &__singular_fct_exact_in_plane_rotation, */
  /* 			     du_dx); */

  /*   return du_dx; */
  /* } */
}

#endif
