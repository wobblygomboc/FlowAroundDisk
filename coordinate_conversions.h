#ifndef OOMPH_COORDINATE_CONVERSIONS
#define OOMPH_COORDINATE_CONVERSIONS

#include "additional_maths.h"
#include "warped_disk_with_torus_geom_obj.h"

namespace CoordinateConversions
{
  const double PI = MathematicalConstants::Pi;

  // QUEHACERES change this to something more generic at some point
  // pointer to a geometric object which provides things like
  // boundary triad vectors and their derivatives
  CylindricallyWarpedCircularDiskWithAnnularInternalBoundary* disk_geom_obj_pt = nullptr;
    
  // ==========================================================================
  // Function which returns a cartesian unit vector in the i direction
  // ==========================================================================
  Vector<double> cartesian_unit_vector(const unsigned& i)
  {
    Vector<double> e(3, 0.0);
    e[i] = 1.0;

    return e;
  }
  
  // ==========================================================================
  // Converts velocities given in arbitrary Lagrangian coordinates to
  // global Cartesian (x,y,z) coordinates
  // ==========================================================================
  void lagrangian_to_eulerian_velocity(const LagrangianCoordinates& lagr_coords,
  				       const Vector<double>& u_lagrangian,
  				       Vector<double>& u_cartesian)
  {
    // Lagrangian triad unit vectors
    Vector<double> a1(3, 0.0);
    Vector<double> a2(3, 0.0);
    Vector<double> a3(3, 0.0);
    
    Vector<double> xi = lagr_coords.vector();
    
    // get 'em
    disk_geom_obj_pt->basis_vectors(lagr_coords, a1, a2, a3);

    // make a vector of the basis vectors for easier indexing
    Vector<Vector<double>> a(3);
    a[0] = a1;
    a[1] = a2;
    a[2] = a3;
    
    // make enough space
    u_cartesian.resize(3, 0.0);
    u_cartesian.initialise(0.0);
    
    // compute the Eulerian coordinates of this point for the derivatives
    // ---------------------------------

    // loop over the Eulerian coordinates
    for(unsigned i=0; i<3; i++)
    {
      // Cartesian unit vector in the ith direction
      Vector<double> e_hat = cartesian_unit_vector(i);

      // loop over the curvilinear coordinates
      for(unsigned j=0; j<3; j++)
      {
	// dot product loop
	for(unsigned k=0; k<3; k++)
	{
	  u_cartesian[i] += e_hat[k] * a[j][k] * u_lagrangian[j];
	}
      }
    }
  }
  
  // ==========================================================================
  // Conversion from the Lagrangian edge coordinates to the global Cartesians
  // ==========================================================================  
  void lagrangian_to_eulerian_coordinates(const LagrangianCoordinates& lagr_coords,
					  Vector<double>& x)
  {    
    // get the Eulerian coordinates of the corresponding disk position
    disk_geom_obj_pt->position(lagr_coords, x);

    // Lagrangian triad unit vectors
    Vector<double> a1(3, 0.0);
    Vector<double> a2(3, 0.0);
    Vector<double> a3(3, 0.0);

    // get 'em
    disk_geom_obj_pt->basis_vectors(lagr_coords, a1, a2, a3);

    // now we move in a normal direction away from the surface by an amount xi_3
    for(unsigned i=0; i<3; i++)
    {
      x[i] += lagr_coords.xi3 * a3[i];
    }
  }
  
  // ==========================================================================
  // Compute the residual vector for an Eulerian to Lagrangian coordinate
  // transformation (with the correct prototype for a black-box FD Newton solve
  // ==========================================================================
  void eulerian_to_lagrangian_coord_residual(const Vector<double>& x,
					     const Vector<double>& unknowns,
					     Vector<double>& residuals)
  {
    // shorthand
    const double xi3 = unknowns[2];

    // reset the residual vector
    residuals.resize(3, 0.0);
    residuals.initialise(0.0);

    // basis vectors at this Lagrangian point on the surface
    Vector<double> a1(3, 0.0);
    Vector<double> a2(3, 0.0);
    Vector<double> a3(3, 0.0);

    // Eulerian position vector of the surface point corresponding to
    // these Lagrangian coords
    Vector<double> r(3, 0.0);

    LagrangianCoordinates lagr_coords(unknowns);
    
    disk_geom_obj_pt->position(lagr_coords, r);
    
    // get the surface vectors from the geometric object at these surface coords
    disk_geom_obj_pt->
      basis_vectors(lagr_coords, a1, a2, a3);
        
    for(unsigned i=0; i<3; i++)
    {
      // (p - r)\cdot a_1 = 0
      residuals[0] += (x[i] - r[i]) * a1[i];

      // (p - r)\cdot a_2 = 0
      residuals[1] += (x[i] - r[i]) * a2[i];

      // (p - r)\cdot a_3 - xi_3 = 0
      residuals[2] += (x[i] - r[i]) * a3[i];
    }
    
    residuals[2] -= xi3;
  }

  // ==========================================================================
  // Perform a black-box finite-diff Newton solve to convert from eulerian
  // to Lagrangian edge coordinates
  // ==========================================================================
  void eulerian_to_lagrangian_coordinates(const Vector<double>& x,
					  LagrangianCoordinates& lagr_coords)
  {
    // get the rotation matrix associated with the geom object
    const DenseMatrix<double> rotation_matrix = disk_geom_obj_pt->rotation_matrix();
    
    // compute the inverse of the rotation matrix
    DenseMatrix<double> inv_rotation_matrix(3, 3, 0.0);
    matrix_inverse(rotation_matrix, inv_rotation_matrix);

    // rotate the Eulerian position back as if the disk was tangent to the x-y plane
    // in order to get a good starting guess for the Lagrangian coordinates
    Vector<double> x_inv_rot(3, 0.0);
    for(unsigned i=0; i<3; i++)
    {
      for(unsigned j=0; j<3; j++)
      {
	x_inv_rot[i] += inv_rotation_matrix(i,j) * x[j];
      }
    }
    
    // the (flat) cylindrical radius of this point to use as an initial guess
    double r0 = sqrt(x_inv_rot[0]*x_inv_rot[0] + x_inv_rot[1]*x_inv_rot[1]);
  
    Vector<double> unknowns(3, 0.0);

    // starting guesses are the Lagrangian coordinates for a flat disk
    unknowns[0] = r0;                                // xi_1
    unknowns[1] = atan2(x_inv_rot[1], x_inv_rot[0]); // xi_2
    unknowns[2] = x_inv_rot[2];                      // xi_3

    // keep track of the last good values we had
    Vector<double> previous_converged_unknowns = unknowns;
      
    // back up the actual disk object
    CylindricallyWarpedCircularDiskWithAnnularInternalBoundary* warped_disk_backup_pt =
      disk_geom_obj_pt;
  
    // disk parameters we're keeping fixed
    const double target_r_curvature = warped_disk_backup_pt->radius_of_curvature();
    const double r_torus = warped_disk_backup_pt->h_annulus();
    
    // number of sequential attempts without convergence before
    // we give up and decide something more fundamental has gone wrong
    const unsigned max_attempts = 4;

    // tolerance on floating-point comparisons
    const double tol = 1e-6;

    double current_r_curvature = target_r_curvature;
    double d_r = target_r_curvature;
    unsigned count = 0;
      
    while (true)
    {
      // make a temporary copy to work with
      auto temp_disk_pt = std::make_unique<
	CylindricallyWarpedCircularDiskWithAnnularInternalBoundary>
	(r_torus, current_r_curvature, rotation_matrix);

      // slightly hacky, but grab the and assign the raw pointer, since the
      // black-box residual fct uses the raw disk_geom_obj_pt pointer
      disk_geom_obj_pt = temp_disk_pt.get();
	
      // hit it with Newton's method with a finite-diff'd Jacobian
      try
      {
	BlackBoxFDNewtonSolver::black_box_fd_newton_solve(
	  &eulerian_to_lagrangian_coord_residual, x, unknowns);

	// if we've got a converged solution at the target epsilon, we're done
	if(abs(target_r_curvature - current_r_curvature) < tol)
	  break;
	else
	{
	  // reset the count
	  count = 0;
	  
	  // otherwise increase the current value and go again
	  current_r_curvature += d_r;

	  // store this solution
	  previous_converged_unknowns = unknowns;
	}
      }
      catch(OomphLibError& error) // didn't converge, initial guess was too poor
      {
	// I've caught the error so shut up!
	error.disable_error_message();
	
	// increment the unconverged solve counter
	count++;

	// if we've hit the limit of attempts, then bail out
	if(count == max_attempts)
	  break;
	
	// otherwise increase the current radius of curvature and go again
	d_r *= 2;
	
	current_r_curvature += d_r;

	// reset the initial guess to the last converged solution
	unknowns = previous_converged_unknowns;
      }
    }
    
    if(count == max_attempts)
    {
      std::ostringstream error_message;
      error_message << "Couldn't find zeta for the bulk point ("
		    << x[0] << ", " << x[1] << ", " << x[2] << ")\n\n";

      throw OomphLibError(error_message.str(),
			  OOMPH_CURRENT_FUNCTION,
			  OOMPH_EXCEPTION_LOCATION);
    }

    // did we get a solution with xi1 < 0?
    if(unknowns[0] < 0)
    {  
      // reset the raw disk pointer
      disk_geom_obj_pt = warped_disk_backup_pt;
	
      std::ostringstream error_message;
      error_message << "Lagrangian coordinate solve returned a negative radial coordinate\n";

      throw OomphLibError(error_message.str(),
			  OOMPH_CURRENT_FUNCTION,
			  OOMPH_EXCEPTION_LOCATION);
    }
    
    // interpret the solve and convert 'azimuthal' angle to the appropriate range
    lagr_coords.xi1 = unknowns[0];
    lagr_coords.xi2 = map_angle_to_range_0_to_2pi(unknowns[1]);
    lagr_coords.xi3 = unknowns[2];

    // reset the raw disk pointer
    disk_geom_obj_pt = warped_disk_backup_pt;
  }
  
  // ==========================================================================
  /// \short Function which computes the derivatives of the Eulerian coordinates
  /// w.r.t. the curviliear coordinates
  // ==========================================================================
  DenseMatrix<double> dx_dxi(const LagrangianCoordinates& lagr_coords)
  {
    // basis vectors
    Vector<double> a1(3, 0.0),
      a2 = a1,
      a3 = a1;

    // derivatives of basis vectors w.r.t. curvilinear coords
    DenseMatrix<double> da1_dxi(3, 2, 0.0),
      da2_dxi = da1_dxi,
      da3_dxi = da1_dxi;

    // derivatives of the disk position vector
    DenseMatrix<double> dr_dxi(3, 2, 0.0);
    
    // make vector for easier indexing
    Vector<double> xi(3, 0.0);
    xi[0] = lagr_coords.xi1;
    xi[1] = lagr_coords.xi2;
    xi[2] = lagr_coords.xi3;
    
    // get stuff from the geom object
    disk_geom_obj_pt->basis_vectors(lagr_coords, a1, a2, a3);
    disk_geom_obj_pt->da_dxi(lagr_coords, da1_dxi, da2_dxi, da3_dxi);
    disk_geom_obj_pt->dposition_dxi(lagr_coords, dr_dxi);

    // make vector for easier indexing
    Vector<DenseMatrix<double>> da_dxi(3);
    da_dxi[0] = da1_dxi;
    da_dxi[1] = da2_dxi;
    da_dxi[2] = da3_dxi;
    
    // output matrix
    DenseMatrix<double> dxdxi(3, 3, 0.0);

    // -------------------------------------
    // compute the derivatives of the Cartesian coordinates w.r.t. the
    // curvilinear coordinates, i.e. the derivatives of:
    //
    // x(xi1,xi2,xi3) = r(xi1,xi2) + xi3 * a_3 (xi1,xi2)
    //
    
    // loop over the Cartesian coordinates
    for(unsigned i=0; i<3; i++)
    {
      // loop over the curvilinear coordinates
      for(unsigned j=0; j<3; j++)
      {
	if(j == 2)
	{
	  dxdxi(i,j) += a3[i];
	}
	else
	{
	  dxdxi(i,j) += dr_dxi(i,j) + lagr_coords.xi3 * da3_dxi(i,j);
	}
      }
    } 

    return dxdxi;
  }

  DenseMatrix<double> dxi_dx(const LagrangianCoordinates& lagr_coords)
  {
    // get the derivatives dx_i/dxi_j
    DenseMatrix<double> dxdxi = dx_dxi(lagr_coords);

    // Calculate the determinant of the matrix
    const double det = dxdxi(0,0) * dxdxi(1,1) * dxdxi(2,2) 
      + dxdxi(0,1) * dxdxi(1,2) * dxdxi(2,0) 
      + dxdxi(0,2) * dxdxi(1,0) * dxdxi(2,1) 
      - dxdxi(0,0) * dxdxi(1,2) * dxdxi(2,1) 
      - dxdxi(0,1) * dxdxi(1,0) * dxdxi(2,2) 
      - dxdxi(0,2) * dxdxi(1,1) * dxdxi(2,0); 
    
    DenseMatrix<double> dxidx(3, 3, 0.0);
    
    //Calculate the inverse of the Lagrangian derivatives of the Eulerian coords
    dxidx(0,0) =  (dxdxi(1,1) * dxdxi(2,2) 
		   - dxdxi(1,2) * dxdxi(2,1)) / det; 
    dxidx(0,1) = -(dxdxi(0,1) * dxdxi(2,2) 
		   - dxdxi(0,2) * dxdxi(2,1)) / det; 
    dxidx(0,2) =  (dxdxi(0,1) * dxdxi(1,2) 
		   - dxdxi(0,2) * dxdxi(1,1)) / det; 
    dxidx(1,0) = -(dxdxi(1,0) * dxdxi(2,2) 
		   - dxdxi(1,2) * dxdxi(2,0)) / det; 
    dxidx(1,1) =  (dxdxi(0,0) * dxdxi(2,2) 
		   - dxdxi(0,2) * dxdxi(2,0)) / det; 
    dxidx(1,2) = -(dxdxi(0,0) * dxdxi(1,2) 
		   - dxdxi(0,2) * dxdxi(1,0)) / det; 
    dxidx(2,0) =  (dxdxi(1,0) * dxdxi(2,1) 
		   - dxdxi(1,1) * dxdxi(2,0)) / det; 
    dxidx(2,1) = -(dxdxi(0,0) * dxdxi(2,1) 
		   - dxdxi(0,1) * dxdxi(2,0)) / det; 
    dxidx(2,2) =  (dxdxi(0,0) * dxdxi(1,1) 
		   - dxdxi(0,1) * dxdxi(1,0)) / det;
    
    return dxidx;
  }
  
  // ==========================================================================
  /// Function which computes dzeta/dx
  // ==========================================================================
  Vector<double> dzeta_dx(const LagrangianCoordinates& lagr_coords)
  {
    // get the derivatives dxi_i/dx_j
    DenseMatrix<double> dxidx = dxi_dx(lagr_coords);
    
    // now extract the zeta (= xi_2) derivatives
    Vector<double> dzetadx(3, 0.0);

    for(unsigned i=0; i<3; i++)
      dzetadx[i] = dxidx(1,i);
    
    return dzetadx;
  }
}

#endif
