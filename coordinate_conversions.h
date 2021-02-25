#ifndef OOMPH_COORDINATE_CONVERSIONS
#define OOMPH_COORDINATE_CONVERSIONS

#include "additional_maths.h"

namespace CoordinateConversions
{
  const double PI = MathematicalConstants::Pi;

  // QUEHACERES change this to something more generic at some point
  // pointer to a geometric object which provides things like
  // boundary triad vectors and their derivatives
  WarpedCircularDiskWithAnnularInternalBoundary* disk_geom_obj_pt = 0x0;
    
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
  // wrapper function, so we can easily switch over all the calls to the disk-like
  // GeomObject once we have the actual disk elements in there
  // ==========================================================================
  void boundary_triad(const double& zeta,
		      Vector<double>& x_disk_edge,
		      Vector<double>& tangent,
		      Vector<double>& normal,
		      Vector<double>& binormal)
  {
    // get the unit normals from the disk-like geometric object at this zeta
    disk_geom_obj_pt->
      boundary_triad(0, zeta, x_disk_edge, tangent, normal, binormal);
  }

  // ==========================================================================
  // wrapper function, so we can easily switch over all the calls to the disk-like
  // GeomObject once we have the actual disk elements in there
  // ==========================================================================
  void dboundary_triad_dzeta(double& zeta,
			     Vector<double>& dr_disk_edge_dzeta,
			     Vector<double>& dtangent_dzeta,
			     Vector<double>& dnormal_dzeta,
			     Vector<double>& dbinormal_dzeta)
  {
    disk_geom_obj_pt->dboundary_triad_dzeta(0, zeta, dr_disk_edge_dzeta,
					    dtangent_dzeta, dnormal_dzeta,
					    dbinormal_dzeta);
  }
  
  // ==========================================================================
  // wrapper function to get the radius of the disk at a given zeta
  // ==========================================================================
  double r_disk(const double& zeta)
  {
    // QUEHACERES for now, it's always 1, will need to change once
    // we've actually got the disk elements in there
    return 1.0;
  }
  
  // ==========================================================================
  // get the Lagrangian triad unit vectors \hat\rho, \hat\zeta and \hat\phi
  // ==========================================================================
  void lagrangian_triad_and_x_disk_edge(const EdgeCoordinates& edge_coords,
					Vector<double>& rho_hat,
					Vector<double>& zeta_hat,
					Vector<double>& phi_hat,
					Vector<double>& x_disk_edge)
  {
    // make enough space
    rho_hat.resize(3);
    zeta_hat.resize(3);
    phi_hat.resize(3);
    x_disk_edge.resize(3);

    // boundary triad (functions of \zeta) in Cartesian coordinates
    Vector<double> tangent(3, 0.0);
    Vector<double> normal(3, 0.0);
    Vector<double> binormal(3, 0.0);
    
    // get the triad unit vectors and disk edge at this zeta
    boundary_triad(edge_coords.zeta, x_disk_edge, tangent, normal, binormal);
    
    // unit vector in the rho direction    
    for(unsigned i=0; i<3; i++)
      rho_hat[i] = normal[i] * cos(edge_coords.phi) + binormal[i] * sin(edge_coords.phi);

    // the GeomObj::boundary_triad() gives back unit vectors, so \hat\zeta is
    // just the tangent vector
    zeta_hat = tangent;
    
    // the - sign in front of binormal component has been cancelled out by the
    // angle flip, since cos(pi-x) = -cos(x)
    for(unsigned i=0; i<3; i++)
      phi_hat[i] = normal[i] * sin(edge_coords.phi) + binormal[i] * cos(edge_coords.phi);
  }

  
  // ==========================================================================
  // Converts velocities given in arbitrary Lagrangian coordinates to
  // global Cartesian (x,y,z) coordinates
  // ==========================================================================
  void lagrangian_to_eulerian_velocity(const EdgeCoordinates& edge_coords,
				       const Vector<Vector<double> >& lagr_unit_vec,
				       const Vector<double>& u_lagrangian,				       
				       Vector<double>& u_cartesian)
  {
    // Matrix of Cartesian unit vectors
    Vector<Vector<double> > e_hat(3);
    e_hat[0] = cartesian_unit_vector(0);
    e_hat[1] = cartesian_unit_vector(1);
    e_hat[2] = cartesian_unit_vector(2);

    // assemble the transformation matrix U_{ij} = e_hat_{ik} xi_{jk}
    DenseMatrix<double> U(3,3, 0.0);
    for(unsigned i=0; i<3; i++)
    {
      for(unsigned j=0; j<3; j++)
      {	
	// do the dot product
	for(unsigned k=0; k<3; k++)
	  U(i,j) += e_hat[i][k] * lagr_unit_vec[j][k]; // xi[j][k];
      }
    }

    // zero everything out
    u_cartesian.resize(3);
    std::fill(u_cartesian.begin(), u_cartesian.end(), 0.0);
    
    // do the conversion u_i = U_ij u'_j
    for(unsigned i=0; i<3; i++)
    {
      for(unsigned j=0; j<3; j++)
      {
	u_cartesian[i] += U(i,j) * u_lagrangian[j];
      }
    }
  }

  // ==========================================================================
  // Conversion from the Lagrangian edge coordinates to the global Cartesians
  // ==========================================================================
  void lagrangian_to_eulerian_coordinates(const EdgeCoordinates& edge_coords,
					  Vector<double>& x)
  {
    // Lagrangian triad unit vectors
    Vector<double> rho_hat(3, 0.0);
    Vector<double> zeta_hat(3, 0.0);
    Vector<double> phi_hat(3, 0.0);

    Vector<double> x_disk_edge(3, 0.0);
    
    // get 'em
    lagrangian_triad_and_x_disk_edge(edge_coords, rho_hat, zeta_hat,
				     phi_hat, x_disk_edge);

    // make enough space
    x.resize(3);
    
    // compute the Eulerian coordinates of this point for the derivatives
    for(unsigned i=0; i<3; i++)
      x[i] = x_disk_edge[i] + edge_coords.rho * rho_hat[i];
  }

  // ==========================================================================
  // Compute the residual vector for an Eulerian to Lagrangian coordinate
  // transformation (with the correct prototype for a black-box FD Newton solve
  // ==========================================================================
  void eulerian_to_lagrangian_coord_residual(const Vector<double>& x,
					     const Vector<double>& unknowns,
					     Vector<double>& residuals)
  {
    residuals.resize(3, 0.0);

    Vector<double> r_disk_edge(3);
    Vector<double> tangent(3);
    Vector<double> binormal(3);
    Vector<double> normal(3);

    double rho  = unknowns[0];
    double zeta = unknowns[1];
    double phi  = unknowns[2];
  
    // get the unit normal from the disk-like geometric object at this zeta
    disk_geom_obj_pt->
      boundary_triad(0, zeta, r_disk_edge, tangent, normal, binormal);

    for(unsigned i=0; i<3; i++)
    {    
      residuals[i] = r_disk_edge[i] +
	rho * (cos(phi) * normal[i] + sin(phi)*binormal[i]) - x[i];
    }
  }

  // ==========================================================================
  // Perform a black-box finite-diff Newton solve to convert from eulerian
  // to Lagrangian edge coordinates
  // ==========================================================================
  void eulerian_to_lagrangian_coordinates(const Vector<double>& x,
					  EdgeCoordinates& edge_coords)
  {
    // the cylindrical radius of this point
    double r0 = sqrt(x[0]*x[0] + x[1]*x[1]) - r_disk(edge_coords.zeta);
  
    Vector<double> unknowns(3, 0.0);
  
    // starting guesses are the Lagrangian coordinates for a flat disk
    unknowns[0] = sqrt(pow(r0, 2) + pow(x[2], 2)); // rho
    unknowns[1] = atan2pi(x[1], x[0]);             // zeta
    unknowns[2] = atan2(x[2], r0);                 // phi

    // keep track of the last good values we had
    Vector<double> previous_converged_unknowns = unknowns;
      
    // back up the actual disk object
    WarpedCircularDiskWithAnnularInternalBoundary* warped_disk_backup_pt =
      disk_geom_obj_pt;

    // disk parameters we're keeping fixed
    const double target_epsilon = warped_disk_backup_pt->epsilon();
    const double r_torus        = warped_disk_backup_pt->h_annulus();
    const unsigned n            = warped_disk_backup_pt->n();

    // number of sequential attempts without convergence before
    // we give up and decide something more fundamental has gone wrong
    const unsigned max_attempts = 4;

    // tolerance on floating-point comparisons
    const double tol = 1e-6;
    
    double current_epsilon = target_epsilon;
    double d_eps = target_epsilon;
    unsigned count = 0;
      
    while (true)
    {
      // make a temporary copy to work with
      std::unique_ptr<WarpedCircularDiskWithAnnularInternalBoundary> temp_disk_pt =
	std::make_unique<WarpedCircularDiskWithAnnularInternalBoundary>(r_torus, 
									current_epsilon,
									n);

      // slightly hacky, but grab the and assign the raw pointer, since the
      // black-box residual fct uses the raw disk_geom_obj_pt pointer
      disk_geom_obj_pt = temp_disk_pt.get();
	
      // hit it with Newton's method with a finite-diff'd Jacobian
      try
      {
	BlackBoxFDNewtonSolver::black_box_fd_newton_solve(
	  &eulerian_to_lagrangian_coord_residual, x, unknowns);

	// if we've got a converged solution at the target epsilon, we're done
	if(abs(target_epsilon - current_epsilon) < tol)
	  break;
	else
	{
	  // reset the count
	  count = 0;
	  
	  // otherwise increase the current value and go again
	  current_epsilon += d_eps;

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
	
	// otherwise reduce the current value and go again
	d_eps *= 0.5;
	
	current_epsilon -= d_eps;

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

    // did we get a solution with rho < 0?
    if(unknowns[0] < 0)
    {
      // make a temporary copy to work with
      std::unique_ptr<WarpedCircularDiskWithAnnularInternalBoundary> temp_disk_pt =
	std::make_unique<WarpedCircularDiskWithAnnularInternalBoundary>(r_torus, 
									current_epsilon,
									n);

      // slightly hacky, but grab the and assign the raw pointer, since the
      // black-box residual fct uses the raw disk_geom_obj_pt pointer
      disk_geom_obj_pt = temp_disk_pt.get();
      
      // if so, then presumably flipping the angle by pi gives the same
      // solution with a positive rho - lets check
      unknowns[0] = -unknowns[0];
      unknowns[2] += MathematicalConstants::Pi;

      BlackBoxFDNewtonSolver::black_box_fd_newton_solve(
	  &eulerian_to_lagrangian_coord_residual, x, unknowns);

      if(unknowns[0] < 0)
      {
	// reset the raw disk pointer
	disk_geom_obj_pt = warped_disk_backup_pt;
	
	std::ostringstream error_message;
	error_message << "Something weird has happened - the initial Lagrangian "
		      << "coordinate solve returned a negative radial coordinate, "
		      << "so flipping the sign and reflecting the elevation angle "
		      << "should have resulted in a valid solution, but it didn't "
		      << "converge.\n\n";

	throw OomphLibError(error_message.str(),
			    OOMPH_CURRENT_FUNCTION,
			    OOMPH_EXCEPTION_LOCATION);
      }
    }
    
    // interpret the solve and convert angles to the appropriate range;
    // we're using [0,2pi] for zeta but [-pi:pi] for phi as the +/- jump
    // is a useful identifier for the top/bottom surface of the plate
    edge_coords.rho  = unknowns[0];
    edge_coords.zeta = map_angle_to_range_0_to_2pi(unknowns[1]);
    edge_coords.phi  = map_angle_to_range_plus_minus_pi(unknowns[2]);

    // reset the raw disk pointer
    disk_geom_obj_pt = warped_disk_backup_pt;
  }

  // ==========================================================================
  //
  // ==========================================================================
  void dx_dzeta_residual(const Vector<double>& edge_coords,
			 const Vector<double>& dx_dzeta,
			 Vector<double>& residuals)
  { 
    // interpret
    double rho  = edge_coords[0];
    double zeta = edge_coords[1];
    double phi  = edge_coords[2];
          
    // get the derivatives of the triad unit vectors
    Vector<double> dr_disk_edge_dzeta(3, 0);
    Vector<double> dtangent_dzeta(3, 0);
    Vector<double> dnormal_dzeta(3, 0);
    Vector<double> dbinormal_dzeta(3, 0);
      
    disk_geom_obj_pt->
      dboundary_triad_dzeta(0, zeta, dr_disk_edge_dzeta,
			    dtangent_dzeta, dnormal_dzeta, dbinormal_dzeta);
    
    residuals.resize(3);

    // residual equation is the derivative of the equation for the normal
    // distance of the bulk point from the plane at this zeta, i.e.
    // d/dx ((x-r)*tangent)
    for(unsigned i=0; i<residuals.size(); i++)
    {                 
      residuals[i] = dr_disk_edge_dzeta[i] +
	rho * (cos(phi) * dnormal_dzeta[i] + sin(phi)*dbinormal_dzeta[i]) - dx_dzeta[i];
    }
  }
  
  // ==========================================================================
  /// Function which computes dx/dzeta
  // ==========================================================================
  DenseMatrix<double> dx_dxi(const EdgeCoordinates& edge_coords)
  {
    // interpret
    double rho  = edge_coords.rho;
    double zeta = edge_coords.zeta;
    double phi  = edge_coords.phi;

    // boundary triad (functions of \zeta) in Cartesian coordinates
    Vector<double> x_disk_edge(3, 0.0);
    Vector<double> tangent(3, 0.0);
    Vector<double> normal(3, 0.0);
    Vector<double> binormal(3, 0.0);
        
    // get the triad unit vectors and disk edge at this zeta
    boundary_triad(edge_coords.zeta, x_disk_edge, tangent, normal, binormal);
    
    // get the derivatives of the triad unit vectors
    Vector<double> dr_disk_edge_dzeta(3, 0);
    Vector<double> dtangent_dzeta(3, 0);
    Vector<double> dnormal_dzeta(3, 0);
    Vector<double> dbinormal_dzeta(3, 0);     
    
    dboundary_triad_dzeta(zeta, dr_disk_edge_dzeta, dtangent_dzeta,
			  dnormal_dzeta, dbinormal_dzeta);
    
    DenseMatrix<double> dxdxi(3, 3, 0.0);

    // derivatives dx_i/dxi_j is the derivative of the equation for the cartesian
    // position, x_i = r_disk_i + rho(cos(phi)*s_i + sin(phi)*n_i)
    for(unsigned i=0; i<3; i++)
    {
      // dx_i/drho
      dxdxi(i,0) = cos(phi) * normal[i] + sin(phi)*binormal[i];
      
      // dx_i/dzeta
      dxdxi(i,1) = dr_disk_edge_dzeta[i] +
	rho * (cos(phi) * dnormal_dzeta[i] + sin(phi)*dbinormal_dzeta[i]);

      // dx_i/dphi
      dxdxi(i,2) = rho * (-sin(phi) * normal[i] + cos(phi)*binormal[i]);      
    } 

    return dxdxi;
  }

  DenseMatrix<double> dxi_dx(const EdgeCoordinates& edge_coords)
  {
    // get the derivatives dx_i/dxi_j
    DenseMatrix<double> dxdxi = dx_dxi(edge_coords);

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
  Vector<double> dzeta_dx(const EdgeCoordinates& edge_coords)
  {
    // QUEHACERES delete
    /* const double tol = 1e-8; */

    /* // get the derivatives dx_i/dxi_j */
    /* DenseMatrix<double> dxdxi = dx_dxi(edge_coords); */

    /* // Calculate the determinant of the matrix */
    /* const double det = dxdxi(0,0) * dxdxi(1,1) * dxdxi(2,2)  */
    /*   + dxdxi(0,1) * dxdxi(1,2) * dxdxi(2,0)  */
    /*   + dxdxi(0,2) * dxdxi(1,0) * dxdxi(2,1)  */
    /*   - dxdxi(0,0) * dxdxi(1,2) * dxdxi(2,1)  */
    /*   - dxdxi(0,1) * dxdxi(1,0) * dxdxi(2,2)  */
    /*   - dxdxi(0,2) * dxdxi(1,1) * dxdxi(2,0);  */

    // get the derivatives dxi_i/dx_j
    DenseMatrix<double> dxidx = dxi_dx(edge_coords);
    
    // now extract the zeta (= xi_2) derivatives
    Vector<double> dzetadx(3, 0.0);

    for(unsigned i=0; i<3; i++)
      dzetadx[i] = dxidx(1,i);
    
    return dzetadx;
  }
}

#endif
