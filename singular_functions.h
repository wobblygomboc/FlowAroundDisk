#ifndef OOMPH_SINGULAR_FUNCTIONS_HEADER
#define OOMPH_SINGULAR_FUNCTIONS_HEADER

// the analytic solution for flow around an edge
#include "moffatt_solution.h"

namespace SingularFunctions
{
  double PI = MathematicalConstants::Pi;
  
  void dzeta_dx_residual(const Vector<double>& parameters,
			 const Vector<double>& unknowns,
			 Vector<double>& residuals)
  {    
    // interpret parameters
    // ---------------------

    // spatial coordinates of the bulk point
    Vector<double> x(3,0);
    x[0] = parameters[0];
    x[1] = parameters[1];
    x[2] = parameters[2];

    double zeta = parameters[3];
          
    Vector<double> dzeta_dx(3,0);
    dzeta_dx[0] = unknowns[0];
    dzeta_dx[1] = unknowns[1];
    dzeta_dx[2] = unknowns[2];

    double n = Global_Parameters::n;
    double epsilon = Global_Parameters::Epsilon;
    
    double w          =       epsilon * cos(n*zeta);
    double dw_dzeta   =    -n*epsilon * sin(n*zeta);

    Vector<double> r_disk_edge(3);
    Vector<double> tangent(3);
    Vector<double> surface_normal(3);
    Vector<double> normal(3);
    
    // get the triad vectors from the disk-like geometric object at this zeta
    Global_Parameters::Warped_disk_with_boundary_pt->
      boundary_triad(0, zeta, r_disk_edge, tangent, normal, surface_normal);
    
    Vector<double> dtangent_dzeta(3, 0);
    Vector<double> dnormal_dzeta(3, 0);
    Vector<double> dbinormal_dzeta(3, 0);
      
    Global_Parameters::Warped_disk_with_boundary_pt->
      dboundary_triad_dzeta(0, zeta, dtangent_dzeta, dnormal_dzeta, dbinormal_dzeta);
    
    residuals.resize(3);

    // residual equation is the derivative of the equation for the normal
    // distance of the bulk point from the plane at this zeta, i.e.
    // d/dx ((x-r)*tangent)
    for(unsigned i=0; i<residuals.size(); i++)
    {                 
      residuals[i] = (::Analytic_Functions::delta(0,i) + sin(zeta)*dzeta_dx[i])*tangent[0]
	+ (::Analytic_Functions::delta(1,i) - cos(zeta)*dzeta_dx[i])*tangent[1]
	+ (::Analytic_Functions::delta(2,i) - dw_dzeta*dzeta_dx[i])*tangent[2]
	+ (x[0] - cos(zeta))*dtangent_dzeta[0]*dzeta_dx[i]
	+ (x[1] - sin(zeta))*dtangent_dzeta[1]*dzeta_dx[i]
	+ (x[2] - w)*dtangent_dzeta[2]*dzeta_dx[i];
    }
  }

  // ==========================================================================

  // get the Lagrangian triad unit vectors \hat\rho, \hat\zeta and \hat\phi
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
    mVector tangent(3);
    mVector binormal(3);
    mVector normal(3);
    unsigned b_dummy = 0;
    
    // get the unit normals from the disk-like geometric object at this zeta
    Global_Parameters::Warped_disk_with_boundary_pt->
      boundary_triad(b_dummy, edge_coords.zeta, x_disk_edge, tangent,
		     normal, binormal);
    
    // unit vector in the rho direction    
    for(unsigned i=0; i<3; i++)
      rho_hat[i] = -normal[i] * cos(edge_coords.phi) + binormal[i] * sin(edge_coords.phi);

    // the GeomObj::boundary_triad() gives back unit vectors, so \hat\zeta is
    // just the tangent vector
    zeta_hat = tangent;
    
    // the - sign in front of binormal component has been cancelled out by the
    // angle flip, since cos(pi-x) = -cos(x)
    for(unsigned i=0; i<3; i++)
      phi_hat[i] = normal[i] * sin(edge_coords.phi) + binormal[i] * cos(edge_coords.phi);
  }
  
  // Converts velocities given in (\rho,\zeta,\phi) components to
  // global Cartesian (x,y,z) coordinates
  void lagrangian_to_eulerian_velocity(const EdgeCoordinates& edge_coords,
				       const Vector<double>& u_lagrangian,
				       Vector<double>& u_cartesian)
  {
    Vector<double> rho_hat;
    Vector<double> zeta_hat;
    Vector<double> phi_hat;
    Vector<double> x_disk_edge;
    
    // get the Lagrangian unit vectors at this point
    lagrangian_triad_and_x_disk_edge(edge_coords, rho_hat, zeta_hat,
				     phi_hat, x_disk_edge);
       
    // Matrix of Cartesian unit vectors
    Vector<mVector> e_hat(3);
    e_hat[0] = mVector::e_x();
    e_hat[1] = mVector::e_y();
    e_hat[2] = mVector::e_z();

    // Matrix of Lagrangian unit vectors
    Vector<mVector> xi(3);
    xi[0] = rho_hat;
    xi[1] = zeta_hat;
    xi[2] = phi_hat;

    // assemble the transformation matrix U_{ij} = e_hat_{ik} xi_{jk}
    DenseMatrix<double> U(3,3, 0.0);
    for(unsigned i=0; i<3; i++)
    {
      for(unsigned j=0; j<3; j++)
      {	
	// do the dot product
	for(unsigned k=0; k<3; k++)
	  U(i,j) += e_hat[i][k] * xi[j][k];
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
 

  /// \short Function to convert derivatives of Lagrangian velocity components
  /// w.r.t. Lagrangian coordinates to Cartesian (Eulerian) velocity gradient tensor
  void lagrangian_to_eulerian_velocity_gradient(const EdgeCoordinates& edge_coords,
						const Vector<double>& u_lagrangian,
						const DenseMatrix<double>& u_lagrangian_derivs,
						DenseMatrix<double>& du_dx)
  {
    // shorthands
    double rho  = edge_coords.rho;
    double zeta = edge_coords.zeta;
    double phi  = edge_coords.phi;

    // Lagrangian triad unit vectors
    Vector<double> rho_hat(3, 0.0);
    Vector<double> zeta_hat(3, 0.0);
    Vector<double> phi_hat(3, 0.0);

    Vector<double> x_disk_edge_dummy(3, 0.0);
    
    // get 'em
    lagrangian_triad_and_x_disk_edge(edge_coords, rho_hat, zeta_hat,
				     phi_hat, x_disk_edge_dummy);

    // compute the rho vector as the magnitude rho times the unit vector
    Vector<double> rho_vector(3, 0.0);
    for(unsigned i=0; i<3; i++)
      rho_vector[i] = rho * rho_hat[i];

    // get the Cartesian coordinates of this point
    Vector<double> x(3, 0.0);
    lagrangian_to_eulerian_coordinates(edge_coords, x);
    
    // Eulerian derivatives of the boundary triad vectors
    DenseMatrix<double> dtangent_dx(3,3,0.0);
    DenseMatrix<double> dnormal_dx(3,3,0.0);
    DenseMatrix<double> dbinormal_dx(3,3,0.0);	

    Vector<double> dtangent_dzeta(3, 0);
    Vector<double> dnormal_dzeta(3, 0);
    Vector<double> dbinormal_dzeta(3, 0);

    // get the triad derivatives
    Global_Parameters::Warped_disk_with_boundary_pt->
      dboundary_triad_dzeta(0, edge_coords.zeta, dtangent_dzeta,
			    dnormal_dzeta, dbinormal_dzeta);
   
    // starting guess for dzeta_dx is the dzeta_dx for a flat disk
    Vector<double> unknowns(3, 0);
    unknowns[0] = -sin(zeta) / (sqrt(x[0]*x[0] + x[1]*x[1]));
    unknowns[1] =  cos(zeta) / (sqrt(x[0]*x[0] + x[1]*x[1]));
    unknowns[2] = 0;
      
    Vector<double> parameters(4);
    parameters[0] = x[0];
    parameters[1] = x[1];
    parameters[2] = x[2];
    parameters[3] = zeta;
      
    // do the solve to get the boundary zeta
    try
    {
      BlackBoxFDNewtonSolver::black_box_fd_newton_solve(
	&dzeta_dx_residual, parameters, unknowns);
    }
    catch(const std::exception e)
    {
      std::ostringstream error_message;
      error_message << "Couldn't find dzeta_dx for the bulk point ("
		    << x[0] << ", " << x[1] << ", " << x[2] << ")\n\n";

      throw OomphLibError(error_message.str(),
			  OOMPH_CURRENT_FUNCTION,
			  OOMPH_EXCEPTION_LOCATION);
    }

    Vector<double> dzeta_dx(3,0);
      
    // interpret the solve
    dzeta_dx[0] = unknowns[0];
    dzeta_dx[1] = unknowns[1];
    dzeta_dx[2] = unknowns[2];
      
    for(unsigned i=0; i<3; i++)
    {
      for(unsigned j=0; j<3; j++)
      {
	dtangent_dx(i,j)  = dtangent_dzeta[i]  * dzeta_dx[j];
	dnormal_dx(i,j)   = dnormal_dzeta[i]   * dzeta_dx[j];
	dbinormal_dx(i,j) = dbinormal_dzeta[i] * dzeta_dx[j];
      }
    }
      
    // polar derivatives w.r.t. (rho, zeta, phi)
    Vector<double> du_lagrangian_drho(3,0);
    du_lagrangian_drho[0] = u_lagrangian_derivs(0,0);   // du_rho/drho
    du_lagrangian_drho[1] = u_lagrangian_derivs(1,0);   // du_zeta/drho
    du_lagrangian_drho[2] = u_lagrangian_derivs(2,0);   // du_phi/drho

    // tangent derivatives 
    Vector<double> du_lagrangian_dzeta(3,0);
    du_lagrangian_dzeta[0] = u_lagrangian_derivs(0,1);   // du_rho/dzeta
    du_lagrangian_dzeta[1] = u_lagrangian_derivs(1,1);   // du_zeta/dzeta
    du_lagrangian_dzeta[2] = u_lagrangian_derivs(2,1);   // du_phi/dzeta
    
    Vector<double> du_lagrangian_dphi(3,0);      
    du_lagrangian_dphi[0] = u_lagrangian_derivs(0,2);  // du_rho/dphi
    du_lagrangian_dphi[0] = u_lagrangian_derivs(1,2);  // du_zeta/dphi
    du_lagrangian_dphi[2] = u_lagrangian_derivs(2,2);  // du_phi/dphi

    // derivatives of the Lagrangian unit vectors w.r.t. Cartesian directions
    DenseMatrix<double> drho_hat_dx(3,3,0);
    DenseMatrix<double> dzeta_hat_dx(3,3,0);    // QUEHACERES use this!
    DenseMatrix<double> dphi_hat_dx(3,3,0);
      
    double k = Global_Parameters::n;
    double dw_dzeta = -Global_Parameters::Epsilon * k * sin(k * zeta);
	
    DenseMatrix<double> drho_vector_dx(3,3,0);

    // ------------------------------------------------
    // QUEHACERES this is only true for a disk with x-y radius 1 defined by w=w(r,zeta),
    // won't be true once it starts deforming - need to get interpolated
    // dw_dzeta and dzeta_dx derivatives from the disk elements
    
    // drho_x_dx
    drho_vector_dx(0,0) = 1 + sin(zeta)*dzeta_dx[0];
    // drho_x_dy
    drho_vector_dx(0,1) = 0 + sin(zeta)*dzeta_dx[1];
    // drho_x_dz
    drho_vector_dx(0,2) = 0 + sin(zeta)*dzeta_dx[2];
    // drho_y_dx
    drho_vector_dx(1,0) = 0 - cos(zeta)*dzeta_dx[0];
    // drho_y_dy
    drho_vector_dx(1,1) = 1 - cos(zeta)*dzeta_dx[1];
    // drho_y_dz
    drho_vector_dx(1,2) = 0 - cos(zeta)*dzeta_dx[2];
    // drho_z_dx
    drho_vector_dx(2,0) = 0 - dw_dzeta * dzeta_dx[0];
    // drho_z_dy
    drho_vector_dx(2,1) = 0 - dw_dzeta * dzeta_dx[1];
    // drho_z_dz
    drho_vector_dx(2,2) = 1 - dw_dzeta * dzeta_dx[2];

    // boundary triad (functions of \zeta) in Cartesian coordinates
    mVector tangent(3);
    mVector binormal(3);
    mVector normal(3);
    unsigned b_dummy = 0;
    
    // get the unit triad vectors from the disk-like geometric object at this zeta
    Global_Parameters::Warped_disk_with_boundary_pt->
      boundary_triad(b_dummy, edge_coords.zeta, x_disk_edge_dummy, tangent,
		     normal, binormal);

    // coordinates of this point in the n-s coordinate system
    double n = 0;
    double s = 0;

    // derivatives w.r.t. the global Cartesian system
    mVector ds_dx(3,0);
    mVector dn_dx(3,0);
    
    for(unsigned i=0; i<3; i++)
    {
      n += rho_vector[i] * binormal[i];
      s += rho_vector[i] * normal[i];
	
      for(unsigned j=0; j<3; j++)
      {
	ds_dx[i] += drho_vector_dx(j,i)*normal[j]   + rho_vector[j]*dnormal_dx(j,i);
	dn_dx[i] += drho_vector_dx(j,i)*binormal[j] + rho_vector[j]*dbinormal_dx(j,i);
      }
    }
      
    // Cartesian derivatives of the Moffatt coordinates
    mVector drho_dx(3,0);
    mVector dphi_dx(3,0);

    for(unsigned i=0; i<3; i++)
    {
      drho_dx[i] = (s*ds_dx[i] + n*dn_dx[i]) / rho;
      dphi_dx[i] = (n*ds_dx[i] - s*dn_dx[i]) / (rho*rho);
    } 
      
    for(unsigned i=0; i<3; i++)
    {
      for(unsigned j=0; j<3; j++)
      {
	drho_hat_dx(i,j) = dphi_dx[j]*phi_hat[i] - cos(phi)*dnormal_dx(i,j) +  
	  sin(phi)*dbinormal_dx(i,j);

	dphi_hat_dx(i,j) = -dphi_dx[j]*rho_hat[i] + sin(phi)*dnormal_dx(i,j) + 
	  cos(phi)*dbinormal_dx(i,j);

	// zeta is the coordinate in the tangent direction, so these are the same
	dzeta_hat_dx(i,j) = dtangent_dx(i,j);
      }
    }

    // ----------------------------------------------------------------
      
    // the Cartesian derivatives of the Moffatt vectors
    Vector<DenseMatrix<double> > dxi_dx(3);
    for(unsigned k=0; k<3; k++)
    {
      dxi_dx[k].resize(3,3);
    }
    for(unsigned i=0; i<3; i++)
    {
      for(unsigned j=0; j<3; j++)
      {
	dxi_dx[0](i,j) = drho_hat_dx(i,j);
	dxi_dx[1](i,j) = dzeta_hat_dx(i,j);
	dxi_dx[2](i,j) = dphi_hat_dx(i,j);
      }
    }

    // Matrix of Cartesian unit vectors
    Vector<mVector> e_hat(3);
    e_hat[0] = mVector::e_x();
    e_hat[1] = mVector::e_y();
    e_hat[2] = mVector::e_z();

    // Matrix of Lagrangian unit vectors
    Vector<mVector> xi(3);
    xi[0] = rho_hat;
    xi[1] = zeta_hat;
    xi[2] = phi_hat;
    
    // loop over the Cartesian components (row)
    for(unsigned i=0; i<3; i++)
    {
      // loop over the Cartesian components (column)
      for(unsigned j=0; j<3; j++)
      {	  
	// loop over the Moffatt vectors
	for(unsigned k=0; k<3; k++)
	{
	  // do the dot product
	  for(unsigned l=0; l<3; l++)
	  {
	    du_dx(i,j) += e_hat[i][l] * dxi_dx[k](l,j) * u_lagrangian[k] +
	      e_hat[i][l]*xi[k][l] * (du_lagrangian_drho[k]  * drho_dx[j] +
				      du_lagrangian_dzeta[k] * dzeta_dx[j] + 
				      du_lagrangian_dphi[k]  * dphi_dx[j]);
	  }
	}
      }
    }
  }
  
  
  // //////////////////////////////////////////////////////////////////////////
  // Main function to compute the singular function and gradient, independent
  // of broadside or in-plane modes - simply takes two constants A and B and
  // computes the Moffatt solution and it's Cartesian velocity gradients.
  // is_lower_disk_element useful for outputting the solution, as it allows for
  // getting the pressure jump correct across the plate
  void singular_fct_and_gradient_moffatt(const EdgeCoordinates& edge_coords,
					 const double& A, const double& B,
					 Vector<double>& u_cartesian,
					 DenseMatrix<double>& du_dx)
  {
    // how far to shift the radius to compute "infinity" when the actual
    // radius is zero
    const double dr = 1e-5;
    
    // tolerance for some floating-point comparisons
    const double tol = 1e-8;

    // shorthands
    double rho  = edge_coords.rho;
    double zeta = edge_coords.zeta;
    double phi  = edge_coords.phi;

    // catch the case where we're sat exactly at the singularity
    if(rho < tol/10.0)
      rho = dr;
    
    // polar derivatives of polar velocity components,
    // i.e. dur_dr, dur_dphi, duphi_dr, duphi_dphi    
    DenseMatrix<double> u_moffatt_derivatives(2,2);

    // get the 2D polar Moffat solution (third component is pressure)
    mVector u_moffatt(3);
    moffatt_solution(rho, phi, A, B, u_moffatt, u_moffatt_derivatives);
	
    // ----------------
    // now use the outer normal to convert the rzp velocity into Cartesians

    u_cartesian.resize(3);
    u_cartesian.initialise(0.0);
    
    // cartesian derivatives
    du_dx.resize(3,3);
    du_dx.initialise(0.0);

    mVector u_lagrangian(3);
    u_lagrangian[0] = u_moffatt[0];  // u_rho
    u_lagrangian[1] = 0;             // u_zeta
    u_lagrangian[2] = u_moffatt[1];  // u_phi
 
    // do the conversion
    lagrangian_to_eulerian_velocity(edge_coords, u_lagrangian, u_cartesian);

    // and add the pressure, which is a scalar so no vector conversions
    u_cartesian.push_back(u_moffatt[2]);
      
    // now get the Eulerian velocity gradient tensor
    // -------------------------------------

    DenseMatrix<double> u_lagrangian_derivs(3,3, 0.0);

    u_lagrangian_derivs(0,0) = u_moffatt_derivatives(0,0);  // du_rho/drho
    u_lagrangian_derivs(0,1) = 0.0;                       // du_rho/dzeta
    u_lagrangian_derivs(0,2) = u_moffatt_derivatives(0,1);  // du_rho/dphi

    // u_zeta = 0 so no du_zeta derivatives
      
    u_lagrangian_derivs(2,0) = u_moffatt_derivatives(1,0);  // du_phi/drho
    u_lagrangian_derivs(2,1) = 0.0;                       // du_phi/dzeta
    u_lagrangian_derivs(2,2) = u_moffatt_derivatives(1,1);  // du_phi/dphi
	
    lagrangian_to_eulerian_velocity_gradient(edge_coords,
					     u_lagrangian,
					     u_lagrangian_derivs,
					     du_dx);    
  }

  void gradient_of_singular_fct_moffatt_finite_diff(const EdgeCoordinates& edge_coords,
						    const double& A, const double& B,
						    DenseMatrix<double>& du_dx)
  {
    // QUEHACERES for debug, do a finite diff (this only works for flat disk)
    // disable but keep code alive as we'll have to modify the derivative stuff
    // when we've got disk elements providing interpolated derivatives so this will
    // be a useful to check

    // coordinate increment for finite-difference gradient
    const double fd_dx = 1e-6;

    // get the Cartesian coords of this point
    Vector<double> x(3,0.0);
    lagrangian_to_eulerian_coordinates(edge_coords, x);

    DenseMatrix<double> du_dx_dummy(3,3,0.0);
      
    // get the solution at the current point
    Vector<double> u0(4, 0.0);
    singular_fct_and_gradient_moffatt(edge_coords, A, B, u0, du_dx_dummy);
      
    for(unsigned j=0; j<3; j++)
    {
      // coordinate vectors with increments applied to their components
      mVector x_plus_dx = x;

      // add the increments in each direction
      x_plus_dx[j] += fd_dx;

      double zeta_dx = atan2pi(x_plus_dx[1], x_plus_dx[0]);
      mVector x_disk_edge_dx(3, 0.0);
      mVector tangent_dx(3, 0.0);
      mVector normal_dx(3, 0.0);
      mVector binormal_dx(3, 0.0);
      unsigned b_dummy = 0;
	
      // get the unit normal from the disk-like geometric object at this zeta
      Global_Parameters::Warped_disk_with_boundary_pt->
	boundary_triad(b_dummy, zeta_dx, x_disk_edge_dx, tangent_dx,
		       normal_dx, binormal_dx);

      // convert back to polars
      mVector rho_vector_dx = x_plus_dx - x_disk_edge_dx;
      double rho_dx = rho_vector_dx.magnitude();
	  
      double phi_dx = atan2pi(rho_vector_dx*binormal_dx, -rho_vector_dx*normal_dx);

      // get the solution at the incremented position
      Vector<double> u_plus_dx(4, 0.0);
      singular_fct_and_gradient_moffatt(edge_coords, A, B, u_plus_dx, du_dx_dummy);
	
      // now do the finite diff to get the velocity gradients
      for(unsigned i=0; i<3; i++)
	du_dx(i,j) = (u_plus_dx[i] - u0[i]) / fd_dx;
    }
  }
  
  void singular_fct_and_gradient_broadside(const EdgeCoordinates& edge_coords,
					   Vector<double>& u,
					   DenseMatrix<double>& du_dx)
  {
    // parameters for broadside motion
    // (minus sign accounts for flow driven by disk motion rather
    // than the far-field conditions used in the Moffatt derivation)
    double A = 0;
    double B = -Global_Parameters::u_disk_rigid_body[2];

    // forward
    singular_fct_and_gradient_moffatt(edge_coords, A, B, u, du_dx);

    // QUEHACERES for debug - generalise this @@@@@@@@@@@
    // add the pressure zero-level to the pressure
    u[3] += Global_Parameters::p0;
  }
  
  void singular_fct_and_gradient_in_plane(const EdgeCoordinates& edge_coords,
					  Vector<double>& u,
					  DenseMatrix<double>& du_dx)
  {
    // parameters for in-plane motion
    // (minus sign accounts for flow driven by disk motion rather
    // than the far-field conditions used in the Moffatt derivation)
    double A = -Global_Parameters::u_disk_rigid_body[0];
    double B = 0; 

    // forward
    singular_fct_and_gradient_moffatt(edge_coords, A, B, u, du_dx);

    // QUEHACERES if this works, think! then write a good comment and generalise
    u[1] = -u[1];

    // QUEHACERES @@@@@@@@@@@@@@@@@@
    u[3] += edge_coords.p0;
  }

  // Asymptotic singular solution for the rotation of a disk about it's axis of
  // symmetry, i.e. +ve z-axis rotation
  Vector<double> singular_fct_in_plane_rotation(const EdgeCoordinates& edge_coords)
  {
    Vector<double> u_lagr(4, 0.0);

    // this solution has u_rho = u_z = p = 0

    // u_zeta expanded to O(\rho^{1/2})
    u_lagr[1] = 1 - (4./PI) * sqrt(1 + cos(edge_coords.phi)) * edge_coords.rho;
    
    // scale by translational velocity (since this function is to be used
    // in combination with in-plane translation)
    u_lagr[1] *= Global_Parameters::u_disk_rigid_body[0];
      
    Vector<double> u_cartesian(3, 0.0);
    lagrangian_to_eulerian_velocity(edge_coords, u_lagr, u_cartesian);

    // and add the pressure
    u_cartesian.push_back(u_lagr[3]);
    
    return u_cartesian;
  }

  DenseMatrix<double>
    gradient_of_singular_fct_in_plane_rotation(const EdgeCoordinates& edge_coords)
  {
    // get the in-plane rotation solution in Lagrangian coords
    Vector<double> u_lagrangian = singular_fct_in_plane_rotation(edge_coords);

    // derivatives of the Lagrangian velocity components w.r.t. the
    // Lagrangian coordinates
    DenseMatrix<double> u_lagrangian_derivs(3, 3, 0.0);

    // only non-zero Lagrangian derivatives are du_zeta/drho and du_zeta/dphi

    // du_zeta/drho
    u_lagrangian_derivs(1,0) = -2 * sqrt(1 + cos(edge_coords.phi)) /
      (PI * sqrt(edge_coords.rho));

    // du_zeta/dphi
    u_lagrangian_derivs(1,2) = 2 * sqrt(edge_coords.rho) * sin(edge_coords.phi) /
      (PI * sqrt(1 + cos(edge_coords.phi)) );

    // Eulerian derivatives
    DenseMatrix<double> du_dx(3,3, 0.0);
    
    // now convert to Eulerian velocity gradient tensor
    lagrangian_to_eulerian_velocity_gradient(edge_coords, u_lagrangian,
					     u_lagrangian_derivs, du_dx);

    return du_dx;
  }
   
  DenseMatrix<double> gradient_of_singular_fct_broadside(const EdgeCoordinates& edge_coords)
  {
    // dummy solution vector
    Vector<double> u;

    // velocity gradient tensor
    DenseMatrix<double> du_dx;
    
    // forward
    singular_fct_and_gradient_broadside(edge_coords, u, du_dx);

    return du_dx;
  }
  
  DenseMatrix<double> gradient_of_singular_fct_in_plane(const EdgeCoordinates& edge_coords)
  {
    // dummy solution vector
    Vector<double> u;

    // velocity gradient tensor
    DenseMatrix<double> du_dx;
    
    // forward
    singular_fct_and_gradient_in_plane(edge_coords, u, du_dx);

    return du_dx;
  }
  Vector<double> singular_fct_broadside(const EdgeCoordinates& edge_coords)
  {
    // create a dummy gradient tensor
    DenseMatrix<double> du_dx;

    // solution vector
    Vector<double> u; 
    
    // forward 
    singular_fct_and_gradient_broadside(edge_coords, u, du_dx);
    
    return u;
  }

  Vector<double> singular_fct_in_plane(const EdgeCoordinates& edge_coords)
  {    
    // create a dummy gradient tensor
    DenseMatrix<double> du_dx;

    // solution vector
    Vector<double> u;
    
    // forward
    singular_fct_and_gradient_in_plane(edge_coords, u, du_dx);
    
    return u;
  }

  // @@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
  // for debug, functions which match the typedef for a singular function
  // to allow us to subtract the exact solution rather than the Moffatt solution
  
  Vector<double> wrapper_to_exact_solution_broadside(const EdgeCoordinates& edge_coords)
  {
    // convert edge coords to Cartesian
    Vector<double> x(3, 0.0);
    lagrangian_to_eulerian_coordinates(edge_coords, x);

    // broadside motion
    Vector<double> u_disk(3, 0.0);
    u_disk[2] = 1.0;

    Vector<double> omega_disk(3, 0.0);
    
    Vector<double> u(4, 0.0);
    FlatDiskExactSolutions::total_exact_solution(x, u_disk, omega_disk, u);
    
    return u;
  }

  DenseMatrix<double> wrapper_to_exact_velocity_gradient_broadside(const EdgeCoordinates& edge_coords)
  {
    // convert edge coords to Cartesian
    Vector<double> x(3, 0.0);
    lagrangian_to_eulerian_coordinates(edge_coords, x);

    // broadside motion
    Vector<double> u_disk(3, 0.0);
    u_disk[2] = 1.0;

    Vector<double> omega_disk(3, 0.0);
    
    DenseMatrix<double> du_dx(3, 3, 0.0);
    FlatDiskExactSolutions::total_exact_velocity_gradient(x, u_disk, omega_disk, du_dx);
    
    return du_dx;
  }
  
  Vector<double> wrapper_to_exact_solution_in_plane(const EdgeCoordinates& edge_coords)
  {
    // convert edge coords to Cartesian
    Vector<double> x(3, 0.0);
    lagrangian_to_eulerian_coordinates(edge_coords, x);

    // in-plane motion, uz = 0
    Vector<double> u_disk(3, 0.0);
    u_disk[0] = Global_Parameters::u_disk_rigid_body[0];
    u_disk[1] = Global_Parameters::u_disk_rigid_body[1];
    
    Vector<double> omega_disk(3, 0.0);
    
    Vector<double> u(4, 0.0);
    FlatDiskExactSolutions::total_exact_solution(x, u_disk, omega_disk, u);
    
    return u;
  }

  DenseMatrix<double> wrapper_to_exact_velocity_gradient_in_plane(const EdgeCoordinates& edge_coords)
  {
    // convert edge coords to Cartesian
    Vector<double> x(3, 0.0);
    lagrangian_to_eulerian_coordinates(edge_coords, x);

    // in-plane motion
    Vector<double> u_disk(3, 0.0);
    u_disk[0] = Global_Parameters::u_disk_rigid_body[0];
    u_disk[1] = Global_Parameters::u_disk_rigid_body[1];

    Vector<double> omega_disk(3, 0.0);
    
    DenseMatrix<double> du_dx(3, 3, 0.0);
    FlatDiskExactSolutions::total_exact_velocity_gradient(x, u_disk, omega_disk, du_dx);
    
    return du_dx;
  }
  
}

#endif
