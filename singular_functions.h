#ifndef OOMPH_SINGULAR_FUNCTIONS_HEADER
#define OOMPH_SINGULAR_FUNCTIONS_HEADER

// the analytic solution for flow around an edge
#include "moffatt_solution.h"

namespace SingularFunctions
{  
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
  
  // //////////////////////////////////////////////////////////////////////////
  // Main function to compute the singular function and gradient, independent
  // of broadside or in-plane modes - simply takes two constants A and B and
  // computes the Moffatt solution and it's Cartesian velocity gradients.
  // is_lower_disk_element useful for outputting the solution, as it allows for
  // getting the pressure jump correct across the plate
  void singular_fct_and_gradient(const EdgeCoordinates& edge_coords,
				 const double& A, const double& B,
				 Vector<double>& u_cartesian,
				 DenseMatrix<double>& du_dx)
  {
    const double infinity = 1003;

    // tolerance for some floating-point comparisons
    const double tol = 1e-8;

    double rho  = edge_coords.rho;
    double zeta = edge_coords.zeta;
    double phi  = edge_coords.phi;

    double b_dummy = 0;
    mVector x_disk_edge(3);
    mVector tangent(3);
    mVector binormal(3);
    mVector normal(3);
    
    // get the unit normal from the disk-like geometric object at this zeta
    Global_Parameters::Warped_disk_with_boundary_pt->
      boundary_triad(b_dummy, zeta, x_disk_edge, tangent,
		     normal, binormal);
        
    // unit vector in the rho direction
    mVector rho_hat = -normal * cos(phi) + binormal * sin(phi);

    // the - sign in front of binormal component has been cancelled out by the angle
    // flip, since cos(pi-x) = -cos(x)
    mVector phi_hat = normal * sin(phi) + binormal * cos(phi);

    // compute the rho vector, the vector from the edge of the disk at this
    // zeta to the point in question
    mVector rho_vector = rho_hat * rho;

    // compute the Eulerian coordinates of this point for the derivatives
    mVector x = x_disk_edge + rho_vector;    
    
    // polar derivatives of polar velocity components,
    // i.e. dur_dr, dur_dphi, duphi_dr, duphi_dphi
    DenseMatrix<double> u_polar_derivatives(2,2);

    // get the 2D polar Moffat solution (third component is pressure)
    mVector u_polar(3);
    moffatt_solution(rho, phi, A, B, u_polar, u_polar_derivatives);
	
    // ----------------
    // now use the outer normal to convert the rzp velocity into Cartesians

    u_cartesian.resize(3);

    // cartesian derivatives
    du_dx.resize(3,3);
    
    Vector<mVector> e_hat(3);
    e_hat[0] = mVector::e_x();
    e_hat[1] = mVector::e_y();
    e_hat[2] = mVector::e_z();

    Vector<mVector> xi(3);
    xi[0] = rho_hat;
    xi[1] = tangent;
    xi[2] = phi_hat;
    
    TransformationMatrix U(3,3);
    for(unsigned i=0; i<3; i++)
    {
      for(unsigned j=0; j<3; j++)
      {
	U(i,j) = e_hat[i] * xi[j];
      }
    }

    mVector u_moffat(3);
    u_moffat[0] = u_polar[0];
    u_moffat[1] = 0; 
    u_moffat[2] = u_polar[1];
 
    // do the conversion
    if(rho > tol)
    {
      // convert polar velocities to Cartesian
      u_cartesian = U*u_moffat;
    
      DenseMatrix<double> dtangent_dx(3,3,0.0);
      DenseMatrix<double> dnormal_dx(3,3,0.0);
      DenseMatrix<double> dbinormal_dx(3,3,0.0);	

      Vector<double> dtangent_dzeta(3, 0);
      Vector<double> dnormal_dzeta(3, 0);
      Vector<double> dbinormal_dzeta(3, 0);
      
      Global_Parameters::Warped_disk_with_boundary_pt->
	dboundary_triad_dzeta(0, zeta, dtangent_dzeta, dnormal_dzeta, dbinormal_dzeta);
   
      // starting guess for dzeta_dx is the dzeta_dx for a flat disk
      Vector<double> unknowns(3, 0);
      unknowns[0] = -sin(zeta) / (sqrt(x[0]*x[0]+x[1]*x[1]));
      unknowns[1] =  cos(zeta) / (sqrt(x[0]*x[0]+x[1]*x[1]));
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
      
      // polar derivatives w.r.t. (rho, t, phi)
      Vector<double> du_moffatt_drho(3,0);
      du_moffatt_drho[0] = u_polar_derivatives(0,0);   // dur_dr
      du_moffatt_drho[2] = u_polar_derivatives(1,0);   // duphi_dr

      // QUEHACERES tangent derivatives 
      
      Vector<double> du_moffatt_dphi(3,0);      
      du_moffatt_dphi[0] = u_polar_derivatives(0,1);  // dur_dphi
      du_moffatt_dphi[2] = u_polar_derivatives(1,1);  // duphi_dphi
      
      DenseMatrix<double> drho_hat_dx(3,3,0);
      DenseMatrix<double> dphi_hat_dx(3,3,0);

      // coordinates of this point in the n-s coordinate system
      double n = 0;
      double s = 0;

      // derivatives w.r.t. the global Cartesian system
      mVector ds_dx(3,0);
      mVector dn_dx(3,0);
      
      double k = Global_Parameters::n;
      double dw_dzeta = -Global_Parameters::Epsilon * k * sin(k * zeta);
	
      DenseMatrix<double> drho_vector_dx(3,3,0);

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
	  dxi_dx[1](i,j) = dtangent_dx(i,j);
	  dxi_dx[2](i,j) = dphi_hat_dx(i,j);
	}
      }
            
      // loop over the Cartesian components (row)
      for(unsigned i=0; i<3; i++)
      {
	// loop over the Cartesian components (column)
      	for(unsigned j=0; j<3; j++)
      	{
      	  du_dx(i,j) = 0;
	  
	  // loop over the Moffatt vectors
      	  for(unsigned k=0; k<3; k++)
      	  {
	    // do the dot product
      	    for(unsigned l=0; l<3; l++)
      	    {
      	      du_dx(i,j) += e_hat[i][l] * dxi_dx[k](l,j) * u_moffat[k] +
      		e_hat[i][l]*xi[k][l] * (du_moffatt_drho[k] * drho_dx[j] + du_moffatt_dphi[k]*dphi_dx[j]);
      	    }
      	  }
      	}
      }
      
      // and finally pressure, which is a scalar so no vector conversions
      u_cartesian.push_back(u_polar[2]);
    }
    else
    {
      // zero from no-slip BCs on disk
      u_cartesian[0] = 0;
      u_cartesian[1] = 0;
      u_cartesian[2] = 0;

      // infinite pressure at the edge
      u_cartesian.push_back(infinity);
    }
  }
  
  void singular_fct_and_gradient_broadside(const EdgeCoordinates& edge_coords,
					   Vector<double>& u,
					   DenseMatrix<double>& du_dx)
  {
    // parameters for broadside motion
    double A = 0;
    double B = Global_Parameters::U_far_field_broadside; // ;1;

    // forward
    singular_fct_and_gradient(edge_coords, A, B, u, du_dx);

    
  }
  
  void singular_fct_and_gradient_in_plane(const EdgeCoordinates& edge_coords,
					  Vector<double>& u,
					  DenseMatrix<double>& du_dx)
  {
    // parameters for in-plane motion
    double A = Global_Parameters::U_far_field_in_plane; // 1;
    double B = 0;

    // forward
    singular_fct_and_gradient(edge_coords, A, B, u, du_dx);
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


  double acot(const double& x)
  {
    // arccot(x) = arctan(1/x)
    return atan2(1,x);
  }

  // from c++11 complex header
  template <typename _Tp>
    std::complex<_Tp> acosh(const std::complex<_Tp>& __z)
  {
    // Kahan's formula.
    return _Tp(2.0) * std::log(std::sqrt(_Tp(0.5) * (__z + _Tp(1.0)))
  			       + std::sqrt(_Tp(0.5) * (__z - _Tp(1.0))));
  }

  double csch(const double& x)
  {
    return 1.0/sinh(x);
  }

  double Sech(const double& x)
  {
    return 1/cosh(x);
  }
  
  void gupta_solution_and_gradient(const Vector<double>& x,
				   Vector<double>& u_cartesian,
				   DenseMatrix<double>& du_dx,
				   const bool& is_lower_disk_element = false)
  {
    const double infinity = 103;

    
    // make enough space
    u_cartesian.resize(4,0);
    du_dx.resize(3,3,0);
        
    // cylindrical coordinates
    double r   = sqrt(x[0]*x[0] + x[1]*x[1]);
    double phi = atan2(x[1],x[0]);
    double z   = x[2];
    
    double Pi = MathematicalConstants::Pi;
    
    // Gupta plate velocity (set to unity, scaling will be done externally)
    double V = 1;

    // plate radius
    double a = 1;

    // mass?
    double M = 1;
    
    double p0 = 0;
    
    double tol = 1e-8;

    // if this is a point sitting on the plate but the flag has been
    // specified to say it's a lower disk element, then set the
    // z coordinate to a very small negative number to handle the
    // branch cut in the pressure above/below the disk
    if(is_lower_disk_element && abs(z) < tol)
    {
      z = -tol;
    }
    
    Vector<double> r_edge(3, 0.0);
    r_edge[0] = cos(phi);
    r_edge[1] = sin(phi);

    Vector<double> rho_vec(3, 0.0);
    for(unsigned i=0; i<3; i++)
    {
      rho_vec[i] = x[i] - r_edge[i];
    }
    
    double rho = sqrt(pow(rho_vec[0],2) + pow(rho_vec[1],2) + pow(rho_vec[2],2));
    
    // catch the case of the very edge of the disk
    tol = 1e-8;
    if(rho < tol)
    {
      u_cartesian[2] = 0.0;
      u_cartesian[3] = infinity;
      return;
    }
    
    // --------------------------------------------------------
    // solution from Al Maskari report;
    // mu and nu are oblate spheroidal coordinates
    
    std::complex<double> arg(r/a, z/a);
    
    double mu = std::real(acosh(arg));
    double nu = std::imag(acosh(arg));

    double ur = (2*V*sin(nu)*cos(nu)*pow(sinh(mu),2)) /
      (Pi*cosh(mu)*(pow(sin(nu),2) + pow(sinh(mu),2) ));

    double uz = (2*V*acot(sinh(mu)))/Pi + (2*V*pow(sin(nu),2)*sinh(mu)) /
      (Pi*(pow(sin(nu),2) + pow(sinh(mu),2) ) );

    double p = p0 + 4*M*V*sin(nu) / (Pi*a*(pow(sinh(mu),2) + pow(sin(nu),2)) );
        
    // convert to Cartesians (flow is axisymmetric so no azimuthal component)
    u_cartesian[0] = ur * cos(phi);
    u_cartesian[1] = ur * sin(phi);
    u_cartesian[2] = uz;

    u_cartesian[3] = p;
    
    // from mathematica
    
    double dmu_dz = -std::imag(std::complex<double>(1,0)/
			       (a*sqrt(std::complex<double>(-a + r, z)/a)*
				sqrt(std::complex<double>(a + r, z)/a)));
   
    double dnu_dz = std::real(std::complex<double>(1,0)/
			      (a*sqrt(std::complex<double>(-a + r, z)/a)*
			       sqrt(std::complex<double>(a + r, z)/a)));

    double dmu_dr = std::real(std::complex<double>(1,0)/
			      (a*sqrt(std::complex<double>(-a + r, z)/a)*
			       sqrt(std::complex<double>(a + r, z)/a)));
    double dnu_dr = std::imag(std::complex<double>(1,0)/
			      (a*sqrt(std::complex<double>(-a + r,z)/a)*
			       sqrt(std::complex<double>(a + r, z)/a)));
    
    // dux_dx
    du_dx(0,0) = (V*(-(r*Power(Cos(phi),2)*Power(Sech(mu),2)*
		       (Sin(4*nu)*(5*Sinh(mu) + Sinh(3*mu)) + 
			Sin(2*nu)*(-7*Sinh(3*mu) + Sinh(5*mu)))*dmu_dr)/8. + 
		     Sinh(mu)*Tanh(mu)*((-Cos(2*nu) + Cosh(2*mu))*Power(Sin(phi),2)*Sin(2*nu) + 
					2*r*Power(Cos(phi),2)*(-1 + Cos(2*nu)*Cosh(2*mu))*dnu_dr)))/
      (2.*Pi*r*Power(Power(Sin(nu),2) + Power(Sinh(mu),2),2));

    // dux_dy
    du_dx(0,1) = (V*Cos(phi)*Sin(phi)*(4*r*Cos(nu)*Sin(nu)*
				       (-Power(Sinh(mu),3) + Sech(mu)*Power(Sin(nu),2)*Tanh(mu) + 
					Sinh(mu)*(Power(Sin(nu),2) + Power(Tanh(mu),2)))*dmu_dr + 
				       2*Sinh(mu)*Tanh(mu)*(-(Sin(2*nu)*(Power(Sin(nu),2) + Power(Sinh(mu),2))) + 
							    r*(-1 + Cos(2*nu)*Cosh(2*mu))*dnu_dr)))/
      (2.*Pi*r*Power(Power(Sin(nu),2) + Power(Sinh(mu),2),2));
    
    // dux_dz
    du_dx(0,2) = (V*Cos(phi)*(2*Cos(nu)*Sin(nu)*(-Power(Sinh(mu),3) + 
						 Sech(mu)*Power(Sin(nu),2)*Tanh(mu) + 
						 Sinh(mu)*(Power(Sin(nu),2) + Power(Tanh(mu),2)))*dmu_dz + 
			      (-1 + Cos(2*nu)*Cosh(2*mu))*Sinh(mu)*Tanh(mu)*dnu_dz))/
      (Pi*Power(Power(Sin(nu),2) + Power(Sinh(mu),2),2));

    // duy_dx
    du_dx(1,0) = (V*Cos(phi)*Sin(phi)*(4*r*Cos(nu)*Sin(nu)*
				       (-Power(Sinh(mu),3) + Sech(mu)*Power(Sin(nu),2)*Tanh(mu) + 
					Sinh(mu)*(Power(Sin(nu),2) + Power(Tanh(mu),2)))*dmu_dr + 
				       2*Sinh(mu)*Tanh(mu)*(-(Sin(2*nu)*(Power(Sin(nu),2) + Power(Sinh(mu),2))) + 
							    r*(-1 + Cos(2*nu)*Cosh(2*mu))*dnu_dr)))/
      (2.*Pi*r*Power(Power(Sin(nu),2) + Power(Sinh(mu),2),2));
    
    //duy_dy
    du_dx(1,1) = (V*(-(r*Power(Sech(mu),2)*Power(Sin(phi),2)*
		       (Sin(4*nu)*(5*Sinh(mu) + Sinh(3*mu)) + 
			Sin(2*nu)*(-7*Sinh(3*mu) + Sinh(5*mu)))*dmu_dr)/8. + 
		     Sinh(mu)*Tanh(mu)*(Power(Cos(phi),2)*(-Cos(2*nu) + Cosh(2*mu))*Sin(2*nu) + 
					2*r*(-1 + Cos(2*nu)*Cosh(2*mu))*Power(Sin(phi),2)*dnu_dr)))/
      (2.*Pi*r*Power(Power(Sin(nu),2) + Power(Sinh(mu),2),2));
    
    // duy_dz
    du_dx(1,2) = (V*Sin(phi)*(2*Cos(nu)*Sin(nu)*(-Power(Sinh(mu),3) + 
						 Sech(mu)*Power(Sin(nu),2)*Tanh(mu) + 
						 Sinh(mu)*(Power(Sin(nu),2) + Power(Tanh(mu),2)))*dmu_dz + 
			      (-1 + Cos(2*nu)*Cosh(2*mu))*Sinh(mu)*Tanh(mu)*dnu_dz))/
      (Pi*Power(Power(Sin(nu),2) + Power(Sinh(mu),2),2));

    // duz_dx
    du_dx(2,0) = (2*V*Cos(phi)*Sinh(mu)*Tanh(mu)*((Power(Sin(nu),2)*(-3 + Power(Sin(nu),2)) - 
						   (1 + Power(Sin(nu),2))*Power(Sinh(mu),2))*dmu_dr + 
						  Cosh(mu)*Sin(2*nu)*Sinh(mu)*dnu_dr))/
      (Pi*Power(Power(Sin(nu),2) + Power(Sinh(mu),2),2));
    
    // duz_dy
    du_dx(2,1) = (2*V*Sin(phi)*Sinh(mu)*Tanh(mu)*((Power(Sin(nu),2)*(-3 + Power(Sin(nu),2)) - 
						   (1 + Power(Sin(nu),2))*Power(Sinh(mu),2))*dmu_dr + 
						  Cosh(mu)*Sin(2*nu)*Sinh(mu)*dnu_dr))/
      (Pi*Power(Power(Sin(nu),2) + Power(Sinh(mu),2),2));
    
    // duz_dz
    du_dx(2,2) = (-(V*Sinh(mu)*((5 + Cos(2*nu))*Power(Sin(nu),2) + 
				2*(1 + Power(Sin(nu),2))*Power(Sinh(mu),2))*Tanh(mu)*dmu_dz) + 
		  2*V*Sin(2*nu)*Power(Sinh(mu),3)*dnu_dz)/
      (Pi*Power(Power(Sin(nu),2) + Power(Sinh(mu),2),2));
  }

  void in_plane_solution_and_gradient(const Vector<double>& x,
				      Vector<double>& u,
				      DenseMatrix<double>& du_dx)
  {
    // make enough space
    u.resize(3, 0.0);

    // get the Cartesian velocity
    u = SherwoodInPlaneSolution::CartesianComponents::u(x);

    // add the pressure
    u.push_back(SherwoodInPlaneSolution::CartesianComponents::pressure(x));

    // get the gradients
    du_dx.resize(3,3, 0.0);
    du_dx = SherwoodInPlaneSolution::CartesianComponents::dudx(x);
  }
  
  void poiseuille_solution_and_gradient(const Vector<double>& x,
  					Vector<double>& u,
  					DenseMatrix<double>& du_dx)
  {
    double h = Global_Parameters::Box_half_height * 2;
    double w = Global_Parameters::Box_half_width  * 2;

    double xc = x[0] + w/2.0;
    double z  = x[2] + h/2.0;
      
    u.resize(4);

    // Poiseuille flow: ux = z(h-z)
    u[0] = z * (h - z);
    u[1] = 0;
    u[2] = 0;

    // constant pressure gradient across box    
    u[3] = 2.0 * (w - xc);
    
    du_dx.resize(3,3,0);

    // only non-zero velocity gradient is dux_dz
    du_dx(0,2) = h - 2*z; 
  }

  // sum the total contributions of all singular functions and their gradients
  void total_singular_solution_and_gradient(const Vector<double>& x,
					    Vector<double>& u,
					    DenseMatrix<double>& du_dx)
  {
    // make sure we've got enough space
    u.resize(4, 0.0);
    du_dx.resize(3,3, 0.0);

    // temporary vector and matrix to hold the individual contributions
    Vector<double> u_temp(4, 0.0);
    DenseMatrix<double> du_dx_temp(3,3, 0.0);
    
    // far-field velocities
    double U_x = Global_Parameters::U_far_field_in_plane;
    double U_z = Global_Parameters::U_far_field_broadside;
    
    if (U_x != 0)
    {
      in_plane_solution_and_gradient(x, u_temp, du_dx_temp);

      // scale and add to total (minus sign because of moving frame)
      for(unsigned i=0; i<4; i++)
	u[i] -= u_temp[i] * U_x;

      // now shift to the rest frame of the disk for consistency with Moffatt solution
      u[0] += U_x;
      
      for(unsigned i=0; i<3; i++)
      {
	for(unsigned j=0; j<3; j++)
	{
	  du_dx(i,j) -= du_dx_temp(i,j) * U_x;
	}
      }      
    }

    // broadside contribution
    if (U_z != 0)
    {
      gupta_solution_and_gradient(x, u_temp, du_dx_temp);

      // scale and add to total (minus sign because of moving frame)
      for(unsigned i=0; i<4; i++)
	u[i] -= u_temp[i] * U_z;

      // now shift to the rest frame of the disk for consistency with Moffatt solution
      u[2] += U_z;
      
      for(unsigned i=0; i<3; i++)
      {
	for(unsigned j=0; j<3; j++)
	{
	  du_dx(i,j) -= du_dx_temp(i,j) * U_z;
	}
      }      
    }
    
    if (Global_Parameters::Omega_y != 0)
    {
      // QUEHACERES write this function!
      // y_rotation_solution_and_gradient(x, u_temp, du_dx_temp);

      // QUEHACERES bring this back once the above function is written
      // // scale and add to total
      // for(unsigned i=0; i<4; i++)
      // 	u[i] += u_y_rotation[i] * Global_Parameters::Omega_y; // no length because disk radius a=1

      // for(unsigned i=0; i<3; i++)
      // {
      // 	for(unsigned j=0; j<3; j++)
      // 	{
      // 	  du_dx(i,j) += du_dx_y_rotation(i,j) * Global_Parameters::Omega_y;
      // 	}
      // }
    }
    
    if (Global_Parameters::Omega_z != 0)
    {
      // QUEHACERES write this function!
      // z_rotation_solution_and_gradient(x, u_temp, du_dx_temp);

      // QUEHACERES bring this back once the above function is written
      // // scale and add to total
      // for(unsigned i=0; i<4; i++)
      // 	u[i] += u_z_rotation[i] * Global_Parameters::Omega_z; // no length because disk radius a=1

      // for(unsigned i=0; i<3; i++)
      // {
      // 	for(unsigned j=0; j<3; j++)
      // 	{
      // 	  du_dx(i,j) += du_dx_z_rotation(i,j) * Global_Parameters::Omega_z;
      // 	}
      // }
    }
  }
}


#endif
