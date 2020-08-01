#ifndef OOMPH_EXACT_SOLUTIONS_FINITE_DISK_HEADER
#define OOMPH_EXACT_SOLUTIONS_FINITE_DISK_HEADER

namespace FlatDiskExactSolutions
{  
  // shorthand
  const double pi = MathematicalConstants::Pi;

  const double infinity = 1003;

  // how far to shift the radius to compute "infinity"
  double dr = -1e-5;
    
  // inverse cotangent function
  double acot(const double& x)
  {
    // catch the singular case 
    if(x == 0)
      return pi/2.0;
    else
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

  // hyperbolic cosecant function
  double csch(const double& x)
  {
    return 1.0/sinh(x);
  }

  // hyperbolic secant function
  double Sech(const double& x)
  {
    return 1/cosh(x);
  }

  // function to compute the (lambda,zeta) oblate spheroidal coordinates
  void oblate_spheroidal_coordinates(const double& r, const double& z,
				     double& lambda, double& zeta)
  {
    std::complex<double> arg(r, z);
    
    double mu = std::real(acosh(arg));
    double nu = std::imag(acosh(arg));

    lambda = sinh(mu);
    zeta   = sin(nu);
  }

  Vector<double> velocity_cylindrical_to_cartesian(const double& u_r,
						   const double& u_theta,
						   const double& u_z,
						   const double& theta)
  {
    Vector<double> u(3,0);

    // convert to Cartesians
    u[0] = u_r*cos(theta) - u_theta*sin(theta);
    u[1] = u_r*sin(theta) + u_theta*cos(theta);
    u[2] = u_z;

    return u;
  }

  Vector<double> shift_x_away_from_disk_edge(const Vector<double>& x)
  {
    // tolerance for being on the edge of the disk
    const double tol = 1e-8;

    // temporary vector to hold the adjusted coordinates
    Vector<double> x_shifted = x;
    
    // azimuthal angle
    double theta = atan2(x[1],x[0]);
    
    // Cartesian coordinates of edge of disk at this zeta
    Vector<double> r_edge(3, 0.0);
    r_edge[0] = cos(theta);
    r_edge[1] = sin(theta);

    // vector from edge of disk to this point
    Vector<double> rho_vec(3, 0.0);
    for(unsigned i=0; i<3; i++)
    {
      rho_vec[i] = x[i] - r_edge[i];
    }
    
    // distance from edge of disk to this point
    double rho = sqrt(pow(rho_vec[0],2) + pow(rho_vec[1],2) + pow(rho_vec[2],2));
    
    // is it on the disk edge?
    if(rho < tol)
    {
      // increment the cylindrical (x-y) radius by dr
      x_shifted[0] += dr * cos(theta);
      x_shifted[1] += dr * sin(theta);
    }
    
    return x_shifted; 
  }

  /* struct RigidBodyMotion */
  /* { */
  /* RigidBodyMotion() : U_x(0), U_z(0), Omega_z(0), Omega_minus_y(0) {} */
    
  /*   double U_x; */
  /*   double U_z; */
  /*   double Omega_z; */
  /*   double Omega_minus_y; */
  /* }; */
  
  // //////////////////////////////////////////////////////////////////////////
  // //////////////////////////////////////////////////////////////////////////
  // //////////////////////////////////////////////////////////////////////////
  // //////////////////////////////////////////////////////////////////////////

  
  // ==========================================================================
  /// \short Exact solutions given in terms of oblate spheroidal coordinates
  /// (Tanzosh & Stone (1996))
  // ==========================================================================

  /// solution for broadside motion, i.e. translation along the axis of
  /// rotational symmetry (z-axis)
  void broadside_translation_solution(const Vector<double>& x_const,
				      Vector<double>& u_cartesian)
  {
    // check if this point is exactly on the edge of the disk,
    // and move it slightly away if so
    Vector<double> x = shift_x_away_from_disk_edge(x_const);
    
    // cylindrical coordinates
    double r   = sqrt(x[0]*x[0] + x[1]*x[1]);
    double theta = atan2(x[1],x[0]);

    // the solution is only valid for the half-space z>=0,
    // so flip the z coordinate if we're in the lower half-space
    // and apply the appropriate symmetries after
    bool in_lower_half_space = (x[2] < 0);
    const double z = in_lower_half_space ? -x[2] : x[2];
    
    double Pi = MathematicalConstants::Pi;
    
    // Gupta plate velocity (set to unity, scaling will be done externally)
    double V = 1;

    // plate radius
    double a = 1;

    // mass?
    double M = 1;
    
    double p0 = 0;
        
    // --------------------------------------------------------
    // solution from Tanzosh & Stone (1996);
    // mu and nu are oblate spheroidal coordinates
    double lambda = 0;
    double zeta   = 0;

    // get 'em
    oblate_spheroidal_coordinates(r, z, lambda, zeta);

    // QUEHACERES delete?
    std::complex<double> arg(r, z);
    double mu = std::real(acosh(arg));
    double nu = std::imag(acosh(arg));

    double u_r = (2.0/Pi)*sqrt((1 - zeta*zeta)/(1+lambda*lambda))*
      (lambda*lambda*zeta) / (lambda*lambda + zeta*zeta);

    // flow is axisymmetric so no azimuthal component
    double u_theta = 0;
	
    double u_z = (2/Pi)*(acot(lambda) + zeta*zeta*lambda / (zeta*zeta + lambda*lambda) );

    double p = (4.0/Pi)*zeta / (lambda*lambda + zeta*zeta);

    // apply the appropriate symmetries for points in the lower half-space
    if(in_lower_half_space)
    {
      // u_z is symmetric about z=0, u_r, and p are anti-symmetric
      u_r = -u_r;
      p   = -p;
    }
    
    // convert to Cartesians
    u_cartesian = velocity_cylindrical_to_cartesian(u_r, u_theta, u_z, theta);

    // and add pressure to solution vector
    u_cartesian.push_back(p);

   
  }

  // //////////////////////////////////////////////////////////////////////////
  // //////////////////////////////////////////////////////////////////////////
  
  // velocity gradients for broadside translation along z-axis
  void broadside_translation_velocity_gradient(const Vector<double>& x_const,
					       DenseMatrix<double>& du_dx,
					       const bool& is_lower_disk_element = false)
  {
    // check if this point is exactly on the edge of the disk,
    // and move it slightly away if so
    Vector<double> x = shift_x_away_from_disk_edge(x_const);
    
    // ========================================================================
    // Derivatives
    // ========================================================================

    // cylindrical coords
    double r     = sqrt(x[0]*x[0] + x[1]*x[1]);
    double theta = atan2(x[1], x[0]);
    double z     = x[2];

    // disk radius
    double a = 1.0;

    // disk velocity (will be scaled later)
    double V = 1.0;
    
    // mu and nu are oblate spheroidal coordinates
    double lambda = 0;
    double zeta   = 0;

    // get 'em
    oblate_spheroidal_coordinates(r, z, lambda, zeta);

    // QUEHACERES delete?
    std::complex<double> arg(r, z);
    double mu = std::real(acosh(arg));
    double nu = std::imag(acosh(arg));
    
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

    // make enough space
    du_dx.resize(3,3,0);
	
    // dux_dx
    du_dx(0,0) = (V*(-(r*Power(Cos(theta),2)*Power(Sech(mu),2)*
		       (Sin(4*nu)*(5*lambda + Sinh(3*mu)) + 
			Sin(2*nu)*(-7*Sinh(3*mu) + Sinh(5*mu)))*dmu_dr)/8. + 
		     lambda*Tanh(mu)*((-Cos(2*nu) + Cosh(2*mu))*Power(Sin(theta),2)*Sin(2*nu) + 
					2*r*Power(Cos(theta),2)*(-1 + Cos(2*nu)*Cosh(2*mu))*dnu_dr)))/
      (2.*pi*r*Power(Power(zeta,2) + Power(lambda,2),2));

    // dux_dy
    du_dx(0,1) = (V*Cos(theta)*Sin(theta)*(4*r*Cos(nu)*zeta*
				       (-Power(lambda,3) + Sech(mu)*Power(zeta,2)*Tanh(mu) + 
					lambda*(Power(zeta,2) + Power(Tanh(mu),2)))*dmu_dr + 
				       2*lambda*Tanh(mu)*(-(Sin(2*nu)*(Power(zeta,2) + Power(lambda,2))) + 
							    r*(-1 + Cos(2*nu)*Cosh(2*mu))*dnu_dr)))/
      (2.*pi*r*Power(Power(zeta,2) + Power(lambda,2),2));
    
    // dux_dz
    du_dx(0,2) = (V*Cos(theta)*(2*Cos(nu)*zeta*(-Power(lambda,3) + 
						 Sech(mu)*Power(zeta,2)*Tanh(mu) + 
						 lambda*(Power(zeta,2) + Power(Tanh(mu),2)))*dmu_dz + 
			      (-1 + Cos(2*nu)*Cosh(2*mu))*lambda*Tanh(mu)*dnu_dz))/
      (pi*Power(Power(zeta,2) + Power(lambda,2),2));

    // duy_dx
    du_dx(1,0) = (V*Cos(theta)*Sin(theta)*(4*r*Cos(nu)*zeta*
				       (-Power(lambda,3) + Sech(mu)*Power(zeta,2)*Tanh(mu) + 
					lambda*(Power(zeta,2) + Power(Tanh(mu),2)))*dmu_dr + 
				       2*lambda*Tanh(mu)*(-(Sin(2*nu)*(Power(zeta,2) + Power(lambda,2))) + 
							    r*(-1 + Cos(2*nu)*Cosh(2*mu))*dnu_dr)))/
      (2.*pi*r*Power(Power(zeta,2) + Power(lambda,2),2));
    
    //duy_dy
    du_dx(1,1) = (V*(-(r*Power(Sech(mu),2)*Power(Sin(theta),2)*
		       (Sin(4*nu)*(5*lambda + Sinh(3*mu)) + 
			Sin(2*nu)*(-7*Sinh(3*mu) + Sinh(5*mu)))*dmu_dr)/8. + 
		     lambda*Tanh(mu)*(Power(Cos(theta),2)*(-Cos(2*nu) + Cosh(2*mu))*Sin(2*nu) + 
					2*r*(-1 + Cos(2*nu)*Cosh(2*mu))*Power(Sin(theta),2)*dnu_dr)))/
      (2.*pi*r*Power(Power(zeta,2) + Power(lambda,2),2));
    
    // duy_dz
    du_dx(1,2) = (V*Sin(theta)*(2*Cos(nu)*zeta*(-Power(lambda,3) + 
						 Sech(mu)*Power(zeta,2)*Tanh(mu) + 
						 lambda*(Power(zeta,2) + Power(Tanh(mu),2)))*dmu_dz + 
			      (-1 + Cos(2*nu)*Cosh(2*mu))*lambda*Tanh(mu)*dnu_dz))/
      (pi*Power(Power(zeta,2) + Power(lambda,2),2));

    // duz_dx
    du_dx(2,0) = (2*V*Cos(theta)*lambda*Tanh(mu)*((Power(zeta,2)*(-3 + Power(zeta,2)) - 
						   (1 + Power(zeta,2))*Power(lambda,2))*dmu_dr + 
						  Cosh(mu)*Sin(2*nu)*lambda*dnu_dr))/
      (pi*Power(Power(zeta,2) + Power(lambda,2),2));
    
    // duz_dy
    du_dx(2,1) = (2*V*Sin(theta)*lambda*Tanh(mu)*((Power(zeta,2)*(-3 + Power(zeta,2)) - 
						   (1 + Power(zeta,2))*Power(lambda,2))*dmu_dr + 
						  Cosh(mu)*Sin(2*nu)*lambda*dnu_dr))/
      (pi*Power(Power(zeta,2) + Power(lambda,2),2));
    
    // duz_dz
    du_dx(2,2) = (-(V*lambda*((5 + Cos(2*nu))*Power(zeta,2) + 
				2*(1 + Power(zeta,2))*Power(lambda,2))*Tanh(mu)*dmu_dz) + 
		  2*V*Sin(2*nu)*Power(lambda,3)*dnu_dz)/
      (pi*Power(Power(zeta,2) + Power(lambda,2),2));
  }

  // //////////////////////////////////////////////////////////////////////////
  // //////////////////////////////////////////////////////////////////////////
  
  /// Exact solution for in-plane motion along x-axis.
  void in_plane_translation_solution(const Vector<double>& x_const,
				     Vector<double>& u)
  {
    // check if this point is exactly on the edge of the disk,
    // and move it slightly away if so
    Vector<double> x = shift_x_away_from_disk_edge(x_const);
    
    // x-y radius
    const double r = sqrt(x[0]*x[0] + x[1]*x[1]);

    // azimuthal angle
    const double theta = atan2(x[1], x[0]);

    // the solution is only valid for the half-space z>=0,
    // so flip the z coordinate if we're in the lower half-space
    // and apply the appropriate symmetries after
    bool in_lower_half_space = (x[2] < 0);
    const double z = in_lower_half_space ? -x[2] : x[2];
      
    // oblate spheroidal coordinates
    double lambda = 0;
    double zeta   = 0;
    
    // get 'em
    oblate_spheroidal_coordinates(r, z, lambda, zeta);
    
    // using notation of Tanzosh & Stone (1996), theta is the azimuthal angle
    double u_r = (2.0 * cos(theta) / (3.0*pi)) * (
      3*acot(lambda) - (lambda*zeta*zeta)/(lambda*lambda + zeta*zeta) +
      (pow(lambda,3) * (1-zeta*zeta))/((1+lambda*lambda)*(lambda*lambda + zeta*zeta)) );
    
    double u_theta = (2.0 * sin(theta) / (3.0*pi)) * (
      -3*acot(lambda) + (lambda*zeta*zeta)/(lambda*lambda + zeta*zeta) +
      (pow(lambda,3) * (1-zeta*zeta))/((1+lambda*lambda)*(lambda*lambda + zeta*zeta)) );;
    
    double u_z = cos(theta) * (4.0/(3.0*pi)) *
      sqrt((1 - zeta*zeta)/(1+lambda*lambda)) * (lambda*lambda*zeta) /
      (lambda*lambda + zeta*zeta);
    
    double p = (8.0*cos(theta)/(3.0*pi)) *
      sqrt((1 - zeta*zeta)/(1+lambda*lambda)) * lambda / (lambda*lambda + zeta*zeta);

    // apply the appropriate symmetries for points in the lower half-space
    if(in_lower_half_space)
    {
      // u_r, u_theta and p are symmetric about z=0, u_z is anti-symmetric
      u_z = -u_z;
    }
    
    // convert
    u = velocity_cylindrical_to_cartesian(u_r, u_theta, u_z, theta);

    // and add pressure to solution vector
    u.push_back(p);
  }

  // //////////////////////////////////////////////////////////////////////////
  // //////////////////////////////////////////////////////////////////////////
  
  /// \short Exact solution for in-plane rotation (about z-axis).
  /// No need to check for the point being exactly on the edge, since the
  /// solution is regular for this mode. 
  void in_plane_rotation_solution(const Vector<double>& x,
				  Vector<double>& u)
  { 
    // x-y radius
    const double r = sqrt(x[0]*x[0] + x[1]*x[1]);

    // azimuthal angle
    const double theta = atan2(x[1], x[0]);

    // the solution is only valid for the half-space z>=0,
    // so flip the z coordinate if we're in the lower half-space
    // and apply the appropriate symmetries after
    bool in_lower_half_space = (x[2] < 0);
    const double z = in_lower_half_space ? -x[2] : x[2];
    
    // oblate spheroidal coordinates
    double lambda = 0;
    double zeta   = 0;
    
    // get 'em
    oblate_spheroidal_coordinates(r, z, lambda, zeta);

    // by definition, for a pure in-plane rotation we have
    // only azimuthal motion
    double u_r = 0;    
    double u_z = 0;
    double p   = 0;

    double u_theta = (2.0/pi) * sqrt((1 - zeta*zeta)/(1+lambda*lambda)) *
      ((1+lambda*lambda)*acot(lambda) - lambda);

    // solution is symmetric, no need to flip any components
    
    // convert
    u = velocity_cylindrical_to_cartesian(u_r, u_theta, u_z, theta);
    
    // and add pressure to solution vector
    u.push_back(p);
  }

  // //////////////////////////////////////////////////////////////////////////
  // //////////////////////////////////////////////////////////////////////////

  /// \short Exact solution for out-of-plane rotation (about -y axis).
  void out_of_plane_rotation_solution(const Vector<double>& x_const,
				      Vector<double>& u)
  {
    // check if this point is exactly on the edge of the disk,
    // and move it slightly away if so
    Vector<double> x = shift_x_away_from_disk_edge(x_const);
    
    // x-y radius
    const double r = sqrt(x[0]*x[0] + x[1]*x[1]);

    // azimuthal angle
    const double theta = atan2(x[1], x[0]);

    // the solution is only valid for the half-space z>=0,
    // so flip the z coordinate if we're in the lower half-space
    // and apply the appropriate symmetries after
    bool in_lower_half_space = (x[2] < 0);
    const double z = in_lower_half_space ? -x[2] : x[2];
      
    // oblate spheroidal coordinates
    double lambda = 0;
    double zeta   = 0;
    
    // get 'em
    oblate_spheroidal_coordinates(r, z, lambda, zeta);

    double u_r = (2*cos(theta)/pi) * lambda*zeta *
      ( -acot(lambda) + (lambda / (lambda*lambda + zeta*zeta)) *
	(1 + (1 - zeta*zeta) / (1 + lambda*lambda) ) );
    
    double u_theta = (2*sin(theta)/pi) * lambda*zeta *
      ( acot(lambda) - lambda / (lambda*lambda + zeta*zeta) );
    
    double u_z = (2*cos(theta)/pi) * sqrt((1 - zeta*zeta)/(1+lambda*lambda)) *
      ( (1+lambda*lambda)*acot(lambda) -
	lambda*(lambda*lambda - zeta*zeta) / (lambda*lambda + zeta*zeta) );
    
    double p = (8*cos(theta)/pi) * sqrt((1 - zeta*zeta)/(1+lambda*lambda)) *
      lambda / (lambda*lambda + zeta*zeta);

    if(in_lower_half_space)
    {
      // u_z is symmetric about z=0, u_r, and p are anti-symmetric
      u_r = -u_r;
      p   = -p;
    }
    
    // convert
    u = velocity_cylindrical_to_cartesian(u_r, u_theta, u_z, theta);
    
    // and add pressure to solution vector
    u.push_back(p);
  }

  // //////////////////////////////////////////////////////////////////////////
  // //////////////////////////////////////////////////////////////////////////

  /* void poiseuille_solution_and_gradient(const Vector<double>& x, */
  /* 					Vector<double>& u, */
  /* 					DenseMatrix<double>& du_dx) */
  /* { */
  /*   double h = Global_Parameters::Box_half_height * 2; */
  /*   double w = Global_Parameters::Box_half_width  * 2; */

  /*   double xc = x[0] + w/2.0; */
  /*   double z  = x[2] + h/2.0; */
      
  /*   u.resize(4); */

  /*   // Poiseuille flow: ux = z(h-z) */
  /*   u[0] = z * (h - z); */
  /*   u[1] = 0; */
  /*   u[2] = 0; */

  /*   // constant pressure gradient across box     */
  /*   u[3] = 2.0 * (w - xc); */
    
  /*   du_dx.resize(3,3,0); */

  /*   // only non-zero velocity gradient is dux_dz */
  /*   du_dx(0,2) = h - 2*z;  */
  /* } */

  // //////////////////////////////////////////////////////////////////////////
  // //////////////////////////////////////////////////////////////////////////
  
  // sum the total contributions of the linear combination of exact solutions
  void total_exact_solution(const Vector<double>& x,
			    // ### QUEHACERES delete   /* const RigidBodyMotion& rigid_body_motion, */
			    const Vector<double> u_disk,
			    const Vector<double> omega_disk,
			    Vector<double>& u)
  {
    // make sure we've got enough space
    u.resize(4, 0.0);

    // zero out the solution, since we're going to be adding to it below
    std::fill(u.begin(), u.end(), 0.0);

    // interpret the linear and angular velocities
    const double U_in_plane    =  u_disk[0];
    const double U_broadside   =  u_disk[2];
    const double Omega_minus_y = -omega_disk[1];
    const double Omega_z       =  omega_disk[2];
      
    // temporary vector and matrix to hold the individual contributions
    Vector<double> u_temp(4, 0.0);
  
    if (U_in_plane != 0)
    {
      // get the exact in plane velocity and pressure
      in_plane_translation_solution(x, u_temp);      
      
      // scale and add to total // ### QUEHACERES delete: (minus sign because of moving frame)
      for(unsigned i=0; i<4; i++)
	u[i] += u_temp[i] * U_in_plane;

      /* // now shift to the rest frame of the disk for consistency with Moffatt solution */
      /* u[0] += U_x; */

      // QUEHACERES not doing gradients for now
      /* for(unsigned i=0; i<3; i++) */
      /* { */
      /* 	for(unsigned j=0; j<3; j++) */
      /* 	{ */
      /* 	  du_dx(i,j) -= du_dx_temp(i,j) * U_x; */
      /* 	} */
      /* }       */
    }

    // broadside contribution
    if (U_broadside != 0)
    {
      broadside_translation_solution(x, u_temp);

      // scale and add to total  // ### QUEHACERES delete (minus sign because of moving frame)
      for(unsigned i=0; i<4; i++)
	u[i] += u_temp[i] * U_broadside;

      /* // now shift to the rest frame of the disk for consistency with Moffatt solution */
      /* u[2] += U_broadside; */
    }
    
    if (Omega_minus_y != 0)
    {
      out_of_plane_rotation_solution(x, u_temp);

      // radius and azimuth
      const double r     = sqrt(x[0]*x[0] + x[1]*x[1]);
      const double theta = atan2(x[1],x[0]);
      
      // scale and add to total
      for(unsigned i=0; i<4; i++)
      	u[i] += u_temp[i] * Omega_minus_y; // no length because disk radius a=1

      /* // now shift to the frame of the disk for consistency with Moffatt solution */
      /* // \bm u = \bm \Omega \ctimes \bm r */
      /* u[0] += -x[2] * Omega_minus_y; // u_x = -r_z * Omega_minus_y */
      /* u[2] += -x[0] * Omega_minus_y; // u_z = -r_x * Omega_minus_y */
    }
    
    if (Omega_z != 0)
    {
      in_plane_rotation_solution(x, u_temp);

      // radius and azimuth
      const double r     = sqrt(x[0]*x[0] + x[1]*x[1]);
      const double theta = atan2(x[1],x[0]);
      
      // scale and add to total
      for(unsigned i=0; i<4; i++)
      	u[i] += u_temp[i] * Omega_z; // no length because disk radius a=1

      // ### QUEHACERES delete
      /* // now shift to the frame of the disk for consistency with Moffatt solution */
      /* // \bm u = \bm \Omega \ctimes \bm r */
      /* u[0] += -x[1] * Omega_z; // u_x = -r_y * Omega_z */
      /* u[1] +=  x[0] * Omega_z; // u_y = +r_x * Omega_z */
    }
  }

  // @@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
  // @@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
  // @@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
  // QUEHACERES adding back in for debug of finite diff dudx
  // @@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
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
  // @@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
  // @@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
  
  
  // sum the total contributions of all singular functions
  void total_exact_velocity_gradient(const Vector<double>& x,				     
				     const Vector<double> u_disk,
				     const Vector<double> omega_disk,
				     DenseMatrix<double>& du_dx)
  {
    // coordinate increment for finite-difference gradient
    const double fd_dx = 1e-6;

    // make enough space
    du_dx.resize(3,3, 0.0);
      
    // get the solution at the current point
    Vector<double> u0(4,0);
    total_exact_solution(x, u_disk, omega_disk, u0);

    for(unsigned j=0; j<3; j++)
    {
      // coordinate vectors with increments applied to their components
      Vector<double> x_plus_dx = x;

      // add the increments in each direction
      x_plus_dx[j] += fd_dx;
      
      // get the solution at the incremented coordinates
      Vector<double> u_plus_dx(4, 0.0);
      total_exact_solution(x_plus_dx, u_disk, omega_disk, u_plus_dx);

      // now do the finite diff to get the velocity gradients
      for(unsigned i=0; i<3; i++)
	du_dx(i,j) = (u_plus_dx[i] - u0[i]) / fd_dx;
    }
  }
  
} // end of FlatDiskExactSolutions namespace


#endif
