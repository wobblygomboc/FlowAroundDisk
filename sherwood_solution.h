#ifndef OOMPH_SHERWOOD_SOLUTION_HEADER
#define OOMPH_SHERWOOD_SOLUTION_HEADER

// get the Bessel functions from the Boost library
#include <boost/math/special_functions/bessel.hpp>

// solutions from Sherwood (2012) for a disk of radius 1 in the x-y plane,
// translating with velocity 1 in the x-direction
namespace SherwoodInPlaneSolution
{
  // function pointer for an integrand
  typedef double(*IntegrandFct) (const double& int_var, const Vector<double>& args);

  const double pi = MathematicalConstants::Pi;
  
  double sinc(const double& x)
  {
    if (x == 0)
      return 1;
    else
      return sin(x) / x;
  }
  
  // Generic integration via trapezium rule
  double trapezium_integral(const IntegrandFct& fct, const Vector<double>& args,
			    const double& lower_lim,
			    const double& upper_lim,
			    const unsigned& n) 
  { 
    // Grid spacing 
    double h = (upper_lim - lower_lim) / double(n);
  
    // Computing sum of first and last terms     
    double integral = fct(lower_lim, args) + fct(upper_lim, args); 
  
    // Adding middle terms in above formula 
    for (unsigned i = 1; i < n; i++)
    {
      double xi = lower_lim + double(i)*h;      
      integral += 2.0 * fct(xi, args); 
    }
    
    return integral * (h/2.0);
  }

  // Sherwood integrands for in-plane motion
  namespace Integrand
  {
    // Shorthand - return the Bessel function of the first kind J_n(x)
    double besselJ(const unsigned& n, const double& x)
    {
      return boost::math::cyl_bessel_j(n,x);
    }

    // x-y radial velocity component integrand
    double ur_integrand(const double& k, const Vector<double>& args)
    {
      // interpret the arguments
      double r     = args[0];
      double theta = args[1];

      // this solution is only valid for z>0 and is symmetric about z=0,
      // so if z is negative then give it a flip      
      double z = args[2] > 0 ? args[2] : -args[2];
      
      double integrand = (2*cos(theta)/(3*pi)) *
	((sinc(k) + z * sin(k)) * exp(-k*z) *(besselJ(2,k*r) -
					      besselJ(0, k*r)) +
	 4*sinc(k) * exp(-k*z) * besselJ(0, k*r) );
    
      return integrand;
    }

    // azimuthal velocity component integrand
    double uzeta_integrand(const double& k, const Vector<double>& args)
    {
      // interpret the arguments
      double r     = args[0];
      double theta = args[1];
      
      // this solution is only valid for z>0 and is symmetric about z=0,
      // so if z is negative then give it a flip      
      double z = args[2] > 0 ? args[2] : -args[2];
      
      double integrand = 2*sin(theta)/(3*pi) *
      ( (sinc(k) + z * sin(k)) * exp(-k*z) * (besselJ(2, k*r) +
					      besselJ(0, k*r)) -
	4*(sinc(k)*exp(-k*z)*besselJ(0, k*r)) );
    
      return integrand;
    }

    // vertical velocity component integrand
    double uz_integrand(const double& k, const Vector<double>& args)
    {
      // interpret the arguments
      double r     = args[0];
      double theta = args[1];

      int z_sign = args[2] / abs(args[2]);

      // this solution is only valid for z>0 and is antisymmetric about z=0,
      // so if z is negative then give it a flip (and then flip integrand)
      double z = z_sign * args[2];

      double integrand = (4*z*cos(theta) / (3.0*pi)) *
	(sin(k) * exp(-k*z) * besselJ(1, k*r));

      // solution is antisymmetric
      integrand *= z_sign;
      
      return integrand;      
    }

    // pressure integrand
    double p_integrand(const double& k, const Vector<double>& args)
    {
      // interpret the arguments
      double r     = args[0];
      double theta = args[1];

      // this solution is only valid for z>0 and is symmetric about z=0,
      // so if z is negative then give it a flip      
      double z = args[2] > 0 ? args[2] : -args[2];

      double integrand = (8*cos(theta) / (3.0*pi)) *
	(sin(k) * exp(-k*z) * besselJ(1, k*r));
    
      return integrand;
    }
  } // end of Integrands namespace

  namespace CylindricalComponents
  {
    // number of trapezium rule steps
    const unsigned nstep = 500;

    // max value of k to represent infinity
    const double k_max   = 20;
    
    // radial velocity component
    double ur(const double& r, const double& zeta, const double& z)
    {
      // package the cylindrical components into an argument vector
      Vector<double> args(3,0);
      args[0] = r;
      args[1] = zeta;
      args[2] = z;

      // do the integral
      return trapezium_integral(&Integrand::ur_integrand,
				args, 0, k_max, nstep);
    }

    // azimuthal velocity component
    double uzeta(const double& r, const double& zeta, const double& z)
    {
      // package the cylindrical components into an argument vector
      Vector<double> args(3,0);
      args[0] = r;
      args[1] = zeta;
      args[2] = z;

      // do the integral
      return trapezium_integral(&Integrand::uzeta_integrand,
				args, 0, k_max, nstep);
    }

    // vertical velocity component
    double uz(const double& r, const double& zeta, const double& z)
    {
      // package the cylindrical components into an argument vector
      Vector<double> args(3,0);
      args[0] = r;
      args[1] = zeta;
      args[2] = z;

      // do the integral
      return trapezium_integral(&Integrand::uz_integrand,
				args, 0, k_max, nstep);
    }

    double pressure(const double& r, const double& zeta, const double& z)
    {
      // package the cylindrical components into an argument vector
      Vector<double> args(3,0);
      args[0] = r;
      args[1] = zeta;
      args[2] = z;

      // do the integral
      return trapezium_integral(&Integrand::p_integrand,
				args, 0, k_max, nstep);
    }
  } // end of CylindricalComponents namespace

  namespace CartesianComponents
  {
    // Cartesian velocity as a function of the Cartesian coordinates
    Vector<double> u(const Vector<double>& x)
    {
      // get cylindrical coordinates
      double r    = sqrt(x[0]*x[0] + x[1]*x[1]);
      double zeta = atan2(x[1],x[0]);        
      double z    = x[2];
      
      // catesian output velocity
      Vector<double> u_cartesian(3, 0.0);

      // get the cylindrical components
      double ur    = CylindricalComponents::ur(r, zeta, z);
      double uzeta = CylindricalComponents::uzeta(r, zeta, z);
      double uz    = CylindricalComponents::uz(r, zeta, z);

      // convert from cylindrical components to Cartesian
      u_cartesian[0] = (ur * cos(zeta)) - (uzeta * sin(zeta));
      u_cartesian[1] = (ur * sin(zeta)) + (uzeta * cos(zeta));
      u_cartesian[2] = uz;

      return u_cartesian;
    }

    // pressure (scalar, so just get it straight from the cylindrical coordinates)
    double pressure(const Vector<double>& x)
    {
      // get cylindrical coordinates
      double r    = sqrt(x[0]*x[0] + x[1]*x[1]);
      double zeta = atan2(x[1],x[0]);        
      double z    = x[2];
      
      return CylindricalComponents::pressure(r, zeta, z);
    }

    DenseMatrix<double> dudx(const Vector<double>& x)
    {
      // increment in the coordinate for finite-diff
      const double fd_dx = 1e-6;

      // get the solution at the current point
      Vector<double> u0 = u(x);

      // coordinate vectors with increments applied to their components
      Vector<double> x_plus_dx = x;
      Vector<double> x_plus_dy = x;
      Vector<double> x_plus_dz = x;

      x_plus_dx[0] += fd_dx;
      x_plus_dy[1] += fd_dx;	    
      x_plus_dz[2] += fd_dx;
      
      // get the solution at the incremented coordinates
      Vector<double> u_plus_dx = u(x_plus_dx);
      Vector<double> u_plus_dy = u(x_plus_dy);
      Vector<double> u_plus_dz = u(x_plus_dz);

      DenseMatrix<double> du_dx(3, 3, 0.0);

      // now do the finite diff to get the velocity gradients

      // dux_dxi
      du_dx(0,0) = (u_plus_dx[0] - u0[0]) / fd_dx; // dux_dx
      du_dx(0,1) = (u_plus_dy[0] - u0[0]) / fd_dx; // dux_dy
      du_dx(0,2) = (u_plus_dz[0] - u0[0]) / fd_dx; // dux_dz

      // duy_dxi
      du_dx(1,0) = (u_plus_dx[1] - u0[1]) / fd_dx; // duy_dx
      du_dx(1,1) = (u_plus_dy[1] - u0[1]) / fd_dx; // duy_dy
      du_dx(1,2) = (u_plus_dz[1] - u0[1]) / fd_dx; // duy_dz

      // duz_dxi
      du_dx(2,0) = (u_plus_dx[2] - u0[2]) / fd_dx; // duz_dx
      du_dx(2,1) = (u_plus_dy[2] - u0[2]) / fd_dx; // duz_dy
      du_dx(2,2) = (u_plus_dz[2] - u0[2]) / fd_dx; // duz_dz

      return du_dx;
    }
      
  } // end of CartesianComponents namespace
} // end of SherwoodInPlaneSolution
#endif
