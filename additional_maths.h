#ifndef OOMPH_ADDITIONAL_MATHS_HEADER
#define OOMPH_ADDITIONAL_MATHS_HEADER

//Generic routines
#include "generic.h"

#include "mathematica_definitions.h"


inline double atan2pi(const double& y, const double& x)
{
  // Polar angle
  double theta = std::atan2(y,x);

  // prevent atan2 negative angle fuckery that causes a discontinuity at theta=pi
  if (y < 0.0)
  {
    theta += 2.0 * MathematicalConstants::Pi;
  }

  return theta;
}

template <typename T>
inline int sgn(T val) {
  return (T(0) < val) - (val < T(0));
}


inline double delta(const unsigned& i, const unsigned& j)
{
  return static_cast<double>(i == j);
}

// inverse cotangent function
inline double acot(const double& x)
{
  // catch the limiting case 
  if(x == 0)
    return MathematicalConstants::Pi / 2.0;
  else
    // arccot(x) = arctan(1/x)
    return atan2(1, x);
}

// secant, for mathematica output
inline double Sec(const double& theta)
{
  return 1.0/cos(theta);
}

// hyperbolic cosecant function
inline double csch(const double& x)
{
  return 1.0/sinh(x);
}

// hyperbolic secant function
inline double Sech(const double& x)
{
  return 1/cosh(x);
}

// function to map any angle (+ve or -ve) to [0,2pi)
inline double map_angle_to_range_0_to_2pi(const double& signed_angle)
{
  // shorthand
  const double two_pi = 2.0 * MathematicalConstants::Pi;

  // effectively implements Python modulo division
  // (C++ would return a negative angle if we did -angle % 2pi)    
  return signed_angle - two_pi * floor(signed_angle/two_pi);    
}

inline double map_angle_to_range_plus_minus_pi(const double& angle)
{
  // shorthand
  const double pi = MathematicalConstants::Pi;
  
  const double angle_zero_to_2pi = map_angle_to_range_0_to_2pi(angle);

  return angle_zero_to_2pi > pi ? angle_zero_to_2pi - 2*pi : angle_zero_to_2pi;
}

#if __cplusplus == 199711L
// from c++11 complex header
template <typename _Tp>
std::complex<_Tp> acosh(const std::complex<_Tp>& __z)
{
  // Kahan's formula.
  return _Tp(2.0) * std::log(std::sqrt(_Tp(0.5) * (__z + _Tp(1.0)))
			     + std::sqrt(_Tp(0.5) * (__z - _Tp(1.0))));
}

#endif

/// Rotation matrix for right-handed rotation of theta about e_x
DenseMatrix<double> R_x(const double& theta)
{  
  DenseMatrix<double> R(3,3,0.0);
  R(0,0) = 1;
  R(1,1) = cos(theta);
  R(1,2) = -sin(theta);
  R(2,1) = sin(theta);
  R(2,2) = cos(theta);

  return R;
}

/// Rotation matrix for right-handed rotation of theta about e_y
DenseMatrix<double> R_y(const double& theta)
{  
  DenseMatrix<double> R(3, 3, 0.0);
  R(0,0) = cos(theta);
  R(0,2) = sin(theta);
  R(1,1) = 1;
  R(2,0) = -sin(theta);
  R(2,2) = cos(theta);

  return R;
}

/// Rotation matrix for right-handed rotation of theta about e_z
DenseMatrix<double> R_z(const double& theta)
{  
  DenseMatrix<double> R(3, 3, 0.0);
  R(0,0) = cos(theta);
  R(0,1) = -sin(theta);
  R(1,0) = sin(theta);
  R(1,1) = cos(theta);
  R(2,2) = 1;

  return R;
}

DenseMatrix<double> identity_matrix()
{
  DenseMatrix<double> M(3, 3, 0.0);

  for(unsigned i=0; i<3; i++)
    M(i,i) = 1.0;

  return M;
}

struct OblateSpheroidalCoordinates
{
  // default constructor that initialises everything to zero
OblateSpheroidalCoordinates() : lambda(0.0), zeta(0.0), theta(0.0), r(0.0), z(0.0) {}
  
  // optional constructor with doubles
OblateSpheroidalCoordinates(const double& r_, const double& theta_, const double& z_) :
  lambda(0.0), zeta(0.0), theta(theta_), r(r_), z(z_)
    {
      std::complex<double> arg(r_, z_);
      
      // compute the intermediate / alternative oblate spheroidal coords
      double mu = std::real(acosh(arg));
      double nu = std::imag(acosh(arg));
      
      // and convert
      lambda = sinh(mu);
      zeta   = sin(nu);
    }
  
  // optional constructor for 2D coordinates from cylindricals
OblateSpheroidalCoordinates(const double& r_, const double& z_) :
  OblateSpheroidalCoordinates(r_, 0.0, z_) { }
			    
  // copy constructor
OblateSpheroidalCoordinates(const OblateSpheroidalCoordinates& oblate_sph_coords) :
  lambda(oblate_sph_coords.lambda), zeta(oblate_sph_coords.zeta),
    theta(oblate_sph_coords.theta), r(oblate_sph_coords.r), z(oblate_sph_coords.z) {}
  
  // assignment operator
  void operator=(const OblateSpheroidalCoordinates& oblate_sph_coords)
    {
      // assign the oblate spheroidal coordinates
      lambda = oblate_sph_coords.lambda;
      zeta   = oblate_sph_coords.zeta;
      theta  = oblate_sph_coords.theta;

      // and compute the cylindrical components
      r = (lambda*lambda + 1) * (1 - zeta*zeta);
      z = lambda * zeta;
    }

    // allows for effective signed zero, as a cleaner way of keeping track
    // of whether we're sat on the upper or lower surface of the disk
    int sign_of_z() const
    {
      return std::signbit(z) ? -1 : 1;
    }

  // returns matrix of the derivatives of the oblate spheroidal coordinates
  // (lambda, zeta, theta) w.r.t. the cylindrical coordinates (r,theta,z)
  // Matrix is: d(lambda,zeta,theta)/d(r,theta,z) (note order of theta)
  // Computed by finite difference
  void doblate_spheroidal_dcylindrical_by_fd(DenseMatrix<double>& dobl_sph_dcyl)
    {
      const int sign = sign_of_z();
      
      // finite-diff step in the cylindrical coordinate
      // (in the appropriate direction to avoid crossing the x-y plane)
      const double fd_step = 1e-8 * (double)(sign);

      // make space
      dobl_sph_dcyl.resize(3, 3, 0.0);

      // take the finite diff steps in the cyclindrical coordinates
      OblateSpheroidalCoordinates obl_sph_plus_dr(r+fd_step, z);
      OblateSpheroidalCoordinates obl_sph_plus_dz(r, z+fd_step);

      // N.B. lambda is non-negative, so no direction correction
      
      // dlambda/dr
      dobl_sph_dcyl(0,0) = /* (double)(sign) * */ (obl_sph_plus_dr.lambda - lambda)/fd_step;

      // dlambda/dz
      dobl_sph_dcyl(0,2) = /* (double)(sign) *  */(obl_sph_plus_dz.lambda - lambda)/fd_step;
      
      // dzeta/dr
      dobl_sph_dcyl(1,0) = (double)(sign) * (obl_sph_plus_dr.zeta - zeta)/fd_step;

      // dzeta/dz
      dobl_sph_dcyl(1,2) = (double)(sign) * (obl_sph_plus_dz.zeta - zeta)/fd_step;

      // oblate spheroidal coords are orthogonal, so only nonzero theta
      // derivative is theta
      
      //dtheta/dtheta
      dobl_sph_dcyl(2,1) = 1.0;
    }
  
  // returns matrix of the derivatives of the oblate spheroidal coordinates
  // (lambda, zeta, theta) w.r.t. the cylindrical coordinates (r,theta,z)
  void doblate_spheroidal_dcylindrical(DenseMatrix<double>& dobl_sph_dcyl)
    {
      //* / QUEHACERES redirect for now */
      /* doblate_spheroidal_dcylindrical_by_fd(dobl_sph_dcyl); */
      /* return; */
      
      dobl_sph_dcyl.resize(3, 3, 0.0);

      // --------------------
      // dlambda/dr
      dobl_sph_dcyl(0,0) =
	(Cosh(Log(Power(r,2) + Power(z,2) + 
		  Sqrt(Power(-1 + r,2) + Power(z,2))*
		  Sqrt(Power(1 + r,2) + Power(z,2)) + 
		  2*Power(Power(-1 + r,2) + Power(z,2),0.25)*
		  Power(Power(1 + r,2) + Power(z,2),0.25)*
		  (r*Cos((ArcTan(-1 + r,z) + ArcTan(1 + r,z))/2.) + 
		   z*Sin((ArcTan(-1 + r,z) + ArcTan(1 + r,z))/2.)))/2.)*
	 (r*Power(Power(-1 + r,2) + Power(z,2),0.25)*
	  Power(Power(1 + r,2) + Power(z,2),0.25)*
	  (-1 + Power(r,2) + Power(z,2) + 
	   Sqrt(Power(-1 + r,2) + Power(z,2))*Sqrt(Power(1 + r,2) + Power(z,2)))
	  + (1 + 2*Power(r,4) + Power(z,2) + Power(r,2)*(-3 + 2*Power(z,2)))*
	  Cos((ArcTan(-1 + r,z) + ArcTan(1 + r,z))/2.) + 
	  2*r*z*(Power(r,2) + Power(z,2))*
	  Sin((ArcTan(-1 + r,z) + ArcTan(1 + r,z))/2.)))/
	(Power(Power(-1 + r,2) + Power(z,2),0.75)*
	 Power(Power(1 + r,2) + Power(z,2),0.75)*
	 (Power(r,2) + Power(z,2) + 
	  Sqrt(Power(-1 + r,2) + Power(z,2))*Sqrt(Power(1 + r,2) + Power(z,2)) + 
	  2*Power(Power(-1 + r,2) + Power(z,2),0.25)*
	  Power(Power(1 + r,2) + Power(z,2),0.25)*
	  (r*Cos((ArcTan(-1 + r,z) + ArcTan(1 + r,z))/2.) + 
	   z*Sin((ArcTan(-1 + r,z) + ArcTan(1 + r,z))/2.))));
  
      // --------------------
      // dlambda/dz
      dobl_sph_dcyl(0,2) = (Cosh(Log(Power(r,2) + Power(z,2) + 
				     Sqrt(Power(-1 + r,2) + Power(z,2))*
				     Sqrt(Power(1 + r,2) + Power(z,2)) + 
				     2*Power(Power(-1 + r,2) + Power(z,2),0.25)*
				     Power(Power(1 + r,2) + Power(z,2),0.25)*
				     (r*Cos((ArcTan(-1 + r,z) + ArcTan(1 + r,z))/2.) + 
				      z*Sin((ArcTan(-1 + r,z) + ArcTan(1 + r,z))/2.)))/2.)*
			    (z*Power(Power(-1 + r,2) + Power(z,2),0.25)*
			     Power(Power(1 + r,2) + Power(z,2),0.25)*
			     (1 + Power(r,2) + Power(z,2) + 
			      Sqrt(Power(-1 + r,2) + Power(z,2))*Sqrt(Power(1 + r,2) + Power(z,2)))
			     + 2*r*z*(Power(r,2) + Power(z,2))*
			     Cos((ArcTan(-1 + r,z) + ArcTan(1 + r,z))/2.) + 
			     (1 - Power(r,2) + (3 + 2*Power(r,2))*Power(z,2) + 2*Power(z,4))*
			     Sin((ArcTan(-1 + r,z) + ArcTan(1 + r,z))/2.)))/
	(Power(Power(-1 + r,2) + Power(z,2),0.75)*
	 Power(Power(1 + r,2) + Power(z,2),0.75)*
	 (Power(r,2) + Power(z,2) + 
	  Sqrt(Power(-1 + r,2) + Power(z,2))*Sqrt(Power(1 + r,2) + Power(z,2)) + 
	  2*Power(Power(-1 + r,2) + Power(z,2),0.25)*
	  Power(Power(1 + r,2) + Power(z,2),0.25)*
	  (r*Cos((ArcTan(-1 + r,z) + ArcTan(1 + r,z))/2.) + 
	   z*Sin((ArcTan(-1 + r,z) + ArcTan(1 + r,z))/2.))));

      // --------------------
      // dzeta/dr
      dobl_sph_dcyl(1,0) = (-(r*z*Power(Power(-1 + r,2) + Power(z,2),0.25)*
			      Power(Power(1 + r,2) + Power(z,2),0.25)*
			      (1 + 2*Power(r,2) + 2*Power(z,2) + 
			       2*Sqrt(Power(-1 + r,2) + Power(z,2))*
			       Sqrt(Power(1 + r,2) + Power(z,2)))) - 
			    z*(3*Power(r,4) + (1 + Power(z,2))*
			       (1 + Power(z,2) + Sqrt(Power(-1 + r,2) + Power(z,2))*
				Sqrt(Power(1 + r,2) + Power(z,2))) + 
			       Power(r,2)*(-2 + 4*Power(z,2) + 
					   Sqrt(Power(-1 + r,2) + Power(z,2))*Sqrt(Power(1 + r,2) + Power(z,2))
				 ))*Cos((ArcTan(-1 + r,z) + ArcTan(1 + r,z))/2.) + 
			    r*(-1 + Power(r,2) - (3 + 2*Power(r,2))*Power(z,2) - 2*Power(z,4))*
			    Sin((ArcTan(-1 + r,z) + ArcTan(1 + r,z))/2.))/
	(Power(Power(-1 + r,2) + Power(z,2),0.75)*
	 Power(Power(1 + r,2) + Power(z,2),0.75)*
	 Power(Power(r + Power(Power(-1 + r,2) + Power(z,2),0.25)*
		     Power(Power(1 + r,2) + Power(z,2),0.25)*
		     Cos((ArcTan(-1 + r,z) + ArcTan(1 + r,z))/2.),2) + 
	       Power(z + Power(Power(-1 + r,2) + Power(z,2),0.25)*
		     Power(Power(1 + r,2) + Power(z,2),0.25)*
		     Sin((ArcTan(-1 + r,z) + ArcTan(1 + r,z))/2.),2),1.5));
  
      // --------------------
      // dzeta_dz
      dobl_sph_dcyl(1,2) = (4*Power(r,6) + (1 + Power(z,2))*
			    (-1 - Power(z,2) + Sqrt(Power(-1 + r,2) + Power(z,2))*
			     Sqrt(Power(1 + r,2) + Power(z,2))) + 
			    Power(r,4)*(-9 + 8*Power(z,2) + 
					4*Sqrt(Power(-1 + r,2) + Power(z,2))*Sqrt(Power(1 + r,2) + Power(z,2)))
			    + Power(r,2)*(6 + 4*Power(z,4) - 
					  5*Sqrt(Power(-1 + r,2) + Power(z,2))*
					  Sqrt(Power(1 + r,2) + Power(z,2)) + 
					  Power(z,2)*(6 + 4*Sqrt(Power(-1 + r,2) + Power(z,2))*
						      Sqrt(Power(1 + r,2) + Power(z,2)))) + 
			    2*r*Power(Power(-1 + r,2) + Power(z,2),0.25)*
			    Power(Power(1 + r,2) + Power(z,2),0.25)*
			    ((2 + 3*Power(r,4) - Sqrt(Power(-1 + r,2) + Power(z,2))*
			      Sqrt(Power(1 + r,2) + Power(z,2)) + 
			      Power(z,2)*(3 + Power(z,2) + 
					  Sqrt(Power(-1 + r,2) + Power(z,2))*
					  Sqrt(Power(1 + r,2) + Power(z,2))) + 
			      Power(r,2)*(-5 + 4*Power(z,2) + 
					  Sqrt(Power(-1 + r,2) + Power(z,2))*
					  Sqrt(Power(1 + r,2) + Power(z,2))))*
			     Cos((ArcTan(-1 + r,z) + ArcTan(1 + r,z))/2.) + 
			     2*r*z*(Power(r,2) + Power(z,2))*
			     Sin((ArcTan(-1 + r,z) + ArcTan(1 + r,z))/2.)))/
	(2.*(Power(-1 + r,2) + Power(z,2))*(Power(1 + r,2) + Power(z,2))*
	 Power(Power(r + Power(Power(-1 + r,2) + Power(z,2),0.25)*
		     Power(Power(1 + r,2) + Power(z,2),0.25)*
		     Cos((ArcTan(-1 + r,z) + ArcTan(1 + r,z))/2.),2) + 
	       Power(z + Power(Power(-1 + r,2) + Power(z,2),0.25)*
		     Power(Power(1 + r,2) + Power(z,2),0.25)*
		     Sin((ArcTan(-1 + r,z) + ArcTan(1 + r,z))/2.),2),1.5));

      // dtheta/dtheta
      dobl_sph_dcyl(2,1) = 1.0;
    }

  
  // the oblate spheroidal coordinates
  double lambda;
  double zeta;
  double theta;
  
  double r;
  double z;
};

// Curvilinear coordinates which parameterise the surface of the disk;
// xi1 is the 'radial' coordinate, xi2 is 'azimuthal' and xi3 is the surface normal.
// 
// N.B.: Deliberately not overloading the [] operator so that the user is
// forced to type .xi1 etc. to prevent any confusion with global Cartesian coordinates
struct LagrangianCoordinates
{
  // default constructor that initialises everything to zero
LagrangianCoordinates() : xi1(0.0), xi2(0.0), xi3(0.0), ncoord(3) {}

  // optional constructor with doubles
LagrangianCoordinates(const double& xi1_, const double& xi2_, const double& xi3_) :
  xi1(xi1_), xi2(xi2_), xi3(xi3_), ncoord(3) {}

  LagrangianCoordinates(const Vector<double>& xi)
    {
      if(xi.size() != 3)
      {
	ostringstream error_message;
	error_message << "Error, wrong number of Lagrangian coordinates; expecting 3, given "
		      << xi.size();
	throw OomphLibError(error_message.str(),
			    OOMPH_CURRENT_FUNCTION,
			    OOMPH_EXCEPTION_LOCATION);
      }

      xi1 = xi[0];
      xi2 = xi[1];
      xi3 = xi[2];

      ncoord = 3;
    }
  
  // copy constructor
LagrangianCoordinates(const LagrangianCoordinates& lagr_coords) :
  xi1(lagr_coords.xi1), xi2(lagr_coords.xi2), xi3(lagr_coords.xi3), ncoord(3) {}

  // assignment operator
  void operator=(const LagrangianCoordinates& lagr_coords)
    {
      xi1 = lagr_coords.xi1;
      xi2 = lagr_coords.xi2;
      xi3 = lagr_coords.xi3;

      ncoord = 3;
    }

    // helper to get the three components as a vector
    Vector<double> vector() const
    {
      Vector<double> vec(3, 0.0);
      vec[0] = xi1;
      vec[1] = xi2;
      vec[2] = xi3;

      return vec;
    }
  
  // the curvilinear coordinates
  double xi1;
  double xi2;
  double xi3;

  // separately keep track of the sign of xi3; this allows for a cleaner
  // representation of points on the underside of the disk without having
  // to subtract small amounts - effectively allows for a signed zero
  int sign_of_xi3() const
    {
      return std::signbit(xi3) ? -1 : 1;
    }
  
  // number of coordinates
  unsigned ncoord;
};

namespace Analytic_Functions
{
  // reflect an angle the z-axis
  inline double reflect_angle_wrt_z_axis(const double& theta)
  {
    // flip by subtracting from pi, then make map back to [0,2pi)
    return map_angle_to_range_0_to_2pi(MathematicalConstants::Pi - theta);
  }

  inline DenseMatrix<double> strain_rate(const DenseMatrix<double>& du_dx)
  {
    DenseMatrix<double> strain_rate(3,3,0);
    
    for(unsigned i=0; i<3; i++)
    {
      for(unsigned j=0; j<3; j++)
      {
	strain_rate(i,j) = 0.5*(du_dx(i,j) + du_dx(j,i));
      }
    }

    return strain_rate;
  }
  
  /// \short Newtonian stress tensor
  inline DenseMatrix<double> stress(const DenseMatrix<double>& strain_rate,
				    const double& p)
  {
    const unsigned dim = 3;
      
    // stress tensor
    DenseMatrix<double> stress(dim,dim);
    
    for(unsigned i=0; i<dim; i++)
    {
      for(unsigned j=0; j<dim; j++)
      {
	// Newtonian constitutive relation
	stress(i,j) = -p*delta(i,j) + 2.0*strain_rate(i,j);
      }
    }

    return stress;
  }
  
} // end of Analytic_Functions namespace

inline void matrix_inverse(const DenseMatrix<double>& M, DenseMatrix<double>& M_inv)
{
  // Calculate the determinant of the matrix
  const double det = M(0,0) * M(1,1) * M(2,2) 
    + M(0,1) * M(1,2) * M(2,0) 
    + M(0,2) * M(1,0) * M(2,1) 
    - M(0,0) * M(1,2) * M(2,1) 
    - M(0,1) * M(1,0) * M(2,2) 
    - M(0,2) * M(1,1) * M(2,0);
    
  //Calculate the inverse of the Lagrangian derivatives of the Eulerian coords
  M_inv(0,0) =  (M(1,1) * M(2,2) 
		 - M(1,2) * M(2,1)) / det; 
  M_inv(0,1) = -(M(0,1) * M(2,2) 
		 - M(0,2) * M(2,1)) / det; 
  M_inv(0,2) =  (M(0,1) * M(1,2) 
		 - M(0,2) * M(1,1)) / det; 
  M_inv(1,0) = -(M(1,0) * M(2,2) 
		 - M(1,2) * M(2,0)) / det; 
  M_inv(1,1) =  (M(0,0) * M(2,2) 
		 - M(0,2) * M(2,0)) / det; 
  M_inv(1,2) = -(M(0,0) * M(1,2) 
		 - M(0,2) * M(1,0)) / det; 
  M_inv(2,0) =  (M(1,0) * M(2,1) 
		 - M(1,1) * M(2,0)) / det; 
  M_inv(2,1) = -(M(0,0) * M(2,1) 
		 - M(0,1) * M(2,0)) / det; 
  M_inv(2,2) =  (M(0,0) * M(1,1) 
		 - M(0,1) * M(1,0)) / det;  
}

#endif
