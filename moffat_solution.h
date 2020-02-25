// =============================================================================
// Christian Vaquero-Stainer
// 17/12/19
//
// The solution for the flow around a semi-infinite plate at rest at theta=0,2pi
// (see Panton)
//
// A and B are constants determined by the far-field flow which determine the
// angle of separation of the flow.
//
// A = 0 corresponds to normal incident flow with attachment streamlines,
// B=0 corresponds to zero angle of attack, i.e. attachment streamline at
// theta=pi. 
// The angle of separation is a function of the ratio A/B 
// =============================================================================
#ifndef OOMPH_MOFFAT_SOLUTION_HEADER
#define OOMPH_MOFFAT_SOLUTION_HEADER

namespace Analytic_Functions
{
  // solution in the original coordinate frame, i.e. a semi-infinite plate at x>0
  void moffat_solution(const double& r, const double& theta,
		       const double& A, const double& B,
		       Vector<double>& u_polar, DenseMatrix<double>& grad_u_polar)
  {
    // make sure we've gone enough space
    u_polar.resize(3);
    
    // radial velocity component
    u_polar[0] = sqrt(r) * (+ A * (-(3.0/2.0)*sin((3.0/2.0)*theta) + 0.5*sin(0.5*theta))
			    - B * ((3.0/2.0)*cos((3.0/2.0)*theta) - (3.0/2.0)*cos(0.5*theta)) );

    // azimuthal velocity component
    u_polar[1] = -(3.0/2.0)*sqrt(r) * (+ A * (cos((3.0/2.0)*theta) - cos(0.5*theta))
				       - B * (sin((3.0/2.0)*theta) - 3*sin(0.5*theta)) );
    
    // analytic result for p
    u_polar[2] = -(2.0/sqrt(r)) * ( A * sin(theta/2.0) + 3.0*B*cos(theta/2.0) );
    
    // ------------------------------
    // matrix of velocity derivatives in polar coordinates

    // Du_r/Dr
    grad_u_polar(0,0) = -(sin(theta/2.)*(A + 3*A*cos(theta) - 3*B*sin(theta)))/(2.*sqrt(r));
    
    // Du_r/Dtheta
    grad_u_polar(0,1) = (sqrt(r)*(A*cos(theta/2.) - 9*A*cos((3*theta)/2.)
				  - 3*B*sin(theta/2.) + 9*B*sin((3*theta)/2.)))/4.;
    
    // Du_theta/Dr
    grad_u_polar(1,0) = (3*pow(sin(theta/2.),2)*(A*cos(theta/2.) - B*sin(theta/2.)))/sqrt(r);
    
    // Du_theta/Dtheta
    grad_u_polar(1,1) = (3*sqrt(r)*sin(theta/2.)*(A + 3*A*cos(theta) - 3*B*sin(theta)))/2.;
  }
}

#endif
