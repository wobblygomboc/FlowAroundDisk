//LIC// ====================================================================
//LIC// This file forms part of oomph-lib, the object-oriented, 
//LIC// multi-physics finite-element library, available 
//LIC// at http://www.oomph-lib.org.
//LIC// 
//LIC//    Version 1.0; svn revision $LastChangedRevision: 1307 $
//LIC//
//LIC// $LastChangedDate: 2018-01-18 11:30:14 +0000 (Thu, 18 Jan 2018) $
//LIC// 
//LIC// Copyright (C) 2006-2016 Matthias Heil and Andrew Hazel
//LIC// 
//LIC// This library is free software; you can redistribute it and/or
//LIC// modify it under the terms of the GNU Lesser General Public
//LIC// License as published by the Free Software Foundation; either
//LIC// version 2.1 of the License, or (at your option) any later version.
//LIC// 
//LIC// This library is distributed in the hope that it will be useful,
//LIC// but WITHOUT ANY WARRANTY; without even the implied warranty of
//LIC// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
//LIC// Lesser General Public License for more details.
//LIC// 
//LIC// You should have received a copy of the GNU Lesser General Public
//LIC// License along with this library; if not, write to the Free Software
//LIC// Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA
//LIC// 02110-1301  USA.
//LIC// 
//LIC// The authors may be contacted at oomph-lib@maths.man.ac.uk.
//LIC// 
//LIC//====================================================================

// ============================================================================
// custom defines
// ============================================================================

#define xPRINT_SINGULAR_JACOBIAN

#define xCOMPUTE_INITIAL_JACOBIAN_EIGENFUNCTIONS

// do we want to duplicate the edge nodes (for better output)?
#define DUPLICATE_EDGE_NODES

// wraps the problem.newton_solve() in a try/catch block and prints jacobian
// (but cocks up output of other exceptions)
#define xPRINT_SINGULAR_JACOBIAN

// QUEHACERES for debug - shift the lower plate nodes down by a small amount
// so the the upper and lower parts of the plate are physically separated
#define SHIFT_LOWER_PLATE_NODES_BY_DZ

#ifdef SHIFT_LOWER_PLATE_NODES_BY_DZ
const double lower_plate_z_shift = -0.987e-8;
#endif

// NNODE_1D for the singular line elements - may affect stability of the scheme
#define SINGULAR_ELEMENT_NNODE_1D 3
 
#define xUSE_FD_JACOBIAN

// ============================================================================
// Includes
// ============================================================================

#include <fenv.h>

// needed to reset the FPU control word which Triangle messes around with
#include <fpu_control.h>

//Generic routines
#include "generic.h"

#include "external_src/oomph_superlu_4.3/slu_ddefs.h"

// singular elements
#include "navier_stokes_sing_face_element.h"

#include "warped_disk_with_torus_geom_obj.h"

// The mesh
#include "meshes/triangle_mesh.h"
 
// Get the mesh
#include "meshes/tetgen_mesh.h" 
#include "meshes/refineable_tetgen_mesh.h"
#include "meshes/gmsh_tet_mesh.h"

// Get the faceted surfaces
#include "tetmesh_faceted_surfaces.h"

// analytic solution for in-plane motion based on integrals
#include "sherwood_solution.h"

// include classes for vector and matrix algebra
#include "additional_maths.h"

// definitions to convert CForm output from Mathematica
#include "mathematica_definitions.h"

// Exact solutions for the 4 distinct modes of rigid body motion
// of a flat disk
#include "exact_solutions_finite_disk.h"

// functions to convert coordinates, velocities and gradients between
// the Lagrangian edge coordinates and the global Cartesian system
#include "coordinate_conversions.h"

using namespace oomph;

// ============================================================================
// Globals
// ============================================================================

double TetMeshBase::Tolerance_for_boundary_finding = 1.0e-8;
bool TetMeshBase::Shut_up_about_nonplanar_boundary_nodes = true;

///////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////

//=============================================================
/// Namespace for problem parameters
//=============================================================
namespace Global_Parameters
{
  // QUEHACERES hacky, but makes life easier for the hacky FSI fixed-point loop
  Problem* Problem = nullptr;
  
  // QUEHACERES delete after debug
  double broadside_amplitude_pin          = 0.0;
  double broadside_rotation_amplitude_pin = 0.0;
  double in_plane_amplitude_pin           = 0.0;
  double in_plane_rotation_amplitude_pin  = 0.0;
    
  string output_directory = "RESLT";
  
  /// (Half-)width of the box
  double Box_half_width = 1.5;

  /// (Half)height of the box
  double Box_half_height = 0.5; //1.0;

  // velocity of the whole disk (rigid)
  Vector<double> u_disk_rigid_translation(3, 0.0);

  // angular velocity of whole disk about positive Cartesian axes
  Vector<double> omega_disk(3, 0.0);

  // do we want to enforce no-slip (and on top boundary, zero-traction)?
  // if this is false then we're doing exact flat disk BCs
  bool No_slip_boundaries = false;
  
  // amplitude and wavenumber of the warped disk
  double Epsilon = 0;
  unsigned n = 5;

  // amount of pressure regularisation to add to the functional which we're
  // minimising in the torus region
  double Pressure_regularisation_factor = 0.1;

  double Velocity_regularisation_factor = 1.0;
  
  // Convergence tolerance for the outer FSI Newton solve
  double FSI_convergence_tol = 1e-6;

  // Finite-diff step size for outer FSI solve (problem is linear, so
  // let's default to a nice big step to avoid noise)
  double FSI_FD_step = 1e-2;
  
  // QUEHACERES pressure zero-level for debug (this is constant for broadside, will
  // need to vary (probably?) for in-plane)
  double p0 = 0;

  double p_offset = 0;
  
  // cross-sectional radius of the torus
  double R_torus = 0.2;

  // the number of line segments making up half the perimeter of the disk
  unsigned Half_nsegment_disk = 15;

  // number of vertices on the cross-sectional circles of the torus
  unsigned Nvertex_torus = 8; 

  // number of elements in the singular line mesh which holds the
  // singular amplitudes.
  // Default: same number of elements as on the edge of the disk
  unsigned Nsingular_line_element = Half_nsegment_disk * 4 - 1;

  // relative density difference between the fluid and the elastic sheet,
  // i.e. (rho_s - rho_f) / rho_f
  double Relative_density_difference = 0.0;

  // Balance of (buoyancy-reduced) gravity vs viscosity
  double Archimedes_number = 0.0;

  // Non-dimensional mass of the disk
  double Mass = 1.0;
  
  unsigned Nplot_for_bulk = 5;

  // how much to shift the edge radius for outputting "infinite" pressure
  // (this will be updated based on the average element size in the torus)
  double Drho_for_infinity = 1e-5;
  
  // size of the radially aligned disks used for outputting the solution at points
  // around the disk edge
  // double disk_on_disk_radius = 0.2;

  // zero is default, i.e. maximum possible element volume which respects
  // torus boundary discretisation
  double Target_element_volume_in_torus_region = 0;

  double Initial_global_element_volume = 1.0;
  
  // fake amplitudes to impose for debug
  double Imposed_singular_amplitude_broadside = 0;
  double Imposed_singular_amplitude_in_plane  = 0;
  
  // split corner elements which have all their nodes pinned on the outer boundary
  bool Split_corner_elements = true; 

  // do the disk-on-disk plots (expensive so default off)
  bool Do_disk_on_disk_plots = false;
  
  // store the warped disk object so that we can use it to get
  // surface normals
  WarpedCircularDiskWithAnnularInternalBoundary* Warped_disk_with_boundary_pt = nullptr;

  MeshAsGeomObject* mesh_as_geom_object_pt = nullptr;

  bool Only_subtract_first_singular_term = false;

  // brutal global finite-diff solver
  bool Use_fd_lu_solver = false;

  // --------------------------------------------------------------------------

  // function to compute the body force, which is the reduced gravity
  // (gravity scaled by the relative difference in densities between the
  // fluid and solid)
  void reduced_gravity(const Vector<double>& x, Vector<double>& body_force)
  {
    body_force.resize(3, 0.0);

    body_force[0] = 0.0;
    body_force[1] = 0.0;
    body_force[2] = -Relative_density_difference;
  }
  
  // compute the total disk velocity u = u_trans + u_rot
  Vector<double> disk_velocity(const Vector<double>& x)
  {
    // total velocity
    Vector<double> u_disk(3, 0.0);

    // velocity associated with rigid-body rotation 
    Vector<double> u_rotational(3, 0.0);
	
    // u = \omega \ctimes r
    u_rotational[0] = omega_disk[1] * x[2] - omega_disk[2] * x[1];
    u_rotational[1] = omega_disk[2] * x[0] - omega_disk[0] * x[2];
    u_rotational[2] = omega_disk[0] * x[1] - omega_disk[1] * x[0];

    // total velocity is the sum of the velocities associated with
    // rigid body rotation and translation
    u_disk[0] = u_rotational[0] + u_disk_rigid_translation[0];
    u_disk[1] = u_rotational[1] + u_disk_rigid_translation[1];
    u_disk[2] = u_rotational[2] + u_disk_rigid_translation[2];
    
    return u_disk;
  }
  
  Vector<double> compute_singular_amplitudes_from_disk_velocity(const double& zeta)
  {
  
    // broadside speed is the z-velocity component
    double u_broadside =  u_disk_rigid_translation[2];
    
    // in-plane speed is the magnitude of the x-y velocity
    double u_in_plane = sqrt( pow(u_disk_rigid_translation[0], 2) +
			      pow(u_disk_rigid_translation[1], 2) );

    // azimuthal angle of the in-plane translation vector
    double zeta_translation = atan2pi(u_disk_rigid_translation[1],
				      u_disk_rigid_translation[0]);
  
    // azimuthal angle of the out-of-plane rotation vector
    double zeta_rotation = atan2pi(omega_disk[1],
				   omega_disk[0]);

    // magnitude of the out-of-plane rotation vector
    double omega_out_of_plane = sqrt(pow(omega_disk[0], 2) +
				     pow(omega_disk[1], 2));


    // broadside amplitude is the broadside translation plus the
    // modulated contribution from out-of-plane rotations (with the factor of 2
    // which appears for the out-of-plane solution to be locally equivalent
    // to broadside motion)
    double c_broadside = u_broadside + 2 * omega_out_of_plane * sin(zeta - zeta_rotation);

    // in-plane amplitude is the in-plane speed modulated by
    // a in-phase function of the angle of this point relative to the
    // translation vector
    double c_in_plane = u_in_plane * cos(zeta - zeta_rotation);
  
    // in-plane rotation amplitude is the in-plane speed modulated by
    // a function pi/2 out of phase with the angle of this point relative to the
    // translation vector
    double c_in_plane_rotation = u_in_plane * sin(zeta - zeta_translation);

    // now package them up for return
    Vector<double> amplitudes(3, 0.0);

    // if we're doing the exact solution, then we don't want any modulation
    if(CommandLineArgs::command_line_flag_has_been_set("--subtract_exact_solution"))
    {
      c_broadside = u_broadside > 0.0 ? 1 : 0;
      c_in_plane = u_in_plane > 0.0 ? 1 : 0;
      c_in_plane_rotation = 0; // QUEHACERES not doing this for the time being
    }

    amplitudes[0] = c_broadside;
    amplitudes[1] = c_in_plane;
    amplitudes[2] = c_in_plane_rotation;
  
    return amplitudes;
  }
  
  // hacky - the asymptotic singular functions to subtract
#include "singular_functions.h"
}

///////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////
namespace Analytic_Functions
{
  // derivative of the functional which we're minimising, w.r.t. the solution
  Vector<double> dfunctional_du(const Vector<double>& u)
  {
    Vector<double> dpi_du(u.size(), 0.0);

    // // QUEHACERES just use zero for now
    // functional \Pi = sum_j 1/2|u_j|^2,
    // so d\Pi/d u_i = |u_i|
    for(unsigned i=0; i<u.size(); i++)
      dpi_du[i] = (u[i]); // abs(u[i]); 

    return dpi_du;
  }
  
  // exact solution for the flat disk, handling linear combinations of
  // in-plane and broadside motion
  void exact_solution_flat_disk(const Vector<double>& x,
				Vector<double>& u_exact)
  {
    u_exact.resize(4, 0.0);

    // forward
    FlatDiskExactSolutions::total_exact_solution(x,
						 Global_Parameters::u_disk_rigid_translation,
						 Global_Parameters::omega_disk,
						 u_exact);
  }

  // 'exact' velocity gradients from linear combination of rigid-body modes
  // (computed via finite-diff of exact solution)
  void exact_velocity_gradient_flat_disk(const Vector<double>& x,
					 DenseMatrix<double>& dudx)
  {
    dudx.resize(3, 3, 0.0);
   
    // forward
    FlatDiskExactSolutions::total_exact_velocity_gradient(x,
    							  Global_Parameters::u_disk_rigid_translation,
    							  Global_Parameters::omega_disk,
    							  dudx);    
  }

  // generic function to compute traction;
  // get the velocity gradient from the selected analytic solution,
  // compute the strain rate, the stress and then the traction
  void prescribed_traction(const double& t,
			   const Vector<double>& x,
			   const Vector<double>& outer_unit_normal,
			   Vector<double>& traction)
  {
    // make sure we've got enough space
    traction.resize(3,0);

    // zero out the traction
    std::fill(traction.begin(), traction.end(), 0.0);

    // if we're doing the no-slip problem, then we're enforcing
    // zero traction on the top surface, so leave the traction empty
    if(Global_Parameters::No_slip_boundaries)
      return;
     
    Vector<double> u(4,0);
    DenseMatrix<double> du_dx(3,3,0);

    // get the linear combination of all modes
    exact_solution_flat_disk(x, u);

    // get the total velocity gradients
    exact_velocity_gradient_flat_disk(x, du_dx);
    
    // interpret the total pressure from the total solution
    double p = u[3] + Global_Parameters::p_offset;       

    // compute the strain rate
    DenseMatrix<double> total_strain_rate = strain_rate(du_dx);
    
    // compute the stress from the strain rate
    DenseMatrix<double> total_stress = stress(total_strain_rate, p);

    // compute the traction t_i = tau_ij n_j
    for(unsigned i=0; i<3; i++)
    {
      for(unsigned j=0; j<3; j++)
      {
	traction[i] += total_stress(i,j) * outer_unit_normal[j];
      }
    }    
  }
}

///////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////

//=============================================================================
// face elements, for attaching to boundaries to measure surface areas and
// to output values (to check boundary conditions)
//=============================================================================
template <class ELEMENT>
class RigidDiskFaceElement : public virtual FaceGeometry<ELEMENT>, 
				public virtual FaceElement
{
 
public:

  ///Constructor, which takes a "bulk" element and the value of the index
  ///and its limit
  RigidDiskFaceElement(FiniteElement* const& element_pt, 
			  const int& face_index) : 
    FaceGeometry<ELEMENT>(), FaceElement()
    { 
      //Attach the geometrical information to the element. N.B. This function
      //also assigns nbulk_value from the required_nvalue of the bulk element
      element_pt->build_face_element(face_index, this);
      
      //Set the dimension from the dimension of the first node
      Dim = node_pt(0)->ndim();

      // flip the sign of the normal, since outer_unit_normal() gives the normal
      // which points out of the bulk element we're attached to, but since
      // we're attached to the disk we want to point into the bulk element,
      // i.e. if we're on the top surface of the disk the bulk element is above
      // and so the default is for the normal to point downwards, but the upper
      // disk normal points upwards
      normal_sign() = -normal_sign();
    }

  RigidDiskFaceElement(const RigidDiskFaceElement&)
    {
      BrokenCopy::broken_copy("RigidDiskFaceElement");
    }

  void operator=(const RigidDiskFaceElement&)
    {
      BrokenCopy::broken_assign("RigidDiskFaceElement");
    }
  
  /// \short Output function:  
  ///  Output x, y, u, v, p, theta at Nplot^DIM plot points
  // Start of output function
  void output(std::ostream &outfile, const unsigned &nplot)
    {

      // Get pointer to assocated bulk element
      ELEMENT* bulk_el_pt = dynamic_cast<ELEMENT*>(bulk_element_pt());
   
      //Local coordinates
      Vector<double> s(Dim-1);

      Vector<double> s_bulk(Dim);
   
      //Eulerian coordinates
      Vector<double> x(Dim, 0.0);

      // Velocity from bulk element
      Vector<double> velocity(Dim);
      
      // Tecplot header info
      outfile << this->tecplot_zone_string(nplot);

      // Loop over plot points
      unsigned num_plot_points=this->nplot_points(nplot);
      for (unsigned iplot=0; iplot<num_plot_points; iplot++)
      {
	// Get local coordinates of plot point
	this->get_s_plot(iplot,nplot,s);
     
	this->get_local_coordinate_in_bulk(s, s_bulk);
     
	//Get x position from bulk
	bulk_el_pt->interpolated_x(s_bulk, x);
	
	bulk_el_pt->interpolated_u_nst(s_bulk,velocity);
	
	// output Eulerian coordinates
	for(unsigned i=0; i<Dim; i++)
	{
	  outfile << x[i] << " ";
	}

	// output fluid velocities
	for(unsigned i=0; i<Dim; i++)
	{
	  outfile << velocity[i] << " ";
	}
     
	// Output the fluid pressure
	outfile << bulk_el_pt->interpolated_p_nst(s_bulk) << std::endl;
      }

      this->write_tecplot_zone_footer(outfile,nplot);
      
    } //End of output function

  void interpolated_x(const Vector<double>& s, Vector<double>& x)
    {
      // number of dimensions in the bulk
      unsigned dim_bulk = Dim + 1;
      
      // local coordinates in bulk element
      Vector<double> s_bulk(dim_bulk);
      
      // Get pointer to assocated bulk element
      ELEMENT* bulk_el_pt = dynamic_cast<ELEMENT*>(bulk_element_pt());

      this->get_local_coordinate_in_bulk(s, s_bulk);
      
      //Get x position from bulk
      bulk_el_pt->interpolated_x(s_bulk, x);      
    }

  void interpolated_u_nst(const Vector<double>& s, Vector<double>& u) const
    {
      // number of dimensions in the bulk
      unsigned dim_bulk = Dim + 1;
      
      // local coordinates in bulk element
      Vector<double> s_bulk(dim_bulk);
      
      // Get pointer to assocated bulk element
      ELEMENT* bulk_el_pt = dynamic_cast<ELEMENT*>(bulk_element_pt());

      this->get_local_coordinate_in_bulk(s, s_bulk);
      
      //Get velocity from bulk
      bulk_el_pt->interpolated_u_nst(s_bulk, u);      
    }

  double interpolated_p_nst(const Vector<double>& s) const
    {
      // number of dimensions in the bulk
      unsigned dim_bulk = Dim + 1;
      
      // local coordinates in bulk element
      Vector<double> s_bulk(dim_bulk);
      
      // Get pointer to assocated bulk element
      ELEMENT* bulk_el_pt = dynamic_cast<ELEMENT*>(bulk_element_pt());

      this->get_local_coordinate_in_bulk(s, s_bulk);
      
      //Get pressure from bulk
      return bulk_el_pt->interpolated_p_nst(s_bulk);
    }

  void strain_rate(const Vector<double>& s, DenseMatrix<double>& _strain_rate) const
    {
      // number of dimensions in the bulk
      unsigned dim_bulk = Dim + 1;
      
      // local coordinates in bulk element
      Vector<double> s_bulk(dim_bulk);
      
      // Get pointer to assocated bulk element
      ELEMENT* bulk_el_pt = dynamic_cast<ELEMENT*>(bulk_element_pt());

      this->get_local_coordinate_in_bulk(s, s_bulk);
      
      //Get pressure from bulk
      return bulk_el_pt->strain_rate(s_bulk, _strain_rate);
    }
  
  Vector<double> get_contribution_to_normal_stress() const
    {
      //Find out how many nodes there are
      const unsigned n_node = this->nnode();

      // shorthand
      const unsigned Dim = this->Dim;
    
      //Set up memory for the shape and test functions
      Shape psi(n_node), test(n_node);
     
      //Set the value of Nintpt
      const unsigned n_intpt = this->integral_pt()->nweight();
     
      //Set the Vector to hold local coordinates
      Vector<double> s(Dim-1);

      // Saves result of integration
      Vector<double> total_force(Dim, 0.0);
      
      //Loop over the integration points
      for(unsigned ipt=0; ipt<n_intpt; ipt++)
      {       
	//Assign values of s
	for(unsigned i=0; i<(Dim-1); i++)
	{
	  s[i] = this->integral_pt()->knot(ipt,i);
	}

	//Get the integral weight
	double w = this->integral_pt()->weight(ipt);
       
	//Find the shape and test functions and return the Jacobian
	//of the mapping
	double J = this->J_eulerian(s); 
       
	//Premultiply the weights and the Jacobian
	double W = w*J;
      	  
	// Compute outer unit normal at the specified local coordinate
	// to compute scaled and unscaled flux of singular solution
	Vector<double> unit_normal(Dim);
	this->outer_unit_normal(s, unit_normal);

	// strain rate tensor
	DenseMatrix<double> strain_rate_fe(Dim, Dim, 0.0);

	strain_rate(s, strain_rate_fe);

	// get the pressure
	double p = interpolated_p_nst(s);
	  
	// get contribution to the total FE stress on this element
	DenseMatrix<double> stress(Dim, Dim);

	stress = Analytic_Functions::stress(strain_rate_fe, p);

	// now add the weighted traction to the total force
	for(unsigned i=0; i<Dim; i++)
	{
	  for(unsigned j=0; j<Dim; j++)
	  {
	    total_force[i] += stress(i,j) * unit_normal[j] * W;
	  }
	}
      } // end loop over integration points

      return total_force;
    }
   
      
private:
  
  unsigned Dim;
};

///////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////

//=========================================================================
/// Class that solves Navier-Stokes flow around a 2D disk
//=========================================================================
template <class ELEMENT>
class FlowAroundDiskProblem : public Problem
{
  
public:

  /// Constructor
  FlowAroundDiskProblem();
  
  /// Destructor
  ~FlowAroundDiskProblem()
    {
      //Delete the objects
      unsigned nh = Inner_boundary_pt.size();
      for(unsigned h=0; h<nh; h++)
      {
	delete Inner_boundary_pt[h];
      }
      delete Outer_boundary_pt;
    }
  
  /// Actions before adapt (empty)
  void actions_before_adapt()
    {}

  /// Totally new mesh; build elements and apply boundary conditions
  void actions_after_adapt()
    {
      // Complete problem setup
      complete_problem_setup();
    }
 
  /// Update the problem specs before solve: (empty)
  void actions_before_newton_solve(){}

  /// Update the problem specs before solve (empty)
  void actions_after_newton_solve(){}
 
  /// Doc the solution
  void doc_solution(const unsigned& nplot);

  DocInfo& doc_info()
    {
      return Doc_info;
    }

  // function to directly impose the singular amplitude and bypass the 
  // proper calculation
  void impose_fake_singular_amplitude(const bool& impose_zero_amplitude = false);

  /// Assign nodal values to be the exact singular solution for debug
  void set_values_to_singular_solution(const bool& broadside = true);

  /// \short set the nodal values to the exact non-singular part of the solution
  void set_values_to_exact_non_singular_solution() const;
  
  // function to validate the singular stress by assigning the singular solution
  // to each node, then computing the velocity gradients via finite difference and
  // comparing these to the analytic gradients.
  void validate_singular_stress(const bool& broadside = true);

  // for debug
  void validate_exact_solution_divergence() const;

  /// Helper to output the pin status of eaech nodal dof in each submesh
  void output_submesh_dof_pin_status() const;

  /// \short Helper to pin the dofs in the singular line mesh corresponding to
  /// a particular singular function
  void pin_singular_function(const unsigned& sing_fct_id, const double& val=0.0);

  // QUEHACERES for debug: compute the condition number of the Jacobian
  // matrix
  double compute_jacobian_condition_number();

  void compute_and_assign_smallest_eigensolution(const double& threshold=1e-12);

  void hacky_fixed_point_fsi_solve();

  void fsi_residual_sweep();
  
  /// QUEHACERES good comment here
  void get_total_disk_force_and_moment_residuals(
    const Vector<double>& centre_of_rotation,
    Vector<double>& residuals);

  /// Helper function to apply boundary conditions
  void apply_boundary_conditions();
  
private:
 
  /// Apply BCs and make elements functional
  void complete_problem_setup();
 
  void delete_face_elements()
    {
      // Loop over the flux elements
      unsigned n_element = Traction_boundary_condition_mesh_pt->nelement();
      for(unsigned e=0;e<n_element;e++)
      {
	// Kill
	delete Traction_boundary_condition_mesh_pt->element_pt(e);
      }
   
      // Wipe the mesh
      Traction_boundary_condition_mesh_pt->flush_element_and_node_storage();
      
      // Loop over the bc elements
      n_element = Face_mesh_for_bc_pt->nelement();
      for(unsigned e=0; e<n_element; e++)
      {
      	// Kill
      	delete Face_mesh_for_bc_pt->element_pt(e);
      }
   
      // Wipe the mesh
      Face_mesh_for_bc_pt->flush_element_and_node_storage();

      // Loop over the integral face elements
      n_element = Face_mesh_for_singularity_integral_pt->nelement();
      for(unsigned e=0; e<n_element; e++)
      {
	delete Face_mesh_for_singularity_integral_pt->element_pt(e);
      }
      
      Face_mesh_for_singularity_integral_pt->flush_element_and_node_storage();

      // delete stress jump elements
      n_element = Face_mesh_for_stress_jump_pt->nelement();
      for(unsigned e=0; e<n_element; e++)
      {
	delete Face_mesh_for_stress_jump_pt->element_pt(e);
      }
      Face_mesh_for_stress_jump_pt->flush_element_and_node_storage();

      // delete the singular line elements
      for(unsigned e=0; e<Singular_fct_element_mesh_pt->nelement(); e++)
      {
	// the destructor of the line element deletes the newly created nodes
	delete Singular_fct_element_mesh_pt->element_pt(e);
      }
      Singular_fct_element_mesh_pt->flush_element_and_node_storage();
      
      for(unsigned e=0; e<Singular_fct_element_mesh_lower_pt->nelement(); e++)
      {
	// the destructor of the line element deletes the newly created nodes
	delete Singular_fct_element_mesh_lower_pt->element_pt(e);
      }      
      Singular_fct_element_mesh_lower_pt->flush_element_and_node_storage();
      
      for(unsigned e=0; e<Singular_fct_element_mesh_upper_pt->nelement(); e++)
      {
	// the destructor of the line element deletes the newly created nodes
	delete Singular_fct_element_mesh_upper_pt->element_pt(e);
      }      
      Singular_fct_element_mesh_upper_pt->flush_element_and_node_storage();
    }
  
  /// \short Helper function to create the face elements needed to:
  /// - impose Dirichlet BCs
  /// - impose the additional traction required from the augmented region to
  ///   the surrounding bulk elements
  /// - compute the reciprocity integral to determine the singular amplitude
  void create_face_elements();
  
  // function to populate the vectors Elements_on_upper[lower]_disk_surface_pt
  // with the elements which have at least one node on the disk.
  // (These will be used when the plate nodes are duplicated, so that we
  // know which elements need to be told that their node pointers have changed).
  void identify_elements_on_upper_and_lower_disk_sufaces();

  void duplicate_plate_nodes_and_add_boundaries();

  /// \short Function to populate a map which maps each torus boundary
  /// element in the augmented region to it's corresponding boundary element
  /// in the bulk region
  void build_map_from_augmented_to_bulk_region_torus_boundary_elements();

  // function which takes a boundary element in the augmented region
  // and finds the corresponding boundary element in the bulk region
  void find_corresponding_element_on_bulk_side_of_augmented_boundary(
    ELEMENT* const& augmented_elem_pt,
    const int& face_index,
    ELEMENT*& corresponding_elem_pt) const;
  
  /// \short function to populate a map which maps each node in the mesh to
  /// a set of all the elements it is associated with
  void generate_node_to_element_map();

  /// \short function to setup the edge coordinate system (\rho, \zeta, \phi) for
  /// elements in the torus region. If use_plot_points is true, then the
  /// coordinates are computed for the elements plot points, otherwise it is
  /// for their knots (integration points).
  void setup_edge_coordinates_and_singular_element(
    FiniteElement* const& elem_pt,
    Vector<EdgeCoordinates>& edge_coordinates_at_knot_pt,
    Vector<std::pair<GeomObject*, Vector<double> > >& line_element_and_local_coordinate,
    const bool& use_plot_points = false) const;

  /// \short function to compute the edge coordinates and corresponding
  /// element and local coordinates in the singular line mesh
  /// (which computes the appropriate amplitude) from Eulerian coordinates
  void get_edge_coordinates_and_singular_element(
    FiniteElement* const& elem_pt,
    const Vector<double>& x,
    EdgeCoordinates& edge_coordinates_at_point,
    std::pair<GeomObject*, Vector<double> >& line_element_and_local_coordinate) const;
  
  /// \short Helper to generate the 1D line mesh which runs around the outside of the
  /// disk and provides the amplitude of the singularity as a function of the
  /// boundary coordinate, i.e. c(zeta_bound).
  void create_one_d_singular_element_mesh(const unsigned& nsingular_el);
  
  /// Setup disk on disk plots
  void setup_disk_on_disk_plots();

  // --------------------------------------------------------------------------
  // Meshes
  // --------------------------------------------------------------------------

  /// Bulk mesh
  RefineableTetgenMesh<ELEMENT>* Bulk_mesh_pt;

  /// \short Face element mesh which imposes the necessary traction
  /// onto the bulk elements on the boundary of the augmented region
  Mesh* Face_mesh_for_stress_jump_pt;
  
  /// \short Meshes of face elements used to compute the amplitudes of the singular
  /// functions (one mesh per singular function)
  Mesh* Face_mesh_for_singularity_integral_pt;

  /// \short Mesh of face
  Mesh* Traction_boundary_condition_mesh_pt;   

  /// Meshes containing the line elements which store the amplitudes of the
  /// singular functions
  Mesh* Singular_fct_element_mesh_upper_pt;
  Mesh* Singular_fct_element_mesh_lower_pt;

  /// Mesh which combines the upper and lower halves, so we don't have to repeat code
  Mesh* Singular_fct_element_mesh_pt;

  /// Mesh of face elements which impose Dirichlet boundary conditions in augmented region
  Mesh* Face_mesh_for_bc_pt;
  
  /// Mesh of elements within the torus region for the computation of Z2
  RefineableTetgenMesh<ELEMENT>* Torus_region_mesh_pt;
   
  /// \short Enumeration for IDs of FaceElements (used to figure out
  /// who's added what additional nodal data...)
  enum{ bla_hierher, Stress_jump_el_id=1, BC_el_id=1, Lambda_hat_hat_id=2 };

  // IDs to identify each singular function
  enum {Sing_fct_id_broadside=100,
	Sing_fct_id_broadside_rotation,
	Sing_fct_id_in_plane,
	Sing_fct_id_in_plane_rotation};

  // --------------------------------------------------------------------------
  
  /// Mesh as geom object representation of mesh  
  MeshAsGeomObject* Mesh_as_geom_object_pt;
  
  // Create the mesh as Geom Object
  MeshAsGeomObject* Face_mesh_as_geom_object_pt;

  // ### QUEHACERES delete
  // MeshAsGeomObject* Singular_line_mesh_upper_as_geom_object_pt;
  // MeshAsGeomObject* Singular_line_mesh_lower_as_geom_object_pt;
  MeshAsGeomObject* Singular_line_mesh_as_geom_object_pt;
  
  /// Storage for the outer boundary object
  TetMeshFacetedClosedSurface* Outer_boundary_pt;

  /// Inner boundary
  Vector<TetMeshFacetedSurface*> Inner_boundary_pt;
  
  /// First boundary ID for outer boundary
  unsigned First_boundary_id_for_outer_boundary;

  /// ID of the top boundary, used for when we want to do a prescribed traction problem
  unsigned Top_outer_boundary_id;

  /// ID of right boundary, i.e. with normal n = (1,0,0)
  unsigned Outer_traction_boundary_id;
  
  // Disk with torus round the edges
  //--------------------------------

  /// (zero-based) Region ID for torus around edge of warped disk
  unsigned Torus_region_id;

  /// First boundary ID for lower disk surface that is surrounded by torus
  unsigned First_lower_disk_boundary_id;
 
  /// Last boundary ID for lower disk surface that is surrounded by torus
  unsigned Last_lower_disk_boundary_id;

  /// First boundary ID for the upper disk surface
  unsigned First_upper_disk_boundary_id;

  /// Last boundary ID for the upper disk surface
  unsigned Last_upper_disk_boundary_id;
  
  /// First boundary ID for torus surrounding edge of disk
  unsigned First_torus_boundary_id;

  /// Last boundary ID for torus surrounding edge of disk
  unsigned Last_torus_boundary_id;
 
  /// \short Storage for one-based boundary IDs for boundaries on disk within
  ///  the torus region
  Vector<unsigned> One_based_boundary_id_for_disk_within_torus;

  /// \short Storage for one-based boundary IDs for boundaries on disk 
  /// outside the torus region
  Vector<unsigned> One_based_boundary_id_for_disk_outside_torus;

  Vector<unsigned> Boundary_id_for_upper_disk_within_torus;
  Vector<unsigned> Boundary_id_for_upper_disk_outside_torus;

  /// Combined list of boundary IDs (upper and lower disk) within the torus
  Vector<double> Disk_boundary_ids_in_torus;
  
  /// \short vectors to hold pointers to the elements which have at least one
  /// node on the disk surface
  std::set<ELEMENT*> Elements_on_upper_disk_surface_pt;
  std::set<ELEMENT*> Elements_on_lower_disk_surface_pt;

  // for debug, list of elements which are not considered "boundary" elements
  // as they do not have a face on the disk boundaries, but do have at least
  // one node on the disk
  std::set<ELEMENT*> Nonboundary_elements_with_node_on_upper_disk_surface_pt;
  std::set<ELEMENT*> Nonboundary_elements_with_node_on_lower_disk_surface_pt;

  /// \short a map which takes a node on a disk boundary and returns a set of
  /// the elements on the upper surface which share this node, and the index
  /// of this node within that element
  std::map<Node*, std::set<std::pair<ELEMENT*, unsigned> > >
  Disk_node_to_upper_disk_element_and_index_map;

  /// \short a map which gives a set of all the elements which are
  /// associated with a given node
  std::map<Node*, std::set<ELEMENT*> > Node_to_element_map;

  /// Plot points along a radial line
  Vector<std::pair<GeomObject*,Vector<double> > > Radial_sample_point_pt;

  /// \short map which takes each of the nodes in the augmented region
  /// and gives the edge coordinates (\rho, \zeta, \phi).
  std::map<Node*, Vector<double> >* Node_to_edge_coordinates_map_pt;

  /// \short Map of nodes duplicated by the stress-jump face elements;
  /// map is old -> new node
  std::map<Node*, Node*> Stress_jump_duplicate_node_map;

  /// \short Map which returns the corresponding bulk region element which
  /// shares a face with a given augmented region boundary element
  std::map<ELEMENT*, std::map<unsigned, ELEMENT*> > Torus_boundary_aug_to_non_aug_by_face_index_map;
  
  // QUEHACERES taking out the vectorisation for now
  // ///This is vectorised
  // /// since the nodes which have zero y-coordinate can have a zeta of either
  // /// 0 or 2pi, and the interpolation will need to take this into account.
  

  // -----------------------------------------------------
  // disk on disk output stuff
  
  /// Flag to force update on geom object representations 
  bool Geom_objects_are_out_of_date;
  
  /// The Line Visualiser.
  LineVisualiser* LV_pt;

  /// \short Number of "disks on disk" around the edge where solution is
  /// to be visualised
  unsigned Ndisk_on_disk_plot;

  /// \short Number of azimuthal plot points in "disks on disk" plots 
  /// around the edge where solution is to be visualised
  unsigned Nphi_disk_on_disk_plot;

  /// \short Number of radial plot points in "disks on disk" plots 
  /// around the edge where solution is to be visualised
  unsigned Nrho_disk_on_disk_plot;

  /// Disk_on_disk_plot_point[k_disk][i_rho][i_phi]
  Vector<Vector<Vector<std::pair<
			 Vector<double>,std::pair<GeomObject*,Vector<double> > > > > >
  Disk_on_disk_plot_point;
 
  /// Sanity check: Exact bounded volume
  double Exact_bounded_volume;

  /// Number of singular functions we will subtract
  unsigned Nsingular_function;
  
  /// Number of dimensions in the problem
  unsigned Dim;
  
  DocInfo Doc_info;

  unsigned Fsi_iteration_counter;
};



//========================================================================
/// Constructor
//========================================================================
template <class ELEMENT>
FlowAroundDiskProblem<ELEMENT>::FlowAroundDiskProblem() :
  Bulk_mesh_pt(nullptr), Face_mesh_for_stress_jump_pt(nullptr),
  Face_mesh_for_singularity_integral_pt(nullptr), Traction_boundary_condition_mesh_pt(nullptr),
  Singular_fct_element_mesh_upper_pt(nullptr), Singular_fct_element_mesh_lower_pt(nullptr),
  Singular_fct_element_mesh_pt(nullptr), Face_mesh_for_bc_pt(nullptr), Torus_region_mesh_pt(nullptr),
  Mesh_as_geom_object_pt(nullptr), Face_mesh_as_geom_object_pt(nullptr),
  Singular_line_mesh_as_geom_object_pt(nullptr), Outer_boundary_pt(nullptr)
{ 
#ifdef OOMPH_HAS_MPI
  std::cout << "This code has been compiled with mpi support \n " 
	    << "and is running on " << this->communicator_pt()->nproc() 
	    << " processors. " << std::endl;
#else
  std::cout << "This code has NOT been compiled mpi support" 
	    << std::endl;
#endif

  // Set the output directory
  Doc_info.set_directory(Global_Parameters::output_directory);

  // set the number of dimensions
  Dim = 3;

  // initialise the FSI iteration counter to zero
  Fsi_iteration_counter = 0;
    
  // OUTER BOUNDARY
  //===============

  // Start boundary IDs for outer boundary from some crazy offset
  // (just for testing). By default the one-based boundary IDs go from
  // 1 to 6; let's start from 1001.
  unsigned outer_boundary_id_offset = 1000;

  //Make the outer boundary object
  Outer_boundary_pt = new CubicTetMeshFacetedSurface(
    Global_Parameters::Box_half_width,
    Global_Parameters::Box_half_height,
    outer_boundary_id_offset);

  // // Look, we can visualise the faceted surface!
  // Outer_boundary_pt->output("outer_faceted_surface.dat");

  // First oomph-lib (zero-based!) boundary ID for outer boundary
  First_boundary_id_for_outer_boundary = outer_boundary_id_offset;
 
  // For sanity check:
  Exact_bounded_volume = 
    2.0*Global_Parameters::Box_half_width*
    2.0*Global_Parameters::Box_half_width*
    2.0*Global_Parameters::Box_half_height;

  // INTERNAL BOUNDARIES
  //====================

  // A warped disk surrounded by a torus
  //------------------------------------

  // Warped disk with specified amplitude and wavenumber for warping
  
  // Thickness of annular region on disk = radius of torus surrounding the
  // edge
  double h_annulus = Global_Parameters::R_torus;
  Global_Parameters::Warped_disk_with_boundary_pt = 
    new WarpedCircularDiskWithAnnularInternalBoundary(h_annulus,
						      Global_Parameters::Epsilon,
						      Global_Parameters::n);

  // set the pointer for coordinate conversions
  CoordinateConversions::disk_geom_obj_pt = Global_Parameters::Warped_disk_with_boundary_pt;
    
  // Enumerate the boundaries making up the disk starting with this
  // one-based ID
  unsigned first_one_based_disk_with_torus_boundary_id = 9001;
    // last_one_based_boundary_for_disk_id+200;
 
  // These get returned
  unsigned last_one_based_disk_with_torus_boundary_id = 0;
  unsigned first_one_based_torus_boundary_id = 0;
  unsigned last_one_based_torus_boundary_id = 0;
 
  // One-based region ID for torus
  unsigned one_based_torus_region_id = 4;

  // Build disk with torus around the edge
  DiskWithTorusAroundEdgeTetMeshFacetedSurface* disk_with_torus_pt = 
    new DiskWithTorusAroundEdgeTetMeshFacetedSurface(
      Global_Parameters::Warped_disk_with_boundary_pt,
      Global_Parameters::Half_nsegment_disk,
      Global_Parameters::R_torus,
      Global_Parameters::Nvertex_torus,
      first_one_based_disk_with_torus_boundary_id,
      one_based_torus_region_id, 
      last_one_based_disk_with_torus_boundary_id,
      first_one_based_torus_boundary_id,
      last_one_based_torus_boundary_id,
      One_based_boundary_id_for_disk_within_torus,
      One_based_boundary_id_for_disk_outside_torus);

  // QUEHACERES this seems to be needed to reset FPU which Triangle
  // messes around with!
  {
    fpu_control_t cw = (_FPU_DEFAULT & ~_FPU_EXTENDED) | _FPU_DOUBLE;
    _FPU_SETCW(cw);
  }
  
  /// Keep track of (zero-based) IDs
  Torus_region_id = one_based_torus_region_id-1; 
  First_lower_disk_boundary_id = 
    first_one_based_disk_with_torus_boundary_id-1;
  Last_lower_disk_boundary_id = 
    last_one_based_disk_with_torus_boundary_id-1;
  First_torus_boundary_id = first_one_based_torus_boundary_id-1;
  Last_torus_boundary_id = last_one_based_torus_boundary_id-1;

  // @@@@@@@@@@@@@@@@@@@@@@@@@@@@
  // QUEHACERES debug

  oomph_info << "\nFirst_lower_disk_boundary_id:                      "
	     << First_lower_disk_boundary_id << "\n"
	     << "Last_lower_disk_boundary_id:                       "
	     << Last_lower_disk_boundary_id << "\n"
	     << "One_based_boundary_id_for_disk_within_torus[0]:    "
	     << One_based_boundary_id_for_disk_within_torus[0] << "\n"
	     << "One_based_boundary_id_for_disk_within_torus[end]:  "
	     << One_based_boundary_id_for_disk_within_torus
    [One_based_boundary_id_for_disk_within_torus.size()-1] << "\n"
	     << "One_based_boundary_id_for_disk_outside_torus[0]:   "
	     << One_based_boundary_id_for_disk_outside_torus[0] << "\n"
	     << "One_based_boundary_id_for_disk_outside_torus[end]: "
	     << One_based_boundary_id_for_disk_outside_torus
    [One_based_boundary_id_for_disk_outside_torus.size()-1] << "\n\n"
	     << "First_torus_boundary_id:                           "
	     << First_torus_boundary_id << "\n"
	     << "Last_torus_boundary_id:                            "
	     << Last_torus_boundary_id << "\n\n";

  // @@@@@@@@@@@@@@@@@@@@@@@@@@@@
  
  // // Look, we can visualise the faceted surface!
  // disk_with_torus_pt->output("warped_disk_with_torus_faceted_surface.dat");
 
  // Add as inner boundary for mesh
  Inner_boundary_pt.push_back(disk_with_torus_pt);

  // Build the mesh
  //--------------- 

  // Not needed, of course, but here to test out the handling
  // of timesteppers
  add_time_stepper_pt(new Steady<1>);

  if(Global_Parameters::Split_corner_elements)
    oomph_info << "\nSplitting corner elements to avoid locking\n\n";
  
  Vector<double> target_element_volume_in_region(1);
  target_element_volume_in_region[0] =
    Global_Parameters::Target_element_volume_in_torus_region;
    
  bool use_attributes = false;
 
  Bulk_mesh_pt =
    new RefineableTetgenMesh<ELEMENT>(Outer_boundary_pt,
				      Inner_boundary_pt,
				      Global_Parameters::Initial_global_element_volume,			      
				      this->time_stepper_pt(),
				      use_attributes,
				      Global_Parameters::Split_corner_elements,
				      &target_element_volume_in_region);
  
  // Problem is linear so we don't need to transfer the solution to the
  // new mesh; we keep it on for self-test purposes...
  Bulk_mesh_pt->disable_projection();

  /// Mesh as geom object representation of mesh
  Mesh_as_geom_object_pt = nullptr;

  // Make new geom object
  Mesh_as_geom_object_pt = new MeshAsGeomObject(Bulk_mesh_pt);

  // hacky - set the global pointer to the GeomObj representation of the mesh
  Global_Parameters::mesh_as_geom_object_pt = Mesh_as_geom_object_pt;
  
  // Number of "disks on disk" around the edge where solution is
  // to be visualised
  Ndisk_on_disk_plot = 4;

  // Number of azimuthal plot points in "disks on disk" plots 
  // around the edge where solution is to be visualised
  // 9 should give us radial lines at 45deg intervals
  Nphi_disk_on_disk_plot = 9;// 30;
 
  // Number of radial plot points in "disks on disk" plots 
  // around the edge where solution is to be visualised
  Nrho_disk_on_disk_plot = 20; // 50;
 
  // Geom object representations need to be (re)built
  Geom_objects_are_out_of_date = true;

  // Set error estimator for bulk mesh
  Z2ErrorEstimator* error_estimator_pt = new Z2ErrorEstimator;
  Bulk_mesh_pt->spatial_error_estimator_pt() = error_estimator_pt;

  // Set targets for spatial adaptivity
  Bulk_mesh_pt->max_permitted_error() = 0.0005; 
  Bulk_mesh_pt->min_permitted_error() = 0.00001;

  // *Before* we duplicate any nodes / generally mess with the mesh, build
  // a map of elements on the torus boundary from inside->outside
  // --------------------------------------------------

  build_map_from_augmented_to_bulk_region_torus_boundary_elements();
  
  // --------------------------------------------------
  
  // populate the vectors which contain pointers to the elements which
  // have nodes on the disk and are identified as being on the upper or
  // lower disk surface. Also populates the face mesh
  identify_elements_on_upper_and_lower_disk_sufaces();

  // Duplicate plate nodes and add upper boundaries
  duplicate_plate_nodes_and_add_boundaries();

  // set the number of singular functions to subtract (for output purposes)
  if(CommandLineArgs::command_line_flag_has_been_set("--subtract_exact_solution"))
    Nsingular_function = 2;
  else
    Nsingular_function = 4;

  // Create face elements that compute contribution to amplitude residual
  //---------------------------------------------------------------------
  Face_mesh_for_singularity_integral_pt = new Mesh;
    
  // Create face elements which handle the jump in stress on the torus boundary
  //-----------------------------------    
  Face_mesh_for_stress_jump_pt = new Mesh;
   
  // Create face elements for imposition of BC
  Face_mesh_for_bc_pt = new Mesh;

  // Traction boundary condition 
  Traction_boundary_condition_mesh_pt = new Mesh;
  
  // Build the face elements
  create_face_elements();

  // initialise
  Singular_fct_element_mesh_upper_pt = new Mesh;
  Singular_fct_element_mesh_lower_pt = new Mesh;

  Singular_fct_element_mesh_pt = new Mesh;
  
  // make the line mesh of elements that sit around the outer edge of the disk
  // and provide the singular functions as functions of the edge coordinates
  // unsigned nsingular_line_element = Global_Parameters::Half_nsegment_disk * 4;
  
  create_one_d_singular_element_mesh(Global_Parameters::Nsingular_line_element);

  oomph_info << "Number of singular amplitude elements: "
	     << Singular_fct_element_mesh_pt->nelement() << "\n";
  
  Singular_line_mesh_as_geom_object_pt =
    new MeshAsGeomObject(Singular_fct_element_mesh_pt);
    
  // Add 'em to mesh
  add_sub_mesh(Bulk_mesh_pt);
  add_sub_mesh(Singular_fct_element_mesh_pt);
  add_sub_mesh(Traction_boundary_condition_mesh_pt);
  add_sub_mesh(Face_mesh_for_bc_pt);
  add_sub_mesh(Face_mesh_for_stress_jump_pt);
  
  build_global_mesh();

  Torus_region_mesh_pt = new RefineableTetgenMesh<ELEMENT>;
  
  // create a Z2 error estimator for this region mesh
  Z2ErrorEstimator* torus_region_error_estimator_pt = new Z2ErrorEstimator;
  Torus_region_mesh_pt->spatial_error_estimator_pt() = torus_region_error_estimator_pt;

  // Complete problem setup
  complete_problem_setup();

    // buckle up
  if(Global_Parameters::Use_fd_lu_solver)
  {
    linear_solver_pt() = new FD_LU;
  }
  else
  {
#ifdef OOMPH_HAS_MUMPS
    oomph_info << "Using MUMPS linear solver\n\n";
  
    MumpsSolver* mumps_linear_solver_pt = new MumpsSolver;

    mumps_linear_solver_pt->enable_suppress_warning_about_MPI_COMM_WORLD();
  
    // set it
    linear_solver_pt() = mumps_linear_solver_pt;
#endif
  }
  
  // moved assign eqns to main body since we may still be pinning tingz
  // after the constructor
}

// ////////////////////////////////////////////////////////////////////////////
// ////////////////////////////////////////////////////////////////////////////
// ////////////////////////////////////////////////////////////////////////////
// ////////////////////////////////////////////////////////////////////////////

//========================================================================
/// \short Function to populate a map which maps each torus boundary
/// element in the augmented region to it's corresponding boundary element
/// in the bulk region
//========================================================================
template <class ELEMENT>
void FlowAroundDiskProblem<ELEMENT>::
build_map_from_augmented_to_bulk_region_torus_boundary_elements()
{
  // generate the node to element look-up
  generate_node_to_element_map();
  
  unsigned region_id = Torus_region_id;
  for (unsigned b = First_torus_boundary_id;
       b <= Last_torus_boundary_id; b++)
  {
    unsigned nel = Bulk_mesh_pt->nboundary_element_in_region(b, region_id);
    for (unsigned e=0; e<nel; e++)
    {
      ELEMENT* augmented_elem_pt = dynamic_cast<ELEMENT*>(
	Bulk_mesh_pt->boundary_element_in_region_pt(b, region_id, e));
        
      // What is the index of the face of the bulk element at the boundary
      int face_index = Bulk_mesh_pt->
	face_index_at_boundary_in_region(b, region_id, e);
      
      // before we build the stress jump element (which duplicates nodes),
      // find the corresponding element which sit on this boundary and shares
      // this face, but in the bulk region (so we can output things like
      // traction on both sides of the region boundary)
      // ---------------------------------------------------------------
      ELEMENT* corresponding_elem_pt = nullptr;
      find_corresponding_element_on_bulk_side_of_augmented_boundary(augmented_elem_pt,
								    face_index,
								    corresponding_elem_pt);

      // get a temperary copy of any existing face index->non-aug element map
      std::map<unsigned, ELEMENT*> face_index_to_non_aug_elem_map =
	Torus_boundary_aug_to_non_aug_by_face_index_map[augmented_elem_pt];

      // add this new one to the copy of the existing map
      face_index_to_non_aug_elem_map[face_index] = corresponding_elem_pt;
      
      // now put the amended map back in 
      Torus_boundary_aug_to_non_aug_by_face_index_map[augmented_elem_pt] =
	face_index_to_non_aug_elem_map;
    }
  }
}

//========================================================================
/// \short Function to populate a map which maps each node in the mesh to
/// a set of all the elements it is associated with
//========================================================================
template <class ELEMENT>
void FlowAroundDiskProblem<ELEMENT>::generate_node_to_element_map()
{
  Node_to_element_map.clear();
  
  // populate the lookup which gives all the elements associated with a given
  // node. We loop over all the elements in the mesh, and for each, we then
  // loop over all the nodes and insert the element pointer into the lookup for
  // that node. 
  unsigned nel_total = Bulk_mesh_pt->nelement();
  for(unsigned i=0; i<nel_total; i++)
  {
    ELEMENT* el_pt = dynamic_cast<ELEMENT*>(Bulk_mesh_pt->element_pt(i));
    
    unsigned nnode = el_pt->nnode();
    for(unsigned j=0; j<nnode; j++)
    {
      Node* node_pt = el_pt->node_pt(j);
      Node_to_element_map[node_pt].insert(el_pt);
    }
  }  
}

//========================================================================
/// \short Function to populate the vectors
/// Elements_on_upper[lower]_disk_surface_pt
/// with the elements which have at least one node on the disk, and also a
/// map which maps each disk node to its associated upper disk elements
/// and the nodal index of the node within each of these elements. 
/// (These will be used when the plate nodes are duplicated, so that we
/// know which elements need to be told that their node pointers have changed).
//========================================================================
template <class ELEMENT>
void FlowAroundDiskProblem<ELEMENT>::identify_elements_on_upper_and_lower_disk_sufaces()
{
  oomph_info << "\nIdentifying upper and lower disk elements...\n\n";

  double t_start = TimingHelpers::timer();

  // get the node to element look-up
  generate_node_to_element_map();

  // ---------------------------------------------------------------------------
  
  // clear vectors (for mesh adaption)
  Elements_on_upper_disk_surface_pt.clear();
  Elements_on_lower_disk_surface_pt.clear();

  // Step 1:
  // populate the vectors of elements which specify the elements which are on
  // the upper and lower surfaces of the disk. The strategy is to loop over
  // the boundary elements, attach face elements, then check the sign of the
  // outer unit normal to see if the element is above or below the disk.ex
  // This will work provided the initial conditions are such that no part of
  // the disk is vertical or beyond (curled over). This also doesn't catch
  // the elements which have nodes or edges on the disk but no faces; these
  // will be found afterwards by brute force
  
  for (unsigned ibound = First_lower_disk_boundary_id;
        ibound <= Last_lower_disk_boundary_id; ibound++)
  { 
    unsigned nel = Bulk_mesh_pt->nboundary_element(ibound);
    for(unsigned e=0; e<nel; e++)
    { 
      ELEMENT* el_pt =
	dynamic_cast<ELEMENT*>(Bulk_mesh_pt->boundary_element_pt(ibound,e));

      // What is the index of the face of the bulk element at the boundary
      int face_index = Bulk_mesh_pt->face_index_at_boundary(ibound,e);

      // Build the corresponding face element
      auto surface_element_pt =
	std::make_unique<RigidDiskFaceElement<ELEMENT>>(el_pt, face_index);

      Vector<double> unit_normal(3, 0.0);
      surface_element_pt->outer_unit_normal(0, unit_normal);

      // if the z component of the surface normal is positive,
      // then the element is sat on the top surface of the disk
      if(unit_normal[2] > 0.0)
      {
	// add to our list of upper surface elements
	Elements_on_upper_disk_surface_pt.insert(el_pt);

	// also add entries for all of its nodes which are on the boundary,
	// with the corresponding index of the node.
	for(unsigned j=0; j<el_pt->nnode(); j++)
	{
	  Node* node_pt = el_pt->node_pt(j);
	  
	  if(node_pt->is_on_boundary(ibound))
	  {
	    // add this boundary node to the map of nodes to upper plate elements
	    // and their associated indices
	    std::pair<ELEMENT*, unsigned> entry(el_pt, j);
	    Disk_node_to_upper_disk_element_and_index_map[node_pt].insert(entry);
	  }
	}
      }
      else 
      {
	// otherwise, it must be below
	Elements_on_lower_disk_surface_pt.insert(el_pt);
      }
    }
  }
  
  // QUEHACERES move this to after so that we output all of them including the touching elems
  // =================================================================
  // QUEHACERES debug - output the elements with faces on the disk
  // =================================================================
  
  // just output the vertices for now
  unsigned nplot = 2;
  char filename[500];
  ofstream some_file;
  
  sprintf(filename, "%s/elements_on_upper_disk_surface.dat", Doc_info.directory().c_str());
  some_file.open(filename);
  
  for(ELEMENT* el_pt : Elements_on_upper_disk_surface_pt)
  {
    el_pt->output(some_file, nplot);
  }

  some_file.close();

  sprintf(filename, "%s/elements_on_lower_disk_surface.dat", Doc_info.directory().c_str());
  some_file.open(filename);
  
  for(ELEMENT* el_pt : Elements_on_lower_disk_surface_pt)
  {
    el_pt->output(some_file, nplot);  
  }

  some_file.close();

  oomph_info << "number of elements with faces on the upper disk surface: "
	     << Elements_on_upper_disk_surface_pt.size() << "\n"
	     << "number of elements with faces on the lower disk surface: "
	     << Elements_on_lower_disk_surface_pt.size() << "\n\n";
  
  // =================================================================
  // now the fiddly bit - identify the elements which have nodes on the disk
  // but no faces. We have to do this by brute force - the algorithm is as
  // follows:
  // - we loop over all the nodes on the disk boundaries, and for each node Nj;
  // - we loop over the vector of upper plate elements to find an element
  //   which contains this node;
  // - we attach a face element to this element on the face which touches the
  //   disk;
  // - we compute the outer surface normal at the node of interest;
  // - we then loop over all the elements in the set of elements associated with
  //   this node, and for each;
  // - we check if we've already identified it as a boundary element.
  //   If not, then we have an element which has a node on the boundary but
  //   no face, and we have an outer unit normal for the disks upper surface
  //   at this nodal location.
  // - We theb check at the distance to the other nodes in this element in the
  //   direction of the outer unit normal.
  //   We'll do a majority vote, so if most of the other nodes have a positive
  //   distance in the normal direction, then this element is above the disk,
  //   otherwise its below the disk. The reason to majority vote is that the
  //   element may have an edge on the disk, i.e. 3 TaylorHood nodes,
  //   and the plate may be curved in the normal direction along this edge,
  //   giving those nodes a positive distance in the normal direction. 
  // =================================================================

  // get the dimensionality of this problem
  unsigned dim = dynamic_cast<FiniteElement*>(Bulk_mesh_pt->element_pt(0))->
    node_pt(0)->ndim();

  // QUEHACERES debug, store any weird elements here so we can output them 
  Vector<ELEMENT*> dodgy_upper_element_pt;
  Vector<ELEMENT*> dodgy_lower_element_pt;

  Vector<Vector<double> > dodgy_upper_element_nodal_distance;
  Vector<Vector<double> > dodgy_lower_element_nodal_distance;
  
  for (unsigned ibound = First_lower_disk_boundary_id;
       ibound <= Last_lower_disk_boundary_id; ibound++)
  { 
    
    unsigned nboundary_node = Bulk_mesh_pt->nboundary_node(ibound);
    for(unsigned n=0; n<nboundary_node; n++)
    {
      // get a pointer to this boundary node
      Node* node_of_interest_pt = Bulk_mesh_pt->boundary_node_pt(ibound, n);
      
      // compute the x-y radius of this point
      double nodal_radius = 0;
      for(unsigned i=0; i<2; i++)
      {
      	nodal_radius += pow(node_of_interest_pt->x(i), 2);
      }
      nodal_radius = sqrt(nodal_radius);

      // don't want to add elements which only touch the edge if we're not
      // duplicating the edge nodes
#ifndef DUPLICATE_EDGE_NODES
      if (abs(nodal_radius - 1) < 1e-6)
	continue;
#endif
      
      // vector to store the outer unit normal for the upper surface at this node
      Vector<double> unit_normal(3, 0.0);
           
      // get the number of elements on this boundary
      unsigned nel = Bulk_mesh_pt->nboundary_element(ibound);

      for(unsigned e=0; e<nel; e++)
      {
	// grab a pointer to this boundary element
	ELEMENT* el_pt =
	  dynamic_cast<ELEMENT*>(Bulk_mesh_pt->boundary_element_pt(ibound, e));

	// now check if this element is in our list of upper surface elements
	typename std::set<ELEMENT*>::iterator it_upper =
	  std::find(Elements_on_upper_disk_surface_pt.begin(),
		    Elements_on_upper_disk_surface_pt.end(), el_pt);

	// and for a double check if this element is in our list of
	// lower surface elements
	typename std::set<ELEMENT*>::iterator it_lower =
	  std::find(Elements_on_lower_disk_surface_pt.begin(),
		    Elements_on_lower_disk_surface_pt.end(), el_pt);
	
	
	if(it_upper != Elements_on_upper_disk_surface_pt.end())
	{
	  // we found it, lets attach a face element to it to get the
	  // outer unit normal
	  
	  // What is the index of the face of the bulk element at the boundary
	  int face_index = Bulk_mesh_pt->face_index_at_boundary(ibound, e);

	  // Build the corresponding face element
	  auto surface_element_pt =
	    std::make_unique<RigidDiskFaceElement<ELEMENT>>(el_pt, face_index);

	  // find the index which corresponds to this node
	  int nodal_index = -1;
	  for(unsigned j=0; j<surface_element_pt->nnode(); j++)
	  {
	    if(surface_element_pt->node_pt(j) == node_of_interest_pt)
	    {
	      nodal_index = j;
	      break;
	    }
	  }

	  // check if we found it
	  if(nodal_index == -1)
	  {	    
	    // if we didn't, lets look in the next element on this boundary
	    continue;
	  }
	  
	  // get the local coordinates of this node
	  Vector<double> s_node(dim-1, 0.0);
	  surface_element_pt->local_coordinate_of_node(nodal_index, s_node); 

	  // get the outer unit normal
	  surface_element_pt->outer_unit_normal(s_node, unit_normal);
	  	  
	  // now we've found the element and got the outer unit normal, we
	  // can stop searching the other elements on this boundary
	  break;
	}
	else if(it_lower != Elements_on_lower_disk_surface_pt.end())
	{
	  // element is on the lower surface
	  // QUEHACERES do something useful - or not?	  
	}
	else
	{
	  std::ostringstream error_message;
	  
	  // weird shit has happened, shout then die
	  error_message << "Weird shit has happened: Found an element which "
			<< "is apparently on a disk boundary but isn't in "
			<< "either of our lists of upper and lower disk elements\n";

	  throw OomphLibError(error_message.str(),
                              OOMPH_CURRENT_FUNCTION,
                              OOMPH_EXCEPTION_LOCATION);
	}
      }
            
      // now we use our map to get a set of all the elements associated with the
      // node of interest, and for each, we check if we've already added it to
      // our list

      // get the set of elements associated with this node
      std::set<ELEMENT*> element_set = Node_to_element_map.at(node_of_interest_pt);
      
      for(ELEMENT* el_pt : element_set)
      {
	// now check if this element is already in our lists of
	// upper and lower elements
	typename std::set<ELEMENT*>::iterator it_upper =
	  std::find(Elements_on_upper_disk_surface_pt.begin(),
		    Elements_on_upper_disk_surface_pt.end(), el_pt);

	typename std::set<ELEMENT*>::iterator it_lower =
	  std::find(Elements_on_lower_disk_surface_pt.begin(),
		    Elements_on_lower_disk_surface_pt.end(), el_pt);

	// if we've already found it then skip it.
	if(it_upper != Elements_on_upper_disk_surface_pt.end() ||
	   it_lower != Elements_on_lower_disk_surface_pt.end() )
	{
	  continue;
	}
	
	unsigned nnode = el_pt->nnode();
	 
	// Now we need to figure out if this element is on the upper
	// or lower disk surface.
	// 
	// We loop over the nodes of this element and get the
	// distance of each from the jth node in the direction of the
	// outer unit normal to the upper surface.
	
	unsigned npositive = 0;
	unsigned nnegative = 0;

	Vector<double> nodal_distance;
	nodal_distance.clear();
	
	for(unsigned k=0; k<nnode; k++)
	{
	  Node* node_k_pt = el_pt->node_pt(k);

	  // if we're at the same node then store the index for our map later,
	  // and skip the processing since the distance will be zero
	  if (node_k_pt == node_of_interest_pt)	  
	    continue;
	  
	  // now do a dot product to get the distance of this node
	  // from the node of interest in the normal direction,
	  // i.e. \delta \bm x\cdot \bm n
	      
	  double normal_distance = 0;	      
	      
	  for(unsigned i=0; i<dim; i++)
	  {
	    double dx = node_k_pt->x(i) - node_of_interest_pt->x(i);
	    normal_distance += dx * unit_normal[i];
	  }

	  nodal_distance.push_back(normal_distance);

	  // add to the counters
	  if(normal_distance >= 0)
	    npositive++;
	  else
	    nnegative++;
	}

	bool element_is_on_upper_surface = npositive > nnegative;
		  
	if(element_is_on_upper_surface)
	{
	  // a non-boundary element can only have at most one edge touching the
	  // disk surface (if it had a face it would be considered a normal
	  // boundary element), so for a Taylor Hood element (NNODE_1D=3) on
	  // the upper surface, this means no more than 2 nodes can have a
	  // positive distance from the boundary node of interest in the outer
	  // normal direction... if this has happened then something weird has
	  // gone wrong, throw an error
	  if(nnegative > 2)
	  {
	    // add this element to our list of dodgy ones so we can output later
	    dodgy_upper_element_pt.push_back(el_pt);

	    // add the vector of nodal distances to our upper vector
	    dodgy_upper_element_nodal_distance.push_back(nodal_distance);
	  }
	  
	  // Let's add this element to our list of upper elements
	  Elements_on_upper_disk_surface_pt.insert(el_pt);
	  
	  // add to our debug vector too
	  Nonboundary_elements_with_node_on_upper_disk_surface_pt.insert(el_pt);
	  
	  // no breaking, because this non-boundary element may have multiple
	  // nodes on the boundary
	}
	else
	{
	  // a non-boundary element can only have at most one edge touching the
	  // disk surface (if it had a face it would be considered a normal
	  // boundary element), so for a Taylor Hood element (NNODE_1D=3) on
	  // the lower surface, this means no more than 2 nodes can have a
	  // positive distance from the boundary node of interest in the outer
	  // normal direction... if this has happened then something weird has
	  // gone wrong, throw an error
	  if(npositive > 2)
	  {	    
	    // add this element to our list of dodgy ones so we can output later
	    dodgy_lower_element_pt.push_back(el_pt);

	    // add the vector of nodal distances to our upper vector
	    dodgy_upper_element_nodal_distance.push_back(nodal_distance);
	  }

	  // Let's add this element to our list of lower elements
	  Elements_on_lower_disk_surface_pt.insert(el_pt);

	  // add to our debug vector too
	  Nonboundary_elements_with_node_on_lower_disk_surface_pt.insert(el_pt);
	  
	  // no breaking, because this non-boundary element may have multiple
	  // nodes on the boundary
	}

	// we keep going here and don't break if we've found a
	// non-"boundary" element containing this boundary node
	// because there may be multiple elements who share this node
	// who are also not "boundary" elements
      }
     
    } // end loop over boundary nodes
    
  } // end loop over boundaries

  // we can now fill in the entry into our map which gives the set of
  // elements and associated nodal indices for each boundary node
  for(ELEMENT* el_pt : Nonboundary_elements_with_node_on_upper_disk_surface_pt)
  {
    for(unsigned j=0; j<el_pt->nnode(); j++)
    {
      Node* node_pt = el_pt->node_pt(j);
      
      for(unsigned b=First_lower_disk_boundary_id; b<=Last_lower_disk_boundary_id; b++)
      {
	if(node_pt->is_on_boundary(b))
	{
	  Disk_node_to_upper_disk_element_and_index_map[node_pt].insert(
	    std::pair<ELEMENT*, unsigned>(el_pt, j));
	}
      }
    }
  }

  // -------------------------------------------------------------
  // output the dodgy stuff
  // -------------------------------------------------------------
  if(!dodgy_upper_element_pt.empty())
  { 
    sprintf(filename, "%s/dodgy_upper_elements.dat", Doc_info.directory().c_str());
    some_file.open(filename);
   
    for(ELEMENT* dodgy_el_pt : dodgy_upper_element_pt)
    {
      dodgy_el_pt->output(some_file, nplot);
    }
    some_file.close();

    sprintf(filename, "%s/dodgy_upper_element_nodal_distances.dat", Doc_info.directory().c_str());
    some_file.open(filename);
   
    for(Vector<Vector<double> >::iterator dodgy_it = dodgy_upper_element_nodal_distance.begin();
	dodgy_it != dodgy_upper_element_nodal_distance.end(); dodgy_it++)
    {
      // get the vector of distances
      Vector<double> distances = *dodgy_it;

      for(Vector<double>::iterator it = distances.begin(); 
	    it != distances.end(); it++)
      {
	some_file << *it << std::endl;
      }

      some_file << std::endl;
    }
    some_file.close();
  }
  
  if(!dodgy_lower_element_pt.empty())
  {
    sprintf(filename, "%s/dodgy_lower_elements.dat", Doc_info.directory().c_str());
    some_file.open(filename);
    
  
    for(ELEMENT* dodgy_el_pt : dodgy_lower_element_pt)
    {
      dodgy_el_pt->output(some_file, nplot);
    }
    some_file.close();

    sprintf(filename, "%s/dodgy_lower_element_nodal_distances.dat", Doc_info.directory().c_str());
    some_file.open(filename);
   
    for(Vector<double> distances : dodgy_lower_element_nodal_distance)
    {
      for(double distance : distances)
      {
	some_file << distance << std::endl;
      }

      some_file << std::endl;
    }
    some_file.close();
  }

  // and print warning message
  if(!dodgy_upper_element_pt.empty() || !dodgy_lower_element_pt.empty())
  {
        oomph_info << "Found some weird elements which are apparently not proper "
	       << "boundary elements \n(i.e. they don't have a whole face on "
	       << "the disk) but which seem to have more \nthan 3 nodes on/very "
	       << "close to its surface. This is probably because the surface \n"
	       << "normal has been computed at one node per element but we are "
	       << "using quadratic \nshape functions and so the normal will "
	       << "rotate within each element for a curved \nsurface.\n"
	       << "  These have been output to: "
	       << Doc_info.directory() << "/dodgy_upper[lower]_elements.dat, and the "
	       << "distances \nof each node from \nthe surface in the surface normal "
	       << "direction has been output to: \n"
	       << Doc_info.directory() << "/dodgy_upper[lower]_element_nodal_distances.dat - "
	       << "you may want to check them!\n\n";
  }
  
  double t_end = TimingHelpers::timer();
  oomph_info << "Identification time: " << t_end - t_start << "s\n";
}


// ////////////////////////////////////////////////////////////////////////////
// ////////////////////////////////////////////////////////////////////////////
// ////////////////////////////////////////////////////////////////////////////
// ////////////////////////////////////////////////////////////////////////////
// ////////////////////////////////////////////////////////////////////////////


//========================================================================
/// \short Function to duplicate the disk nodes to allow for a pressure
/// jump across the plate. New nodes are added to new boundaries for the
/// upper plate surface, and also put back onto any other boundaries that
/// the original nodes were on (except for lower plate boundaries).
//========================================================================
template <class ELEMENT>
void FlowAroundDiskProblem<ELEMENT>::duplicate_plate_nodes_and_add_boundaries()
{
  oomph_info << "\nDuplicating plate nodes and adding upper disk boundaries...\n\n";

  std::set<Node*> lower_disk_nodes_set;
  for(unsigned b=First_lower_disk_boundary_id;
      b<=Last_lower_disk_boundary_id; b++)
  {
    for(unsigned j=0; j<Bulk_mesh_pt->nboundary_node(b); j++)
    {
      Node* node_pt = Bulk_mesh_pt->boundary_node_pt(b,j);

      // just chuck it in, can't have duplicates in a set
      lower_disk_nodes_set.insert(node_pt);
    }
  }

  double t_start = TimingHelpers::timer();
  
  // number of disk boundaries we currently have (before duplication)
  unsigned ndisk_boundary =
    Last_lower_disk_boundary_id - First_lower_disk_boundary_id + 1;
  
  // update the first and last boundary IDs so we can loop over the
  // upper boundaries later.
  // Outer boundaries are enumered as outer boundaries from 1000 onwards,
  // disk boundaries from 9000, and torus boundaries following on from the
  // disk boundaries, so we'll take the last torus boundary and add another
  // offset just to be sure
  unsigned upper_disk_boundary_offset = 1000;
  
  First_upper_disk_boundary_id =
    Last_torus_boundary_id + upper_disk_boundary_offset;

  Last_upper_disk_boundary_id =
    First_upper_disk_boundary_id + ndisk_boundary - 1;

  // increase the number of boundaries in the mesh to accomodate
  // the extra plate boundaries for the upper surface
  Bulk_mesh_pt->set_nboundary(Last_upper_disk_boundary_id + 1); 

  // map to keep track of nodes we've already duplicated;
  // map is original node -> new node
  std::map<Node*,Node*> existing_duplicate_node_pt;

  // counter to keep track of how many boundaries we've done, used to increment
  // the new boundary IDs from the starting ID
  unsigned boundary_counter = 0;
  unsigned new_boundary_id;
  
  for (unsigned b=First_lower_disk_boundary_id; b<=Last_lower_disk_boundary_id; b++)
  {
    new_boundary_id = First_upper_disk_boundary_id + boundary_counter;
    
    unsigned nnode = Bulk_mesh_pt->nboundary_node(b);
    for (unsigned j=0; j<nnode; j++)
    {
      // pointer to the current node (which will stay on the lower disk
      // surface)
      Node* original_node_pt = Bulk_mesh_pt->boundary_node_pt(b,j);

      // compute the x-y radius
      double x = original_node_pt->x(0);
      double y = original_node_pt->x(1);      
      double r = sqrt(x*x + y*y);

      // tolerance on the outer radius
      double tol = 1e-8;
      
      // is the x-y radius 1 (to within the tolerance)? 
      bool node_is_on_edge_of_disk = (abs(1-r) < tol);
      
      // Look for this original node in the map; if we find it
      // it's already been duplicated earlier
      std::map<Node*,Node*>::iterator existing_duplicate_it =
	existing_duplicate_node_pt.find(original_node_pt);

      // this is going to be the new node for the upper disk surface
      // (either a fresh one we will create, or a previously created upper node)
      BoundaryNode<Node>* new_upper_disk_node_pt;

      bool already_duplicated_this_node =
	existing_duplicate_it != existing_duplicate_node_pt.end();
      
      // if we've already duplicated we don't want to create another new node,
      // but we still need to add it to the any new boundaries
      if(!already_duplicated_this_node)
      { 
	// ----------------------------------------------------------------------
	// Step 1: Duplicate the current node and copy over all its attributes
	// ----------------------------------------------------------------------
      
	// get key attributes from the old node
	unsigned n_dim           = original_node_pt->ndim();
	unsigned n_position_type = original_node_pt->nposition_type();
	unsigned n_value         = original_node_pt->nvalue();
	
	// if this node is on the edge we don't want to duplicate it, so just set
	// the pointer to point to the original node. Otherwise, copy over all the
	// info to the new one

#ifndef DUPLICATE_EDGE_NODES
	if(node_is_on_edge_of_disk)
	  new_upper_disk_node_pt = dynamic_cast<BoundaryNode<Node>*>(original_node_pt);
	else
#endif
	{
	  // create a new node
	  new_upper_disk_node_pt =
	    new BoundaryNode<Node>(this->time_stepper_pt(),
				   n_dim, n_position_type, n_value);
      
	  // get the number of time history values each node has
	  unsigned ntstorage = this->time_stepper_pt()->ntstorage();
	
	  // copy over all the nodal values at each time step
	  for(unsigned t=0; t<ntstorage; t++)
	  {
	    // It has the same coordinates...
	    for (unsigned i=0; i<n_dim; i++)
	    {
	      new_upper_disk_node_pt->x(t,i) = original_node_pt->x(t,i);
	    }
      
	    // ...and the same values
	    for (unsigned i=0; i<n_value; i++)
	    {
	      new_upper_disk_node_pt->set_value(t, i, original_node_pt->value(t,i));
	    }	
	  }
	}
      
	// ----------------------------------------------------------------------
	// Step 2: Tell all the elements on the upper surface about the new node
	//         (the old node now corresponds to the lower surface)
	// ----------------------------------------------------------------------
      
	// get the set containing the elements which share this node, and
	// the corresponding node index of the node within the element
	std::set<std::pair<ELEMENT*, unsigned> > upper_disk_element_set =
	  Disk_node_to_upper_disk_element_and_index_map[original_node_pt];

	typename std::set<std::pair<ELEMENT*, unsigned> >::iterator it;

	// now we loop over each of these elements and update the node
	// pointer to point to the newly created node
	for(it = upper_disk_element_set.begin();
	    it != upper_disk_element_set.end(); it++)
	{
	  ELEMENT* el_pt      = it->first;
	  unsigned node_index = it->second;

	  // switch the pointer to the new node
	  el_pt->node_pt(node_index) = new_upper_disk_node_pt;
	}
      }
      else
      {
	// if we already duplicated the current node, then use the existing
	// duplicate rather than creating another 
	new_upper_disk_node_pt = dynamic_cast<BoundaryNode<Node>*>(
	  existing_duplicate_node_pt[original_node_pt]);
      }
      
      // ----------------------------------------------------------------------
      // Step 3: Add the new upper node to the new upper surface boundary,
      //         add it to all the same boundaries as the original (except the
      //         lower disk surface boundaries), and add it to the bulk mesh.
      //         Also want to remove lower surface nodes from upper surface
      //         boundaries.
      // ----------------------------------------------------------------------

      // tell the new node which boundary it's on
      // new_upper_disk_node_pt->add_to_boundary(new_boundary_id);

      // calling this both calls the node pointers function to tell it it's on
      // the boundary, and also updates the mesh's list of boundary nodes      
      Bulk_mesh_pt->add_boundary_node(new_boundary_id,
				      new_upper_disk_node_pt);
      
      // Get/set boundary coordinates
      if ( original_node_pt->boundary_coordinates_have_been_set_up() )
      {
	// get number of coordinates on the original plate boundary
	unsigned ncoords = original_node_pt->ncoordinates_on_boundary(b);
		
	Vector<double> boundary_zeta(ncoords);

	// get 'em from original plate boundary
	original_node_pt->get_coordinates_on_boundary(b, boundary_zeta);
		
	// set 'em for new plate boundary
	new_upper_disk_node_pt->set_coordinates_on_boundary(new_boundary_id,
							    boundary_zeta);
      }
      else
      {
	// hierher throw? (Doesn't happen at the moment, i.e. 
	// when this diagnostic was finally commented out)
             
	oomph_info << "No boundary coordinates have been set up"
		   << " for new local node " << j
		   << " at : "
		   << original_node_pt->x(0) << " "
		   << original_node_pt->x(1)
		   << std::endl;
      }	 

      // get a (pointer to a) set of the boundaries that the original node is on
      std::set<unsigned>* original_node_boundaries_pt;
      original_node_pt->get_boundaries_pt(original_node_boundaries_pt);

      // grab a local copy that we can iterate over
      std::set<unsigned> original_node_boundaries = *original_node_boundaries_pt;
      
      // loop over these and only add the new node to them if it is a boundary
      // for which there is an upper disk element on the boundary which shares
      // this node but not a lower disk boundary.
      for(unsigned ibound : original_node_boundaries)
      {	
	// is the current boundary one of the lower disk boundaries?
	bool is_lower_disk_boundary = (ibound <= Last_lower_disk_boundary_id) &&
	  (ibound >= First_lower_disk_boundary_id);

	// is the current boundary one of the newly created ones? 
	bool is_new_upper_disk_boundary = (ibound <= Last_upper_disk_boundary_id) &&
	  (ibound >= First_upper_disk_boundary_id);
	
	// We don't want to add this new node to the lower boundaries
	// (unless it's on the edge of the disk, in which case it won't
	// be added anyway as the node hasn't been duplicated so it's already
	// on the lower boundaries). 	
	if (is_lower_disk_boundary
#ifndef DUPLICATE_EDGE_NODES
	    && !node_is_on_edge_of_disk
#endif
	  )
	{	  
	  continue;
	}
	else if(is_new_upper_disk_boundary)
	{
	  // if this is a newly created boundary, then we are presumably on an
	  // edge node, so we can just go ahead and add it to this boundary

	  if(!node_is_on_edge_of_disk)
	  {
	    // something weird has happened, the only way a node can be
	    // on both lower and new upper disk boundaries is if it's an
	    // edge node
	    
	    std::ostringstream error_message;
	    error_message << "Something weird has happened, a node seems "
			  << "to be on both a lower (" << b << ") and upper ("
			  << ibound << ") disk boundary "
			  << "but isn't an edge node\n"
			  << "Node coordinates: "
			  << new_upper_disk_node_pt->x(0) << " "
			  << new_upper_disk_node_pt->x(1) << " "
			  << new_upper_disk_node_pt->x(2) << "\n\n";
	    
	    throw OomphLibError(error_message.str(),
				OOMPH_CURRENT_FUNCTION,
				OOMPH_EXCEPTION_LOCATION);
	  }
	  
	  new_upper_disk_node_pt->add_to_boundary(ibound);

	  // get number of coordinates on the original plate boundary
	  unsigned ncoords = original_node_pt->ncoordinates_on_boundary(b);
	  
	  Vector<double> boundary_zeta(ncoords);
	  
	  // get 'em
	  original_node_pt->get_coordinates_on_boundary(b, boundary_zeta);
	  
	  // set 'em
	  new_upper_disk_node_pt->set_coordinates_on_boundary(ibound, boundary_zeta);
	}
	else //Otherwise we want to add it to all the same boundaries
	{	  
	  
	  if ( original_node_pt->boundary_coordinates_have_been_set_up() )
	  {
	    unsigned nboundary_el = Bulk_mesh_pt->nboundary_element(ibound);
	    bool is_only_upper_disk_boundary = true;
	    	    
	    for(unsigned e=0; e<nboundary_el; e++)
	    {
	      ELEMENT* el_pt = dynamic_cast<ELEMENT*>(
		Bulk_mesh_pt->boundary_element_pt(ibound, e));

	      // look for this boundary element in our list of "non-boundary"
	      // elements touching the lower surface of the disk
	      typename std::set<ELEMENT*>::iterator it = std::find(
		Nonboundary_elements_with_node_on_lower_disk_surface_pt.begin(),
		Nonboundary_elements_with_node_on_lower_disk_surface_pt.end(),
		el_pt);

	      // if we find it, then this boundary isn't exclusive to the upper
	      // plate nodes
	      if(it != Nonboundary_elements_with_node_on_lower_disk_surface_pt.end())
	      {
		is_only_upper_disk_boundary = false;
		break;
	      }

	      
	    }
	    // as a double check, also look for this boundary element in our
	    // list of "non-boundary" elements touching the upper surface
	    if(is_only_upper_disk_boundary)
	    {
	      is_only_upper_disk_boundary = false;
	      
	      for(unsigned e=0; e<nboundary_el; e++)
	      {
		ELEMENT* el_pt = dynamic_cast<ELEMENT*>(
		  Bulk_mesh_pt->boundary_element_pt(ibound, e));
	      
		typename std::set<ELEMENT*>::iterator it = std::find(
		  Nonboundary_elements_with_node_on_upper_disk_surface_pt.begin(),
		  Nonboundary_elements_with_node_on_upper_disk_surface_pt.end(),
		  el_pt);

		if(it != Nonboundary_elements_with_node_on_lower_disk_surface_pt.end())
		{
		  is_only_upper_disk_boundary = true;
		  break;
		}
	      }

	      if(!is_only_upper_disk_boundary)
	      {
		std::ostringstream error_message;
		oomph_info << "Something weird has happened, this disk node "
			   << "seems to be on a boundary (" <<ibound << ") which "
			   << "none of its associated upper or lower disk "
			   << "elements are on.\n"
			   << "Nodal coordinates: "
			   << new_upper_disk_node_pt->x(0) << " "
			   << new_upper_disk_node_pt->x(1) << " "
			   << new_upper_disk_node_pt->x(2) << "\n\n";
		
		throw OomphLibError(error_message.str(),
				    OOMPH_CURRENT_FUNCTION,
				    OOMPH_EXCEPTION_LOCATION);
	      }
	    }
	    
	    if(is_only_upper_disk_boundary)
	    {
	      new_upper_disk_node_pt->add_to_boundary(ibound);

	      // get number of coordinates on the original plate boundary
	      unsigned ncoords = original_node_pt->ncoordinates_on_boundary(b);
	  
	      Vector<double> boundary_zeta(ncoords);
	  
	      // get 'em
	      original_node_pt->get_coordinates_on_boundary(b, boundary_zeta);
	  
	      // set 'em
	      new_upper_disk_node_pt->set_coordinates_on_boundary(ibound, boundary_zeta);

	      // if this is only an upper disk boundary, remove the original
	      // node from it, unless it's an edge node
	      if(!node_is_on_edge_of_disk)
	      	original_node_pt->remove_from_boundary(ibound);	      
	    }
	  }
	  else
	  {
	    // hierher throw? (Doesn't happen at the moment, i.e. 
	    // when this diagnostic was finally commented out)
             
	    oomph_info << "No boundary coordinates have been set up"
		       << " for new local node " << j
		       << " at : "
		       << original_node_pt->x(0) << " "
		       << original_node_pt->x(1)
		       << std::endl;
	  }	
	}
      }

      // add it to our list of duplicates
      existing_duplicate_node_pt[original_node_pt] = new_upper_disk_node_pt;

      // and add the new node to the bulk mesh if we haven't already
      if(!already_duplicated_this_node
#ifndef DUPLICATE_EDGE_NODES
	 && !node_is_on_edge_of_disk
#endif
	)
      {	
	Bulk_mesh_pt->add_node_pt(new_upper_disk_node_pt);
      }
      
    } // end loop over boundary nodes

    boundary_counter++;
    
  } // end loop over boundaries

  std::set<Node*> upper_disk_nodes_set;
  for(unsigned b=First_upper_disk_boundary_id;
      b<=Last_upper_disk_boundary_id; b++)
  {
    for(unsigned j=0; j<Bulk_mesh_pt->nboundary_node(b); j++)
    {
      Node* node_pt = Bulk_mesh_pt->boundary_node_pt(b,j);

      // just chuck it in, can't have duplicates in a set
      upper_disk_nodes_set.insert(node_pt);
    }
  }
  
  unsigned offset = First_upper_disk_boundary_id - First_lower_disk_boundary_id;
 
  // update the list of boundary IDs in/outside the torus region
  unsigned n = One_based_boundary_id_for_disk_within_torus.size();
  Boundary_id_for_upper_disk_within_torus.resize(n);
    
  for(unsigned i=0; i<n; i++)
  {
    // get the new (zero-based) boundary ID
    unsigned new_id = One_based_boundary_id_for_disk_within_torus[i] + offset - 1;

    // and add it to the list
    Boundary_id_for_upper_disk_within_torus[i] = new_id;
  }
  
  n = One_based_boundary_id_for_disk_outside_torus.size();
  Boundary_id_for_upper_disk_outside_torus.resize(n);
  
  for(unsigned i=0; i<n; i++)
  {
    // get the new (zero-based) boundary ID
    unsigned new_id = One_based_boundary_id_for_disk_outside_torus[i] + offset - 1;

    // and add it to the list (with the 1 added)
    Boundary_id_for_upper_disk_outside_torus[i] = new_id;
  }
  
  unsigned nlower_disk_nodes = lower_disk_nodes_set.size();
  unsigned nupper_disk_nodes = upper_disk_nodes_set.size();
  
  oomph_info << "Number of plate nodes before duplication: "
	     << nlower_disk_nodes << "\n";
  oomph_info << "Number of plate nodes after duplication: "
	     << nlower_disk_nodes + nupper_disk_nodes << "\n\n";

  // and finally, update this since we've fiddled the nodes on the plate
  // boundaries. N.B. this doesn't update element-in-region info, so new
  // boundaries aren't "in" the torus region
  Bulk_mesh_pt->setup_boundary_element_info();

  // combine the upper and lower boundary IDs because nodes have been duplicated
  // so need to attach face elements onto the upper and lower surfaces of the disk
  for(unsigned i=0; i<One_based_boundary_id_for_disk_within_torus.size(); i++)
    Disk_boundary_ids_in_torus.push_back(One_based_boundary_id_for_disk_within_torus[i] - 1);

  for(unsigned i=0; i<Boundary_id_for_upper_disk_within_torus.size(); i++)
    Disk_boundary_ids_in_torus.push_back(Boundary_id_for_upper_disk_within_torus[i]);

  oomph_info << "\nFirst_upper_disk_boundary_id: " << First_upper_disk_boundary_id 
	     << "\nLast_upper_disk_boundary_id:  " << Last_upper_disk_boundary_id << "\n\n";
  
  double t_end = TimingHelpers::timer();
  oomph_info << "Time to duplicate nodes and add boundaries: " << t_end - t_start << "s\n\n";

  // regenerate the node to element look-up
  generate_node_to_element_map();

  // ----------------------
  // QUEHACERES check that we don't have any nodes which aren't attached to elements
#ifdef PARANOID
  for(unsigned j=0; j<Bulk_mesh_pt->nnode(); j++)
  {
    Node* node_pt = Bulk_mesh_pt->node_pt(j);

    if(Node_to_element_map.find(node_pt) == Node_to_element_map.end())
    {
      oomph_info << "\n===========\nWARNING: Node " << j << "("
		 << node_pt->x(0) << ", "
		 << node_pt->x(1) << ", "
		 << node_pt->x(2) << ") has no associated elements\n\n";
    }      
  }
#endif

#ifdef SHIFT_LOWER_PLATE_NODES_BY_DZ
  
  // QUEHACERES shift the lower plate nodes down a touch so we can idendify them  
  for(std::map<Node*,Node*>::iterator it=existing_duplicate_node_pt.begin();
      it != existing_duplicate_node_pt.end(); it++)
  {
    it->first->x(2) += lower_plate_z_shift;
  }
#endif
  
  bool first_boundary_without_nodes = true;
  
  unsigned new_boundaries_with_some_nodes = 0;
  unsigned new_boundaries_with_no_nodes = 0;
  
  
  for(unsigned b=First_upper_disk_boundary_id; b<=Last_upper_disk_boundary_id; b++)
  {
    unsigned nnode = Bulk_mesh_pt->nboundary_node(b);

    if (nnode > 0)
      new_boundaries_with_some_nodes++;
    else
    {
      new_boundaries_with_no_nodes++;
      if(first_boundary_without_nodes)
      {
	first_boundary_without_nodes = false;
      }
    }
  }

#ifdef PARANOID
  // For extra comfort,  check the sync between
  // Disk_node_to_upper_disk_element_and_index_map and
  // Elements_on_upper_disk_surface_pt
  
  std::set<ELEMENT*> unique_elements_from_map;

  typename std::map<Node*, std::set<std::pair<ELEMENT*, unsigned> > >::iterator it;
  for(it = Disk_node_to_upper_disk_element_and_index_map.begin();
      it != Disk_node_to_upper_disk_element_and_index_map.end(); it++)
  {
    // get the set
    std::set<std::pair<ELEMENT*, unsigned> > upper_disk_element_set = it->second;

    // iterate over the second and output the nodes
    for(typename std::set<std::pair<ELEMENT*, unsigned> >::iterator set_it =
	  upper_disk_element_set.begin(); set_it != upper_disk_element_set.end(); set_it++)
    {
      ELEMENT* el_pt = set_it->first;

      unique_elements_from_map.insert(el_pt);
    }
  }
  
  if(unique_elements_from_map != Elements_on_upper_disk_surface_pt)
  {
    std::ostringstream error_message;
    error_message << "Error: there is a discrepancy between the set of \n"
		  << "upper disk elements and the map from disk nodes \n"
		  << "to upper disk elements\n\n";
    
    throw OomphLibError(error_message.str(),
			OOMPH_CURRENT_FUNCTION,
			OOMPH_EXCEPTION_LOCATION); 
  }
#endif
  
}

//========================================================================
/// \short function to compute the edge coordinates and corresponding
/// element and local coordinates in the singular line mesh
/// (which computes the appropriate amplitude) from Eulerian coordinates
//========================================================================
template <class ELEMENT>
void FlowAroundDiskProblem<ELEMENT>::get_edge_coordinates_and_singular_element(
  FiniteElement* const& elem_pt,
  const Vector<double>& x,
  EdgeCoordinates& edge_coordinates_at_point,
  std::pair<GeomObject*, Vector<double> >& line_element_and_local_coordinate) const
{
  // tolerance on point away from the x-axis
  const double tol = 1e-8;
  
  // ------------------------------------------
  // Step 1: Compute the (\rho,\zeta,\phi) coordinates
  // ------------------------------------------

  CoordinateConversions::eulerian_to_lagrangian_coordinates(x, edge_coordinates_at_point);
  
  // if this point as an angle of pi but is on the lower side of the disk
  // rather than the upper, set it's angle to -pi to get the pressure jump right.
  // N.B. this will only matter for plot points where the output is done at the
  // nodes; integration points are always within the element
  if(abs(edge_coordinates_at_point.phi - MathematicalConstants::Pi) < tol)
  {
    // try and cast this element just to make sure we don't have a face element
    const ELEMENT* bulk_el_pt = dynamic_cast<const ELEMENT*>(elem_pt);

    if(bulk_el_pt != 0)
    {
      // now search for it in the list of lower disk elements
      if(Elements_on_lower_disk_surface_pt.find(const_cast<ELEMENT*>(bulk_el_pt)) !=
	 Elements_on_lower_disk_surface_pt.end())
      {	  
	edge_coordinates_at_point.phi = -MathematicalConstants::Pi;
      }
    }
  }
  // equally if the angle is -pi but is on the top surface, change it to +pi
  else if(abs(edge_coordinates_at_point.phi + MathematicalConstants::Pi) < tol)
  {
    // try and cast this element just to make sure we don't have a face element
    const ELEMENT* bulk_el_pt = dynamic_cast<const ELEMENT*>(elem_pt);

    if(bulk_el_pt != 0)
    {
      // now search for it in the list of lower disk elements
      if(Elements_on_upper_disk_surface_pt.find(const_cast<ELEMENT*>(bulk_el_pt)) !=
	 Elements_on_upper_disk_surface_pt.end())
      {	  
	edge_coordinates_at_point.phi = MathematicalConstants::Pi;
      }
    }
  }
        
  // ----------------------------------------------------------------------
  // Step 2: Find the line element and corresponding local coordinates
  //         which correspond to this value of zeta, and then save these
  //         values for this point
  // ----------------------------------------------------------------------

  Vector<double> s_line(1, 0.0);
  Vector<double> zeta_vec(1, 0.0);
  zeta_vec[0] = edge_coordinates_at_point.zeta;

  GeomObject* geom_obj_pt = nullptr;

  Singular_line_mesh_as_geom_object_pt->locate_zeta(zeta_vec,
						    geom_obj_pt,
						    s_line);

  // if we still didn't find it then shout and die, because we seem to have a zeta
  // which doesn't exist in our line meshes
  if(geom_obj_pt == 0)
  {
    ostringstream error_message;

    error_message << "Zeta: " << edge_coordinates_at_point.zeta
		  << " not found in the singular line meshes";
    throw OomphLibError(error_message.str(),
			OOMPH_CURRENT_FUNCTION,
			OOMPH_EXCEPTION_LOCATION);
  }
    
  // combine the geometric object representation of the line element pointer
  // and its local coordinate which represent this zeta
  line_element_and_local_coordinate = std::make_pair(geom_obj_pt, s_line);
}

//========================================================================
/// \short function to setup the edge coordinates (\rho, \zeta, \phi) for
/// elements in the torus region. If use_plot_points is true, then the
/// coordinates are computed for the elements plot points, otherwise it is
/// for their knots (integration points).
//========================================================================
template <class ELEMENT>
void FlowAroundDiskProblem<ELEMENT>::setup_edge_coordinates_and_singular_element(
  FiniteElement* const& elem_pt,
  Vector<EdgeCoordinates>& edge_coordinates_at_point,
  Vector<std::pair<GeomObject*, Vector<double> > >& line_element_and_local_coordinate,
  const bool& use_plot_points) const
{  
  unsigned npt = 0;

  if(use_plot_points)
    npt = elem_pt->nplot_points(Global_Parameters::Nplot_for_bulk);
  else
    npt = elem_pt->integral_pt()->nweight();

  // make enough space for the coordinates and element pairs
  edge_coordinates_at_point.resize(npt);
  line_element_and_local_coordinate.resize(npt);
  
  //Loop over the integration/plot points and compute 
  for(unsigned ipt=0; ipt<npt; ipt++)
  {
    // local coords in this bulk element
    Vector<double> s(Dim, 0.0);

    if(use_plot_points)
    {
      // Get local coordinates of plot point
      elem_pt->get_s_plot(ipt, Global_Parameters::Nplot_for_bulk, s);
    }
    else
    {
      //Assign values of s
      for(unsigned i=0; i<Dim; i++) 
      {
	s[i] = elem_pt->integral_pt()->knot(ipt, i);
      }
    }
    
    // get the interpolated Eulerian coordinates for this knot point
    Vector<double> x(Dim, 0.0);
    
    for(unsigned i=0; i<Dim; i++) 
    {
      x[i] = elem_pt->interpolated_x(s,i);
    }

    // now get the edge coordinates and corresponding singular line mesh element
    // for these Eulerian coordinates
    EdgeCoordinates edge_coords;
    std::pair<GeomObject*, Vector<double> > line_elem_and_coord;
    
    get_edge_coordinates_and_singular_element(elem_pt, x,
					      edge_coords,
					      line_elem_and_coord);

    // and add them to the list
    line_element_and_local_coordinate[ipt] = line_elem_and_coord;
    edge_coordinates_at_point[ipt] = edge_coords;
  }    
}


/// \Short Helper function which creates the line mesh of 1D elements which sit
/// on the outer edge of the disk and provide the amplitude of the sinuglar
/// function as a function of the outer boundary coordinate, i.e. c(\zeta).
template <class ELEMENT>
void FlowAroundDiskProblem<ELEMENT>::
create_one_d_singular_element_mesh(const unsigned& nsingular_el)
{
  // build the bastard! QUEHACERES sort out NNODE_1D 
  Singular_fct_element_mesh_pt =
    new CircularLineMesh<ScalableSingularityForNavierStokesLineElement<SINGULAR_ELEMENT_NNODE_1D> >(
      nsingular_el);

  // now set the nodal zeta coordinates for all but the last element
  for(unsigned e=0; e < (nsingular_el - 1); e++)
  {
    dynamic_cast<ScalableSingularityForNavierStokesLineElement<SINGULAR_ELEMENT_NNODE_1D>*>(
      Singular_fct_element_mesh_pt->element_pt(e))->setup_zeta_nodal();
  }

  // and now do the last element, passing the flag which sets the zeta of the
  // last node in the last element equal to 2pi not 0 so that we don't get
  // interpolation errors in the last element
  bool use_zeta_2pi_instead_of_0 = true;
  dynamic_cast<ScalableSingularityForNavierStokesLineElement<SINGULAR_ELEMENT_NNODE_1D>*>(
    Singular_fct_element_mesh_pt->element_pt(nsingular_el - 1))->
    setup_zeta_nodal(use_zeta_2pi_instead_of_0);
 
#ifdef PARANOID

  // sanity check - is the number of nodes in the two meshes the same as the
  // number of vertices requested for the meshing of the disk? (minus two,
  // because the two halves of the mesh overlap at 2 nodes
  oomph_info << "Number of nodes in singular line mesh:    "
	     << Singular_fct_element_mesh_pt->nnode() << "\n"
	     << "Number of nodes on disk edge:             "
	     << Global_Parameters::Half_nsegment_disk * 4 << "\n" << std::endl;
#endif
}



// ///////////////////////////////////////////////////////////////////////
// ///////////////////////////////////////////////////////////////////////
// ///////////////////////////////////////////////////////////////////////

// QUEHACERES delete
//========================================================================
/// Setup disk on disk plots
//========================================================================
template <class ELEMENT>
void FlowAroundDiskProblem<ELEMENT>::setup_disk_on_disk_plots()
{
  if(!Global_Parameters::Do_disk_on_disk_plots)
  {
    oomph_info << "Not making geom object" << std::endl;
    return;
  }
  
  oomph_info << "Starting make geom object" << std::endl;
  double t_start=TimingHelpers::timer();
   
  // Make space for plot points: Disk_on_disk_plot_point[k_disk][i_rho][i_phi]
  Disk_on_disk_plot_point.resize(Ndisk_on_disk_plot);
  for (unsigned i=0; i < Ndisk_on_disk_plot; i++)
  {
    Disk_on_disk_plot_point[i].resize(Nrho_disk_on_disk_plot);
    for (unsigned j=0; j < Nrho_disk_on_disk_plot; j++)
    {
      Disk_on_disk_plot_point[i][j].resize(Nphi_disk_on_disk_plot);
    }
  }

  Vector<double> r_edge(3);
  Vector<double> normal(3);  
  Vector<double> tangent(3); 
  Vector<double> normal_normal(3);   
  Vector<double> x(3);     
  Vector<double> s(3);  
  Vector<double> rho_and_phi(2);
  GeomObject* geom_object_pt = nullptr;
  
  for (unsigned k=0; k < Ndisk_on_disk_plot; k++)
  {
    double theta = double(k) / double(Ndisk_on_disk_plot) *
      2.0*MathematicalConstants::Pi;
    
    Global_Parameters::Warped_disk_with_boundary_pt->
      boundary_triad(0,theta, r_edge, tangent, normal, normal_normal);
    
    for (unsigned i=0; i < Nrho_disk_on_disk_plot; i++) 
    {
      double rho_min = 0.0;
      double rho_max = Global_Parameters::R_torus * 1.2;
      double rho = rho_min + (rho_max - rho_min) * double(i) /
	double(Nrho_disk_on_disk_plot-1);
      
      rho_and_phi[0] = rho;
      
      for (unsigned j=0; j < Nphi_disk_on_disk_plot; j++)
      {
        double phi = double(j) / double(Nphi_disk_on_disk_plot-1) *
	  2.0*MathematicalConstants::Pi;
        rho_and_phi[1] = phi;
	
        x[0] = r_edge[0] + rho*cos(phi)*normal[0] + rho*sin(phi)*normal_normal[0];
        x[1] = r_edge[1] + rho*cos(phi)*normal[1] + rho*sin(phi)*normal_normal[1];
        x[2] = r_edge[2] + rho*cos(phi)*normal[2] + rho*sin(phi)*normal_normal[2];

        Mesh_as_geom_object_pt->locate_zeta(x, geom_object_pt, s);
	
        if (geom_object_pt == 0)
	{
          oomph_info << "Point : " 
                     << x[0] << " " 
                     << x[1] << " " 
                     << x[2] << " "
                     << " not found in setup of disk on disk plots" 
                     << std::endl;
	}        
        Disk_on_disk_plot_point[k][i][j] =
	  std::make_pair(rho_and_phi, std::make_pair(geom_object_pt,s));
      }
    }
  }

  oomph_info << "Completed setup of disk on disk plots. This took " 
             << TimingHelpers::timer()-t_start << " sec"
             << std::endl;
}

//========================================================================
/// Complete problem setup
//========================================================================
template <class ELEMENT>
void FlowAroundDiskProblem<ELEMENT>::complete_problem_setup()
{
  // map which keeps track of the first index of each node associated with
  // the PDE-constraint Lagrange multipliers
  std::map<Node*, unsigned> node_to_first_lm_index_map;
  
  // loop over all the elements and tell them how many singular functions
  // there are (for output purposes only)
  for(unsigned e=0; e<Bulk_mesh_pt->nelement(); e++)
  {
    ELEMENT* elem_pt = dynamic_cast<ELEMENT*>(Bulk_mesh_pt->element_pt(e));

    // for output purposes
    elem_pt->set_nsingular_fct(Nsingular_function);

    // add the Lagrange multipliers which weakly enforce the momentum and
    // continuity equations    
    elem_pt->add_lagrange_multiplier_dofs(node_to_first_lm_index_map);

    // constitutive relation
    elem_pt->stress_fct_pt() = &Analytic_Functions::stress;

    // QUEHACERES removing, effectively absorbing the hydrostatic pressure
    // into the pressure gradient term
    // // body force acting on fluid due to reduced gravity
    // elem_pt->body_force_fct_pt() = &Global_Parameters::reduced_gravity;
  }
  
  // Loop over the elements to set up element-specific
  // things that cannot be handled by constructor

  oomph_info << "Computing Lagrangian coordinates and corresponding singular "
	     << "line elements for integration/plot points in augmented region..."
	     << std::endl;
    
  // Bulk elements in torus region
  unsigned region_id = Torus_region_id;
  unsigned n_el = Bulk_mesh_pt->nregion_element(region_id);
  for (unsigned e=0; e<n_el; e++)
  {
    ELEMENT* torus_region_el_pt = dynamic_cast<ELEMENT*>(
      Bulk_mesh_pt->region_element_pt(region_id, e));

    // pass in the function which computes the derivative of the functional
    // which we're minimising in the augmented region
    torus_region_el_pt->dfunctional_du_fct_pt() =
      &Analytic_Functions::dfunctional_du;

    // set the pressure regularisation factor
    torus_region_el_pt->set_pressure_regularisation_factor(
      Global_Parameters::Pressure_regularisation_factor);

    // QUEHACERES 
    torus_region_el_pt->set_velocity_regularisation_factor(
      Global_Parameters::Velocity_regularisation_factor);
    
    // tell the element it's augmented, so that it uses the PDE-constrained min
    // residuals rather than standard Stokes
    torus_region_el_pt->set_augmented_element();
    
    // the edge coordinates for each of this elements plot points
    Vector<EdgeCoordinates> edge_coordinates_at_plot_point;

    // the edge coordinates for each of this elements knot points
    Vector<EdgeCoordinates> edge_coordinates_at_knot;
    
    // the line element and its local coordinate which correspond to each of
    // this elements plot/knot points
    Vector<std::pair<GeomObject*, Vector<double> > >
      line_element_and_local_coordinates;

    // we're doing bulk elements here, so we only need to know about the
    // the singularity for plotting purposes not the actual solve, so want
    // the rzp coordinates at plot points rather than integration points
    bool use_plot_points = true;

    // -------------------------------
    // first we'll do the plot points
    
    // get the coordinates
    setup_edge_coordinates_and_singular_element(torus_region_el_pt,
						edge_coordinates_at_plot_point,
						line_element_and_local_coordinates,
						use_plot_points);
    
    // now tell the element about them
    torus_region_el_pt->set_edge_coordinates_at_plot_point(edge_coordinates_at_plot_point);	

    torus_region_el_pt->set_line_element_and_local_coordinate_at_plot_point(
    line_element_and_local_coordinates);

    // -------------------------------
    // now the knot points

    use_plot_points = false;

    // get the coordinates
    setup_edge_coordinates_and_singular_element(torus_region_el_pt,
      edge_coordinates_at_knot,
      line_element_and_local_coordinates,
      use_plot_points);
    
    // now tell the element about them
    torus_region_el_pt->set_edge_coordinates_at_knot(edge_coordinates_at_knot);
	
    torus_region_el_pt->set_line_element_and_local_coordinate_at_knot(
    line_element_and_local_coordinates);
  }
  
  // ==========================================================================
  
  if (!CommandLineArgs::command_line_flag_has_been_set("--dont_subtract_singularity"))
  {
    // Stress jump elements
    unsigned n_element = Face_mesh_for_stress_jump_pt->nelement();
    for(unsigned e=0; e<n_element; e++)
    {
      // Upcast from GeneralisedElement to the present element
      NavierStokesWithSingularityStressJumpFaceElement<ELEMENT>* stress_jump_el_pt =
	dynamic_cast<NavierStokesWithSingularityStressJumpFaceElement<ELEMENT>*>(
	  Face_mesh_for_stress_jump_pt->element_pt(e));

      // the edge coordinates for each of this elements knot points
      Vector<EdgeCoordinates> edge_coordinates_at_knot_point;

      // the line element and its local coordinate which correspond to each of
      // this elements knot points
      Vector<std::pair<GeomObject*, Vector<double> > >
	line_element_and_local_coordinate;

      // the singularity contributes to the residuals of this face element, so
      // we need to compute the (\rho,\zeta,\phi) coordinates at the knot points
      // instead of at the plot points
      bool use_plot_points = false;

      // get the coordinates
      setup_edge_coordinates_and_singular_element(stress_jump_el_pt,
						  edge_coordinates_at_knot_point,
						  line_element_and_local_coordinate,
						  use_plot_points);
      // now tell the element about them
      stress_jump_el_pt->set_edge_coordinates_at_knot(edge_coordinates_at_knot_point);
      
      stress_jump_el_pt->set_line_element_and_local_coordinate_at_knot(
	line_element_and_local_coordinate);

      // QUEHACERES for debug
      stress_jump_el_pt->exact_traction_fct_pt() = &Analytic_Functions::prescribed_traction;
    }
  }
  
  // BC elements
  unsigned n_element =  Face_mesh_for_bc_pt->nelement();
  for(unsigned e=0; e<n_element; e++)
  {
    // Upcast from GeneralisedElement to the present element
    NavierStokesWithSingularityBCFaceElement<ELEMENT>* bc_el_pt =
      dynamic_cast<NavierStokesWithSingularityBCFaceElement<ELEMENT>*>(
    Face_mesh_for_bc_pt->element_pt(e));

    // the edge coordinates for each of this elements knot points
    Vector<EdgeCoordinates> edge_coordinates_at_knot_point;

    // the line element and its local coordinate which correspond to each of
    // this elements knot points
    Vector<std::pair<GeomObject*, Vector<double> > >
      line_element_and_local_coordinate;

    // the singularity contributes to the residuals of this face element, so
    // we need to compute the (\rho,\zeta,\phi) coordinates at the knot points
    // instead of at the plot points
    bool use_plot_points = false;

    // get the coordinates
    setup_edge_coordinates_and_singular_element(bc_el_pt,
						edge_coordinates_at_knot_point,
						line_element_and_local_coordinate,
						use_plot_points);
    
    // now tell the element about them
    bc_el_pt->set_edge_coordinates_at_knot(edge_coordinates_at_knot_point);
      
    bc_el_pt->set_line_element_and_local_coordinate_at_knot(
    line_element_and_local_coordinate);

    // set the external data pointers since this face element contributes to
    // the residuals for the singular amplitudes
    bc_el_pt->set_singular_amplitudes_as_external_data();
  }

  // ------------------------------------------------
  // and now loop over the elements in the singular line meshes and tell them
  // about the singular functions
    
  for(unsigned e=0; e<Singular_fct_element_mesh_pt->nelement(); e++)
  {
    ScalableSingularityForNavierStokesLineElement<SINGULAR_ELEMENT_NNODE_1D>* sing_el_pt =
      dynamic_cast<ScalableSingularityForNavierStokesLineElement<SINGULAR_ELEMENT_NNODE_1D>*>(
	Singular_fct_element_mesh_pt->element_pt(e));

    // QUEHACERES exact solution for debug
    if(CommandLineArgs::command_line_flag_has_been_set("--subtract_exact_solution"))
    {
      // Broadside singular function ----------------------------------
      sing_el_pt->add_unscaled_singular_fct_and_gradient_pt(
	&Global_Parameters::SingularFunctions::wrapper_to_exact_solution_broadside,
	&Global_Parameters::SingularFunctions::wrapper_to_exact_velocity_gradient_broadside,
	Sing_fct_id_broadside);

      // In-plane singular function ----------------------------------
      sing_el_pt->add_unscaled_singular_fct_and_gradient_pt(
	&Global_Parameters::SingularFunctions::wrapper_to_exact_solution_in_plane,
	&Global_Parameters::SingularFunctions::wrapper_to_exact_velocity_gradient_in_plane,
	Sing_fct_id_in_plane);

      // QUEHACERES maybe add in-plane rotation here too
    }
    else
    {
      // Broadside singular function ------------------------------------------
      
      // QUEHACERES sticking in all 4 full solutions 11/08

      // QUEHACERES experimental - use full soln not asymptotic one! 19/01 @@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
      sing_el_pt->add_unscaled_singular_fct_and_gradient_pt(
      	&Global_Parameters::SingularFunctions::__singular_fct_exact_broadside_translation,
      	&Global_Parameters::SingularFunctions::__gradient_of_singular_fct_exact_broadside_translation,
      	Sing_fct_id_broadside);

      sing_el_pt->add_unscaled_singular_fct_and_gradient_pt(
      	&Global_Parameters::SingularFunctions::__singular_fct_exact_broadside_rotation,
      	&Global_Parameters::SingularFunctions::__gradient_of_singular_fct_exact_broadside_rotation,
      	Sing_fct_id_broadside_rotation);

      // In-plane singular function -------------------------------------------
      sing_el_pt->add_unscaled_singular_fct_and_gradient_pt(
      	&Global_Parameters::SingularFunctions::__singular_fct_exact_in_plane_translation,
      	&Global_Parameters::SingularFunctions::__gradient_of_singular_fct_exact_in_plane_translation,
      	Sing_fct_id_in_plane);

      sing_el_pt->add_unscaled_singular_fct_and_gradient_pt(
      	&Global_Parameters::SingularFunctions::__singular_fct_exact_in_plane_rotation,
      	&Global_Parameters::SingularFunctions::__gradient_of_singular_fct_exact_in_plane_rotation,
      	Sing_fct_id_in_plane_rotation);
      // @@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
      
      // sing_el_pt->add_unscaled_singular_fct_and_gradient_pt(
      // 	&Global_Parameters::SingularFunctions::singular_fct_exact_asymptotic_broadside,
      // 	&Global_Parameters::SingularFunctions::gradient_of_singular_fct_exact_asymptotic_broadside,
      // 	Sing_fct_id_broadside);

      // sing_el_pt->add_unscaled_singular_fct_and_gradient_pt(
      // 	&Global_Parameters::SingularFunctions::singular_fct_exact_asymptotic_broadside_rotation,
      // 	&Global_Parameters::SingularFunctions::
      // 	gradient_of_singular_fct_exact_asymptotic_broadside_rotation,
      // 	Sing_fct_id_broadside_rotation);

      // // sing_el_pt->add_unscaled_singular_fct_and_gradient_pt(
      // // 	&Global_Parameters::SingularFunctions::singular_fct_exact_asymptotic_broadside_rotation_az,
      // // 	&Global_Parameters::SingularFunctions::
      // // 	gradient_of_singular_fct_exact_asymptotic_broadside_rotation_az,
      // // 	Sing_fct_id_broadside_az);
	
      // // In-plane singular function -------------------------------------------

      // sing_el_pt->add_unscaled_singular_fct_and_gradient_pt(
      // 	&Global_Parameters::SingularFunctions::singular_fct_exact_asymptotic_in_plane_full,
      // 	&Global_Parameters::SingularFunctions::gradient_of_singular_fct_exact_asymptotic_in_plane_full,
      // 	Sing_fct_id_in_plane);
	    
      // // // QUEHACERES using the asymptotically expanded exact solution for now 12/08
      // // sing_el_pt->add_unscaled_singular_fct_and_gradient_pt(
      // // 	&Global_Parameters::SingularFunctions::singular_fct_exact_asymptotic_in_plane,
      // // 	&Global_Parameters::SingularFunctions::gradient_of_singular_fct_exact_asymptotic_in_plane,
      // // 	Sing_fct_id_in_plane);

      // // // u_zeta velocity component for in-plane singular function ----------------
      // // // QUEHACERES if this works, change the enum name for the ID
      // // sing_el_pt->add_unscaled_singular_fct_and_gradient_pt(
      // // 	&Global_Parameters::SingularFunctions::singular_fct_exact_asymptotic_in_plane_zeta,
      // // 	&Global_Parameters::SingularFunctions::gradient_of_singular_fct_exact_asymptotic_in_plane_zeta,
      // // 	Sing_fct_id_in_plane_rotation);

      // // QUEHACERES back in! 15/1
      // // // QUEHACERES taking out the proper rotational bit for the time being
      // // In-plane rotation singular function ----------------------------------
      // sing_el_pt->add_unscaled_singular_fct_and_gradient_pt(
      // 	&Global_Parameters::SingularFunctions::singular_fct_in_plane_rotation,
      // 	&Global_Parameters::SingularFunctions::gradient_of_singular_fct_in_plane_rotation,
      // 	Sing_fct_id_in_plane_rotation);
    }

    // set the function pointer to the function which computes dzeta/dx
    sing_el_pt->dzeta_dx_fct_pt() = &CoordinateConversions::dzeta_dx;
  }

  // add the elements in the torus region to the torus region mesh, so that it
  // can be used to compute the Z2 error
  unsigned nel = Bulk_mesh_pt->nregion_element(Torus_region_id);
  for (unsigned e=0; e<nel; e++)
  {
    ELEMENT* el_pt = dynamic_cast<ELEMENT*>(
    Bulk_mesh_pt->region_element_pt(Torus_region_id, e));
    
    Torus_region_mesh_pt->add_element_pt(el_pt);
  }

  // tell the elements on the lower side of the disk that they're on the lower disk
  for(ELEMENT* el_pt : Elements_on_lower_disk_surface_pt)
  {
    el_pt->set_lower_disk_element();
  }
  
  // Apply bcs  
  apply_boundary_conditions();
}

//==start_of_find_corresponding_element_on_bulk_side_of_augmented_boundary==
/// \short Helper function which takes a boundary element in the augmented 
/// region and finds the corresponding boundary element in the bulk region
//========================================================================
template <class ELEMENT>
void FlowAroundDiskProblem<ELEMENT>::
find_corresponding_element_on_bulk_side_of_augmented_boundary(ELEMENT* const& augmented_elem_pt,
							      const int& face_index,
							      ELEMENT*& corresponding_elem_pt) const
{
  // First, build a temporary face element onto this
  // augmented region bulk element
  // (using TractionElements as a random face element since we only want the nodes)
  auto face_el_pt =
    std::make_unique<NavierStokesTractionElement<ELEMENT>>(augmented_elem_pt, face_index);
  
  // assemble a list of this elements nodes
  Vector<Node*> face_el_nodes;
  for(unsigned j=0; j<face_el_pt->nnode(); j++)
  {
    face_el_nodes.push_back(face_el_pt->node_pt(j));
  }

  // initialise pointer which will store the bulk region element we're looking for
  corresponding_elem_pt = nullptr;

  // now get the elements associated with the first node in this face element
  std::set<ELEMENT*> element_pt_set = Node_to_element_map.at(face_el_nodes[0]);

  // now look through each of these elements and check which one contains all the
  // face nodes
  for(ELEMENT* elem_pt : element_pt_set)
  {
    // don't want to find the same element!
    if(elem_pt == augmented_elem_pt)
      continue;
    
    unsigned matching_node_count = 0;
    for(unsigned k=0; k<elem_pt->nnode(); k++)
    {
      Node* node_pt = elem_pt->node_pt(k);
      
      if(std::find(face_el_nodes.begin(), face_el_nodes.end(), node_pt)
	 != face_el_nodes.end() )
	matching_node_count++;
    }

    // if the count is equal to the number of face nodes then we've found a
    // bulk region bulk element which shares the face we're interested in, so we're done
    if(matching_node_count == face_el_nodes.size())
    {
      corresponding_elem_pt = elem_pt;
      break;
    }
  }

  if(corresponding_elem_pt == nullptr)
  {
    std::ostringstream error_message;
    error_message << "Error: couldn't find a bulk element which shares a face with "
		  << "this augmented region boundary element.\n";
    
    throw OomphLibError(error_message.str(),
			OOMPH_CURRENT_FUNCTION,
			OOMPH_EXCEPTION_LOCATION);
  }
}

//==start_of_create_face_elements=========================================
/// \short Helper function to create face elements needed to:
/// - impose Dirichlet BCs
/// - impose the additional traction required from the augmented region to
///   the surrounding bulk elements
/// - compute the reciprocity integral to determine the singular amplitude
//========================================================================
template <class ELEMENT>
void FlowAroundDiskProblem<ELEMENT>::create_face_elements()
{  
  // Flux jump elements on boundary of torus
  //----------------------------------------
  // NOTE: Since these duplicate nodes, these elements must be
  //----------------------------------------------------------
  //       constructed first!
  //       ------------------
  if(!CommandLineArgs::command_line_flag_has_been_set("--dont_subtract_singularity"))
  {    
    // hierher
    ofstream some_file;
    std::ostringstream filename;
    
    filename << Doc_info.directory() << "/stress_jump_elements.dat";

    some_file.open(filename.str().c_str());
        
    // Where are we?
    unsigned region_id = Torus_region_id;
    for (unsigned b = First_torus_boundary_id;
	 b <= Last_torus_boundary_id; b++)
    {
      unsigned nel = Bulk_mesh_pt->nboundary_element_in_region(b, region_id);
      for (unsigned e=0; e<nel; e++)
      {
	ELEMENT* el_pt = dynamic_cast<ELEMENT*>(
	  Bulk_mesh_pt->boundary_element_in_region_pt(b, region_id, e));
        
	// What is the index of the face of the bulk element at the boundary
	int face_index = Bulk_mesh_pt->
	  face_index_at_boundary_in_region(b, region_id, e);

	// get the non-augmented element corresponding to this augmented element
	// and face index from the map of maps
	ELEMENT* non_aug_el_pt =
	  Torus_boundary_aug_to_non_aug_by_face_index_map.at(el_pt).at(face_index);
	
	  // Torus_boundary_element_augmented_to_bulk_map.at(el_pt);
	
	// Build the corresponding flux jump element
	NavierStokesWithSingularityStressJumpFaceElement<ELEMENT>* 
	  stress_jump_element_pt =
	  new NavierStokesWithSingularityStressJumpFaceElement<ELEMENT>
	  (el_pt, face_index, non_aug_el_pt, Stress_jump_duplicate_node_map,
	   Stress_jump_el_id, Lambda_hat_hat_id);
        
	//Add the flux jump element to the mesh
	Face_mesh_for_stress_jump_pt->add_element_pt(stress_jump_element_pt);

	// hierher
	stress_jump_element_pt->output(some_file);
      }
    }

    // hierher
    some_file.close();
     
    // Now add all duplicated nodes to the meshes and sort out boundary lookups
    filename.str("");
    filename << Doc_info.directory() << "/duplicated_nodes.dat";
    some_file.open(filename.str().c_str());
   
    for (std::map<Node*,Node*>::iterator it = Stress_jump_duplicate_node_map.begin();
	 it != Stress_jump_duplicate_node_map.end(); it++)
    {
      Node* original_node_pt = it->first;
      Node* new_node_pt = it->second;

      // add the new one to the stress jump mesh
      Face_mesh_for_stress_jump_pt->add_node_pt(new_node_pt);

      // find out what boundaries the old node was on in the bulk mesh
      std::set< unsigned >* boundaries_pt;
      original_node_pt->get_boundaries_pt(boundaries_pt);

      // now add the new node to the same boundaries in the bulk mesh
      // so that the boundary_node_pt() lookups still work
      for(std::set<unsigned>::iterator it_bound = boundaries_pt->begin();
	  it_bound != boundaries_pt->end(); it_bound++)
      {
	Bulk_mesh_pt->add_boundary_node(*it_bound, new_node_pt);
      }
      
      some_file << it->second->x(0) << " " 
		<< it->second->x(1) << " " 
		<< it->second->x(2) << " "
		<< std::endl;
    }
    some_file.close();
  
  
    // Now loop over bulk elements in torus region ("torus" around singularity)
    //-------------------------------------------------------------------------
    // and swap over any of their nodes that have been replaced
    //---------------------------------------------------------
    region_id = Torus_region_id;
    unsigned n_el = Bulk_mesh_pt->nregion_element(region_id);
    for (unsigned e=0; e<n_el; e++)
    {
      ELEMENT* bulk_el_pt = dynamic_cast<ELEMENT*>(
	Bulk_mesh_pt->region_element_pt(region_id, e));
     
      // Loop over all nodes and check if they're amongst the replaced
      // ones
      unsigned nnod = bulk_el_pt->nnode();
      for (unsigned j=0; j<nnod; j++)
      {
	Node* nod_pt = bulk_el_pt->node_pt(j);
       
	// Find this original node in the map; if we find it
	// it's already been duplicated
	std::map<Node*,Node*>::iterator it = Stress_jump_duplicate_node_map.find(nod_pt);
	if (it != Stress_jump_duplicate_node_map.end())
	{
	  // Use the existing duplicate node
	  bulk_el_pt->node_pt(j) = it->second;
	}
      }   
    }

    // BC elements live on disk inside torus
    //--------------------------------------
          
    // now add BC face elements to the disk boundaries
    for (unsigned i=0; i<Disk_boundary_ids_in_torus.size(); i++)
    {
      unsigned b = Disk_boundary_ids_in_torus[i];
      unsigned n_element = Bulk_mesh_pt->nboundary_element(b);
      
      // Loop over the bulk elements adjacent to boundary b
      for(unsigned e=0; e<n_element; e++)
      {
	// Get pointer to the bulk element that is adjacent to boundary b
	ELEMENT* bulk_elem_pt = dynamic_cast<ELEMENT*>(
	  Bulk_mesh_pt->boundary_element_pt(b,e));
        
	//Find the index of the face of element e along boundary b 
	int face_index = Bulk_mesh_pt->face_index_at_boundary(b,e);
        
	// Build the corresponding bc element
	NavierStokesWithSingularityBCFaceElement<ELEMENT>* bc_element_pt =
	  new NavierStokesWithSingularityBCFaceElement<ELEMENT>
	  (bulk_elem_pt,face_index, BC_el_id);
        
	//Add the bc element to the surface mesh
	Face_mesh_for_bc_pt->add_element_pt(bc_element_pt);
      }
    }    
  } // end if(! --dont_subtract_singularity)
   
  // ========================================================================
  // Now find the upper outer boundary and stick traction elements on it
  
  for(unsigned ibound = First_boundary_id_for_outer_boundary;
      ibound<First_boundary_id_for_outer_boundary+6; ibound++)
  {
    // get a pointer to the first element on this outer face
    ELEMENT* el_pt = dynamic_cast<ELEMENT*>(
      Bulk_mesh_pt->boundary_element_pt(ibound,0));
     
    // What is the index of the face of the bulk element at the boundary
    int face_index = Bulk_mesh_pt->face_index_at_boundary(ibound,0);

    // Build the corresponding face element
    NavierStokesTractionElement<ELEMENT>* surface_element_pt =
      new NavierStokesTractionElement<ELEMENT>(el_pt, face_index);

    // get the outer unit normal
    Vector<double> outer_unit_normal(3);
    surface_element_pt->outer_unit_normal(0, outer_unit_normal);

    // clean up
    delete surface_element_pt;
    surface_element_pt = nullptr;
    
    // check if we've got the top face, i.e. with n = (0,0,1)
    double tol = 1e-6;

    if(abs(outer_unit_normal[0])   > tol ||
       abs(outer_unit_normal[1])   > tol ||
       abs(outer_unit_normal[2]-1) > tol )
    {     
      continue;
    }

    // now we've found it, save this boundary ID to save us searching for it again    
    Outer_traction_boundary_id = ibound;
           
    unsigned n_element = Bulk_mesh_pt->nboundary_element(ibound);
    for(unsigned e=0; e<n_element; e++)
    {
      //Create Pointer to bulk element adjacent to the boundary
      ELEMENT* bulk_elem_pt = dynamic_cast<ELEMENT*>
	(Bulk_mesh_pt->boundary_element_pt(ibound, e));
         
      //Get Face index of boundary in the bulk element
      int face_index = Bulk_mesh_pt->face_index_at_boundary(ibound,e);
         
      //Create corresponding face element
      NavierStokesPdeConstrainedOptimisationTractionElement<ELEMENT>*
	traction_element_pt =
	new NavierStokesPdeConstrainedOptimisationTractionElement<ELEMENT>
	(bulk_elem_pt, face_index);

      // Set the pointer to the prescribed traction function
      traction_element_pt->traction_fct_pt() =
	&Analytic_Functions::prescribed_traction;

      //Attach it to the mesh
      Traction_boundary_condition_mesh_pt->add_element_pt(traction_element_pt);
    }
  }

} // end of create_face_elements()

//==start_of_apply_bc=====================================================
/// Helper function to apply boundary conditions
//========================================================================
template <class ELEMENT>
void FlowAroundDiskProblem<ELEMENT>::apply_boundary_conditions()
{
  // Identify boundary ids of pinned nodes 
  Vector<unsigned> pinned_boundary_id;
  
  for (unsigned ibound = First_lower_disk_boundary_id;
       ibound <= Last_lower_disk_boundary_id; ibound++)
  {
    pinned_boundary_id.push_back(ibound);
  }

  for (unsigned ibound = First_upper_disk_boundary_id;
       ibound <= Last_upper_disk_boundary_id; ibound++)
  {
    pinned_boundary_id.push_back(ibound);
  }
  
  for (unsigned ibound = First_boundary_id_for_outer_boundary;
       ibound < First_boundary_id_for_outer_boundary+6; ibound++)
  {
    pinned_boundary_id.push_back(ibound);
  }

  // Loop over pinned boundaries
  unsigned num_pin_bnd = pinned_boundary_id.size();
  for (unsigned bnd=0; bnd<num_pin_bnd; bnd++)
  {
    unsigned ibound = pinned_boundary_id[bnd];
    unsigned num_nod = Bulk_mesh_pt->nboundary_node(ibound);
    
    if (num_nod == 0)
    {
      std::ostringstream error_message;
      error_message << "No boundary nodes on boundary " 
  		    << ibound << "! Something's gone wrong!\n";
      throw OomphLibError(error_message.str(),
  			  OOMPH_CURRENT_FUNCTION,
  			  OOMPH_EXCEPTION_LOCATION);
    }

    // if we're doing a traction problem and this is the traction boundary, then
    // we don't want to pin it with Dirchlet conditions
    
    if(ibound == Outer_traction_boundary_id)
      continue;
    
    for (unsigned inod=0; inod<num_nod; inod++)
    {
      // grab a pointer to this boundary node
      Node* node_pt = Bulk_mesh_pt->boundary_node_pt(ibound, inod);

      Vector<double> x(3);
      x[0] = node_pt->x(0);
      x[1] = node_pt->x(1);
      x[2] = node_pt->x(2);

      // get the sum of all analytic solution modes on the boundary
      Vector<double> u(4, 0.0);
      
      if((ibound >= First_boundary_id_for_outer_boundary &&
	 ibound < First_lower_disk_boundary_id) )
      {
	// if we're doing no-slip on the outer boundaries then
	// just leave u_i = 0, otherwise assign the exact
	// flat-disk solution
	if(!Global_Parameters::No_slip_boundaries)
	{
	  Analytic_Functions::exact_solution_flat_disk(x, u);
	}
      }
      else
      {
	// we're on the disk, apply the prescribed motion
	u = Global_Parameters::disk_velocity(x);
      }
      
      // set and pin 'em
      for(unsigned i=0; i<3; i++)
      {
	node_pt->pin(i);
	node_pt->set_value(i, u[i]);
      }
    }
    
    // Now pin the momentum-enforcing Lagrange multipliers on the outer
    // Dirichlet boundaries

    // loop over the bulk elements on this boundary
    unsigned nel = Bulk_mesh_pt->nboundary_element(ibound);
    
    // QUEHACERES debug
    if(nel == 0)
    {
      oomph_info << "No elements on boundary " << bnd << std::endl;
      abort();
    }
    
    for(unsigned e=0; e<nel; e++)
    {
      // get the boundary bulk element
      ELEMENT* bulk_elem_pt = dynamic_cast<ELEMENT*>(
	Bulk_mesh_pt->boundary_element_pt(ibound, e));

      //Get Face index of boundary in the bulk element
      int face_index = Bulk_mesh_pt->face_index_at_boundary(ibound, e);

      // now loop over the boundary nodes to get their bulk node numbers,
      // so we can tell the bulk element to pin the LMs at those nodes
      for(unsigned j=0; j<bulk_elem_pt->nnode_on_face(); j++)
      {
	// get the number of this node in the bulk elements numbering scheme
	unsigned node_number_in_bulk = bulk_elem_pt->get_bulk_node_number(face_index, j);
	
	// now tell the bulk element to pin the pde-enforcing LMs
	bulk_elem_pt->pin_momentum_lagrange_multipliers(node_number_in_bulk);
      }
    }
  }
   
  // Now unpin nodal values where the bc conditions are enforced
  // by Lagrange multiplier to ensure that the sum of fe and singular
  // solution is correct
  unsigned nel = Face_mesh_for_bc_pt->nelement();
  for (unsigned e=0; e<nel; e++)
  {
    // Get element
    NavierStokesWithSingularityBCFaceElement<ELEMENT>* el_pt =
      dynamic_cast<NavierStokesWithSingularityBCFaceElement<ELEMENT>*>(
	Face_mesh_for_bc_pt->element_pt(e));
     
    // number of nodes in this face element
    unsigned nnod = el_pt->nnode();
      
    // matrix to store velocities at each boundary node
    DenseMatrix<double> nodal_boundary_value(nnod, Dim);

    // Unpin the FE part of the solution and pin bulk Lagrange multipliers
    for (unsigned j=0; j<nnod; j++)
    {
      Node* node_pt = el_pt->node_pt(j);

      // Now need to check if we have any lambda-hat-hats (the set of LMs which
      // enforce continuity of the bulk momentum-enforcing LMs across the
      // augmented boundary), and if so we need to pin them, because the
      // momentum-enforcing LMs are pinned on this Dirichlet boundary.
      
      // to find out if this node has them, we'll grab a map of the IDs of the
      // additional nodal values that have been added to this node...
      map<unsigned, unsigned> map_l = *(
	dynamic_cast<BoundaryNodeBase*>(node_pt)->
	index_of_first_value_assigned_by_face_element_pt() );

      // and then attempt to access the lambda-hat-hat ID from the map
      // (map::at() throws an exception if the key isn't found)
      try
      {
	map_l.at(Lambda_hat_hat_id);

	// if map::at() didn't throw an exception then we have lambda-hat-hats
	// at this node, so pin them
	for(unsigned i=0; i<Dim; i++)
	{
	  // pin lambda-hat-hat if we've got two sets
	  el_pt->pin_lagrange_multiplier_at_specified_local_node(j, i, Lambda_hat_hat_id);
	}
      }
      catch(...)
      {
	// if map::at() threw an exception then the lambda-hat-hat ID wasn't in
	// the map of additional nodal values, so we don't have any and don't
	// need to do anything else
      }

      // Now unpin the FE bit of the solution and get the Eulerian coordinates
      // of this node
      Vector<double> x(Dim, 0.0);
      for(unsigned i=0; i<Dim; i++)
      {
	el_pt->unpin_u_fe_at_specified_local_node(j, i);
	x[i] = node_pt->x(i);
      }

      {
	// total velocity is rigid body contribution + rotational contribution
	Vector<double> u_disk = Global_Parameters::disk_velocity(x);
	
	// assign to the matrix of nodal values
	for(unsigned i=0; i<Dim; i++)
	{
	  // total velocity is rigid body contribution + rotational contribution
	  nodal_boundary_value(j,i) = u_disk[i];	  
	}
      }

      // Now pin the Lagrange multipliers that enforce the PDE-constraints,
      // since we're applying Dirichlet BCs here
      // ------------------------------------------------------------------
      
      // get the bulk element this face element is attached to
      ELEMENT* bulk_el_pt = dynamic_cast<ELEMENT*>(el_pt->bulk_element_pt());

      // get the number of this node in the bulk elements numbering scheme
      unsigned node_number_in_bulk = el_pt->bulk_node_number(j);
      
      // now tell the bulk element to pin the pde-enforcing LMs
      bulk_el_pt->pin_momentum_lagrange_multipliers(node_number_in_bulk);   

      // QUEHACERES delete at some point, seems we shouldn't pin \lambda_p
      // // and we need to pin \lambda_p at a single node to set
      // // it's level (like pressure), but should be unpinned everywhere else
      if( (e==0) && (j==1) )
      {
	// oomph_info << "Pinning \lambda_p at ";

	// for(unsigned i=0; i<Dim; i++)
	//   oomph_info << x[i] << " ";
	
	// oomph_info << std::endl;
	
	// bulk_el_pt->pin_lagrange_multiplier_p(node_number_in_bulk, 0.0);      
      }
            
    } // end loop over BC nodes
    
    // Tell the element about these nodal boundary values
    el_pt->set_nodal_boundary_values(nodal_boundary_value);

    
  } // end loop over bc face elements

  // // QUEHACERES try pinning lambda_p in the non-aug region
  // {
  //   ELEMENT* el_pt = dynamic_cast<ELEMENT*>(Bulk_mesh_pt->region_element_pt(0,0));
  //   el_pt->pin_lagrange_multiplier_p(0, 0.0);
  // }
  
  
  // pin different eigenfunctions if requested
  if(CommandLineArgs::command_line_flag_has_been_set("--pin_broadside_amplitude"))
  {
    pin_singular_function(Sing_fct_id_broadside,
			  Global_Parameters::broadside_amplitude_pin);
  }
  if(CommandLineArgs::command_line_flag_has_been_set("--pin_broadside_rotation_amplitude"))
  {
    pin_singular_function(Sing_fct_id_broadside_rotation,
			  Global_Parameters::broadside_rotation_amplitude_pin);
  }
  if(CommandLineArgs::command_line_flag_has_been_set("--pin_in_plane_amplitude"))
  {
    pin_singular_function(Sing_fct_id_in_plane,
			  Global_Parameters::in_plane_amplitude_pin);
  }
  if(CommandLineArgs::command_line_flag_has_been_set("--pin_in_plane_rotation_amplitude"))
  {
    pin_singular_function(Sing_fct_id_in_plane_rotation,
			  Global_Parameters::in_plane_rotation_amplitude_pin);
  }

} // end apply BCs

//== start of output_submesh_pin_status ==================================
/// Function to output the pin status of each nodal dof in each submesh
//========================================================================
template <class ELEMENT>
void FlowAroundDiskProblem<ELEMENT>::output_submesh_dof_pin_status() const
{
  ofstream pin_file;
  std::ostringstream filename;

  // make a vector of pointers to the submeshes
  Vector<Mesh*> mesh_pt;
  mesh_pt.push_back(Torus_region_mesh_pt);
  mesh_pt.push_back(Face_mesh_for_bc_pt);
  mesh_pt.push_back(Face_mesh_for_stress_jump_pt);
  mesh_pt.push_back(Bulk_mesh_pt);
  mesh_pt.push_back(Singular_fct_element_mesh_pt);
  
  // makea vector of names for each sub-mesh
  Vector<std::string> mesh_name;
  mesh_name.push_back("torus_region");
  mesh_name.push_back("bc");
  mesh_name.push_back("stress_jump");
  mesh_name.push_back("bulk_mesh");
  mesh_name.push_back("singular_line_mesh");
  
  // loop over each submesh 
  for(unsigned k=0; k<mesh_pt.size(); k++)
  {
    // clear the filename
    filename.str("");

    // generate a filename specific to the submesh name
    filename << Doc_info.directory() << "/pinned_nodes_"
	     << mesh_name[k] << ".dat";
    
    pin_file.open(filename.str().c_str());

    std::set<Node*> submesh_node_pt;
    
    // now loop over the elements in this submesh and
    // get pointers to nodes that we haven't output yet
    for(unsigned e=0; e<mesh_pt[k]->nelement(); e++)
    {
      // get a pointer to this element
      FiniteElement* el_pt =
	dynamic_cast<FiniteElement*>(mesh_pt[k]->element_pt(e));

      // now loop over its nodes
      for(unsigned j=0; j<el_pt->nnode(); j++)
      {
	Node* node_pt = el_pt->node_pt(j);

	// check we haven't already output this node from a neighbouring element
	if(submesh_node_pt.find(node_pt) == submesh_node_pt.end())
	{
	  // output the Cartesian coordinates of this node
	  for(unsigned i=0; i<3; i++)
	    pin_file << node_pt->x(i) << " ";

	  // now loop over the nodal dofs and output their pin status
	  for(unsigned i=0; i<node_pt->nvalue(); i++)
	    pin_file << node_pt->is_pinned(i) << " ";

	  // and output the node pointer for good measure
	  pin_file << node_pt << std::endl;

	  // add the node pointer to the list of nodes we've outputted
	  submesh_node_pt.insert(node_pt);
	}
      }
    }
  
    pin_file.close();
  }

  // non-augmented region of bulk mesh
  // ---------------------------------------------------------
  {
    // clear the filename
    filename.str("");

    // generate a filename specific to the submesh name
    filename << Doc_info.directory() << "/pinned_nodes_non-augmented_region.dat";
    
    pin_file.open(filename.str().c_str());

    std::set<Node*> submesh_node_pt;

    // non-augmented region ID
    unsigned nonaug_region_id = 0;
  
    // now loop over the elements in the non-augmented region (region 0)
    for(unsigned e=0; e<Bulk_mesh_pt->nregion_element(nonaug_region_id); e++)
    {
      // get a pointer to this element
      FiniteElement* el_pt =
	dynamic_cast<FiniteElement*>(Bulk_mesh_pt->region_element_pt(nonaug_region_id, e));

      // now loop over its nodes
      for(unsigned j=0; j<el_pt->nnode(); j++)
      {
	Node* node_pt = el_pt->node_pt(j);

	// check we haven't already output this node from a neighbouring element
	if(submesh_node_pt.find(node_pt) == submesh_node_pt.end())
	{
	  // output the Cartesian coordinates of this node
	  for(unsigned i=0; i<3; i++)
	    pin_file << node_pt->x(i) << " ";

	  // now loop over the nodal dofs and output their pin status
	  for(unsigned i=0; i<node_pt->nvalue(); i++)
	    pin_file << node_pt->is_pinned(i) << " ";

	  // and output the node pointer for good measure
	  pin_file << node_pt << std::endl;

	  // add the node pointer to the list of nodes we've outputted
	  submesh_node_pt.insert(node_pt);
	}
      }
    }
    
    pin_file.close();
  }
}

//== start of set_values_to_singular_solution ============================
/// Function to assign the singular solution to all nodes of the mesh
//========================================================================
template <class ELEMENT>
void FlowAroundDiskProblem<ELEMENT>::set_values_to_singular_solution(
  const bool& broadside)
{
  std::ofstream outfile;
  std::ostringstream filename;

  filename << Doc_info.directory() << "/divergence_in_torus_region.dat";

  outfile.open(filename.str().c_str());
  
  unsigned nel = Torus_region_mesh_pt->nelement();
  for(unsigned e=0; e<nel; e++)
  {
    ELEMENT* el_pt = dynamic_cast<ELEMENT*>(Torus_region_mesh_pt->element_pt(e));

    // see if this element is in our list of lower disk elements
    bool on_lower_disk = Elements_on_lower_disk_surface_pt.find(el_pt) !=
      Elements_on_lower_disk_surface_pt.end();
    
    // loop over the nodes in this element
    for(unsigned j=0; j<el_pt->nnode(); j++)
    {
      // get a pointer to this node
      Node* node_pt = el_pt->node_pt(j);
      
      // get Cartesian coordinates of this node
      Vector<double> x(3, 0.0);

      for(unsigned i=0; i<3; i++)
	x[i] = node_pt->x(i);

      // now get the Lagrangian edge coordinates
      // -------------------------------------
      
      // zeta for flat disk
      double zeta = atan2pi(x[1], x[0]);

      double b_dummy = 0;
      mVector x_disk_edge(3);
      mVector tangent(3);
      mVector binormal(3);
      mVector normal(3);
    
      // get the unit normal from the disk-like geometric object at this zeta
      Global_Parameters::Warped_disk_with_boundary_pt->
	boundary_triad(b_dummy, zeta, x_disk_edge, tangent,
		       normal, binormal);

      mVector rho_vector = -(x_disk_edge - x);

      double rho = rho_vector.magnitude();
      
      // Moffat angle (minus sign accounts for the reflection of the moffat solution, which assumes
      // the semi-infinite plate is at x>0 not x<0 as we have with this coordinate system
      double phi = atan2pi(rho_vector*binormal, -rho_vector*normal);

      // catch the 2pi case
      double tol = 1e-6;
      if(abs(phi) < tol && on_lower_disk)
      {	  
	  phi = 2 * MathematicalConstants::Pi;		
      }

      // set the edge coordinates
      EdgeCoordinates edge_coords;
      edge_coords.rho  = rho;
      edge_coords.zeta = zeta;
      edge_coords.phi  = phi;

      // now get the (unscaled) singular solutions
      Vector<double> u_broadside =
	Global_Parameters::SingularFunctions::singular_fct_exact_asymptotic_broadside(edge_coords);
      Vector<double> u_in_plane  =
	Global_Parameters::SingularFunctions::singular_fct_exact_asymptotic_in_plane(edge_coords);

      // and set the nodal values
      for(unsigned i=0; i<3; i++)
	node_pt->set_value(i, u_broadside[i] + u_in_plane[i]);

      // and set the pressure if this node stores it
      if(node_pt->nvalue() % 3 != 0)	
	node_pt->set_value(3, u_broadside[3] + u_in_plane[3]);
      
    } // end loop over nodes

    // for debug, now output the divergence
    el_pt->output_divergence(outfile, Global_Parameters::Nplot_for_bulk);
  }

  outfile.close();
}


//== start of set_values_to_exact_non_singular_solution ==================
/// \short Function to assign the exact non-singular part of the solution
/// to all nodes of the mesh
//========================================================================
template <class ELEMENT>
void FlowAroundDiskProblem<ELEMENT>::set_values_to_exact_non_singular_solution() const
{
  oomph_info << "Setting initial values to exact non-singular solution...\n"
	     << std::endl;
  
  // do elements in the bulk region
  unsigned region_id = 0;
  unsigned nel = Bulk_mesh_pt->nregion_element(region_id);
  for(unsigned e=0; e<nel; e++)
  {
    // get a pointer to the current element
    ELEMENT* elem_pt = dynamic_cast<ELEMENT*>(
      Bulk_mesh_pt->region_element_pt(region_id, e));

    unsigned nnode = elem_pt->nnode();
    for(unsigned j=0; j<nnode; j++)
    {
      Node* node_pt = elem_pt->node_pt(j);
      
      // get the nodal coordinates
      Vector<double> x(3, 0.0);     
      node_pt->position(0,x);
      
      // get the total solution at this point (in the bulk, so
      // the FE solution is the total solution)
      Vector<double> u(4, 0.0);
      FlatDiskExactSolutions::total_exact_solution(x, Global_Parameters::u_disk_rigid_translation, 
						   Global_Parameters::omega_disk,
						   u);
      // set the velocities
      for(unsigned i=0; i<3; i++)
	node_pt->set_value(i, u[i]);

      // set the pressure if it's there
      if(j<4)
	node_pt->set_value(3, u[3]);
    }
  }

  // now do elements in the augmented region
  region_id = Torus_region_id;
  nel = Bulk_mesh_pt->nregion_element(region_id);
  for(unsigned e=0; e<nel; e++)
  {
    // get a pointer to the current element
    ELEMENT* elem_pt = dynamic_cast<ELEMENT*>(
      Bulk_mesh_pt->region_element_pt(region_id, e));

    unsigned nnode = elem_pt->nnode();
    for(unsigned j=0; j<nnode; j++)
    {
      Node* node_pt = elem_pt->node_pt(j);
      
      // get the nodal coordinates
      Vector<double> x(3, 0.0);     
      node_pt->position(0,x);
      
      // get the total solution at this point (in the bulk, so
      // the FE solution is the total solution)
      Vector<double> u_total(4, 0.0);
      FlatDiskExactSolutions::total_exact_solution(x, Global_Parameters::u_disk_rigid_translation, 
						   Global_Parameters::omega_disk,
						   u_total);

      EdgeCoordinates edge_coords;
      std::pair<GeomObject*, Vector<double> > line_element_and_local_coord;
      
      // now get the element which computes the scaled singular functions
      get_edge_coordinates_and_singular_element(elem_pt,
						x,
						edge_coords,
						line_element_and_local_coord);

      // cast the GeomObject to a singular line element      
      ScalableSingularityForNavierStokesLineElement<SINGULAR_ELEMENT_NNODE_1D>* sing_el_pt =
	dynamic_cast<ScalableSingularityForNavierStokesLineElement<SINGULAR_ELEMENT_NNODE_1D>*>
	(line_element_and_local_coord.first);

      // local coordinate in the singular element for the zeta of this plot point
      Vector<double> s_singular_el = line_element_and_local_coord.second;

      
      // get the total singular contribution at this point
      Vector<double> u_sing_total =
	sing_el_pt->total_singular_contribution(edge_coords,
						s_singular_el);
      
      // set the velocities to the non-singular bit
      for(unsigned i=0; i<3; i++)
	node_pt->set_value(i, u_total[i] - u_sing_total[i]);

      // set the pressure if it's there (if there are only velocities and
      // corresponding Lagrange multipliers there will be an odd number of values)

      // vertices are enumerated first, so nodes j <= Dim have pressures
      if(j <= Dim)
      {
	// catch the infinity case
	if(edge_coords.rho < 1e-8)
	{
	  double infinity = 10;
	  node_pt->set_value(3, 0.0); // infinity); //  * cos(edge_coords.zeta));
	}
	else
	{
	  double p = u_total[3] - u_sing_total[3];
	  node_pt->set_value(3, p);
	}
      }
    }
  }

}

//== start of validate_stress ============================================
/// Function to validate the singular stress function by assigning the singular velocity field
// to the nodes of the mesh, then computing the "FE" stress via the navier-stokes
// helper functions and comparing the two
//========================================================================
template <class ELEMENT>
void FlowAroundDiskProblem<ELEMENT>::validate_singular_stress(const bool& broadside)
{
  oomph_info << "Don't call this at the moment, validate_singular_stress() "
	     << "needs updating to use EdgeCoordinates\n";

  abort();
  
  oomph_info << "\nValidating singular stress...\n"
	     << "-----------------------------\n" << std::endl;

  double t_start = TimingHelpers::timer();
  
  // assign \hat u_i to the nodal values
  if(broadside)
  {
    set_values_to_singular_solution();
  }
  else
  {
    set_values_to_singular_solution(false);
  }
   

  char filename[500];

  if(broadside)
  {
    sprintf(filename, "%s/error_in_singular_stress_broadside.dat",
	    Doc_info.directory().c_str());
  }
  else
  {
    sprintf(filename, "%s/error_in_singular_stress_in_plane.dat",
	    Doc_info.directory().c_str());
  }
  
  // open the output file to record the error (tecplot format)
  ofstream stress_error_output(filename);
  
  // open the output file to record the error (plain format)
  if(broadside)
  {
    sprintf(filename, "%s/error_in_singular_stress_broadside_plain.dat",
	    Doc_info.directory().c_str());
  }
  else
  {
    sprintf(filename, "%s/error_in_singular_stress_in_plane_plain.dat",
	    Doc_info.directory().c_str());
  }
  
  ofstream stress_error_output_plain(filename);

  // column headers
  stress_error_output << "x,y,z,err_xx,err_xy,err_xz,err_yx,err_yy,err_yz,"
		      << "err_zx,err_zy,err_zz,sing_xx,sing_xy,sing_xz,"
		      << "sing_yx,sing_yy,sing_yz,sing_zx,sing_zy,sing_zz,"
		      <<"fd_xx,fd_xy,fd_xz,fd_yx,fd_yy,fd_yz,fd_zx,fd_zy,fd_zz"
		      << "p_sing";
  
  // number of plot points per side
  unsigned nplot = 2;

  oomph_info << "Computing singular and 'FE' stress...\n";

  // loop over all the elements in the torus region to compute the error in the stress
  const unsigned nel = Bulk_mesh_pt->nregion_element(Torus_region_id);
  
  for(unsigned e=0; e<nel; e++)
  {
    // get a pointer to this element
    ELEMENT* elem_pt = dynamic_cast<ELEMENT*>(Bulk_mesh_pt->region_element_pt(Torus_region_id,e));

    // dimension of this element
    const unsigned dim = elem_pt->dim();
  
    // write the tecplot header for this element
    stress_error_output << elem_pt->tecplot_zone_string(nplot);
    
    //Set the Vector to hold local coordinates
    Vector<double> s(dim);
 
    // Loop over plot points    
    unsigned num_plot_points = elem_pt->nplot_points(nplot);
    for (unsigned iplot=0; iplot < num_plot_points; iplot++)
    {
      // Get local coordinates of plot point
      elem_pt->get_s_plot(iplot, nplot, s);

      // global coordinates
      Vector<double> x(dim, 0.0);
      
      // get interpolated global coordinates
      for(unsigned i=0; i<dim; i++)
      { 
	x[i] = elem_pt->interpolated_x(s,i);
	stress_error_output << x[i] << " ";
	stress_error_output_plain << x[i] << " ";
      }

      // -----------------------------------------
      // singular stuff
      // -----------------------------------------
      
      // singular solution at this knot (don't care about velocity, just need the pressure for the stress)
      Vector<double> u_sing(dim+1, 0.0);

      if(broadside)
      {
		// QUEHACERES needs to be edge coordinates
	// u_sing = Global_Parameters::SingularFunctions::test_singular_function(x);
      }
      else
      {
	// QUEHACERES needs to be edge coordinates
	// u_sing = Analytic_Functions::singular_fct_in_plane(x);
      }
      
      // extract the singular pressure
      double p_sing = u_sing[dim];
	
      // get the singular velocity gradient
      DenseMatrix<double> du_dx_sing(dim, dim, 0.0);

      if(broadside)
      {
	// QUEHACERES needs to be edge coordinates
	// du_dx_sing = Global_Parameters::SingularFunctions::gradient_of_test_singular_function(x);	
      }
      else
      {
	// QUEHACERES needs to be edge coordinates
	// du_dx_sing = Analytic_Functions::gradient_of_singular_fct_in_plane(x);
      }
      
      // compute the singular strain rate
      DenseMatrix<double> strain_rate_sing(dim, dim, 0.0);

      for(unsigned i=0; i<dim; i++)
      {
	for(unsigned j=0; j<dim; j++)
	{
	  strain_rate_sing(i,j) = 0.5*(du_dx_sing(i,j) + du_dx_sing(j,i));
	}
      }

      // get the singular stress
      DenseMatrix<double> stress_sing(dim, dim, 0.0);
      stress_sing = Analytic_Functions::stress(strain_rate_sing, p_sing);

      // -----------------------------------------
      // "FE" stuff
      // -----------------------------------------

      // FE pressure
      double p_fe = elem_pt->interpolated_p_nst(s);
      
      // compute the "FE" strain-rate
      DenseMatrix<double> strain_rate_fe(dim, dim, 0.0);

      elem_pt->strain_rate(s, strain_rate_fe);

      // compute the "FE" stress
      DenseMatrix<double> stress_fe(dim, dim, 0.0);
      stress_fe = Analytic_Functions::stress(strain_rate_fe, p_fe);

      // -----------------------------------------
      // Error
      // -----------------------------------------
	
      // compute the error
      for(unsigned i=0; i<dim; i++)
      {
	for(unsigned j=0; j<dim; j++)
	{
	  // compute the error between the interpolated "FE" stress and the exact singular stress
	  double error = stress_fe(i,j) - stress_sing(i,j);

	  // output it
	  stress_error_output       << error << " ";
	  stress_error_output_plain << error << " ";
	}
      }
           
      // output actual singular stress
      for(unsigned i=0; i<dim; i++)
      {
      	for(unsigned j=0; j<dim; j++)
      	{
      	  stress_error_output       << stress_sing(i,j) << " ";
      	  stress_error_output_plain << stress_sing(i,j) << " ";
      	}
      }
      
      // output actual FD stress
      for(unsigned i=0; i<dim; i++)
      {
      	for(unsigned j=0; j<dim; j++)
      	{
      	  stress_error_output       << stress_fe(i,j) << " ";
      	  stress_error_output_plain << stress_fe(i,j) << " ";
      	}
      }

      stress_error_output       << p_sing;
      stress_error_output_plain << p_sing;
      
      stress_error_output       << std::endl;
      stress_error_output_plain << std::endl;
      
    } // end loop over plot point

    stress_error_output       << std::endl;    
    
    // Write tecplot footer (e.g. FE connectivity lists)
    elem_pt->write_tecplot_zone_footer(stress_error_output, nplot);
  } // end loop over elements
  
  // done, close the output file
  stress_error_output.close();
  stress_error_output_plain.close();
  
  oomph_info << "Finished singular stress validation after "
	     << TimingHelpers::timer() - t_start << "s.\n";
  
}


//==start_of_impose_fake_singular_amplitude===============================
/// Set the singular amplitude to a prescribed value and bypass the proper calculation
//========================================================================
template <class ELEMENT>
void FlowAroundDiskProblem<ELEMENT>::
impose_fake_singular_amplitude(const bool& impose_zero_amplitude)
{
  // tell all the elements in the singular element mesh about the fake
  // amplitude to impose
  for(unsigned e=0; e<Singular_fct_element_mesh_pt->nelement(); e++)
  {    
    // get a pointer to this singular line element in the upper mesh
    ScalableSingularityForNavierStokesLineElement<SINGULAR_ELEMENT_NNODE_1D>* sing_el_pt =
      dynamic_cast<ScalableSingularityForNavierStokesLineElement<SINGULAR_ELEMENT_NNODE_1D>*>(
	Singular_fct_element_mesh_pt->element_pt(e));

    // defaults here are so if we're subtracting the exact solution
    // we get the right amplitudes (the exact in-plane already has the azimuthal/rotational
    // bit built in, so the default amplitude for this should be zero)
    Vector<double> amplitude_broadside(sing_el_pt->nnode(), 1.0);
    
    // Global_Parameters::Imposed_singular_amplitude_broadside);
    Vector<double> amplitude_in_plane(sing_el_pt->nnode(), 1.0);
    
    // Global_Parameters::Imposed_singular_amplitude_in_plane);

    // same amplitude as for the in-plane, but will be modulated pi/2 out of phase
    Vector<double> amplitude_in_plane_rotation(sing_el_pt->nnode(), 0.0);
    // Global_Parameters::Imposed_singular_amplitude_in_plane);
    
    // if we're subtracting the exact solution for debug, the modulation is already
    // taken care of
    // // QUEHACERES taking this out for the zeta=0 rotated exact solution
    // if(!CommandLineArgs::command_line_flag_has_been_set("--subtract_exact_solution"))
    // {
    // now loop over the nodes and modulate the in-plane amplitude
    // by the azimuthal angle of each
    for(unsigned j=0; j<sing_el_pt->nnode(); j++)
    {
      // get the zeta for this node
      double zeta = sing_el_pt->zeta_nodal(j,0,0);

      // compute the amplitudes, initially set to zero
      Vector<double> amplitudes(Nsingular_function, 0.0);

      if(!impose_zero_amplitude)
      {
	amplitudes =
	  Global_Parameters::compute_singular_amplitudes_from_disk_velocity(zeta);
      }

      // interpret the vector
      amplitude_broadside[j]         = amplitudes[0];
      amplitude_in_plane[j]          = amplitudes[1];
      amplitude_in_plane_rotation[j] = amplitudes[2];

      if(CommandLineArgs::command_line_flag_has_been_set("--set_sing_amplitude_broadside"))
	amplitude_broadside[j] = Global_Parameters::Imposed_singular_amplitude_broadside;
    }

    if(!CommandLineArgs::command_line_flag_has_been_set("--subtract_exact_solution"))
    {
      sing_el_pt->impose_singular_fct_amplitude(Sing_fct_id_in_plane_rotation,
						amplitude_in_plane_rotation);
    }

    // pin the broadside singular function dof and set its amplitude to the
    // imposed amplitude
    sing_el_pt->impose_singular_fct_amplitude(Sing_fct_id_broadside,
					      amplitude_broadside);
    
    // pin the in-plane translation and rotation singular function dofs
    // and set their amplitude to the imposed amplitude
    sing_el_pt->impose_singular_fct_amplitude(Sing_fct_id_in_plane,
    					      amplitude_in_plane);
  }  
}

template <class ELEMENT>
void FlowAroundDiskProblem<ELEMENT>::pin_singular_function(const unsigned& sing_fct_id,
							   const double& val)
{
  unsigned nel = Singular_fct_element_mesh_pt->nelement();

  for(unsigned e=0; e<nel; e++)
  {
    // get a pointer to this singular line element in the upper mesh
    ScalableSingularityForNavierStokesLineElement<SINGULAR_ELEMENT_NNODE_1D>* sing_el_pt =
      dynamic_cast<ScalableSingularityForNavierStokesLineElement<SINGULAR_ELEMENT_NNODE_1D>*>(
	Singular_fct_element_mesh_pt->element_pt(e));

    Vector<double> amplitude(sing_el_pt->nnode(), val);
    
    sing_el_pt->impose_singular_fct_amplitude(sing_fct_id, amplitude);
  }
}

//////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////



//========================================================================
/// Doc the solution
//========================================================================
template <class ELEMENT>
void FlowAroundDiskProblem<ELEMENT>::doc_solution(const unsigned& nplot)
{
  oomph_info << "Documenting solution...\n"
	     << "------------------------\n" << std::endl;
  ofstream some_file;
  ofstream some_file2;
  ofstream face_some_file;
  ofstream coarse_some_file;
  char filename[500];
  
  bool do_bulk_output = true;
  if (CommandLineArgs::command_line_flag_has_been_set("--suppress_bulk_output"))
  {
    do_bulk_output = false;
  }
  
  if (CommandLineArgs::command_line_flag_has_been_set("--output_mesh_quality"))
  {
    oomph_info << "Outputting mesh quality..." << std::endl;
    // Doc mesh quality (Ratio of max. edge length to min. height,
    /// so if it's very large it's BAAAAAD)
    sprintf(filename,"%s/mesh_quality%i.dat",
	    Doc_info.directory().c_str(),
	    Doc_info.number()); 
    ofstream quality_file;
    quality_file.open(filename);
    if (do_bulk_output) Bulk_mesh_pt->assess_mesh_quality(quality_file);
    quality_file.close();
  }
  
  // Output bulk elements in torus region
  //-------------------------------------
  double volume_in_torus_region = 0.0;
  
  bool output_torus_region_soln = false;
  if (CommandLineArgs::command_line_flag_has_been_set("--output_soln_in_torus"))
  {
    oomph_info << "Outputting solution in torus region..." << std::endl;
    output_torus_region_soln = true;
  
    sprintf(filename,"%s/soln_in_torus_region%i.dat",Doc_info.directory().c_str(),
	    Doc_info.number());
    some_file.open(filename);
  }
  oomph_info << "Computing torus volume..." << std::endl;
  
  double functional_integral = 0.0;
  
  unsigned region_id = Torus_region_id;
  unsigned n_el = Bulk_mesh_pt->nregion_element(region_id);
  for (unsigned e=0; e<n_el; e++)
  {
    // get a pointer to the torus-region element
    ELEMENT* el_pt =
      dynamic_cast<ELEMENT*>(Bulk_mesh_pt->region_element_pt(region_id, e));
    
    if (output_torus_region_soln) 
    {
      el_pt->output(some_file, nplot);
    }
    volume_in_torus_region += el_pt->size();

    // Compute the integral of the functional over this element and add to total
    functional_integral += el_pt->integral_of_functional();
  }

  if(output_torus_region_soln)
    some_file.close();

  oomph_info << "Average element volume in torus region: "
	     << volume_in_torus_region / n_el << std::endl;

  oomph_info << "Integral of L2 velocity norm in torus region: "
	     << functional_integral << std::endl;
  
  // set the small but finite edge radius to use for outputting "infinite" pressures
  // at the edge of the disk
  if(Doc_info.number() == 0)
    Global_Parameters::Drho_for_infinity = pow(volume_in_torus_region / n_el, 1./3.)/50.;
  
  // *** QUEHACERES for debug
  {
    sprintf(filename,"%s/traction_bc_mesh%i.dat",Doc_info.directory().c_str(),
	    Doc_info.number());
    some_file.open(filename);
    
    n_el = Traction_boundary_condition_mesh_pt->nelement();
    for(unsigned e=0; e<n_el; e++)
    {
      dynamic_cast<NavierStokesTractionElement<ELEMENT>*>(
	Traction_boundary_condition_mesh_pt->element_pt(e))->output(some_file, nplot);
    }

    some_file.close();
  }
  // ***
  
  // --------------------------------------------------------------------------
  // Plot disks around the perimeter of the disk...
  // --------------------------------------------------------------------------
  if(Global_Parameters::Do_disk_on_disk_plots)
  {
    if (Geom_objects_are_out_of_date)
    {
      // Setup disk on disk plots
      setup_disk_on_disk_plots();
    
      // Now they're not...
      Geom_objects_are_out_of_date=false;
    }

    sprintf(filename,"%s/disk_on_disk%i.dat", Doc_info.directory().c_str(),
	    Doc_info.number());
    some_file.open(filename);
        
    Vector<double> x(3);
    
    for (unsigned k=0; k < Ndisk_on_disk_plot; k++)
    {
      some_file << "ZONE I=" << Nphi_disk_on_disk_plot 
		<< ", J=" << Nrho_disk_on_disk_plot << std::endl;
      
      for (unsigned i=0; i < Nrho_disk_on_disk_plot; i++)
      {
	for (unsigned j=0; j<Nphi_disk_on_disk_plot; j++)
	{		
	  (Disk_on_disk_plot_point[k][i][j].second.first)->
	    position(Disk_on_disk_plot_point[k][i][j].second.second,x);
	  
	  double rho = (Disk_on_disk_plot_point[k][i][j].first)[0];
	  double phi = (Disk_on_disk_plot_point[k][i][j].first)[1];

	  // get the interpolated velocity at this point
	  Vector<double> u(3);
	  dynamic_cast<ELEMENT*>(Disk_on_disk_plot_point[k][i][j].second.first)->
	    interpolated_u_nst(Disk_on_disk_plot_point[k][i][j].second.second, u);

	  // get the interpolated pressure at this point
	  double p = dynamic_cast<ELEMENT*>(Disk_on_disk_plot_point[k][i][j].second.first)->
	    interpolated_p_nst(Disk_on_disk_plot_point[k][i][j].second.second);
	  
	  some_file 
	    << x[0] << " " 
	    << x[1] << " " 
	    << x[2] << " " 
	    << u[0] << " "
	    << u[1] << " "
	    << u[2] << " "
	    << p    << " "
	    << rho << " " 
	    << phi << " "
	    // QUEHACERES put moffat solution here
	    // << Global_Parameters::asymptotic_solution(rho,phi) << " "
	    << std::endl;	  
	}
      }
    }
    some_file.close();
    
    for (unsigned k=0; k < Ndisk_on_disk_plot; k++)
    {
      sprintf(filename,"%s/disk_on_disk_gnuplot_disk%i_%i.dat", Doc_info.directory().c_str(),
	      k, Doc_info.number());
      
      some_file.open(filename);
      
      for (unsigned j=0; j<Nphi_disk_on_disk_plot; j++)
      {
	for (unsigned i=0; i < Nrho_disk_on_disk_plot; i++)
	{			
	  (Disk_on_disk_plot_point[k][i][j].second.first)->
	    position(Disk_on_disk_plot_point[k][i][j].second.second,x);
	  
	  double rho = (Disk_on_disk_plot_point[k][i][j].first)[0];
	  double phi = (Disk_on_disk_plot_point[k][i][j].first)[1];

	  // get the interpolated velocity at this point
	  Vector<double> u(3);
	  dynamic_cast<ELEMENT*>(Disk_on_disk_plot_point[k][i][j].second.first)->
	    interpolated_u_nst(Disk_on_disk_plot_point[k][i][j].second.second, u);

	  // get the interpolated pressure at this point
	  double p = dynamic_cast<ELEMENT*>(Disk_on_disk_plot_point[k][i][j].second.first)->
	    interpolated_p_nst(Disk_on_disk_plot_point[k][i][j].second.second);
	  
	  some_file 
	    << x[0] << " " 
	    << x[1] << " " 
	    << x[2] << " " 
	    << u[0] << " "
	    << u[1] << " "
	    << u[2] << " "
	    << p    << " "
	    << rho << " " 
	    << phi << " "
	    // QUEHACERES put moffat solution here
	    // << Global_Parameters::asymptotic_solution(rho,phi) << " "
	    << std::endl;	  
	}
	some_file << "\n\n";
      }
      some_file.close();
    }    
  }
  
  // Attach face elements to boundary of torus
  //------------------------------------------
  sprintf(filename,"%s/face_elements_on_boundary_of_torus%i.dat",
	  Doc_info.directory().c_str(),
	  Doc_info.number());
  some_file.open(filename);
  double torus_surface_area = 0.0;
  region_id = Torus_region_id;

  for (unsigned b = First_torus_boundary_id; b<=Last_torus_boundary_id; b++)
  {
    unsigned nel = Bulk_mesh_pt->nboundary_element_in_region(b,region_id);
    for (unsigned e=0; e<nel; e++)
    {
      FiniteElement* el_pt = 
	Bulk_mesh_pt->boundary_element_in_region_pt(b,region_id,e);
     
      // What is the index of the face of the bulk element at the boundary
      int face_index = Bulk_mesh_pt->
	face_index_at_boundary_in_region(b,region_id,e);

      // Build the corresponding face element
      auto surface_element_pt =
	std::make_unique<RigidDiskFaceElement<ELEMENT>>(el_pt, face_index);
     
      // Get surface area
      torus_surface_area += surface_element_pt->size();
     
      // Output
      surface_element_pt->output(some_file, nplot);
    }
  }
  
  some_file.close();

  oomph_info << "Torus surface area: " <<  torus_surface_area << std::endl;

  // if(Doc_info.number() > 0)
  // {
  //   sprintf(filename, "%s/augmented_boundary_traction_and_lms%i.csv",
  // 	    Doc_info.directory().c_str(), Doc_info.number());

  //   some_file.open(filename);

  //   // write the header
  //   some_file << "x,y,z,t_aug_x,t_aug_y,t_aug_z,t_bulk_x,t_bulk_y,t_bulk_z,"
  // 	      << "lambda_x,lambda_y,lambda_z" << std::endl;
    
  //   unsigned n_element = Face_mesh_for_stress_jump_pt->nelement();
  //   for(unsigned e=0; e<n_element; e++)
  //   {
  //     // Upcast from GeneralisedElement to the present element
  //     NavierStokesWithSingularityStressJumpFaceElement<ELEMENT>* stress_jump_el_pt =
  // 	dynamic_cast<NavierStokesWithSingularityStressJumpFaceElement<ELEMENT>*>(
  // 	  Face_mesh_for_stress_jump_pt->element_pt(e));

  //     stress_jump_el_pt->output_traction_and_lagrange_multipliers(some_file);
  //   }

  //   some_file.close();
    
  // }

  // ==========================================================================
  // compute the total force on the disk and surface area
  // ==========================================================================

  Vector<double> total_force_on_plate(3, 0.0);

  double disk_in_torus_surface_area      = 0.0;
  double disk_outside_torus_surface_area = 0.0;
  double disk_surface_area               = 0.0;
  
  // Attach face elements to part of disk inside torus
  //--------------------------------------------------
  sprintf(filename,"%s/face_elements_on_disk_in_torus%i.dat",
	  Doc_info.directory().c_str(),
	  Doc_info.number());
  some_file.open(filename);
  
  

  // Loop over the BC face elements, which sit on the upper and lower
  // disk within the torus region
  for(unsigned e=0; e<Face_mesh_for_bc_pt->nelement(); e++)
  {
    // get a pointer to the face element
    NavierStokesWithSingularityBCFaceElement<ELEMENT>* bc_el_pt =
      dynamic_cast<NavierStokesWithSingularityBCFaceElement<ELEMENT>*>(
  	Face_mesh_for_bc_pt->element_pt(e));

    // compute the integrated traction on this element
    Vector<double> force_on_element =
      bc_el_pt->get_contribution_to_normal_stress();

    // and add to the total
    for(unsigned i=0; i<3; i++)
      // minus sign accounts for the BC Face element being attached to a
      // bulk element on the other side of the plate
      total_force_on_plate[i] -= force_on_element[i];

    double surface_area = bc_el_pt->size();

    bc_el_pt->output(some_file, nplot);
    
    disk_surface_area += surface_area;
    disk_in_torus_surface_area += surface_area;
  }

  // QUEHACERES delete
  // region_id = Torus_region_id;
  // unsigned nb = One_based_boundary_id_for_disk_within_torus.size();
  // for (unsigned i=0; i<nb; i++)
  // {
  //   unsigned b = One_based_boundary_id_for_disk_within_torus[i]-1;
  //   unsigned nel = Bulk_mesh_pt->nboundary_element_in_region(b,region_id);
  //   for (unsigned e=0; e<nel; e++)
  //   {
  //     FiniteElement* el_pt = 
  // 	Bulk_mesh_pt->boundary_element_in_region_pt(b,region_id,e);
     
  //     // What is the index of the face of the bulk element at the boundary
  //     int face_index = Bulk_mesh_pt->
  // 	face_index_at_boundary_in_region(b,region_id,e);

  //     // Build the corresponding surface power element
  //     RigidDiskFaceElement<ELEMENT>* surface_element_pt =
  //     	new RigidDiskFaceElement<ELEMENT>(el_pt, face_index);
     
  //     // Get surface area
  //     disk_in_torus_surface_area += surface_element_pt->size();
      
  //     // Output
  //     surface_element_pt->output(some_file, nplot);
     
  //     // ...and we're done!
  //     delete surface_element_pt;
  //     surface_element_pt = nullptr;
  //   }
  // }
  some_file.close();
  oomph_info << "\nDisk in torus surface area: "
	     <<  disk_in_torus_surface_area << std::endl;

 
  // Attach face elements to part of disk outside torus
  //--------------------------------------------------
  sprintf(filename,"%s/face_elements_on_disk_outside_torus%i.dat",
	  Doc_info.directory().c_str(),
	  Doc_info.number());
  some_file.open(filename);
  
  region_id = 0; 
  unsigned nb = One_based_boundary_id_for_disk_outside_torus.size();
  for (unsigned i=0; i<nb; i++)
  {
    unsigned b = One_based_boundary_id_for_disk_outside_torus[i]-1;
    unsigned nel = Bulk_mesh_pt->nboundary_element_in_region(b,region_id);
    for (unsigned e=0; e<nel; e++)
    {
      FiniteElement* el_pt = 
	Bulk_mesh_pt->boundary_element_in_region_pt(b,region_id,e);
     
      // What is the index of the face of the bulk element at the boundary
      int face_index = Bulk_mesh_pt->
	face_index_at_boundary_in_region(b,region_id,e);
      
      // Build the corresponding face element
      auto surface_element_pt =
	std::make_unique<RigidDiskFaceElement<ELEMENT>>(el_pt, face_index);
     
      // Get surface area
      disk_outside_torus_surface_area += surface_element_pt->size();
      disk_surface_area += surface_element_pt->size();
      
      // Output      
      surface_element_pt->output(some_file, nplot);

      // compute the integrated traction on this element
      Vector<double> force_on_element =
  	surface_element_pt->get_contribution_to_normal_stress();

      // and add to the total
      for(unsigned i=0; i<3; i++)
	total_force_on_plate[i] += force_on_element[i];
    }
  }
  some_file.close();
  
  oomph_info << "Disk outside torus surface area: "
	     <<  disk_outside_torus_surface_area << std::endl;

  oomph_info << "Total surface area of disk: ";
  printf("%.03fpi", disk_surface_area / MathematicalConstants::Pi);
  oomph_info << std::endl;

  // Attach face elements to lower disk boundaries
  //--------------------------------------------------
  sprintf(filename,"%s/face_elements_on_lower_disk%i.dat",
	  Doc_info.directory().c_str(),
	  Doc_info.number());
  some_file.open(filename);
 
  for (unsigned b=First_lower_disk_boundary_id; b<=Last_lower_disk_boundary_id; b++)
  {
    unsigned nel = Bulk_mesh_pt->nboundary_element(b);
    for (unsigned e=0; e<nel; e++)
    {
      FiniteElement* el_pt = 
	Bulk_mesh_pt->boundary_element_pt(b,e);
     
      // What is the index of the face of the bulk element at the boundary
      int face_index = Bulk_mesh_pt->
	face_index_at_boundary(b,e);

      // Build the corresponding surface element
      auto surface_element_pt =
	std::make_unique<RigidDiskFaceElement<ELEMENT>>(el_pt, face_index);
     
      // Output
      surface_element_pt->output(some_file, nplot);
    }
  }
  some_file.close();

  oomph_info << "Total force on plate: "
	     << total_force_on_plate[0] << ", "
	     << total_force_on_plate[1] << ", "
	     << total_force_on_plate[2] << "\n" << std::endl;
  
  // output the stress jump elements
  if (!CommandLineArgs::command_line_flag_has_been_set("--dont_subtract_singularity"))      
  {   
    std::ostringstream filename;
    filename << Doc_info.directory() << "/stress_jump_elements"
	     << Doc_info.number() << ".dat";

    some_file.open(filename.str().c_str());

    unsigned nel = Face_mesh_for_stress_jump_pt->nelement();
    for (unsigned e=0; e<nel; e++)
    {
      NavierStokesWithSingularityStressJumpFaceElement<ELEMENT>* 
	stress_jump_element_pt 
	= dynamic_cast<NavierStokesWithSingularityStressJumpFaceElement<ELEMENT>*>(
	  Face_mesh_for_stress_jump_pt->element_pt(e));
		
      stress_jump_element_pt->output(some_file);
    }
    
    some_file.close();
  }

  // Output solution showing element outlines
  //-----------------------------------------
  oomph_info << "Outputting coarse solution..." << std::endl;
  sprintf(filename,"%s/coarse_soln%i.vtu",Doc_info.directory().c_str(),
	  Doc_info.number());
  some_file.open(filename);
  if (do_bulk_output) Bulk_mesh_pt->output_paraview(some_file,2);
  some_file.close();

  // Exact solution (only need to output it once)
  if (Doc_info.number() == 0 &&
    !CommandLineArgs::command_line_flag_has_been_set("--dont_output_exact_solution"))
  {
    oomph_info << "Outputting exact solution..." << std::endl;
    sprintf(filename,"%s/exact_soln.vtu",Doc_info.directory().c_str() );
    some_file.open(filename);

    Bulk_mesh_pt->output_fct_paraview(some_file, nplot,
				      Analytic_Functions::exact_solution_flat_disk);
    
    some_file.close();
  }

  // // Get norm of solution
  // //---------------------
  // sprintf(filename,"%s/norm%i.dat",Doc_info.directory().c_str(),
  // 	  Doc_info.number());
  // some_file.open(filename);
  // double norm_soln = 0.0;
  // Bulk_mesh_pt->compute_norm(norm_soln);  
  // some_file << sqrt(norm_soln) << std::endl;
  // oomph_info << "Norm of computed solution: "   << sqrt(norm_soln)  << endl;
  // some_file.close();

  // Get error from exact solution
  // -----------------------------
  
  oomph_info << "Computing error from exact solution..." << std::endl;
  
  // global norm and error
  double norm    = 0.0;
  double v_error = 0.0;
  double p_error = 0.0;

  sprintf(filename,"%s/error%i.dat",Doc_info.directory().c_str(),
      Doc_info.number());
  
  some_file.open(filename);

  // dummy stream to pass in which we won't open - the compute_error function
  // uses knot points rather than plot points, so can't format for paraview
  ofstream dummy_ofstream;
  
  for(unsigned e=0; e<Bulk_mesh_pt->nelement(); e++)
  {
    // get a pointer to this bulk element
    ELEMENT* el_pt = dynamic_cast<ELEMENT*>(Bulk_mesh_pt->element_pt(e));
      
    // elemental errors and norms
    double el_norm    = 0.0;
    double el_v_error = 0.0;
    double el_p_error = 0.0;
    
    //Calculate the elemental errors for each non-halo element
#ifdef OOMPH_HAS_MPI
    if (!(el_pt->is_halo()))
#endif
    {
    el_pt->compute_error(dummy_ofstream, &Analytic_Functions::exact_solution_flat_disk,
			 el_v_error, el_p_error, el_norm);
    }
    
    //Add each elemental error to the global error
    norm    += el_norm;
    v_error += el_v_error;
    p_error += el_p_error;
    
    // do the actual paraview'able output at plot points
    el_pt->output_error_at_plot_points(some_file, nplot,
				       &Analytic_Functions::exact_solution_flat_disk);
  }

  some_file.close();

  oomph_info << "L2 velocity error in total solution: " << v_error << std::endl;
  oomph_info << "L2 pressure error in total solution: " << p_error << std::endl;

  // // output divergence
  // // --------------------------------------------------------------------------

  // sprintf(filename,"%s/divergence%i.dat",Doc_info.directory().c_str(),
  //     Doc_info.number());
  // some_file.open(filename);
  // for(unsigned e=0; e<Bulk_mesh_pt->nelement(); e++)
  // {
  //   ELEMENT* el_pt = dynamic_cast<ELEMENT*>(Bulk_mesh_pt->element_pt(e));
    
  //   el_pt->output_divergence(some_file, Global_Parameters::Nplot_for_bulk);    
  // }

  // some_file.close();
  
  // output the Z2 error
  // ---------------------------------
  
  // // grab the error in each element from the Z2 estimator
  // Vector<double> elementwise_z2_error(Torus_region_mesh_pt->nelement());
  // Mesh* mesh_pt = dynamic_cast<Mesh*>(Torus_region_mesh_pt);

  // // Use actual value without normalisation!
  // Z2ErrorEstimator* z2_pt = dynamic_cast<Z2ErrorEstimator*>(
  //   Torus_region_mesh_pt->spatial_error_estimator_pt());

  // // keep a copy of the original normalisation
  // double backup = z2_pt->reference_flux_norm();

  // // set the normalisation to 1
  // z2_pt->reference_flux_norm() = 1.0;

  // // get the element-wise z2 errors
  // z2_pt->get_element_errors(mesh_pt, elementwise_z2_error);

  // // Reset the normalisation
  // z2_pt->reference_flux_norm() = backup;
   
  // double z2_integral = 0;

  // sprintf(filename,"%s/elementwise_Z2error%i.dat",
  //         Doc_info.directory().c_str(),
  // 	  Doc_info.number());
  // some_file.open(filename);
    
  // // sum the errors to get a global measure    
  // for(unsigned e=0; e<elementwise_z2_error.size(); e++)
  // {
  //   ELEMENT* el_pt = dynamic_cast<ELEMENT*>(Torus_region_mesh_pt->element_pt(e));

  //   unsigned npts = 5;
  //   Vector<double> s(3);
  //   unsigned num_plot_points = el_pt->nplot_points(npts);
  //   some_file << el_pt->tecplot_zone_string(npts);
    
  //   for(unsigned j=0; j<num_plot_points; j++)
  //   {
  //     el_pt->get_s_plot(j, npts, s);
  //     Vector<double> x(3);
  //     x[0] = el_pt->interpolated_x(s,0);
  //     x[1] = el_pt->interpolated_x(s,1);
  //     x[2] = el_pt->interpolated_x(s,2);
  //     some_file << x[0] << " " << x[1] << " " << x[2] << " "
  // 		<< elementwise_z2_error[e] << endl;
  //   }
    
  //   el_pt->write_tecplot_zone_footer(some_file, npts);
    
  //   // add the weighted conribution of this element to the integral
  //   z2_integral += elementwise_z2_error[e] * el_pt->size();
  // }
  // some_file.close();
  
  // sprintf(filename,"%s/torus_region_z2_error%i.dat", Doc_info.directory().c_str(),
  // 	  Doc_info.number());
  // some_file.open(filename);

  // some_file << z2_integral << " " << z2_integral / volume_in_torus_region << std::endl;
  // some_file.close();

  // oomph_info << "Torus region Z2 error: " << z2_integral << std::endl;
  
  if (CommandLineArgs::command_line_flag_has_been_set("--output_jacobian_sparse"))
  {  
    // residual vector and Jacobian matrix
    DoubleVector r;
    CRDoubleMatrix jac;

    // get 'em
    get_jacobian(r,jac);

    oomph_info << "Number of non-zero Jacobian entries: " << jac.nnz() << "\n";
    oomph_info << "Jacobian density: " << std::setprecision(2)
	       << 100.0 * jac.nnz() / (double)(jac.nrow() * jac.nrow())
	       << "%" << std::endl;
    
    sprintf(filename,"%s/residuals%i.dat", Doc_info.directory().c_str(),
	    Doc_info.number());
    some_file.open(filename);

    for(unsigned i=0; i<r.nrow(); i++)
    {
      some_file << i << " " << r[i] << std::endl;
    }

    some_file.close();
    
    sprintf(filename,"%s/jacobian_sparse%i.dat", Doc_info.directory().c_str(),
	    Doc_info.number());
    some_file.open(filename);

    bool output_bottom_right = true;
    unsigned precision = 0;
    jac.sparse_indexed_output(some_file, precision, output_bottom_right);

    some_file.close();

    oomph_info << "Output sparse Jacobian matrix to " << filename << "\n\n";
  }
  
  if (CommandLineArgs::command_line_flag_has_been_set("--output_jacobian_full"))
  {
    // residual vector and Jacobian matrix
    DoubleVector r;
    CRDoubleMatrix jac;

    // get 'em
    get_jacobian(r,jac);

    sprintf(filename,"%s/jacobian%i.dat",Doc_info.directory().c_str(),Doc_info.number());      
    some_file.open(filename);
      
    for(unsigned i=0; i<jac.nrow(); i++)
    {
      for(unsigned j=0; j<jac.ncol(); j++)
      {
	some_file << jac(i,j) << " ";
      }
      some_file << std::endl;
    }
    some_file.close();
    oomph_info << "\nOutput full Jacobian matrix to " << filename << std::endl;
  }
  
  if (CommandLineArgs::command_line_flag_has_been_set("--describe_dofs"))
  {
    sprintf(filename,"%s/describe_dofs.dat", Doc_info.directory().c_str());
    some_file.open(filename);
    describe_dofs(some_file);
    
    some_file.close();

    oomph_info << "Output description of dofs to " << filename << std::endl;
  }
  if (CommandLineArgs::command_line_flag_has_been_set("--describe_nodes"))
  {
    sprintf(filename,"%s/describe_nodes.dat", Doc_info.directory().c_str());
    some_file.open(filename);

    for(unsigned j=0; j<mesh_pt()->nnode(); j++)
    {
      // grab the node
      Node* node_pt = mesh_pt()->node_pt(j);

      // get it's coordinates
      double x = node_pt->x(0);
      double y = node_pt->x(1);
      double z = node_pt->x(2);

      some_file << j << " " << x << " " << y << " " << z << " "
		<< node_pt->nvalue() << " " << node_pt << std::endl;
    }
    
    some_file.close();
    oomph_info << "Output description of nodes to " << filename << std::endl;
  }

  oomph_info << "Outputting extended solution..." << std::endl;
  
  // Plot "extended solution" showing contributions
  sprintf(filename,"%s/extended_soln%i.dat",Doc_info.directory().c_str(),
	  Doc_info.number());
  some_file.open(filename);
  
  unsigned nel = Bulk_mesh_pt->nelement();
  for (unsigned e=0; e<nel; e++)
  {
    // shouldn't change this, since the maps have been setup for a specific number
    // unsigned npts = Global_Parameters::Nplot_for_bulk;
    ELEMENT* el_pt = dynamic_cast<ELEMENT*>(Bulk_mesh_pt->element_pt(e));

    el_pt->output_with_various_contributions(some_file, nplot,
					     Analytic_Functions::exact_solution_flat_disk);
  }

  some_file.close();    

  sprintf(filename, "%s/singular_line_mesh%i.dat",
	  Doc_info.directory().c_str(), Doc_info.number());
  some_file.open(filename);
  Singular_fct_element_mesh_pt->output(filename, nplot);
  some_file.close();

  //Increment counter for solutions 
  Doc_info.number()++;
  
  oomph_info << "Finished documenting solution.\n\n";
} // end of doc



// ////////////////////////////////////////////////////////////////////////////
// ////////////////////////////////////////////////////////////////////////////
// ////////////////////////////////////////////////////////////////////////////
// ////////////////////////////////////////////////////////////////////////////

void validate_finite_diff()
{
  // disk translational velocity
  Vector<double> u_disk(3, 0.0);

  // broadside
  u_disk[2] = 1; 

  // no rotations
  Vector<double> omega_disk(3, 0.0);

  // coordinates
  Vector<double> x(3, 0.0);

  // we'll use the same z coordinate for the time being
  x[2] = 0.5;

  double x_start = -1.5;
  double x_end   =  1.5;

  Vector<double> u_dummy(4, 0.0);

  ofstream some_file;
  some_file.open("dudx_fd_validation.txt");
  
  unsigned n = 50;
  for(unsigned i=0; i<n; i++)
  {
    for(unsigned j=0; j<n; j++)
    {      
      DenseMatrix<double> dudx_exact(3,3, 0.0);
      DenseMatrix<double> dudx_fd(3,3, 0.0);
      
      x[0] = (x_end - x_start) * double(i) / (double(n+1));
      x[1] = (x_end - x_start) * double(j) / (double(n+1));

      // get the exact version
      FlatDiskExactSolutions::gupta_solution_and_gradient(x, u_dummy, dudx_exact);
      
      // get the finite diff'd version
      FlatDiskExactSolutions::total_exact_velocity_gradient(x, u_disk, omega_disk, dudx_fd);

      // output
      for(unsigned k=0; k<3; k++)
	some_file << x[k] << " ";

      for(unsigned a=0; a<3; a++)
      {
	for(unsigned b=0; b<3; b++)
	{
	  some_file << dudx_exact(a,b) - dudx_fd(a,b) << " ";
	}
      }

      some_file << std::endl;
    }
  }

  some_file.close();
}

template <class ELEMENT>
void FlowAroundDiskProblem<ELEMENT>::validate_exact_solution_divergence() const
{
  ostringstream filename;

  filename << Doc_info.directory() << "exact_solution_divergence.txt";

  ofstream some_file;
  some_file.open(filename.str().c_str());

  // number of sample points in each direction
  const unsigned N = 100;

  const double w = Global_Parameters::Box_half_width;
  const double h = Global_Parameters::Box_half_height;
  
  for(unsigned xi=0; xi<N; xi++)
  {
    for(unsigned zi=0; zi<N; zi++)
    {
      // Cartesian coordinates in x-z plane
      Vector<double> x(3, 0.0);

      x[0] = -w + 2 * w * (double(xi)) / (N-1);
      x[2] = -h + 2 * h * (double(zi)) / (N-1);
      
      DenseMatrix<double> du_dx(3, 3, 0.0);
  
      // get the exact velocity derivatives
      FlatDiskExactSolutions::total_exact_velocity_gradient(x,
							    Global_Parameters::u_disk_rigid_translation,
							    Global_Parameters::omega_disk,
							    du_dx);

      // compute the divergence at this point
      double divergence = 0.0;      
      for(unsigned j=0; j<3; j++)
	divergence += du_dx(j,j);

      // output
      some_file << x[0] << " " << x[1] << " " << x[2] << " " << divergence << std::endl;
    }
  }
  
  some_file.close();
}

typedef struct {
    SuperMatrix *L;
    SuperMatrix *U;
    int *perm_c;
    int *perm_r;
} factors_t;

template <class ELEMENT>
void FlowAroundDiskProblem<ELEMENT>::compute_and_assign_smallest_eigensolution(
  const double& threshold)
{
  oomph_info << "Doing eigenshizzle...\n" << std::endl;
  
   // output the jacobian if it's singular
      
  // residual vector and Jacobian matrix
  DoubleVector r;
  CRDoubleMatrix jac;

  get_jacobian(r, jac);

  char filename[500];
  ofstream some_file;
      
  sprintf(filename,"%s/singular_jacobian_sparse%i.dat", doc_info().directory().c_str(),
	  doc_info().number());
  some_file.open(filename);

  // allows for size detection in python
  bool output_bottom_right = true;
  unsigned precision = 0;
  jac.sparse_indexed_output(some_file, precision, output_bottom_right);

  some_file.close();

  // Get eigenvalues/vectors?
  //=========================
  bool do_eigenvalues = true;
  unsigned n = r.nrow();
  if (do_eigenvalues)
  {
    DenseComplexMatrix DenseA(n);
    DenseComplexMatrix DenseM(n);
    for (unsigned i=0; i<n; i++)
    {
      DenseM(i,i) = complex<double>(1.0, 0.0);
      for (unsigned j=0; j<n; j++)
      {
	DenseA(i,j) = complex<double>(jac(i,j), 0.0);
      }
    }
    // DenseA.output(std::cout);
      
      
    // Make eigensolver
    LAPACK_QZ eigen_solver;
      
    // Storage
    Vector<std::complex<double> >eval(n);
    Vector<Vector<std::complex<double> > > evec(n);
    DoubleVector singular_vector;
      
    // This is a little hack to resize the vector without needing
    // to figure out how this annoying distribution pointer thing works.
    get_dofs(singular_vector);
    singular_vector.initialise(0.0);
      
    // Do it
    eigen_solver.find_eigenvalues(DenseA, DenseM, eval, evec);

    unsigned nzero_eigenvalues = 0;
    
    for(unsigned i=0; i<n; i++)
    {
      std::cout << "Eigenvalue " << i << " : " << eval[i] << "\n";
      if (fabs( real(eval[i]) ) < threshold)
      {
	nzero_eigenvalues++;
	
	for(unsigned j=0; j<n; j++)
	{
	  singular_vector[j] = real( evec[i][j] );
	  std::cout << evec[i][j] << ", ";
	}
      }
      // std::cout << std::endl;
    }
    sprintf(filename, "%s/sing_eigen_vect.dat", doc_info().directory().c_str());
    some_file.open(filename);
	
    for (unsigned i=0; i<n; i++)
    {
      oomph_info << "Singular eigenvector " << i
		 << " " << singular_vector[i] << std::endl;
      some_file << singular_vector[i] << std::endl;
    }
    some_file.close();
    assign_eigenvector_to_dofs(singular_vector);
    doc_solution(5); 

    if(nzero_eigenvalues == 0)
    {
      oomph_info << "Didn't find any eigenvalues below "
		 << threshold <<"!\n" << std::endl;
    }
    else
    {
      oomph_info << "Number of eigenvalues < " << threshold
		 << ": " << nzero_eigenvalues << std::endl;
    }
    
    oomph_info << "Done eigenshizzle\n" << std::endl;
  }
}

template <class ELEMENT>
double FlowAroundDiskProblem<ELEMENT>::compute_jacobian_condition_number()
{
  // residual vector and Jacobian matrix
  DoubleVector r;
  CRDoubleMatrix jac;

  // get 'em
  get_jacobian(r, jac);

  // convert the CRDoubleMatrix Jacobian to a SuperMatrix for SuperLU
  int n = jac.nrow();

  //Number of non-zero entries in the matrix
  int nnz           = jac.nnz();
  double *values    = jac.value();
  int* column_index = jac.column_index();
  int* row_start    = jac.row_start(); 
  
  SuperMatrix A;
  dCreate_CompRow_Matrix(&A, n, n, nnz, values, column_index, row_start,
			 SLU_NC, SLU_D, SLU_GE);
  
  // get the norm of the Jacobian
  double jac_norm = jac.inf_norm();

  oomph_info << "Jacobian inf-norm: " << jac_norm << std::endl;
  
  // create a SuperLU solver 
  std::unique_ptr<SuperLUSolver> linear_solver_pt = std::make_unique<SuperLUSolver>();
   
  // do the LU factorisation
  linear_solver_pt->factorise(&jac);

  // retrive the LU factors

  factors_t* f_factors = (factors_t*) linear_solver_pt->serial_f_factors();

  int info = 0;
  SuperLUStat_t stat;
  
  /* Initialize the statistics variables. */
  StatInit(&stat);
  
  char* norm = (char*)"1";
  
  // get the condition number
  double reciprocal_condition_number = 0.0;
  dgscon(norm, f_factors->L, f_factors->U, jac_norm, &reciprocal_condition_number, &stat, &info);

  return 1.0/reciprocal_condition_number;
}

template <class ELEMENT>
void FlowAroundDiskProblem<ELEMENT>::get_total_disk_force_and_moment_residuals(
  const Vector<double>& centre_of_rotation,
  Vector<double>& residuals)
{
  // save the iteration count, since the blackbox solver is used again during
  // the fluid solve to compute dudx
  unsigned fsi_iter_count = BlackBoxFDNewtonSolver::N_iter_taken;
  
  // and reset the tolerance temporarily
  BlackBoxFDNewtonSolver::Tol = 1e-8;
  BlackBoxFDNewtonSolver::FD_step = 1e-8;
  
  // Integral of forces and torques over the disk
  Vector<double> total_force(3, 0.0);
  Vector<double> total_moment(3, 0.0);

  // list of the non-augmented disk boundary IDs
  Vector<unsigned> non_aug_disk_boundary_ids;
  
  // add the lower (One-based!) boundary IDs
  for(unsigned ibound : One_based_boundary_id_for_disk_outside_torus)
    non_aug_disk_boundary_ids.push_back(ibound - 1);

  // add the upper (zero-based) boundary IDs
  for(unsigned ibound : Boundary_id_for_upper_disk_outside_torus)
    non_aug_disk_boundary_ids.push_back(ibound);
  
  // Part 1/2: loop over the elments outside the torus
  for(unsigned ibound : non_aug_disk_boundary_ids)  
  {    
    unsigned nel = Bulk_mesh_pt->nboundary_element(ibound);
    for(unsigned e=0; e<nel; e++)
    { 
      ELEMENT* el_pt =
  	dynamic_cast<ELEMENT*>(Bulk_mesh_pt->boundary_element_pt(ibound, e));

      // What is the index of the face of the bulk element at the boundary
      int face_index = Bulk_mesh_pt->face_index_at_boundary(ibound, e);

      // Build the corresponding face element
      // (N.B. we're going for the WithSingularity version despite this being
      // the non-augmented region simply to reuse the machinery for
      // integrating forces and torques)
      auto surface_element_pt =
      	std::make_unique<NavierStokesWithSingularityFaceElement<ELEMENT>>
	(el_pt, face_index);

      // integrate the forces and moments on the element
      Vector<double> total_force_on_element(Dim, 0.0);
      Vector<double> total_moment_on_element(Dim, 0.0);
      
      surface_element_pt->integrated_force_and_moment(centre_of_rotation,
						      total_force_on_element,
						      total_moment_on_element);
      // add to totals
      for(unsigned i=0; i<3; i++)
      {
	// QUEHACERES experimental, changing to -ve as the
	// normal direction will be wrong
	total_force[i]  -= total_force_on_element[i];
	total_moment[i] -= total_moment_on_element[i];
      }
    }
  }
 
  // Part 2/2: loop over the BC face elements, which sit on the upper and lower
  //           disk within the torus region
  for(unsigned e=0; e<Face_mesh_for_bc_pt->nelement(); e++)
  {
    // get a pointer to the face element
    NavierStokesWithSingularityBCFaceElement<ELEMENT>* surface_element_pt =
      dynamic_cast<NavierStokesWithSingularityBCFaceElement<ELEMENT>*>(
  	Face_mesh_for_bc_pt->element_pt(e));

    // integrate the forces and moments on the element
    Vector<double> total_force_on_element(Dim, 0.0);
    Vector<double> total_moment_on_element(Dim, 0.0);
    
    surface_element_pt->integrated_force_and_moment(centre_of_rotation,
						    total_force_on_element,
						    total_moment_on_element);
    // add to totals
    for(unsigned i=0; i<3; i++)
    {
      total_force[i]  += total_force_on_element[i];
      total_moment[i] += total_moment_on_element[i];
    }    
  }

  total_force[2] -= Global_Parameters::Mass * sqrt(Global_Parameters::Archimedes_number);
    
  // // get the body force 
  // Vector<double> x_dummy(3, 0.0);
  // Vector<double> body_force(3, 0.0);
  // Global_Parameters::reduced_gravity(x_dummy, body_force);
  
   // for(unsigned i=0; i<3; i++)
   //   total_force[i] -= body_force[i] * 2 * MathematicalConstants::Pi;
  
  // fill in residuals - for equilibrum we require the sum of the forces to be zero
  // and the sum of the torques to be zero, so these form our 6 residuals
  // QUEHACERES changing to just 3 translations
  residuals.resize(1, 0.0); // 6

  residuals[0] = total_force[2];
  
  // for(unsigned i=0; i<3; i++)
  // {
  //   residuals[i]   = total_force[i];
  //   // residuals[3+i] = total_moment[i]; 
  // }

  // reset
  BlackBoxFDNewtonSolver::N_iter_taken = fsi_iter_count;
  BlackBoxFDNewtonSolver::Tol = Global_Parameters::FSI_convergence_tol;
  BlackBoxFDNewtonSolver::FD_step = Global_Parameters::FSI_FD_step;
  
  oomph_info << "FSI iteration " << BlackBoxFDNewtonSolver::N_iter_taken
	     << "\n-----------------\n";

  for(double res : residuals)
    oomph_info << res << "\n";

  oomph_info << std::endl;
}


void hacky_fixed_point_fsi_residual(
  const Vector<double>& centre_of_rotation,
  const Vector<double>& translation_and_rotation_rates,
  Vector<double>& residuals)
{
  // update the disks rigid body parameters with the current guess
  // Global_Parameters::u_disk_rigid_translation[0] = translation_and_rotation_rates[0];
  // Global_Parameters::u_disk_rigid_translation[1] = translation_and_rotation_rates[1];
  Global_Parameters::u_disk_rigid_translation[2] = translation_and_rotation_rates[0];

  // QUEHACERES
  // Global_Parameters::omega_disk[0] = translation_and_rotation_rates[3];
  // Global_Parameters::omega_disk[1] = translation_and_rotation_rates[4];
  // Global_Parameters::omega_disk[2] = translation_and_rotation_rates[5];

  typedef ProjectableTaylorHoodElement<
    TNavierStokesWithSingularityPdeConstrainedMinElement<3> > ELEMENT;
  
  FlowAroundDiskProblem <ELEMENT>* flow_problem =
    dynamic_cast<FlowAroundDiskProblem<ELEMENT>*>(Global_Parameters::Problem);
  
  // apply new BCs, i.e. updated place velocities above
  flow_problem->apply_boundary_conditions();

  // save the iteration count, since the blackbox solver is used again during
  // the fluid solve to compute dudx
  unsigned fsi_iter_count = BlackBoxFDNewtonSolver::N_iter_taken;

  // and reset the tolerance temporarily
  BlackBoxFDNewtonSolver::Tol = 1e-8;
  BlackBoxFDNewtonSolver::FD_step = 1e-8;
  
  // do the solve  
  flow_problem->newton_solve();

  BlackBoxFDNewtonSolver::N_iter_taken = fsi_iter_count;
  BlackBoxFDNewtonSolver::Tol = Global_Parameters::FSI_convergence_tol;
  BlackBoxFDNewtonSolver::FD_step = Global_Parameters::FSI_FD_step;
  
  // get the residuals
  flow_problem->
    get_total_disk_force_and_moment_residuals(centre_of_rotation, residuals);  
}

template <class ELEMENT>
void FlowAroundDiskProblem<ELEMENT>::hacky_fixed_point_fsi_solve()
{
  oomph_info << "\nDoing FSI fixed-point solve (tol = "
	     << Global_Parameters::FSI_convergence_tol 
	     << ")...\n\n" << std::endl;

  // QUEHACERES make things faster for debug
  BlackBoxFDNewtonSolver::Tol = Global_Parameters::FSI_convergence_tol;

  // avoids issues with noise
  BlackBoxFDNewtonSolver::FD_step = Global_Parameters::FSI_FD_step;
  
  // vector of unknowns - 3 linear translational velocities and
  // 3 rotation rates which characterise full 3D rigid body motion
  // [0] = u_x
  // [1] = u_y
  // [2] = u_z
  // [3] = omega_x
  // [4] = omega_y
  // [5] = omega_z
  Vector<double> unknown_translation_and_rotation_rates(1, 0.0); // 6

  // initial guess
  unknown_translation_and_rotation_rates[0] = Global_Parameters::u_disk_rigid_translation[2];
    
  // Centre of rotation is the parameter for this rigid problem
  Vector<double> centre_of_rotation(3, 0.0);

  // solve the bastard!
  try
  {    
    BlackBoxFDNewtonSolver::black_box_fd_newton_solve(
      &hacky_fixed_point_fsi_residual,
      centre_of_rotation, unknown_translation_and_rotation_rates);

    oomph_info << "\n===========================================\n"
	       << "FSI fixed-point iteration scheme converged\n"
	       << "===========================================\n\n"
	       // << "u_x     = " << unknown_translation_and_rotation_rates[0] << "\n"
	       // << "u_y     = " << unknown_translation_and_rotation_rates[1] << "\n"
	       << "u_z     = " << unknown_translation_and_rotation_rates[0] << "\n"
	       // << "omega_x = " << unknown_translation_and_rotation_rates[3] << "\n"
	       // << "omega_y = " << unknown_translation_and_rotation_rates[4] << "\n"
	       // << "omega_z = " << unknown_translation_and_rotation_rates[5] << "\n"
	       << "===========================================\n" << std::endl;
  }
  catch (OomphLibError& error)
  {
    throw OomphLibError("FSI fixed-point iteration scheme didn't converge\n",
			OOMPH_CURRENT_FUNCTION,
			OOMPH_EXCEPTION_LOCATION);
  }
}

template <class ELEMENT>
void FlowAroundDiskProblem<ELEMENT>::fsi_residual_sweep()
{
  const unsigned n=10;

  Vector<double> u_z;
  Vector<double> res;
  
  for(unsigned i=0; i<n; i++)
  {
    // set u_z
    double uz_i = -(double)(i) / (n-1);

    u_z.push_back(uz_i);
    
    Global_Parameters::u_disk_rigid_translation[2] = uz_i;
    
    apply_boundary_conditions();
    
    newton_solve();

    Vector<double> centre_of_rotation_dummy(3, 0.0);
    Vector<double> residuals;
    get_total_disk_force_and_moment_residuals(centre_of_rotation_dummy,
					      residuals);

    res.push_back(residuals[0]);
  }

  oomph_info << "Uz residuals\n" << "--------------" << std::endl;
  for(unsigned i=0; i<n; i++)
  {
    oomph_info << u_z[i] << " " << res[i] << std::endl;
  }
  oomph_info << "--------------" << std::endl;
}

//========================================================================
/// Driver
//========================================================================
int main(int argc, char* argv[])
{
  // set up the multi-processor interface
  MPI_Helpers::init(argc,argv);
  
  oomph_info << "\n\n=======================   Configuration:   =======================\n"  
	     << "=                                                                =\n";
#ifdef USE_FD_JACOBIAN  
  oomph_info << "= - Using finite-diff jacobian                                   =\n";
#else
  oomph_info << "= - Using analytic jacobian                                      =\n";
#endif
  oomph_info << "= - Using formulation with symmetric Jacobian matrix             =\n";
  oomph_info << "= - All elements are PDE-constrained                             =\n"; 
  oomph_info << "= - Using multiple spatially-varying DoFs for singular amplitude =\n";
  oomph_info << "=    with SINGULAR_ELEMENT_NNODE_1D = "
	     << SINGULAR_ELEMENT_NNODE_1D <<           "                          =\n";
  oomph_info << "=                                                                =\n"
	     << "==================================================================\n"
	     << std::endl;  
  
  // keep track of total program runtime
  double t_start = TimingHelpers::timer();
  
  // Store command line arguments
  CommandLineArgs::setup(argc, argv);

  // ==========================================================================
  // Physical problem parameters
  // ==========================================================================

  CommandLineArgs::specify_command_line_flag("--no_slip_boundaries");

  // do a fixed-point iteration FSI solve, using any velocities
  // / rotation rates specified below as initial conditions
  CommandLineArgs::specify_command_line_flag("--rigid_fsi_solve");
  
  CommandLineArgs::specify_command_line_flag("--fsi_convergence_tol",
					     &Global_Parameters::FSI_convergence_tol);

  CommandLineArgs::specify_command_line_flag("--fsi_residual_sweep");
    
  // rigid body velocity of the plate
  CommandLineArgs::specify_command_line_flag("--velocity_x",
					     &Global_Parameters::u_disk_rigid_translation[0]);
  CommandLineArgs::specify_command_line_flag("--velocity_y",
					     &Global_Parameters::u_disk_rigid_translation[1]);
  CommandLineArgs::specify_command_line_flag("--velocity_z",
					     &Global_Parameters::u_disk_rigid_translation[2]);

   // rigid body angular velocities of the plate
  CommandLineArgs::specify_command_line_flag("--angular_velocity_x",
					     &Global_Parameters::omega_disk[0]);
  CommandLineArgs::specify_command_line_flag("--angular_velocity_y",
					     &Global_Parameters::omega_disk[1]);
  CommandLineArgs::specify_command_line_flag("--angular_velocity_z",
					     &Global_Parameters::omega_disk[2]);

  CommandLineArgs::specify_command_line_flag("--p0",
					     &Global_Parameters::p0);
  
  // half width of the container box (disk radius is 1)
  CommandLineArgs::specify_command_line_flag("--box_half_width",
					     &Global_Parameters::Box_half_width);

  // half length of the container box
  CommandLineArgs::specify_command_line_flag("--box_half_height",
					     &Global_Parameters::Box_half_height);
  
  // amplitude of the warping of the disk
  CommandLineArgs::specify_command_line_flag("--epsilon",
					     &Global_Parameters::Epsilon);

  // wavenumber of the warping of the disk
  CommandLineArgs::specify_command_line_flag("--n",
					     &Global_Parameters::n);

  // Non-dimensional mass of the disk
  CommandLineArgs::specify_command_line_flag("--mass", &Global_Parameters::Mass);

  // Archimedes number 
  CommandLineArgs::specify_command_line_flag("--Ar", &Global_Parameters::Archimedes_number);
  
  // specify the relative difference in densities between the solid and fluid
  CommandLineArgs::specify_command_line_flag("--relative_density_difference",
					     &Global_Parameters::Relative_density_difference);
  
  // ==========================================================================
  // Singularity stuff and pinning
  // ==========================================================================
  
  // just do pure FE
  CommandLineArgs::specify_command_line_flag("--dont_subtract_singularity"); 

  // use the singular function subtraction machinery, but set the amplitudes
  // to zero (for debug)
  CommandLineArgs::specify_command_line_flag("--impose_zero_singular_amplitude");

  // not specifying this does the real problem where c is computed
  CommandLineArgs::specify_command_line_flag("--impose_exact_singular_amplitude");

  // pin different singular functions to zero, for debug
  CommandLineArgs::specify_command_line_flag("--pin_broadside_amplitude",
					     &Global_Parameters::broadside_amplitude_pin);  
  CommandLineArgs::specify_command_line_flag("--pin_broadside_rotation_amplitude",
					     &Global_Parameters::broadside_rotation_amplitude_pin);
  CommandLineArgs::specify_command_line_flag("--pin_in_plane_amplitude",
					     &Global_Parameters::in_plane_amplitude_pin);
  CommandLineArgs::specify_command_line_flag("--pin_in_plane_rotation_amplitude",
					     &Global_Parameters::in_plane_rotation_amplitude_pin);
  
  // subtract the exact solution instead of the Moffatt solution
  CommandLineArgs::specify_command_line_flag("--subtract_exact_solution");

  // subtract the source and body force terms which arise from the
  // (2 term) asymptotic expansion of the exact solution not exactly satisfying the
  // Stokes equations
  CommandLineArgs::specify_command_line_flag("--subtract_source_and_body_force");

  // Only subtract the leading order asymptotic term as the singular function;
  // default is to subtract two asymptotic terms for the singular functions
  CommandLineArgs::specify_command_line_flag("--only_subtract_first_asymptotic_term");
    
  // Value of the fake singular amplitude of the broadside mode
  // to impose on the solution for debug
  CommandLineArgs::specify_command_line_flag(
    "--set_sing_amplitude_broadside",
    &Global_Parameters::Imposed_singular_amplitude_broadside);
  
  // Value of the fake singular amplitude of the in-plane mode
  // to impose on the solution for debug
  CommandLineArgs::specify_command_line_flag(
    "--set_sing_amplitude_in_plane",    
    &Global_Parameters::Imposed_singular_amplitude_in_plane);

  // QUEHACERES for debug, initial residual should be zero
  CommandLineArgs::specify_command_line_flag(
    "--set_initial_conditions_to_non_singular_solution");
    
  // Modulate the singular amplitude by cos(zeta)
  CommandLineArgs::specify_command_line_flag("--cosine_amplitude");

  // specify the penalty parameter applied to the pressure term in the
  // functional which is minimised  
  CommandLineArgs::specify_command_line_flag("--pressure_regularisation_factor",
					     &Global_Parameters::Pressure_regularisation_factor);

  CommandLineArgs::specify_command_line_flag("--velocity_regularisation_factor",
					     &Global_Parameters::Velocity_regularisation_factor);
  
  // ==========================================================================
  // Meshing parameters
  // ==========================================================================
  
  // Cross-sectional radius of the torus region
  CommandLineArgs::specify_command_line_flag("--r_torus",
					     &Global_Parameters::R_torus);

  // number of vertices on the cross-sectional circles of the torus
  CommandLineArgs::specify_command_line_flag("--nvertex_torus",
					     &Global_Parameters::Nvertex_torus);

  // Half the number of segments making up the perimeter of the disk
  CommandLineArgs::specify_command_line_flag("--half_nsegment_disk",
					     &Global_Parameters::Half_nsegment_disk);

  CommandLineArgs::specify_command_line_flag("--nsingular_line_element",
					     &Global_Parameters::Nsingular_line_element);

  CommandLineArgs::specify_command_line_flag("--global_element_volume",
					     &Global_Parameters::Initial_global_element_volume);
    
  CommandLineArgs::specify_command_line_flag("--target_element_volume_in_torus", 
					     &Global_Parameters::Target_element_volume_in_torus_region);
  
  // prevent the splitting of corner elements
  CommandLineArgs::specify_command_line_flag("--dont_split_corner_elements");

  // ==========================================================================
  // Output parameters
  // ==========================================================================
  
  // get output directory
  CommandLineArgs::specify_command_line_flag("--dir", &Global_Parameters::output_directory);

    // suppress the (expensive!) bulk element output
  CommandLineArgs::specify_command_line_flag("--suppress_bulk_output");

  // don't output the exact solution
  CommandLineArgs::specify_command_line_flag("--dont_output_exact_solution");
  
  // do we want to output the disk-on-disk azimuthal plots? (also expensive)
  CommandLineArgs::specify_command_line_flag("--disk_on_disk_plots");
  
  // extra output to describe equation numbers
  CommandLineArgs::specify_command_line_flag("--describe_dofs");
  CommandLineArgs::specify_command_line_flag("--describe_nodes");
  
  // number of plot points per side in the bulk elements
  CommandLineArgs::specify_command_line_flag("--nplot", &Global_Parameters::Nplot_for_bulk);

  // set all the torus region elements to the singular solution (very slow!!)
  CommandLineArgs::specify_command_line_flag(
    "--set_initial_conditions_to_singular_solution");

  // doc the solution before any solves
  CommandLineArgs::specify_command_line_flag("--doc_initial_conditions");

  // separately output the solution in the augmented region
  CommandLineArgs::specify_command_line_flag("--output_soln_in_torus");


  // ==========================================================================
  // Debug options
  // ==========================================================================
  
  CommandLineArgs::specify_command_line_flag("--validate_singular_stress_broadside");
  CommandLineArgs::specify_command_line_flag("--validate_singular_stress_in_plane");

  CommandLineArgs::specify_command_line_flag("--output_jacobian_full");
  CommandLineArgs::specify_command_line_flag("--output_jacobian_sparse");
  CommandLineArgs::specify_command_line_flag("--output_initial_jacobian");
  
  CommandLineArgs::specify_command_line_flag("--validate_exact_solution_divergence");

  CommandLineArgs::specify_command_line_flag("--compute_jacobian_condition_number");

  double threshold_for_zero_eigenval = 1e-12;
  CommandLineArgs::specify_command_line_flag("--do_eigenshizzle", &threshold_for_zero_eigenval);
  
  CommandLineArgs::specify_command_line_flag("--use_fd_lu_solver");

  // ==========================================================================
  // Solvers
  // ==========================================================================
  CommandLineArgs::specify_command_line_flag("--use_mumps_solver");
  CommandLineArgs::specify_command_line_flag("--use_hypre_solver");
  CommandLineArgs::specify_command_line_flag("--use_gmres_solver");
  CommandLineArgs::specify_command_line_flag("--use_bicgstab_solver");
  
  // ==========================================================================
  // **************************************************************************
  // ==========================================================================
  
  // // Note that this can make tetgen die!
  // feenableexcept(FE_INVALID | FE_DIVBYZERO | FE_OVERFLOW | FE_UNDERFLOW);

  // Shut up prefix
  oomph_info.output_modifier_pt() = &default_output_modifier;

  // Parse command line
  CommandLineArgs::parse_and_assign();
  
  // Doc what has actually been specified on the command line
  CommandLineArgs::doc_specified_flags();

  oomph_info << "\nPressure regularisation factor: "
	     << Global_Parameters::Pressure_regularisation_factor << "\n" << std::endl;
  
  if (CommandLineArgs::command_line_flag_has_been_set("--disk_on_disk_plots"))
    Global_Parameters::Do_disk_on_disk_plots = true;
    
  if (CommandLineArgs::command_line_flag_has_been_set("--dont_split_corner_elements"))
  {
    Global_Parameters::Split_corner_elements = false;
  }
    
  if(CommandLineArgs::command_line_flag_has_been_set("--use_fd_lu_solver"))
  {
    Global_Parameters::Use_fd_lu_solver = true;
    oomph_info << "Get comfy, we're using the FD_LU solver...\n" << std::endl;
  }

  if(CommandLineArgs::command_line_flag_has_been_set("--no_slip_boundaries"))
  {
    Global_Parameters::No_slip_boundaries = true;
  }
  
  // FlowAroundDiskProblem <ProjectableTaylorHoodElement<
  //   TNavierStokesElementWithSingularity<3,3> > > problem; // TTaylorHoodElement<3>::NNODE_1D>
  FlowAroundDiskProblem <ProjectableTaylorHoodElement<
			   TNavierStokesWithSingularityPdeConstrainedMinElement<3> > > problem;

  EdgeCoordinates edge_coords(1e-3, 0, -MathematicalConstants::Pi);

  DenseMatrix<double> grad_sing = Global_Parameters::SingularFunctions::
    __gradient_of_singular_fct_exact_broadside_translation(edge_coords);

  for(unsigned i=0; i<3; i++)
  {
    for(unsigned j=0; j<3; j++)
    {
      oomph_info << grad_sing(i,j) << " ";
    }
    oomph_info << std::endl;
  }

  // return 0;
  
  // hacky, set the global problem pointer
  Global_Parameters::Problem = &problem;
  
  if(!CommandLineArgs::command_line_flag_has_been_set("--dont_subtract_singularity") &&
     (CommandLineArgs::command_line_flag_has_been_set("--impose_exact_singular_amplitude") ||
      CommandLineArgs::command_line_flag_has_been_set("--impose_zero_singular_amplitude") ||
      CommandLineArgs::command_line_flag_has_been_set("--set_sing_amplitude_broadside") ||
      CommandLineArgs::command_line_flag_has_been_set("--set_sing_amplitude_in_plane")))
  {
    bool impose_zero_amplitude =
      CommandLineArgs::command_line_flag_has_been_set("--impose_zero_singular_amplitude");
    
    problem.impose_fake_singular_amplitude(impose_zero_amplitude);
  }
      
  if(CommandLineArgs::command_line_flag_has_been_set("--only_subtract_first_asymptotic_term"))
  {
    Global_Parameters::Only_subtract_first_singular_term = true;
  }

  // for debug
  if (CommandLineArgs::command_line_flag_has_been_set(
  	"--set_initial_conditions_to_singular_solution") )    
  {
    oomph_info << "Setting initial conditions to singular solution..." << std::endl;
    
    problem.set_values_to_singular_solution();
  }

  if (CommandLineArgs::command_line_flag_has_been_set(
	"--set_initial_conditions_to_non_singular_solution") )    
  {
    problem.set_values_to_exact_non_singular_solution();
  }
  
  if (CommandLineArgs::command_line_flag_has_been_set(
	"--validate_singular_stress_broadside") )
  {
    problem.validate_singular_stress();
    problem.doc_solution(2);
    exit(0);
  }
  else if(CommandLineArgs::command_line_flag_has_been_set(
	"--validate_singular_stress_in_plane") )
  {
    problem.validate_singular_stress(false);
    problem.doc_solution(2);
    exit(0);
  }

  if(CommandLineArgs::command_line_flag_has_been_set("--validate_exact_solution_divergence"))
  {
    oomph_info << "\n\nValidating divergence of exact solution...\n" << std::endl;
    problem.validate_exact_solution_divergence();
    exit(0);
  }

  // for debug - output the pin-status of each dof in each submesh
  problem.output_submesh_dof_pin_status();

  // ==========================================================================
  // Setup equation numbering scheme - no more fiddling around / pinning after this!
  // ==========================================================================
  oomph_info << "--------------------------\n";
  oomph_info << "Number of equations: " << problem.assign_eqn_numbers() << std::endl; 
  oomph_info << "--------------------------\n" << std::endl;

  if(CommandLineArgs::command_line_flag_has_been_set("--compute_jacobian_condition_number"))
  {
    double condition_number = problem.compute_jacobian_condition_number();

    oomph_info << "Jacobian condition number: " << condition_number << std::endl;
    exit(0);
  }

  if(CommandLineArgs::command_line_flag_has_been_set("--do_eigenshizzle"))
  {
    problem.compute_and_assign_smallest_eigensolution(threshold_for_zero_eigenval);

    exit(0);
  }
  
  // QUEHACERES for debug
  problem.newton_solver_tolerance() = 5e-8;

  problem.max_residuals() = 100;
  
  // Number of output points per edge
  unsigned nplot = Global_Parameters::Nplot_for_bulk;

  if(CommandLineArgs::command_line_flag_has_been_set("--doc_initial_conditions"))
  {
    //Output initial guess
    problem.doc_solution(nplot);
  } 
  else
  {
    problem.doc_info().number()++;
  }

  if(CommandLineArgs::command_line_flag_has_been_set("--fsi_residual_sweep"))
  {
    problem.fsi_residual_sweep();

    return 0;
  }
  else if(CommandLineArgs::command_line_flag_has_been_set("--rigid_fsi_solve"))
  {
    problem.hacky_fixed_point_fsi_solve();
    problem.doc_solution(nplot);

    return 0;
  }
  
  unsigned max_adapt = 0; 
  for (unsigned i=0; i<=max_adapt; i++)
  {
#ifndef COMPUTE_INITIAL_JACOBIAN_EIGENFUNCTIONS
    
#ifdef PRINT_SINGULAR_JACOBIAN
    try
    {
#endif
      // QUEHACERES debug
      if (CommandLineArgs::command_line_flag_has_been_set("--output_initial_jacobian"))
      {
      	oomph_info << "outputting jacobian before solve...\n" << std::endl;
	
      	// residual vector and Jacobian matrix
      	DoubleVector r;
      	CRDoubleMatrix jac;

      	problem.get_jacobian(r,jac);
 
      	char filename[500];
      	ofstream some_file;
      
      	sprintf(filename,"%s/jacobian_sparse%i.dat", problem.doc_info().directory().c_str(),
      		problem.doc_info().number());
      	some_file.open(filename);
      
      	jac.sparse_indexed_output(some_file);

      	some_file.close();
      }
      
      // Solve the bastard!
      problem.newton_solve();

#endif
#ifdef PRINT_SINGULAR_JACOBIAN
#ifndef COMPUTE_INITIAL_JACOBIAN_EIGENFUNCTIONS
    }
    catch (const std::exception& ex)
    {
#endif	
      // QUEHACERES
      problem.compute_and_assign_smallest_eigensolution();
      
#ifndef COMPUTE_INITIAL_JACOBIAN_EIGENFUNCTIONS    
    }
#endif
#endif // end PRINT_SINGULAR_JACOBIAN
    
    //Output solution
    problem.doc_solution(nplot);

    if (i != max_adapt)
    {
      problem.adapt();
    }
  }

  oomph_info << "Done, total runtime: " << TimingHelpers::timer() - t_start << "s\n\n";
  
  // Shut down oomph-lib's MPI
  MPI_Helpers::finalize(); 
  
  return 0;
}

  
