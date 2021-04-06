#ifndef OOMPH_FLOW_AROUND_DISK_PROBLEM_HEADER
#define OOMPH_FLOW_AROUND_DISK_PROBLEM_HEADER

// NNODE_1D for the singular line elements - may affect stability of the scheme
#define SINGULAR_ELEMENT_NNODE_1D 3

//Generic routines
#include "generic.h"

// include classes for vector and matrix algebra
#include "additional_maths.h"

// singular elements
#include "navier_stokes_sing_face_element.h"


// the cylindrically warped circular disk GeomObject
#include "warped_disk_with_torus_geom_obj.h"

// singular elements
#include "navier_stokes_sing_face_element.h"

// the cylindrically warped circular disk GeomObject
#include "warped_disk_with_torus_geom_obj.h"

// definitions to convert CForm output from Mathematica
#include "mathematica_definitions.h"

// Exact solutions for the 4 distinct modes of rigid body motion
// of a flat disk
#include "exact_solutions_finite_disk.h"

// the unscaled singular functions which are subtracted in the augmented region
#include "singular_functions.h"

#include "chebyshev_gauss_integration.h"
#include "meshes/annular_mesh.h"

// functions to convert coordinates, velocities and gradients between
// the Lagrangian edge coordinates and the global Cartesian system
#include "coordinate_conversions.h"

// 2D disk elements, which are really acting as a placeholder for
// proper elastic shell elements
#include "rigid_disk_elements.h"

// The mesh
#include "meshes/triangle_mesh.h"
 
// Get the mesh
#include "meshes/tetgen_mesh.h" 
#include "meshes/refineable_tetgen_mesh.h"
#include "meshes/gmsh_tet_mesh.h"

// Get the faceted surfaces
#include "tetmesh_faceted_surfaces.h"

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

  // radius of the cylinder around which the disk is wrapped.
  // Default is ~inf, i.e. flat disk
  double Radius_of_curvature = 1e3;
    
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

  // number of rectangular elements in the xi1 ('radial') and xi2 ('azimuthal')
  // directions for computing the singular contribution to the drag
  unsigned Ndrag_element_xi1 = 5;
  unsigned Ndrag_element_xi2 = 20;
  
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
  CylindricallyWarpedCircularDiskWithAnnularInternalBoundary*
  Warped_disk_with_boundary_pt = nullptr;

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

  Vector<double> compute_singular_amplitudes_from_disk_velocity(
    const double& zeta)
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
    double c_broadside = u_broadside + // 2 * omega_out_of_plane * sin(zeta - zeta_rotation);
      2.0 * omega_out_of_plane * (cos(zeta) + sin(zeta));

    double c_broadside_rotation = -omega_out_of_plane * sin(zeta) / sqrt(2);
      
    // in-plane amplitude is the in-plane speed modulated by
    // an in-phase function of the angle of this point relative to the
    // translation vector
    double c_in_plane = u_in_plane * cos(zeta - zeta_rotation);
  
    // in-plane rotation amplitude is the in-plane speed modulated by
    // a function pi/2 out of phase with the angle of this point relative to the
    // translation vector
    double c_in_plane_rotation = u_in_plane * (sin(zeta - zeta_translation) -
					       cos(zeta - zeta_translation)) +
      omega_disk[2];

    // now package them up for return
    Vector<double> amplitudes(4, 0.0);

    // if we're doing the exact solution, then we don't want any modulation
    if(CommandLineArgs::command_line_flag_has_been_set("--subtract_exact_solution"))
    {
      c_broadside = u_broadside > 0.0 ? 1 : 0;
      c_in_plane = u_in_plane > 0.0 ? 1 : 0;
      c_in_plane_rotation = 0; // QUEHACERES not doing this for the time being
    }

    amplitudes[0] = c_broadside;
    amplitudes[1] = c_broadside_rotation;
    amplitudes[2] = c_in_plane;
    amplitudes[3] = c_in_plane_rotation;
  
    return amplitudes;
  }
}

namespace Analytic_Functions
{
  // derivative of the functional which we're minimising, w.r.t. the solution
  Vector<double> dfunctional_du(const Vector<double>& u)
  {
    Vector<double> dpi_du(u.size(), 0.0);

    // functional \Pi = sum_j 1/2|u_j|^2, so d\Pi/d u_i = u_i
    for(unsigned i=0; i<u.size(); i++)
      dpi_du[i] = (u[i]);

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

  // wrapper function to get the outer unit normal as a function
  // of the Lagrangian disk coordinates 
  void outer_unit_normal(const LagrangianCoordinates& lagr_coords,
			 Vector<double>& unit_normal)
  {
    // basis vectors (a3 is the outer unit normal)
    Vector<double> a1_dummy, a2_dummy;

    // get 'em from the disk GeomObject
    Global_Parameters::Warped_disk_with_boundary_pt->basis_vectors(lagr_coords,
								   a1_dummy,
								   a2_dummy,
								   unit_normal);

    // now set the direction based on the sign of xi3
    // (xi3 = -0.0 indicates a point on the lower side of the disk)
    for(double& ni : unit_normal)
      ni *= lagr_coords.sign_of_xi3();
  }
}

// ////////////////////////////////////////////////////////////////////////////
// ////////////////////////////////////////////////////////////////////////////
// ////////////////////////////////////////////////////////////////////////////


//=========================================================================
/// Class that solves Navier-Stokes flow around a 2D disk
//=========================================================================
template <class ELEMENT>
class FlowAroundDiskProblem : public Problem
{
  
public:
  
  static const unsigned NNODE_1D = ELEMENT::_NNODE_1D_;

  // shorthand for the disk elements which compute the singular drag. We'll
  // use the same integration order as the normal elements
  typedef QSingularDragIntegralDiskElement<2, NNODE_1D> SingularDragElement;
  
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
  void set_values_to_singular_solution() const;

  /// \short set the nodal values to the exact non-singular part of the solution
  void set_values_to_exact_non_singular_solution() const;
  
  // function to validate the singular stress by assigning the singular solution
  // to each node, then computing the velocity gradients via finite difference and
  // comparing these to the analytic gradients.
  void validate_singular_stress() const;

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
  /// for their knots (integration points). If use_undeformed_coords is true, then
  /// the Lagrangian coordinates are just the standard cylindrical coordinates of
  /// each point (rather than performing the conversion from Eulerian->Lagrangian)
  void setup_lagr_coordinates_and_singular_element(
    FiniteElement* const& elem_pt,
    Vector<LagrangianCoordinates>& lagr_coordinates_at_knot_pt,
    Vector<std::pair<GeomObject*, Vector<double> > >& line_element_and_local_coordinate,
    const bool& use_plot_points = false,
    const bool& use_undeformed_coords = false) const;

  /// \short function to compute the edge coordinates and corresponding
  /// element and local coordinates in the singular line mesh
  /// (which computes the appropriate amplitude) from Eulerian coordinates
  void get_lagr_coordinates_and_singular_element(
    FiniteElement* const& elem_pt,
    const Vector<double>& x,
    LagrangianCoordinates& lagr_coordinates_at_point,
    std::pair<GeomObject*, Vector<double> >& line_element_and_local_coordinate,
    const bool& use_undeformed_coords = false) const;

  /// \short check that we don't have any inconsistencies between the sign of the
  /// elevation angle phi and the direction of the unit normal for upper/lower surfaces
  Vector<std::pair<NavierStokesWithSingularityBCFaceElement<ELEMENT>*, unsigned>>
    sanity_check_lagr_coords_on_disk() const;    

  void output_augmented_elements_eulerian_and_lagr_coords() const;
  
  /// \short Helper to generate the 1D line mesh which runs around the outside of the
  /// disk and provides the amplitude of the singularity as a function of the
  /// boundary coordinate, i.e. c(zeta_bound).
  void create_one_d_singular_element_mesh(const unsigned& nsingular_el);

  /// \short create a mesh of quad elements with a mixed Gauss/Chebyshev-Gauss
  /// integration scheme for integrating the singular contribution to the drag
  void create_singular_drag_integration_mesh(const double& r_torus,
					     const unsigned& nel_xi1,
					     const unsigned& nel_xi2);
      
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

  /// \short Mesh of quad elements with mixed gauss/chebyshev-gauss integration 
  /// scheme for integrating the singular contribution to the drag
  Mesh* Singular_drag_integration_mesh_pt;
  Mesh* Singular_drag_integration_mesh_upper_pt;
  Mesh* Singular_drag_integration_mesh_lower_pt;

  // pointer to the GeomObject which characterises the disk
  CylindricallyWarpedCircularDiskWithAnnularInternalBoundary* Warped_disk_with_boundary_pt;
  
  /// \short Enumeration for IDs of FaceElements (used to figure out
  /// who's added what additional nodal data...)
  /// N.B. we pass the same ID for the stress-jump and BC IDs, since this is
  /// really just one continuous field and so at the intersection of the
  /// disk and torus boundary we don't want separate fields being interpolated
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

  // lists of the bulk elements and their corresponding face indices which touch
  // the upper and lower surfaces of the disk
  Vector<std::pair<ELEMENT*, unsigned>> Bulk_element_and_face_index_for_upper_disk_in_torus;
  Vector<std::pair<ELEMENT*, unsigned>> Bulk_element_and_face_index_for_lower_disk_in_torus;
  
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
  std::map<Node*, Vector<double> >* Node_to_lagr_coordinates_map_pt;

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

  bool Have_output_exact_soln;
  
  /// Sanity check: Exact bounded volume
  double Exact_bounded_volume;

  /// Number of singular functions we will subtract
  unsigned Nsingular_function;
  
  /// Number of dimensions in the problem
  unsigned Dim;
  
  DocInfo Doc_info;

  unsigned Fsi_iteration_counter;
};

#endif
