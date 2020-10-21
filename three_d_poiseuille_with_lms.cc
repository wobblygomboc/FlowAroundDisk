
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

#include<fenv.h>

// needed to reset the FPU control word which Triangle messes around with
#include <fpu_control.h>

//Generic routines
#include "generic.h"

// ============================================================================
// custom defines
// ============================================================================

#define USE_FD_JACOBIAN

#define PRINT_SINGULAR_JACOBIAN

// do we want to duplicate the edge nodes (for better output)?
#define xDUPLICATE_EDGE_NODES

// Tetgen or Gmsh
#define DO_TETGEN

// wraps the problem.newton_solve() in a try/catch block and prints jacobian
// (but cocks up output of other exceptions)
#define xPRINT_SINGULAR_JACOBIAN

// QUEHACERES for debug - shift the lower plate nodes down by a small amount
// so the the upper and lower parts of the plate are physically separated
#define SHIFT_LOWER_PLATE_NODES_BY_DZ

#ifdef SHIFT_LOWER_PLATE_NODES_BY_DZ
const double lower_plate_z_shift = -0.987e-8;
#endif

// ============================================================================
// singular elements
#include "navier_stokes_sing_face_element.h"

// The mesh
#include "meshes/triangle_mesh.h"
 
// Get the mesh
#include "meshes/tetgen_mesh.h" 
#include "meshes/refineable_tetgen_mesh.h"

// Get the faceted surfaces
#include "tetmesh_faceted_surfaces.h"

// include classes for vector and matrix algebra
#include "additional_maths.h"

using namespace oomph;

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
  string output_directory = "RESLT";
  
  /// (Half-)width of the box
  double Box_half_width = 0.25;

  /// (Half)height of the box
  double Box_half_height = 0.5;

  // pressure gradient
  double G = 1.0;
  
  unsigned Nplot_for_bulk = 5;
  
  // split corner elements which have all their nodes pinned on the outer boundary
  bool Split_corner_elements = true; 

  MeshAsGeomObject* mesh_as_geom_object_pt;

  double Target_element_volume = 0;

  // are we doing PDE-constrained minimisation or standard FEM?
  bool Do_pde_constrained_problem = true;

  bool Use_fd_lu_solver = false;
}

///////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////
namespace Analytic_Functions
{
  // derivative of the functional which we're minimising, w.r.t. the solution
  Vector<double> dfunctional_du(const Vector<double>& u)
  {
    Vector<double> dpi_du(u.size(), 0.0);

    // QUEHACERES leave at zero for the time being
    // // functional \Pi = 1/2 sum_j u_j^2,
    // // so d\Pi/d u_i = u_i
    // for(unsigned i=0; i<u.size(); i++)
    //   dpi_du[i] = u[i];

    return dpi_du;
  }

  // derivative of the functional which we're minimising, w.r.t. the solution
  double dfunctional_du(const Vector<double>& u, const unsigned& i)
  {
    return u[i];
  }

  // derivative of dPi/du_i with respect to U_jk, i.e.
  // w.r.t the jth velocity component at the kth node
  double d2pi_du_dnodal_unknown(const Vector<Vector<double> > U_ij,
				const unsigned& i, const unsigned& j,
				const unsigned& k, const Shape& shape)
  {
    return (i == j) ? shape[k] : 0.0;
  }
  
  // exact solution for the flat disk, handling linear combinations of
  // in-plane and broadside motion
  void exact_solution_poiseuille(const Vector<double>& x,
				Vector<double>& u_exact)
  {
    u_exact.resize(x.size()+1, 0.0);
    
    // y coordinate as measured from the lower wall
    double y = x[1] + Global_Parameters::Box_half_width;

    // z coord as measured from the in-flow boundary
    double z = x[2] + Global_Parameters::Box_half_height;
      
    // total width
    double w = 2.0 * Global_Parameters::Box_half_width;
    
    // planar Poiseuille flow profile
    u_exact[2] = (Global_Parameters::G / 2.0) * y * (w - y);

    // pressure - applying a constant gradient across the box
    u_exact[3] = 1.0 - z * Global_Parameters::G;    
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
class NavierStokesPdeContraintBoundaryElement : public virtual FaceGeometry<ELEMENT>, 
						public virtual FaceElement
{
 
public:

  ///Constructor, which takes a "bulk" element and the value of the index
  ///and its limit
  NavierStokesPdeContraintBoundaryElement(FiniteElement* const& element_pt, 
					  const int& face_index) : 
    FaceGeometry<ELEMENT>(), FaceElement()
    { 
      //Attach the geometrical information to the element. N.B. This function
      //also assigns nbulk_value from the required_nvalue of the bulk element
      element_pt->build_face_element(face_index, this);
      
      //Set the dimension from the nodal dimension 
      Dim = this->nodal_dimension();
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
	
	bulk_el_pt->interpolated_u_nst(s_bulk, velocity);

	// QUEHACERES debug @@@@@@@@@@
	unsigned nnode = bulk_el_pt->nnode();	  
	Vector<Node*> node_pt_list (nnode);
	{	  
	  for(unsigned j=0; j<nnode; j++)
	  {
	    node_pt_list[j] = bulk_el_pt->node_pt(j);
	  }
	}
	// @@@@@@@@@@@@@@@@@@@@@@@@@@	
	
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
  
  void get_traction(const Vector<double>& s, Vector<double>& traction) const
    {
      //Find out how many nodes there are
      const unsigned n_node = this->nnode();
      
      //Set up memory for the shape and test functions
      Shape psi(n_node), test(n_node);
     
      // shorthand
      ELEMENT* bulk_el_pt = dynamic_cast<ELEMENT*>(this->bulk_element_pt());

      // Compute outer unit normal at the specified local coordinate
      Vector<double> unit_normal(Dim, 0.0);
      this->outer_unit_normal(s, unit_normal);

      // strain rate tensor
      DenseMatrix<double> strain_rate_fe(Dim, Dim, 0.0);
     
      strain_rate(s, strain_rate_fe);

      // get the pressure
      double p = interpolated_p_nst(s);
	  
      // get contribution to the total FE stress on this element
      DenseMatrix<double> stress(Dim, Dim);

      stress = Analytic_Functions::stress(strain_rate_fe, p);

      // zero out the traction to be sure
      traction.resize(Dim, 0.0);
      std::fill(traction.begin(), traction.end(), 0.0);
      
      // now compute the traction t_i = \tau_{ij} n_j
      for(unsigned i=0; i<Dim; i++)
      {
	for(unsigned j=0; j<Dim; j++)
	{	      
	  traction[i] += stress(i,j) * unit_normal[j];
	}
      }
    }
   
  inline void fill_in_contribution_to_residuals(Vector<double>& residuals)
    {
      // shorthand
      ELEMENT* bulk_el_pt = dynamic_cast<ELEMENT*>(this->bulk_element_pt());
    
      //Find out how many nodes there are
      const unsigned n_node = this->nnode();
      const unsigned n_node_bulk = bulk_el_pt->nnode();
    
      // shorthand
      const unsigned Dim = this->Dim;

      // Number of pressure dofs is the number of vertices,
      // which in a lower dimensional face element is the dimensionality of the problem
      const unsigned n_pres = Dim;

      //Set up memory for the shape and test functions and their derivatives
      Shape psi(n_node), test(n_node);

      Shape psi_bulk(n_node_bulk);
	  
      DShape dpsidx(n_node_bulk, Dim);
    
      //Set the value of Nintpt
      const unsigned n_intpt = this->integral_pt()->nweight();
     
      //Set the Vector to hold local coordinates
      Vector<double> s(Dim-1);
    
      //Loop over the integration points
      for(unsigned ipt=0; ipt<n_intpt; ipt++)
      {       
	//Assign values of s
	for(unsigned i=0; i<(Dim-1); i++)
	{
	  s[i] = this->integral_pt()->knot(ipt,i);
	}

	// Get the local bulk coordinates of this FaceElement knot
	Vector<double> s_bulk = this->local_coordinate_in_bulk(s);
      
	//Get the integral weight
	double w = this->integral_pt()->weight(ipt);

	// Get the shape and test functions 
	this->shape_at_knot(ipt,psi);

	// Get the Jacobian of the mapping
	double J = J_eulerian_at_knot(ipt);
	
	// get the derivatives of the shape functions from the bulk element      
	bulk_el_pt->dshape_eulerian(s_bulk, psi_bulk, dpsidx);

	//Premultiply the weights and the Jacobian
	double W = w*J;

	// get the traction at this integration point
	Vector<double> traction(Dim, 0.0);
	get_traction(s, traction);

	// Compute outer unit normal at the specified local coordinate
	Vector<double> unit_normal(Dim, 0.0);
	this->outer_unit_normal(s, unit_normal);
      
	// get the momentum-enforcing Lagrange multiplier field from the bulk
	Vector<double> lambda_momentum = bulk_el_pt->interpolated_lambda(s_bulk);

	// ===================================================

	//Loop over the test functions
	for(unsigned l=0; l<n_node; l++)
	{
	  Node* node_pt = this->node_pt(l);
	  for(unsigned d=0; d<Dim; d++)
	  {
	    // Contributions to the bulk Lagrange multiplier equations which
	    // enforce the Stokes momentum PDEs
	    // ------------------------------------------------------------------
	  
	    // Index of the bulk Lagrange multiplier which enforces momentum PDEs 
	    unsigned lambda_momentum_index =
	      bulk_el_pt->index_of_lagrange_multiplier(node_pt, d);

	    int local_eqn_lagr_mom = this->nodal_local_eqn(l, lambda_momentum_index);
	  
	    if(local_eqn_lagr_mom >= 0)
	    {
	      residuals[local_eqn_lagr_mom] += traction[d] * psi[l] * W;
	    }

	    // get the bulk node number corresponding to this face element vertex node
	    const unsigned node_number_in_bulk = this->bulk_node_number(l);
	  
	    int local_eqn_augmented = this->nodal_local_eqn(l, d);
	    if (local_eqn_augmented >= 0)
	    {	      
	      for(unsigned j=0; j<Dim; j++)
	      {	      
		// boundary contribution from momentum-enforcing LMs
		// QUEHACERES write-up Eq. 4.24 term 4
		residuals[local_eqn_augmented] +=
		  (lambda_momentum[d] * unit_normal[j] +
		   lambda_momentum[j] * unit_normal[d]) * dpsidx(node_number_in_bulk, j)*W;
	      }
	    } 
	    
	  } // end loop over dimensions

	} // end loop over nodes / test functions

	// Boundary contributions of momentum LMs to the pressure residual
	// --------------------------------------------------------------------

	// shape functions from bulk
	Shape psip(bulk_el_pt->npres_nst());
      
	// loop over the pressure shape functions
	// N.B. this works because the vertex nodes are enumerated first, so
	// looping over the first DIM (=n_pres) nodes gives the nodes which store
	// the pressure; so the strategy is to take this face element's vertex
	// node numbers, look up their corresponding numbering in the bulk
	// element (via bulk_node_number(k)), and then use this index to get
	// the value of the correct shape function (which we get via the conversion
	// of the local coordinates of this face element integration point to
	// bulk local coordinates)
	for(unsigned k=0; k<n_pres; k++)
	{
	  // get the bulk node number corresponding to this face element vertex node
	  unsigned node_number_in_bulk = this->bulk_node_number(k);

	  // get the pressure basis functions from the bulk element
	  bulk_el_pt->pshape_nst(s_bulk, psip);

	  // get the local equation number for the pressure at this vertex node
	  int local_eqn_p = this->nodal_local_eqn(k, bulk_el_pt->p_index_nst());

	  if(local_eqn_p >= 0)
	  {
	    for(unsigned j=0; j<Dim; j++)
	    {
	      // QUEHACERES write-up Eq. 4.26, \partial\Gamma_C part of term 2
	      residuals[local_eqn_p] -=
		lambda_momentum[j] * unit_normal[j] * psip[node_number_in_bulk] * W;
	    }	  
	  }
	}
      } // end loop over integration points      
    }
  
private:
  unsigned Dim;
};
///////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////



//=========================================================================
/// Class that solves Navier-Stokes flow around a 2D disk
//=========================================================================
template <class ELEMENT>
class PdeConstrainedPoiseuilleProblem : public Problem
{
  
public:

  /// Constructor
  PdeConstrainedPoiseuilleProblem();
  
  /// Destructor
  ~PdeConstrainedPoiseuilleProblem()
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
    
private:
 
  /// Apply BCs and make elements functional
  void complete_problem_setup();

  /// Helper function to apply boundary conditions
  void apply_boundary_conditions();

  /// Helper to output the pin status of eaech nodal dof in each submesh
  void output_submesh_dof_pin_status() const;
  
  /// \short Helper function to create the face elements needed to:
  /// - impose Dirichlet BCs
  /// - impose the additional traction required from the augmented region to
  ///   the surrounding bulk elements
  /// - compute the reciprocity integral to determine the singular amplitude
  void create_face_elements();
  
  // --------------------------------------------------------------------------
  // Meshes
  // --------------------------------------------------------------------------
  /// Bulk mesh
  RefineableTetgenMesh<ELEMENT>* Bulk_mesh_pt;

  /// Mesh of face elements which impose Dirichlet boundary conditions
  Mesh* Face_mesh_for_bc_pt;

  Mesh* Face_mesh_for_boundary_contributions_pt;
  // --------------------------------------------------------------------------

  /// Mesh as geom object representation of mesh  
  MeshAsGeomObject* Mesh_as_geom_object_pt;
  
  // Create the mesh as Geom Object
  MeshAsGeomObject* Face_mesh_as_geom_object_pt;

  /// Storage for the outer boundary object
  TetMeshFacetedClosedSurface* Outer_boundary_pt;

  /// Inner boundary
  Vector<TetMeshFacetedSurface*> Inner_boundary_pt;
  
  /// First boundary ID for outer boundary
  unsigned First_boundary_id_for_outer_boundary;
 
  /// Sanity check: Exact bounded volume
  double Exact_bounded_volume;

  /// Number of singular functions we will subtract
  unsigned Nsingular_function;
  
  /// Number of dimensions in the problem
  unsigned Dim;
  
  DocInfo Doc_info;
};

//========================================================================
/// Constructor
//========================================================================
template <class ELEMENT>
PdeConstrainedPoiseuilleProblem<ELEMENT>::PdeConstrainedPoiseuilleProblem()
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
  Dim = ELEMENT::_DIM_;

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

  // First oomph-lib (zero-based!) boundary ID for outer boundary
  First_boundary_id_for_outer_boundary = outer_boundary_id_offset;
 
  // For sanity check:
  Exact_bounded_volume = 
    2.0*Global_Parameters::Box_half_width*
    2.0*Global_Parameters::Box_half_width*
    2.0*Global_Parameters::Box_half_height;
  
  // QUEHACERES this seems to be needed to reset FPU which Triangle
  // messes around with!
  {
    fpu_control_t cw = (_FPU_DEFAULT & ~_FPU_EXTENDED) | _FPU_DOUBLE;
    _FPU_SETCW(cw);
  }
  
  // Build the mesh
  //--------------- 

  // Not needed, of course, but here to test out the handling
  // of timesteppers
  add_time_stepper_pt(new Steady<1>);

  if(Global_Parameters::Split_corner_elements)
    oomph_info << "\nSplitting corner elements to avoid locking\n\n";
  
  Vector<double> target_element_volume_in_region(1);
  target_element_volume_in_region[0] = Global_Parameters::Target_element_volume;
  
  bool use_attributes = false;
 
  Bulk_mesh_pt =
    new RefineableTetgenMesh<ELEMENT>(Outer_boundary_pt,
					  Inner_boundary_pt,
					  Global_Parameters::Target_element_volume,
					  this->time_stepper_pt(),
					  use_attributes,
					  Global_Parameters::Split_corner_elements,
					  &target_element_volume_in_region);
  
  // Problem is linear so we don't need to transfer the solution to the
  // new mesh; we keep it on for self-test purposes...
  Bulk_mesh_pt->disable_projection();

  // Add sub-meshes
  add_sub_mesh(Bulk_mesh_pt);
   
  // Create face elements for imposition of BC
  Face_mesh_for_bc_pt = new Mesh;

  Face_mesh_for_boundary_contributions_pt = new Mesh;

  // Build the face elements
  create_face_elements();
 
  add_sub_mesh(Face_mesh_for_bc_pt);

  if(Global_Parameters::Do_pde_constrained_problem)
    add_sub_mesh(Face_mesh_for_boundary_contributions_pt);
  
  build_global_mesh();

  // // Complete problem setup
  complete_problem_setup();
  
  // Setup equation numbering scheme
  oomph_info << "--------------------------\n";
  oomph_info << "Number of equations: " << assign_eqn_numbers() << std::endl; 
  oomph_info << "--------------------------\n" << std::endl;

  // for debug - output the pin-status of each dof in each submesh
  output_submesh_dof_pin_status();

  if(Global_Parameters::Use_fd_lu_solver)
  {
    oomph_info << "Get comfy, we're using the FD_LU solver...\n" << std::endl;
    
    // brutal debug
    linear_solver_pt() = new FD_LU;
  }
}


//========================================================================
/// Complete problem setup
//========================================================================
template <class ELEMENT>
void PdeConstrainedPoiseuilleProblem<ELEMENT>::complete_problem_setup()
{
  // map which keeps track of the first index of each node associated with
  // the PDE-constraint Lagrange multipliers
  std::map<Node*, unsigned> node_to_first_lm_index_map;

  // Loop over the elements to set up element-specific
  // things that cannot be handled by constructor

  unsigned n_el = Bulk_mesh_pt->nelement();
  for (unsigned e=0; e<n_el; e++)
  {
    ELEMENT* elem_pt = dynamic_cast<ELEMENT*>(
      Bulk_mesh_pt->element_pt(e));

    // pass the function which computes stress from the
    // velocity gradient and the pressure
    elem_pt->stress_fct_pt() = &Analytic_Functions::stress;

    // // tell the element it's augmented, so that it uses the PDE-constrained min
    // // residuals rather than standard Stokes
    if(Global_Parameters::Do_pde_constrained_problem)
    {
      elem_pt->set_augmented_element();

      // pass in the function which computes the derivative of the functional
      // which we're minimising in the augmented region
      elem_pt->dfunctional_du_fct_pt() =
	&Analytic_Functions::dfunctional_du;
    
      // add the Lagrange multipliers which weakly enforce the momentum and
      // continuity equations    
      elem_pt->add_lagrange_multiplier_dofs(node_to_first_lm_index_map);
    }
  }
  
  // Apply bcs  
  apply_boundary_conditions();
}


//==start_of_create_face_elements=========================================
/// \short Helper function to create face elements needed to:
/// - impose Dirichlet BCs
/// - impose the additional traction required from the augmented region to
///   the surrounding bulk elements
/// - compute the reciprocity integral to determine the singular amplitude
//========================================================================
template <class ELEMENT>
void PdeConstrainedPoiseuilleProblem<ELEMENT>::create_face_elements()
{
  if(!Global_Parameters::Do_pde_constrained_problem)
    return;
  
  for(unsigned ibound = First_boundary_id_for_outer_boundary;
      ibound<First_boundary_id_for_outer_boundary+6; ibound++)
  {
    for(unsigned e=0; e<Bulk_mesh_pt->nboundary_element(ibound); e++)
    {
      // get a pointer to the first element on this outer face
      FiniteElement* el_pt = Bulk_mesh_pt->boundary_element_pt(ibound, e);
     
      // What is the index of the face of the bulk element at the boundary
      int face_index = Bulk_mesh_pt->face_index_at_boundary(ibound, e);

      // Build the corresponding face element
      NavierStokesPdeContraintBoundaryElement<ELEMENT>* boundary_elem_pt =
	new NavierStokesPdeContraintBoundaryElement<ELEMENT>(el_pt, face_index);

      // add the new element to the mesh
      Face_mesh_for_boundary_contributions_pt->add_element_pt(boundary_elem_pt);

      // Build the corresponding bc element
      NavierStokesWithSingularityBCFaceElement<ELEMENT>* bc_element_pt =
	new NavierStokesWithSingularityBCFaceElement<ELEMENT>
	(el_pt, face_index);

      // since we don't have a singular line mesh in this test harness, make sure
      // things don't die looking for singular element pointers
      bc_element_pt->no_throw_if_no_singular_elem_pt();
      
      //Add the bc element to the surface mesh
      Face_mesh_for_bc_pt->add_element_pt(bc_element_pt);
    }
  }  
} // end of create_face_elements()

//==start_of_apply_bc=====================================================
/// Helper function to apply boundary conditions
//========================================================================
template <class ELEMENT>
void PdeConstrainedPoiseuilleProblem<ELEMENT>::apply_boundary_conditions()
{
  // need to pin one node on the outer boundary to define the
  // pressure level
  bool have_pinned_pressure = 0;

  // QUEHACERES 
  bool have_pinned_pressure2 = 0;
  
  unsigned debug_lambda_p_pin_counter = 0;
  
  // loop over the outer boundaries
  for (unsigned ibound = First_boundary_id_for_outer_boundary;
       ibound < First_boundary_id_for_outer_boundary+6; ibound++)
  {
    unsigned nnode = Bulk_mesh_pt->nboundary_node(ibound);

    // loop the nodes on this boundary
    for(unsigned n=0; n<nnode; n++)
    {
      // get a pointer to this boundary node
      Node* node_pt = Bulk_mesh_pt->boundary_node_pt(ibound, n);

      // get it's Eulerian coords
      Vector<double> x(Dim, 0.0);
      for(unsigned i=0; i<Dim; i++)
	x[i] = node_pt->x(i);
      
      // get the exact Poiseuille flow solution at this point
      Vector<double> u_poiseuille(Dim, 0.0);
      Analytic_Functions::exact_solution_poiseuille(x, u_poiseuille);

      for(unsigned i=0; i<Dim; i++)
      {
	node_pt->pin(i);
	node_pt->set_value(i, u_poiseuille[i]);
      }
    }

    // now sort out the pressure
    if(!have_pinned_pressure)
    {
      // get a pointer to the first element on this outer face
      ELEMENT* el_pt = dynamic_cast<ELEMENT*>(Bulk_mesh_pt->boundary_element_pt(ibound, 0));
     
      // What is the index of the face of the bulk element at the boundary
      int face_index = Bulk_mesh_pt->face_index_at_boundary(ibound,0);

      // Build a face element
      NavierStokesTractionElement<ELEMENT>* face_element_pt =
	new NavierStokesTractionElement<ELEMENT>(el_pt, face_index);

      // get the outer unit normal
      Vector<double> outer_unit_normal(Dim, 0.0);
      face_element_pt->outer_unit_normal(0, outer_unit_normal);

      // // check if we've got the in-flow face, i.e. with n = (0,0,-1)
      const double tol = 1e-6;
      if(abs(outer_unit_normal[0])     < tol &&
	 abs(outer_unit_normal[1])     < tol &&
	 abs(outer_unit_normal[2]+1.0) < tol)
      {
	const unsigned face_el_node_number = 0;
	
	// get a pointer to the first boundary node
	Node* node_pt = face_element_pt->node_pt(face_el_node_number);

	// pin the pressure to 0
	node_pt->pin(Dim);
	node_pt->set_value(Dim, 1.0);

	// update the flag
	have_pinned_pressure = true;
	
	oomph_info << "Pinned pressure to 1 at ";

	for(unsigned i=0; i<Dim; i++)
	  oomph_info << node_pt->x(i) << " ";

	oomph_info << std::endl;
	
	// if we're doing a PDE-constrained problem, tell the bulk element
	// to pin the continuity-enforcing LM at the same node
	if(Global_Parameters::Do_pde_constrained_problem)
	{
	  unsigned node_num_in_bulk =
	    face_element_pt->bulk_node_number(face_el_node_number);

	  el_pt->pin_lagrange_multiplier_p(node_num_in_bulk, 0.0);

	  oomph_info << "Pinned lambda_p to " << 0.0 << " at ";
	  
	  for(unsigned i=0; i<Dim; i++)
	    oomph_info << node_pt->x(i) << " ";

	  oomph_info << std::endl;
	}
      }

      // clean up
      delete face_element_pt;
    }
    
    oomph_info << "debug_lambda_p_pin_counter = " << debug_lambda_p_pin_counter << std::endl;
    
    // now pin the bulk LMs on the boundary if we're doing a PDE-constrained problem
    if(Global_Parameters::Do_pde_constrained_problem)
    {
      unsigned nel = Face_mesh_for_boundary_contributions_pt->nelement();

      for(unsigned e=0; e<nel; e++)
      {
	// get this face element
	NavierStokesPdeContraintBoundaryElement<ELEMENT>* boundary_el_pt =
	  dynamic_cast<NavierStokesPdeContraintBoundaryElement<ELEMENT>*>(
	    Face_mesh_for_boundary_contributions_pt->element_pt(e));

	// get the bulk element this face element is attached to
	ELEMENT* bulk_el_pt = dynamic_cast<ELEMENT*>(boundary_el_pt->bulk_element_pt());

	// number of nodes in this face element
	unsigned nnod = boundary_el_pt->nnode();

	for (unsigned j=0; j<nnod; j++)
	{
	  // get the number of this node in the bulk elements numbering scheme
	  unsigned node_number_in_bulk = boundary_el_pt->bulk_node_number(j);
	  
	  // now tell the bulk element to pin the pde-enforcing LMs
	  bulk_el_pt->pin_momentum_lagrange_multipliers(node_number_in_bulk);
	  // bulk_el_pt->pin_pde_lagrange_multipliers(node_number_in_bulk);
	}
      }
    }
  } // end loop over outer boundaries

  if(!have_pinned_pressure)
  {
    oomph_info << "Couldn't pin pressure\n" << std::endl;
    abort();
  }
  
  // QUEHACERES debug
  std::set<Node*> nodes_with_3_lms_set;
  
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
    DenseMatrix<double> nodal_boundary_value(nnod, Dim, 0.0);

    // Unpin the FE part of the solution and pin bulk Lagrange multipliers
    for (unsigned j=0; j<nnod; j++)
    {
      Node* node_pt = el_pt->node_pt(j);

      Vector<double> x(Dim, 0.0);

      // number of dofs at this node
      const unsigned nvalue = node_pt->nvalue();

      // the vertex nodes are indexed first, so the indices of the vertices
      // of a tet are given by: 0 <= j <= Dim-1
      // (where Dim is the problem dim not the face element dim)
      const bool is_vertex_node = j <= Dim-1;

      // Now deal with the subtle case - if there are two sets of (boundary) Lagrange
      // multipliers at this node, they come from those that enforce the boundary
      // conditions, and those that enforce the jump in stress across the
      // boundary of the augmented region. They enforce the same constraints,
      // so one set needs pinning.
      // If this is not a vertex node (i.e. no pressure), then there will be
      // 9 dofs if there is 1 set of boundary LMs (3 velocity components, 3 momentum-
      // enforcing LMs, and 3 BC-enforcing boundary LMs), or 12 if there are 2 sets
      // of boundary LMs. If this is a vertex node, there will be either 11 dofs
      // (3 velocity components + pressure, 3 mom-LMs, 1 continuity LM, 3 BC LMs) or
      // 14 dofs (the above + 3 continuity-enforcing LMs)
      const bool has_two_sets_of_boundary_lagrange_multipliers =
	(is_vertex_node && (nvalue > 11)) || (!is_vertex_node && (nvalue > 9));

      // QUEHACERES debug
      const double r = sqrt(pow(node_pt->x(0),2) + pow(node_pt->x(1),2));
      
      for(unsigned i=0; i<Dim; i++)
      {	
	el_pt->unpin_u_fe_at_specified_local_node(j, i);

	// // *** QUEHACERES pinning to 0
	// el_pt->node_pt(j)->set_value(i, 0.0);
	
	// *** QUEHACERES pinning all of them for debug ***
	if( has_two_sets_of_boundary_lagrange_multipliers || true )
	{
	  // if we've got both continuity and BC LMs, pin the BC ones
	  // el_pt->pin_lagrange_multiplier_at_specified_local_node(j, i, BC_el_id); 
	  // el_pt->pin_lagrange_multiplier_at_specified_local_node(j, i, Stress_jump_el_id);

	  nodes_with_3_lms_set.insert(node_pt);	    
	}

	// QUEHACERES debug
	x[i] = node_pt->x(i);
      }
      
      {	
	Vector<double> u_bc(Dim, 0.0);

	// get the Poiseuille flow solution
	Analytic_Functions::exact_solution_poiseuille(x, u_bc);
	   
	// assign to the matrix of nodal values
	for(unsigned i=0; i<Dim; i++)
	{	  
	  nodal_boundary_value(j,i) = u_bc[i];
	}
      }

      // Now pin the Lagrange multipliers that enforce the PDE-constraints,
      // since we're applying Dirichlet BCs here
      // ------------------------------------------------------------------
      if(Global_Parameters::Do_pde_constrained_problem)
      {
	// get the bulk element this face element is attached to
	ELEMENT* bulk_el_pt = dynamic_cast<ELEMENT*>(el_pt->bulk_element_pt());

	// get the number of this node in the bulk elements numbering scheme
	unsigned node_number_in_bulk = el_pt->bulk_node_number(j);

	// now tell the bulk element to pin the pde-enforcing LMs
	bulk_el_pt->pin_momentum_lagrange_multipliers(node_number_in_bulk);
	// bulk_el_pt->pin_pde_lagrange_multipliers(node_number_in_bulk);
      }
      
    } // end loop over BC nodes
    
    // Tell the element about these nodal boundary values
    el_pt->set_nodal_boundary_values(nodal_boundary_value);

    
  } // end loop over bc face elements

  oomph_info << "LM pin counter = " << nodes_with_3_lms_set.size() << std::endl;
} // end apply BCs

//== start of output_submesh_pin_status ==================================
/// Function to output the pin status of each nodal dof in each submesh
//========================================================================
template <class ELEMENT>
void PdeConstrainedPoiseuilleProblem<ELEMENT>::output_submesh_dof_pin_status() const
{
  ofstream pin_file;
  std::ostringstream filename;

  // make a vector of pointers to the submeshes
  Vector<Mesh*> mesh_pt;
  mesh_pt.push_back(Face_mesh_for_bc_pt);
  mesh_pt.push_back(Bulk_mesh_pt);
  
  // makea vector of names for each sub-mesh
  Vector<std::string> mesh_name;
  mesh_name.push_back("bc");
  mesh_name.push_back("bulk_mesh");
  
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
	  for(unsigned i=0; i<Dim; i++)
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

//////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////



//========================================================================
/// Doc the solution
//========================================================================
template <class ELEMENT>
void PdeConstrainedPoiseuilleProblem<ELEMENT>::doc_solution(const unsigned& nplot)
{
  bool do_bulk_output = true;
  if (CommandLineArgs::command_line_flag_has_been_set("--suppress_bulk_output"))
  {
    do_bulk_output = false;
  }

  ofstream some_file;
  ofstream some_file2;
  ofstream face_some_file;
  ofstream coarse_some_file;
  char filename[100];
  
  // Output solution showing element outlines
  //-----------------------------------------
  oomph_info << "Outputting coarse solution..." << std::endl;
  sprintf(filename,"%s/coarse_soln%i.vtu",Doc_info.directory().c_str(),
	  Doc_info.number());
  some_file.open(filename);
  if (do_bulk_output) Bulk_mesh_pt->output_paraview(some_file,2);
  some_file.close();

  // Output solution showing element outlines
  //-----------------------------------------
  oomph_info << "Outputting solution..." << std::endl;
  sprintf(filename,"%s/soln%i.vtu",Doc_info.directory().c_str(),
	  Doc_info.number());
  some_file.open(filename);
  if (do_bulk_output) Bulk_mesh_pt->output_paraview(some_file, nplot);
  some_file.close();
  
  // Exact solution (only need to output it once)
  if (Doc_info.number() == 0 &&
    !CommandLineArgs::command_line_flag_has_been_set("--dont_output_exact_solution"))
  {
    oomph_info << "Outputting exact solution..." << std::endl;
    sprintf(filename,"%s/exact_soln.vtu",Doc_info.directory().c_str() );
    some_file.open(filename);

    Bulk_mesh_pt->output_fct_paraview(some_file, nplot,
				      Analytic_Functions::exact_solution_poiseuille);
    
    some_file.close();
  }

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
    el_pt->compute_error(dummy_ofstream, &Analytic_Functions::exact_solution_poiseuille,
			 el_v_error, el_p_error, el_norm);
    }
    
    //Add each elemental error to the global error
    norm    += el_norm;
    v_error += el_v_error;
    p_error += el_p_error;
    
    // do the actual paraview'able output at plot points
    el_pt->output_error_at_plot_points(some_file, nplot,
				       &Analytic_Functions::exact_solution_poiseuille);
  }

  some_file.close();

  oomph_info << "L2 velocity error in total solution: " << v_error << std::endl;
  oomph_info << "L2 pressure error in total solution: " << p_error << std::endl;

  
  if (CommandLineArgs::command_line_flag_has_been_set("--output_jacobian_sparse"))
  {  
    // residual vector and Jacobian matrix
    DoubleVector r;
    CRDoubleMatrix jac;

    // get 'em
    get_jacobian(r,jac);
  
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

  // oomph_info << "Outputting extended solution..." << std::endl;
  
  // // Plot "extended solution" showing contributions
  // sprintf(filename,"%s/extended_soln%i.dat",Doc_info.directory().c_str(),
  // 	  Doc_info.number());
  // some_file.open(filename);
  
  // unsigned nel = Bulk_mesh_pt->nelement();
  // for (unsigned e=0; e<nel; e++)
  // {
  //   // shouldn't change this, since the maps have been setup for a specific number
  //   // unsigned npts = Global_Parameters::Nplot_for_bulk;
  //   ELEMENT* el_pt = dynamic_cast<ELEMENT*>(Bulk_mesh_pt->element_pt(e));

  //   el_pt->output_with_various_contributions(some_file, nplot,
  // 					     Analytic_Functions::exact_solution_poiseuille);
  // }

  // some_file.close();

  // ==========================================================================
  // output the PDE-enforcing LMs from the bulk elements
  // QUEHACERES do this at plot points

  if(Global_Parameters::Do_pde_constrained_problem)
  {
    sprintf(filename,"%s/pde_enforcing_lms%i.dat",Doc_info.directory().c_str(),Doc_info.number());      
    some_file.open(filename);
    
    unsigned nel = Bulk_mesh_pt->nelement();
    for(unsigned e=0; e<nel; e++)
    {
      ELEMENT* elem_pt = dynamic_cast<ELEMENT*>(Bulk_mesh_pt->element_pt(e));

      unsigned nnode = elem_pt->nnode();
      for(unsigned n=0; n<nnode; n++)
      { 
	Node* node_pt = elem_pt->node_pt(n);

	// output Eulerian position
	for(unsigned i=0; i<Dim; i++)
	  some_file << node_pt->x(i) << " ";

	// output the momentum LMs
	for(unsigned i=0; i<Dim; i++)
	{
	  unsigned lambda_i_index = elem_pt->index_of_lagrange_multiplier(n, i);
	  some_file << node_pt->value(lambda_i_index) << " ";
	}

	// output continuity LM
	if(n <= Dim)
	{
	  unsigned lambda_p_index = elem_pt->index_of_lagrange_multiplier_p(n);

	  some_file << node_pt->value(lambda_p_index);
	}
	some_file << std::endl;
      
      }
    }

    some_file.close();
  }
  
  sprintf(filename,"%s/global_mesh_nodal_values%i.dat",Doc_info.directory().c_str(),Doc_info.number());      
  some_file.open(filename);

  for(unsigned j=0; j<mesh_pt()->nnode(); j++)
  {
    Node* node_pt = mesh_pt()->node_pt(j);

    some_file << j << " ";
    
    // output Eulerian position
    for(unsigned i=0; i<Dim; i++)
      some_file << node_pt->x(i) << " ";
      
    for(unsigned n=0; n<node_pt->nvalue(); n++)
    {
      some_file << node_pt->value(n) << " ";
    }

    some_file << std::endl;    
  }
  
  some_file.close();
  
  //Increment counter for solutions 
  Doc_info.number()++;
  
  oomph_info << "Finished documenting solution.\n\n";
} // end of doc






//========================================================================
/// Driver
//========================================================================
int main(int argc, char* argv[])
{  
#ifdef USE_FD_JACOBIAN
  oomph_info << "====== Using finite-diff jacobian ======\n";
#else
  oomph_info << "====== Using analytic jacobian ======\n\n";
#endif

  // set up the multi-processor interface
  MPI_Helpers::init(argc,argv);

  // keep track of total program runtime
  double t_start = TimingHelpers::timer();
  
  // Store command line arguments
  CommandLineArgs::setup(argc,argv);
  
  // extra output to describe equation numbers
  CommandLineArgs::specify_command_line_flag("--describe_dofs");
  CommandLineArgs::specify_command_line_flag("--describe_nodes");
  
  // half width of the container box (disk radius is 1)
  CommandLineArgs::specify_command_line_flag("--box_half_width",
					     &Global_Parameters::Box_half_width);

  // half length of the container box
  CommandLineArgs::specify_command_line_flag("--box_half_height",
					     &Global_Parameters::Box_half_height);

  CommandLineArgs::specify_command_line_flag("--element_volume", &Global_Parameters::Target_element_volume);
  
  // prevent the splitting of corner elements
  CommandLineArgs::specify_command_line_flag("--dont_split_corner_elements");

  // get output directory
  CommandLineArgs::specify_command_line_flag("--dir", &Global_Parameters::output_directory);

  // number of plot points per side in the bulk elements
  CommandLineArgs::specify_command_line_flag("--nplot", &Global_Parameters::Nplot_for_bulk);

  // doc the solution before any solves
  CommandLineArgs::specify_command_line_flag("--doc_initial_conditions");

  CommandLineArgs::specify_command_line_flag("--output_jacobian_full");
  CommandLineArgs::specify_command_line_flag("--output_jacobian_sparse");
  CommandLineArgs::specify_command_line_flag("--output_initial_jacobian");

  CommandLineArgs::specify_command_line_flag("--standard_fem");

  CommandLineArgs::specify_command_line_flag("--use_fd_lu");
    
  // Parse command line
  CommandLineArgs::parse_and_assign();
 
  // Doc what has actually been specified on the command line
  CommandLineArgs::doc_specified_flags();
    
  if (CommandLineArgs::command_line_flag_has_been_set("--dont_split_corner_elements"))
  {
    Global_Parameters::Split_corner_elements = false;
  }

  if (CommandLineArgs::command_line_flag_has_been_set("--standard_fem"))
  {
    Global_Parameters::Do_pde_constrained_problem = false;
  }

  if(CommandLineArgs::command_line_flag_has_been_set("--use_fd_lu"))
  {
    Global_Parameters::Use_fd_lu_solver = true;
  }
  // Shut up prefix
  oomph_info.output_modifier_pt() = &default_output_modifier;

  // PdeConstrainedPoiseuilleProblem <ProjectableTaylorHoodElement<
  //   TNavierStokesElementWithSingularity<3,3> > > problem; // TTaylorHoodElement<3>::NNODE_1D>
  PdeConstrainedPoiseuilleProblem <ProjectableTaylorHoodElement<TNavierStokesElementWithSingularity<3> > > problem;
 
  // // QUEHACERES for debug
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
    
  unsigned max_adapt = 0; 
  for (unsigned i=0; i<=max_adapt; i++)
  {
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

      	char filename[100];
      	ofstream some_file;
      
      	sprintf(filename,"%s/jacobian_sparse%i.dat", problem.doc_info().directory().c_str(),
      		problem.doc_info().number());
      	some_file.open(filename);
      
      	jac.sparse_indexed_output(some_file);

      	some_file.close();
      }
      
      // Solve the bastard!
      problem.newton_solve();
      
#ifdef PRINT_SINGULAR_JACOBIAN
    }
    catch (const std::exception& ex)
    {
      // output the jacobian if it's singular
      
      // residual vector and Jacobian matrix
      DoubleVector r;
      CRDoubleMatrix jac;

      problem.get_jacobian(r,jac);

      char filename[100];
      ofstream some_file;
      
      sprintf(filename,"%s/singular_jacobian_sparse%i.dat", problem.doc_info().directory().c_str(),
	      problem.doc_info().number());
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
	oomph_info << "Doing eigenstuff..." << std::endl;

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
	problem.get_dofs(singular_vector);
	singular_vector.initialise(0.0);
      
	// Do it
	eigen_solver.find_eigenvalues(DenseA, DenseM, eval, evec);

	// tolerance on an eigenvalue being zero
	double tol = 1e-12;
	
	for(unsigned i=0; i<n; i++)
	{
	  // std::cout << "Eigenvalue " << i << " : " << eval[i] << "\n";
	  if (fabs( real(eval[i]) ) < tol)
	  {
	    for(unsigned j=0; j<n; j++)
	    {
	      singular_vector[j] = real( evec[i][j] );
	      // std::cout << evec[i][j] << ", ";
	    }
	  }
	  // std::cout << std::endl;
	}
	sprintf(filename, "%s/sing_eigen_vect.dat", problem.doc_info().directory().c_str());
	some_file.open(filename);
	
	for (unsigned i=0; i<n; i++)
	{
	  // oomph_info << "Singular eigenvector " << i
	  // 	     << " " << singular_vector[i] << std::endl;
	  some_file << singular_vector[i] << std::endl;
	}
	some_file.close();
	problem.assign_eigenvector_to_dofs(singular_vector);
	problem.doc_solution(nplot); 
      
	oomph_info << "Done eigenstuff; bailing\n";
	exit(0);
      }
    }
    
#endif
    
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

