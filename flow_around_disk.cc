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

//Generic routines
#include "generic.h"

// QUEHACERES delete
// // Poisson
#include "poisson.h"

// singular elements
#include "navier_stokes_sing_face_element.h"

// The mesh
#include "meshes/triangle_mesh.h"
 
// Get the mesh
#include "meshes/tetgen_mesh.h" 
#include "meshes/refineable_tetgen_mesh.h"
#include "meshes/gmsh_tet_mesh.h"

// Get the faceted surfaces
#include "tetmesh_faceted_surfaces.h"

// Tetgen or Gmsh
#define DO_TETGEN


using namespace oomph;

///////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////

double atan2pi(const double& y, const double& x)
{
  // Polar angle
  double theta = atan2(y,x);

  // prevent atan2 negative angle fuckery that causes a discontinuity at theta=pi
  if (y < 0.0)
  {
    theta += 2.0 * MathematicalConstants::Pi;
  }

  return theta;
}

//=============================================================
/// Namespace for problem parameters
//=============================================================
namespace Global_Parameters
{
  string output_directory = "RESLT";
  
  /// (Half-)width of the box
  double Box_half_width = 1.5;

  /// (Half)height of the box
  double Box_half_length = 1.0;

  /// Specify how to call gmsh from the command line
  std::string Gmsh_command_line_invocation="/home/mheil/gmesh/bin/bin/gmsh";

  // velocity of the whole disk (rigid)
  Vector<double> disk_velocity(3, 0.0);

  // amplitude and wavenumber of the warped disk
  double epsilon = 0.1;
  unsigned n = 5;

  // offset for boundary ID numbering, essentially the newly created upper
  // disk boundaries will be numbered as the corresponding lower disk boundary
  // plus this offset. The offset will be determined during the duplication
  // of plate nodes.
  unsigned upper_disk_boundary_offset = 0;
}

///////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////

//=============================================================================
// face elements, for attaching to boundaries to measure surface areas and
// to output values (to check boundary conditions)
//=============================================================================
template <class ELEMENT>
class NavierStokesFaceElement : public virtual FaceGeometry<ELEMENT>, 
				public virtual FaceElement
{
 
public:

  ///Constructor, which takes a "bulk" element and the value of the index
  ///and its limit
  NavierStokesFaceElement(FiniteElement* const& element_pt, 
			  const int& face_index) : 
    FaceGeometry<ELEMENT>(), FaceElement()
    { 
      //Attach the geometrical information to the element. N.B. This function
      //also assigns nbulk_value from the required_nvalue of the bulk element
      element_pt->build_face_element(face_index, this);
      
      //Set the dimension from the dimension of the first node
      Dim = node_pt(0)->ndim();
    }

  /// \short Output function:  
  ///  Output x, y, u, v, p, theta at Nplot^DIM plot points
  // Start of output function
  void output(std::ostream &outfile, const unsigned &nplot)
    {

      // Get pointer to assocated bulk element
      ELEMENT* bulk_el_pt = dynamic_cast<ELEMENT*>(bulk_element_pt());
   
      // Elemental dimension
      unsigned dim_el = dim();

      unsigned dim_bulk = dim_el + 1;
      
      //Local coordinates
      Vector<double> s(dim_el);

      Vector<double> s_bulk(dim_bulk);
   
      //Eulerian coordinates
      Vector<double> x(dim_bulk, 0.0);

      // Velocity from bulk element
      Vector<double> velocity(dim_bulk);
      
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
	for(unsigned i=0; i<dim_bulk; i++)
	{
	  outfile << x[i] << " ";
	}

	// output fluid velocities
	for(unsigned i=0; i<dim_bulk; i++)
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
      // Elemental dimension
      unsigned dim_el = dim();

      unsigned dim_bulk = dim_el + 1;
      
      // local coordinates in bulk element
      Vector<double> s_bulk(dim_bulk);
      
      // Get pointer to assocated bulk element
      ELEMENT* bulk_el_pt = dynamic_cast<ELEMENT*>(bulk_element_pt());

      this->get_local_coordinate_in_bulk(s, s_bulk);
      
      //Get x position from bulk
      bulk_el_pt->interpolated_x(s_bulk, x);      
    }

  void interpolated_u_nst(const Vector<double>& s, Vector<double>& u)
    {
      // Elemental dimension
      unsigned dim_el = dim();

      unsigned dim_bulk = dim_el + 1;
      
      // local coordinates in bulk element
      Vector<double> s_bulk(dim_bulk);
      
      // Get pointer to assocated bulk element
      ELEMENT* bulk_el_pt = dynamic_cast<ELEMENT*>(bulk_element_pt());

      this->get_local_coordinate_in_bulk(s, s_bulk);
      
      //Get velocity from bulk
      bulk_el_pt->interpolated_u_nst(s_bulk, u);      
    }

  double interpolated_p_nst(const Vector<double>& s)
    {
      // Elemental dimension
      unsigned dim_el = dim();
      
      unsigned dim_bulk = dim_el + 1;
      
      // local coordinates in bulk element
      Vector<double> s_bulk(dim_bulk);
      
      // Get pointer to assocated bulk element
      ELEMENT* bulk_el_pt = dynamic_cast<ELEMENT*>(bulk_element_pt());

      this->get_local_coordinate_in_bulk(s, s_bulk);
      
      //Get pressure from bulk
      return bulk_el_pt->interpolated_p_nst(s_bulk);
    }

private:
  unsigned Dim;
};
///////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////

//=========================================================================
/// Class that solves Navier-Stokes flow around a 2D disk using Gmsh mesh
//=========================================================================
template<class ELEMENT>
class FlowAroundDiskProblem : public Problem
{

public:

  /// Constructor
  FlowAroundDiskProblem();
  
  /// Destructor (empty)
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
  
private:
 
  /// Apply BCs and make elements functional
  void complete_problem_setup();

  /// Helper function to apply boundary conditions
  void apply_boundary_conditions();

  // function to populate the vectors Elements_on_upper[lower]_disk_surface_pt
  // with the elements which have at least one node on the disk.
  // (These will be used when the plate nodes are duplicated, so that we
  // know which elements need to be told that their node pointers have changed).
  void identify_elements_on_upper_and_lower_disk_sufaces();

  void duplicate_plate_nodes_and_add_boundaries();

  /// \short function to populate a map which maps each node in the mesh to
  /// a set of all the elements it is associated with
  void generate_node_to_element_map();
  
#ifdef DO_TETGEN

  /// Bulk mesh
  RefineableTetgenMesh<ELEMENT>* Bulk_mesh_pt;

#else

  /// Bulk mesh
  RefineableGmshTetMesh<ELEMENT>* Bulk_mesh_pt;

#endif

  // Mesh of upper disk elements for output
  Mesh* Face_mesh_for_upper_disk_pt;
  
  // Create the mesh as Geom Object
  MeshAsGeomObject* Face_mesh_as_geom_object_pt;
  
  /// Storage for the outer boundary object
  TetMeshFacetedClosedSurface* Outer_boundary_pt;

  /// Inner boundary
  Vector<TetMeshFacetedSurface*> Inner_boundary_pt;

  /// First boundary ID for outer boundary
  unsigned First_boundary_id_for_outer_boundary;

  
  // Disk with torus round the edges
  //--------------------------------

  // store the warped disk object so that we can use it to get
  // surface normals
  WarpedCircularDiskWithAnnularInternalBoundary* Warped_disk_with_boundary_pt;
  
  /// Region ID for torus around edge of warped disk
  unsigned Torus_region_id;

  /// First boundary ID for disk that is surrounded by torus
  unsigned First_lower_disk_boundary_id;
 
  /// Last boundary ID for disk that is surrounded by torus
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
  
  // Volumes
  //--------

  /// Sanity check: Exact bounded volume
  double Exact_bounded_volume;

  DocInfo Doc_info;
};



//========================================================================
/// Constructor
//========================================================================
template<class ELEMENT>
FlowAroundDiskProblem<ELEMENT>::FlowAroundDiskProblem()
{
  
#ifdef OOMPH_HAS_MPI
  std::cout << "This code has been compiled with mpi support \n " 
	    << "and is running on " << this->communicator_pt()->nproc() 
	    << " processors. " << std::endl;
#else
  std::cout << "This code has NOT been compiled mpi support" 
	    << std::endl;
#endif

  // Output directory
  Doc_info.set_directory(Global_Parameters::output_directory);
  
  // OUTER BOUNDARY
  //===============

  // Start boundary IDs for outer boundary from some crazy offset
  // (just for testing). By default the one-based boundary IDs go from
  // 1 to 6; let's start from 1001.
  unsigned outer_boundary_id_offset = 1000;

  //Make the outer boundary object
  Outer_boundary_pt = new CubicTetMeshFacetedSurface(
    Global_Parameters::Box_half_width,
    Global_Parameters::Box_half_length,
    outer_boundary_id_offset);

  // // Look, we can visualise the faceted surface!
  // Outer_boundary_pt->output("outer_faceted_surface.dat");

  // First oomph-lib (zero-based!) boundary ID for outer boundary
  First_boundary_id_for_outer_boundary = outer_boundary_id_offset;
 
  // For sanity check:
  Exact_bounded_volume = 
    2.0*Global_Parameters::Box_half_width*
    2.0*Global_Parameters::Box_half_width*
    2.0*Global_Parameters::Box_half_length;

 
  // INTERNAL BOUNDARIES
  //====================

  // A warped disk surrounded by a torus
  //------------------------------------
 
  // Radius of torus region
  double r_torus = 0.2;

  // Warped disk with specified amplitude and wavenumber for warping

  // (Half) number of segments used to represent the disk perimeter
  unsigned half_nsegment = 30; // 30; 
  
  // Thickness of annular region on disk = radius of torus surrounding the
  // edge
  double h_annulus = r_torus;
  Warped_disk_with_boundary_pt = 
    new WarpedCircularDiskWithAnnularInternalBoundary(h_annulus,
						      Global_Parameters::epsilon,
						      Global_Parameters::n);
  
  // Number of vertices around perimeter of torus
  unsigned nvertex_torus=10; //20;

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
      Warped_disk_with_boundary_pt,
      half_nsegment,
      r_torus,
      nvertex_torus,
      first_one_based_disk_with_torus_boundary_id,
      one_based_torus_region_id, 
      last_one_based_disk_with_torus_boundary_id,
      first_one_based_torus_boundary_id,
      last_one_based_torus_boundary_id,
      One_based_boundary_id_for_disk_within_torus,
      One_based_boundary_id_for_disk_outside_torus);


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
  
  // Look, we can visualise the faceted surface!
  disk_with_torus_pt->output("warped_disk_with_torus_faceted_surface.dat");
 
  // Add as inner boundary for mesh
  Inner_boundary_pt.push_back(disk_with_torus_pt);

  // Build the mesh
  //--------------- 

  // Initial element volume
  double initial_element_volume = 1.0;

  // Setup parameters for gmsh
  GmshParameters* gmsh_parameters_pt = 
    new GmshParameters(Outer_boundary_pt,
		       Global_Parameters::Gmsh_command_line_invocation);

  // Element volume
  gmsh_parameters_pt->element_volume() = initial_element_volume;


  // Specify inner boundaries
  gmsh_parameters_pt->internal_surface_pt() = Inner_boundary_pt;

  // Filename for file in which target element size is stored
  // (for disk-based operation of gmsh)
  gmsh_parameters_pt->stem_for_filename_gmsh_size_transfer() = 
    "target_size_on_grid";
  gmsh_parameters_pt->counter_for_filename_gmsh_size_transfer() = 0;

  // Problem is linear so we don't need to transfer the solution to the
  // new mesh; we keep it on for self-test purposes...
  // gmsh_parameters_pt->disable_projection();

  // Redirect gmsh on-screen output
  gmsh_parameters_pt->gmsh_onscreen_output_file_name() = 
    "RESLT/gmsh_on_screen_output.dat";

  // Not needed, of course, but here to test out the handling
  // of timesteppers
  add_time_stepper_pt(new Steady<1>);

#ifdef DO_TETGEN

  // And now build it... 
  Bulk_mesh_pt =
    new RefineableTetgenMesh<ELEMENT>(Outer_boundary_pt,
					 Inner_boundary_pt,
					 initial_element_volume,
					 this->time_stepper_pt());
  
  // Problem is linear so we don't need to transfer the solution to the
  // new mesh; we keep it on for self-test purposes...
  Bulk_mesh_pt->disable_projection();

#else

  // And now build it...
  Bulk_mesh_pt = new RefineableGmshTetMesh<ELEMENT>(gmsh_parameters_pt,
						       this->time_stepper_pt());

#endif

  // Set error estimator for bulk mesh
  Z2ErrorEstimator* error_estimator_pt = new Z2ErrorEstimator;
  Bulk_mesh_pt->spatial_error_estimator_pt() = error_estimator_pt;

  // Set targets for spatial adaptivity
  Bulk_mesh_pt->max_permitted_error() = 0.0005; 
  Bulk_mesh_pt->min_permitted_error() = 0.00001;

  // --------------------------------------------------
  char filename[100];
  ofstream some_file;
  sprintf(filename,"%s/boundaries%i.dat",Doc_info.directory().c_str(),
	  Doc_info.number());
  some_file.open(filename);
  Bulk_mesh_pt->output_boundaries(some_file);
  some_file.close();
  // --------------------------------------------------

  Face_mesh_for_upper_disk_pt = new Mesh;
    
  // populate the vectors which contain pointers to the elements which
  // have nodes on the disk and are identified as being on the upper or
  // lower disk surface. Also populates the face mesh
  identify_elements_on_upper_and_lower_disk_sufaces();

  // Duplicate plate nodes and add upper boundaries
  duplicate_plate_nodes_and_add_boundaries();

  // Add sub-meshes
  add_sub_mesh(Bulk_mesh_pt);

  // QUEHACERES
  // add_sub_mesh(Face_mesh_for_upper_disk_pt);
  
  build_global_mesh();
  
  // Complete problem setup
  complete_problem_setup();
 
// #ifdef OOMPH_HAS_HYPRE
  
  // // Create a new Hypre linear solver
  // HypreSolver* hypre_linear_solver_pt = new HypreSolver;
 
  // // Set the linear solver for problem
  // linear_solver_pt() = hypre_linear_solver_pt;
 
  // // Set some solver parameters
  // hypre_linear_solver_pt->max_iter() = 100;
  // hypre_linear_solver_pt->tolerance() = 1e-10;
  // hypre_linear_solver_pt->amg_simple_smoother() = 1;
  // hypre_linear_solver_pt->disable_doc_time();
  // hypre_linear_solver_pt->enable_hypre_error_messages();
  // hypre_linear_solver_pt->amg_print_level() = 0;
  // hypre_linear_solver_pt->krylov_print_level() = 0;
  // hypre_linear_solver_pt->hypre_method() = HypreSolver::BoomerAMG;
   
// #endif

// #ifdef OOMPH_HAS_MUMPS

//   oomph_info << "Using MUMPS linear solver\n\n";
  
//   MumpsSolver* mumps_linear_solver_pt = new MumpsSolver;

//   mumps_linear_solver_pt->enable_suppress_warning_about_MPI_COMM_WORLD();
  
//   // set it
//   linear_solver_pt() = mumps_linear_solver_pt;

// #endif
  
  // Setup equation numbering scheme
  oomph_info <<"Number of equations: " << assign_eqn_numbers() << std::endl; 

  // @@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
  // QUEHACERES debug
  sprintf(filename, "%s/disk_boundary_ids_in_torus.dat", Doc_info.directory().c_str());
  some_file.open(filename);

  for(Vector<unsigned>::iterator it=One_based_boundary_id_for_disk_within_torus.begin();
      it != One_based_boundary_id_for_disk_within_torus.end(); it++)
  {
    unsigned ibound = (*it)-1;
    unsigned nnode = Bulk_mesh_pt->nboundary_node(ibound);

    for(unsigned j=0; j<nnode; j++)
    {
      Node* node_pt = Bulk_mesh_pt->boundary_node_pt(ibound, j);
      for(unsigned i=0; i<3; i++)
      {
	some_file << node_pt->x(i) << " ";
      }
      some_file << std::endl;
    }
    some_file << "\n\n";
  }

  some_file.close();

  sprintf(filename, "%s/disk_boundary_ids_outside_torus.dat", Doc_info.directory().c_str());
  some_file.open(filename);

  for(Vector<unsigned>::iterator it=One_based_boundary_id_for_disk_outside_torus.begin();
      it != One_based_boundary_id_for_disk_outside_torus.end(); it++)
  {
    unsigned ibound = (*it)-1;
    unsigned nnode = Bulk_mesh_pt->nboundary_node(ibound);

    for(unsigned j=0; j<nnode; j++)
    {
      Node* node_pt = Bulk_mesh_pt->boundary_node_pt(ibound, j);
      for(unsigned i=0; i<3; i++)
      {
	some_file << node_pt->x(i) << " ";
      }
      some_file << std::endl;
    }
    some_file <<"\n\n";
  }

  some_file.close();  
  // @@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
}

// ////////////////////////////////////////////////////////////////////////////
// ////////////////////////////////////////////////////////////////////////////
// ////////////////////////////////////////////////////////////////////////////
// ////////////////////////////////////////////////////////////////////////////


//========================================================================
/// \short Function to populate a map which maps each node in the mesh to
/// a set of all the elements it is associated with
//========================================================================
template<class ELEMENT>
void FlowAroundDiskProblem<ELEMENT>::generate_node_to_element_map()
{
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
template<class ELEMENT>
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
      NavierStokesFaceElement<ELEMENT>* surface_element_pt =
      	new NavierStokesFaceElement<ELEMENT>(el_pt, face_index);

      Vector<double> unit_normal(3);
      surface_element_pt->outer_unit_normal(0, unit_normal);

      // if the z component of the boundary surface normal is negative,
      // then the element is sat on the top surface of the disk
      if(unit_normal[2] < 0)
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

	// and finally, add the face element 
	Face_mesh_for_upper_disk_pt->add_element_pt(surface_element_pt);
      }
      else 
      {
	// otherwise, it must be below
	Elements_on_lower_disk_surface_pt.insert(el_pt);

	// clean up (only if this is a lower element, don't want to delete the
	// surface element if it's an upper element since we're keeping them for
	// output
	delete surface_element_pt;
      }
      
    }
  }
  
  // QUEHACERES move this to after so that we output all of them including the touching elems
  // =================================================================
  // QUEHACERES debug - output the elements with faces on the disk
  // =================================================================
  
  // just output the vertices for now
  unsigned nplot = 2;
  char filename[100];
  ofstream some_file;
  
  sprintf(filename, "%s/elements_on_upper_disk_surface.dat", Doc_info.directory().c_str());
  some_file.open(filename);
  
  for(typename std::set<ELEMENT*>::iterator it = Elements_on_upper_disk_surface_pt.begin();
      it != Elements_on_upper_disk_surface_pt.end(); it++)
  {
    (*it)->output(some_file, nplot);
  }

  some_file.close();

  sprintf(filename, "%s/elements_on_lower_disk_surface.dat", Doc_info.directory().c_str());
  some_file.open(filename);
  
  for(typename std::set<ELEMENT*>::iterator it = Elements_on_lower_disk_surface_pt.begin();
      it != Elements_on_lower_disk_surface_pt.end(); it++)
  {
    (*it)->output(some_file, nplot);
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
  //   We'll do a majority vote, so if most of the other nodes have a negative
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
      
      // vector to store the outer unit normal for the upper surface at this node
      Vector<double> unit_normal(3);
           
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
	  NavierStokesFaceElement<ELEMENT>* surface_element_pt =
	    new NavierStokesFaceElement<ELEMENT>(el_pt, face_index);

	  // find the index which corresponds to this node
	  unsigned nodal_index;
	  for(unsigned j=0; j<surface_element_pt->nnode(); j++)
	  {
	    if(surface_element_pt->node_pt(j) == node_of_interest_pt)
	    {
	      nodal_index = j;
	      break;
	    }
	  }
	  
	  // get the outer unit normal
	  surface_element_pt->outer_unit_normal(nodal_index, unit_normal);

	  delete surface_element_pt;

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
      std::set<ELEMENT*> element_set = Node_to_element_map[node_of_interest_pt];
      
      for(typename std::set<ELEMENT*>::iterator element_iter =
	    element_set.begin(); element_iter != element_set.end();
	  element_iter++)
      {
	ELEMENT* el_pt = *element_iter;// dynamic_cast<ELEMENT*>(Bulk_mesh_pt->element_pt(e));

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
	
	unsigned index_of_boundary_node;

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
	  {
	    index_of_boundary_node = k;
	    continue;
	  }

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
	  
	  // add to the counters (N.B. a zero distance would contribute to
	  // the negative count, which is ok since negative means upper surface)
	  if(normal_distance > 0)
	    npositive++;
	  else
	    nnegative++;
	}

	// since the unit normal is in the outward direction, a negative
	// distance in the normal direction indicates the other nodes are
	// on the upper side of the disk (since the outer unit normal we
	// got was for the upper surface).
	bool element_is_on_upper_surface = nnegative > npositive;
		  
	if(element_is_on_upper_surface)
	{
	  // a non-boundary element can only have at most one edge touching the
	  // disk surface (if it had a face it would be considered a normal
	  // boundary element), so for a Taylor Hood element (NNODE_1D=3) on
	  // the upper surface, this means no more than 2 nodes can have a
	  // positive distance from the boundary node of interest in the outer
	  // normal direction... if this has happened then something weird has
	  // gone wrong, throw an error
	  if(npositive > 2)
	  {
	    // add this element to our list of dodgy ones so we can output later
	    dodgy_upper_element_pt.push_back(el_pt);

	    // add the vector of nodal distances to our upper vector
	    dodgy_upper_element_nodal_distance.push_back(nodal_distance);
	  }

	  // QUEHACERES
	  // // we can now fill in the entry into our map which gives the set of
	  // // elements and associated nodal indices for each boundary node
	  // Disk_node_to_upper_disk_element_and_index_map[node_of_interest_pt].insert(
	  //   std::pair<ELEMENT*, unsigned>(el_pt, index_of_boundary_node));
	  
	  
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
	  // negative distance from the boundary node of interest in the outer
	  // normal direction... if this has happened then something weird has
	  // gone wrong, throw an error	  if(npositive > 2)
	  if(nnegative > 2)
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
  for(typename std::set<ELEMENT*>::iterator it =
	Nonboundary_elements_with_node_on_upper_disk_surface_pt.begin();
      it != Nonboundary_elements_with_node_on_upper_disk_surface_pt.end(); it++)
  {
    ELEMENT* el_pt = *it;

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
  // =================================================================
  // QUEHACERES debug - now output the elements with nodes but
  // no faces on the disk
  // =================================================================
    
  sprintf(filename, "%s/elements_touching_upper_disk_surface.dat", Doc_info.directory().c_str());
  some_file.open(filename);
  
  for(typename std::set<ELEMENT*>::iterator it =
	Nonboundary_elements_with_node_on_upper_disk_surface_pt.begin();
      it != Nonboundary_elements_with_node_on_upper_disk_surface_pt.end(); it++)
  {
    (*it)->output(some_file, nplot);
  }

  some_file.close();

  sprintf(filename, "%s/elements_touching_lower_disk_surface.dat", Doc_info.directory().c_str());
  some_file.open(filename);
  
  for(typename std::set<ELEMENT*>::iterator it =
	Nonboundary_elements_with_node_on_lower_disk_surface_pt.begin();
      it != Nonboundary_elements_with_node_on_lower_disk_surface_pt.end(); it++)
  {
    (*it)->output(some_file, nplot);
  }

  some_file.close();

  // -------------------------------------------------------------
  // output the dodgy stuff
  // -------------------------------------------------------------
  if(!dodgy_upper_element_pt.empty())
  {

    
    sprintf(filename, "%s/dodgy_upper_elements.dat", Doc_info.directory().c_str());
    some_file.open(filename);
   
    for(typename Vector<ELEMENT*>::iterator dodgy_it = dodgy_upper_element_pt.begin();
	dodgy_it != dodgy_upper_element_pt.end(); dodgy_it++)
    {
      (*dodgy_it)->output(some_file, nplot);
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
    
  
    for(typename Vector<ELEMENT*>::iterator dodgy_it = dodgy_lower_element_pt.begin();
	dodgy_it != dodgy_lower_element_pt.end(); dodgy_it++)
    {
      (*dodgy_it)->output(some_file, nplot);
    }
    some_file.close();

    sprintf(filename, "%s/dodgy_lower_element_nodal_distances.dat", Doc_info.directory().c_str());
    some_file.open(filename);
   
    for(Vector<Vector<double> >::iterator dodgy_it = dodgy_lower_element_nodal_distance.begin();
	dodgy_it != dodgy_lower_element_nodal_distance.end(); dodgy_it++)
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
template<class ELEMENT>
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
	if(node_is_on_edge_of_disk)
	  new_upper_disk_node_pt = dynamic_cast<BoundaryNode<Node>*>(original_node_pt);
	else
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
      for(std::set<unsigned>::iterator boundary_it =
	    original_node_boundaries.begin();
	  boundary_it != original_node_boundaries.end(); boundary_it++)
      {
	unsigned ibound = *boundary_it;
	
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
	if (is_lower_disk_boundary && !node_is_on_edge_of_disk)
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
      if(!already_duplicated_this_node && !node_is_on_edge_of_disk)
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
  
  unsigned nlower_disk_nodes = lower_disk_nodes_set.size();
  unsigned nupper_disk_nodes = upper_disk_nodes_set.size();
  
  oomph_info << "Number of plate nodes before duplication: "
	     << nlower_disk_nodes << "\n";
  oomph_info << "Number of plate nodes after duplication: "
	     << nlower_disk_nodes + nupper_disk_nodes << "\n\n";

  //  exit(1);
  
  // Finally, probably need this since we've fiddled the nodes on the plate boundary
  Bulk_mesh_pt->setup_boundary_element_info();

  // @@@@@@@@@@@@@@@@@@@@@@@@@
  // QUEHACERES debug

  char filename[100];
  ofstream some_file;
  unsigned nplot = 2;
  
  sprintf(filename, "%s/elements_on_duplicated_boundary.dat", Doc_info.directory().c_str());
  some_file.open(filename);
    
  for(unsigned ibound = First_upper_disk_boundary_id;
      ibound <= Last_upper_disk_boundary_id; ibound++)
  {
    unsigned nel = Bulk_mesh_pt->nboundary_element(ibound);
    for(unsigned e=0; e<nel; e++)
    {
      ELEMENT* el_pt = dynamic_cast<ELEMENT*>(Bulk_mesh_pt->boundary_element_pt(ibound,e));
      el_pt->output(some_file, nplot);
    }
  }

  some_file.close();

  oomph_info << "\nFirst_upper_disk_boundary_id: " << First_upper_disk_boundary_id 
	     << "\nLast_upper_disk_boundary_id:  " << Last_upper_disk_boundary_id << "\n\n";
  
  double t_end = TimingHelpers::timer();
  oomph_info << "Time to duplicate nodes and add boundaries: " << t_end - t_start << "s\n\n";

  // regenerate the node to element look-up
  generate_node_to_element_map();

  // ----------------------
  // QUEHACERES check that we don't have any nodes which aren't attached to elements

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

  sprintf(filename, "%s/duplicated_node_numbers.dat",
	  Doc_info.directory().c_str());
  some_file.open(filename);

  for(unsigned j=0; j<Bulk_mesh_pt->nnode(); j++)
  {
    Node* node_pt = Bulk_mesh_pt->node_pt(j);
    
    for(std::map<Node*,Node*>::iterator it=existing_duplicate_node_pt.begin();
	it != existing_duplicate_node_pt.end(); it++)
    {
      if(node_pt == it->second)
      {
	some_file << j << " "
		  << node_pt->x(0) << " "
		  << node_pt->x(1) << " "
		  << node_pt->x(2) << "\n";

	break;
      }
    }
  }
  
  some_file.close();
  
  bool first_boundary_without_nodes = true;
  unsigned id_of_first_boundary_without_nodes = 0;

  // @@@@ QUEHACERES @@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
  
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
	id_of_first_boundary_without_nodes = b;
	first_boundary_without_nodes = false;
      }
    }
  }

  oomph_info << "Number of new boundaries with no nodes: "
	     << new_boundaries_with_no_nodes << "\n"
	     << "Number of new boundaries with some nodes: "
	     << new_boundaries_with_some_nodes << "\n"
	     << "ID of first boundary without nodes:          "
	     << id_of_first_boundary_without_nodes << "\n\n";

  // for debug, let's output the number of elements touching the uppper plate
  // which are on lower boundaries (this should just be the edge nodes).

  sprintf(filename, "%s/upper_element_nodes_on_lower_disk_boundary.dat",
	  Doc_info.directory().c_str());
  some_file.open(filename);
  
  for(typename std::set<ELEMENT*>::iterator el_it = Elements_on_upper_disk_surface_pt.begin();
      el_it != Elements_on_upper_disk_surface_pt.end(); el_it++)
  {
    ELEMENT* el_pt = *el_it;

    for(unsigned j=0; j<el_pt->nnode(); j++)
    {
      Node* node_pt = el_pt->node_pt(j);

      for(unsigned b=First_lower_disk_boundary_id; b<=Last_lower_disk_boundary_id; b++)
      {
	if(node_pt->is_on_boundary(b))
	{
	  for(unsigned i=0; i<3; i++)	    
	    some_file << node_pt->x(i) << " ";
	  
	  some_file << std::endl;
	}
      }
    }
  }
  some_file.close();

  sprintf(filename, "%s/upper_element_nodes_on_upper_disk_boundary.dat",
	  Doc_info.directory().c_str());
  some_file.open(filename);
  
  for(typename std::set<ELEMENT*>::iterator it = Elements_on_upper_disk_surface_pt.begin();
      it != Elements_on_upper_disk_surface_pt.end(); it++)
  {
    ELEMENT* el_pt = *it;

    for(unsigned j=0; j<el_pt->nnode(); j++)
    {
      Node* node_pt = el_pt->node_pt(j);

      for(unsigned b=First_upper_disk_boundary_id; b<=Last_upper_disk_boundary_id; b++)
      {
	if(node_pt->is_on_boundary(b))
	{
	  for(unsigned i=0; i<3; i++)	    
	    some_file << node_pt->x(i) << " ";
	  
	  some_file << std::endl;
	}
      }
    }
  }
  some_file.close();

  sprintf(filename, "%s/upper_boundary_nodes_on_lower_disk_boundary.dat",
	  Doc_info.directory().c_str());
  some_file.open(filename);
  
  for(unsigned b_upper=First_upper_disk_boundary_id; b_upper<=Last_upper_disk_boundary_id; b_upper++)
  {
    for(unsigned j=0; j<Bulk_mesh_pt->nboundary_node(b_upper); j++)
    {
      Node* node_pt = Bulk_mesh_pt->boundary_node_pt(b_upper,j);

      for(unsigned b=First_lower_disk_boundary_id; b<=Last_lower_disk_boundary_id; b++)
      {
	if(node_pt->is_on_boundary(b))
	{
	  for(unsigned i=0; i<3; i++)	    
	    some_file << node_pt->x(i) << " ";
	  
	  some_file << std::endl;
	}
      }
    }
  }
  some_file.close();

  sprintf(filename, "%s/upper_boundary_nodes.dat",
	  Doc_info.directory().c_str());
  some_file.open(filename);
  
  for(unsigned b_upper=First_upper_disk_boundary_id; b_upper<=Last_upper_disk_boundary_id; b_upper++)
  {
    for(unsigned j=0; j<Bulk_mesh_pt->nboundary_node(b_upper); j++)
    {
      Node* node_pt = Bulk_mesh_pt->boundary_node_pt(b_upper,j);

      for(unsigned i=0; i<3; i++)	    
	some_file << node_pt->x(i) << " ";
	  
      some_file << std::endl;
    }
  }
  some_file.close();

  sprintf(filename, "%s/duplicated_nodes_on_upper_boundary.dat",
	  Doc_info.directory().c_str());
  some_file.open(filename);

  for(std::map<Node*,Node*>::iterator it=existing_duplicate_node_pt.begin();
      it != existing_duplicate_node_pt.end(); it++)
  {
    Node* node_pt = it->second;
    
    for(unsigned b=First_upper_disk_boundary_id; b<=Last_upper_disk_boundary_id; b++)
    {
      if(node_pt->is_on_boundary(b))
      {
	for(unsigned i=0; i<3; i++)	    
	  some_file << node_pt->x(i) << " ";
	  
	some_file << std::endl;
      }
    }
  }
  
  some_file.close();

  sprintf(filename, "%s/duplicated_nodes_on_lower_boundary.dat",
	  Doc_info.directory().c_str());
  some_file.open(filename);
  
  for(std::map<Node*,Node*>::iterator it=existing_duplicate_node_pt.begin();
      it != existing_duplicate_node_pt.end(); it++)
  {
    Node* node_pt = it->second;
    
    for(unsigned b=First_lower_disk_boundary_id; b<=Last_lower_disk_boundary_id; b++)
    {
      if(node_pt->is_on_boundary(b))
      {
	for(unsigned i=0; i<3; i++)	    
	  some_file << node_pt->x(i) << " ";
	  
	some_file << std::endl;
      }
    }
  }
  
  some_file.close();

  sprintf(filename, "%s/upper_boundary_nodes_from_map.dat", Doc_info.directory().c_str());
  some_file.open(filename);
    
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
      Node* node_pt = el_pt->node_pt(set_it->second);
      for(unsigned b=First_upper_disk_boundary_id; b<=Last_upper_disk_boundary_id; b++)
      {
	if(node_pt->is_on_boundary(b))
	{
      
	  for(unsigned i=0; i<3; i++)
	  {
	    some_file << node_pt->x(i) << " ";
	  }
	  some_file << std::endl;
	}
      }
    }
  }

  some_file.close();
  
  sprintf(filename, "%s/lower_boundary_nodes_from_map.dat", Doc_info.directory().c_str());
  some_file.open(filename);
    
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
      Node* node_pt = el_pt->node_pt(set_it->second);
      for(unsigned b=First_lower_disk_boundary_id; b<=Last_lower_disk_boundary_id; b++)
      {
	if(node_pt->is_on_boundary(b))
	{
      
	  for(unsigned i=0; i<3; i++)
	  {
	    some_file << node_pt->x(i) << " ";
	  }
	  some_file << std::endl;
	}
      }
    }
  }

  some_file.close();

  sprintf(filename, "%s/upper_element_nodes_from_map.dat", Doc_info.directory().c_str());
  some_file.open(filename);
    
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

      for(unsigned j=0; j<el_pt->nnode(); j++)
      {
	Node* node_pt = el_pt->node_pt(j);
	for(unsigned b=First_lower_disk_boundary_id; b<=Last_lower_disk_boundary_id; b++)
	{
	  if(node_pt->is_on_boundary(b))
	  {
      
	    for(unsigned i=0; i<3; i++)
	    {
	      some_file << node_pt->x(i) << " ";
	    }
	    some_file << std::endl;
	  }
	}
      }
    }
  }

  some_file.close();
  
  // QUEHACERES check the sync between Disk_node_to_upper_disk_element_and_index_map
  // and Elements_on_upper_disk_surface_pt
  
  std::set<ELEMENT*> unique_elements_from_map;
  
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

  if(unique_elements_from_map == Elements_on_upper_disk_surface_pt)
    oomph_info << "\n----------\nsets are equal\n-----------\n";
  else
    oomph_info << "\n----------\nsets are NOT equal\n-----------\n";
  
  // @@@@@@@@@@@@@@@@@@@@@@@
  
  // QUEHACERES delete
  // unsigned n_edge_nodes = 0;
  
  // for(std::map<Node*, Node*>::iterator it = existing_duplicate_node_pt.begin();
  //     it != existing_duplicate_node_pt.end(); it++)
  // {
  //   Node* node_pt = it->second;
    
  //   // compute the x-y radius
  //   double x = node_pt->x(0);
  //   double y = node_pt->x(1);      
  //   double r = sqrt(x*x + y*y);

  //   // tolerance on the outer radius
  //   double tol = 1e-8;
      
  //   bool node_is_on_edge_of_disk = (abs(1-r) < tol);

  //   if(node_is_on_edge_of_disk)
  //   {
  //     n_edge_nodes++;
  //     continue;
  //   }
  //   node_pt->x(2) += 0.01;
  // }

  // oomph_info << "n_edge_nodes: " << n_edge_nodes << "\n\n";
  // @@@@@@@@@@@@@@@@@@@@@@@@
}

//========================================================================
/// Complete problem setup
//========================================================================
template<class ELEMENT>
void FlowAroundDiskProblem<ELEMENT>::complete_problem_setup()
{
  // QUEHACERES for now just use the bulk mesh
  // // Create the mesh as Geom Object
  //Face_mesh_as_geom_object_pt = new MeshAsGeomObject(Face_mesh_for_upper_disk_pt);
  Face_mesh_as_geom_object_pt = new MeshAsGeomObject(Bulk_mesh_pt);
  
  unsigned n_sample = 100;

  Vector<double> x(3, 0.0);
  Radial_sample_point_pt.resize(n_sample);
 
  for (unsigned j=0; j<n_sample; j++)
  {
    // radial line along phi=0
    double phi = 0;
    double r = j*1.0/n_sample;
    
    // height of warped disk
    double z = Global_Parameters::epsilon * r*r*cos(Global_Parameters::n*phi);
    
    Radial_sample_point_pt[j].second.resize(3);
    
    x[0] = r;
    x[1] = 0;
    x[2] = z;
    
    // Get the element and its local coordinates
    Face_mesh_as_geom_object_pt->locate_zeta(x,
					     Radial_sample_point_pt[j].first,
					     Radial_sample_point_pt[j].second);
  }
  
  // Apply bcs  
  apply_boundary_conditions();
}

//==start_of_apply_bc=====================================================
/// Helper function to apply boundary conditions
//========================================================================
template<class ELEMENT>
void FlowAroundDiskProblem<ELEMENT>::apply_boundary_conditions()
{  
  ofstream pin_file;
  pin_file.open("pinned_nodes.dat");

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
    for (unsigned inod=0; inod<num_nod; inod++)
    {
      // grab a pointer to this boundary node
      Node* node_pt = Bulk_mesh_pt->boundary_node_pt(ibound,inod);
      
      node_pt->pin(0);
      node_pt->pin(1);
      node_pt->pin(2);
      
      Vector<double> x(3);
      x[0] = node_pt->x(0);
      x[1] = node_pt->x(1);
      x[2] = node_pt->x(2);

      // set plate velocity
      if (( (ibound >= First_lower_disk_boundary_id) &&
  	  (ibound <= Last_lower_disk_boundary_id) ) ||
	  ( (ibound >= First_upper_disk_boundary_id) &&
	    (ibound <= Last_upper_disk_boundary_id) ) ) 
      {
	
  	node_pt->set_value(0, Global_Parameters::disk_velocity[0]);
  	node_pt->set_value(1, Global_Parameters::disk_velocity[1]);
  	node_pt->set_value(2, Global_Parameters::disk_velocity[2]);
      }
      else
      {
  	// outer boundaries
  	node_pt->set_value(0, 0.0);
  	node_pt->set_value(1, 0.0);
  	node_pt->set_value(2, 0.0);
      }
      
      pin_file << x[0] << " " 
  	       << x[1] << " " 
  	       << x[2] << " " 
  	       << std::endl;
    }
  }
 
  pin_file.close();



} // end set bc

//========================================================================
/// Doc the solution
//========================================================================
template<class ELEMENT>
void FlowAroundDiskProblem<ELEMENT>::doc_solution(const unsigned& nplot)
{ 


  bool do_bulk_output=true;
  if (CommandLineArgs::command_line_flag_has_been_set("--suppress_bulk_output"))
  {
    do_bulk_output=false;
  }

  ofstream some_file;
  ofstream some_file2;
  ofstream face_some_file;
  ofstream coarse_some_file;
  char filename[100];


  // Doc mesh quality (Ratio of max. edge length to min. height,
  /// so if it's very large it's BAAAAAD)
  sprintf(filename,"%s/mesh_quality%i.dat",
	  Doc_info.directory().c_str(),
	  Doc_info.number()); 
  ofstream quality_file;
  quality_file.open(filename);
  if (do_bulk_output) Bulk_mesh_pt->assess_mesh_quality(quality_file);
  quality_file.close();

 
  // Output elements adjacent to outer boundary
  //-------------------------------------------
  sprintf(filename,"%s/elements_next_to_outer_boundary%i.dat",
	  Doc_info.directory().c_str(),
	  Doc_info.number());
  some_file.open(filename);
  for (unsigned ibound = First_boundary_id_for_outer_boundary;
       ibound < First_boundary_id_for_outer_boundary+6; ibound++)
  {
    unsigned n_el = Bulk_mesh_pt->nboundary_element(ibound);
    for (unsigned e=0; e<n_el; e++)
    {
      if (do_bulk_output) 
      {
	Bulk_mesh_pt->boundary_element_pt(ibound,e)->
	  output(some_file,nplot);
      }
    }
  }
  some_file.close();


  // Output boundary coordinates on outer boundary
  //-----------------------------------------------
  sprintf(filename,"%s/boundary_coordinates_outer_boundary%i.dat",
	  Doc_info.directory().c_str(),
	  Doc_info.number());
  some_file.open(filename);
  for (unsigned ibound = First_boundary_id_for_outer_boundary;
       ibound < First_boundary_id_for_outer_boundary+6; ibound++)
  {
    Bulk_mesh_pt->Mesh::template 
      doc_boundary_coordinates<ELEMENT>(ibound,some_file);
  }
  some_file.close();

  // Output boundary coordinates on outer boundary
  //-----------------------------------------------
  unsigned n_b = Bulk_mesh_pt->nboundary();
  oomph_info << "number of boundaries in bulk mesh: " << n_b << std::endl;
  sprintf(filename,"%s/boundary_coordinates%i.dat",
	  Doc_info.directory().c_str(),
	  Doc_info.number());
  some_file.open(filename);
  for (unsigned ibound=0; ibound<n_b; ibound++)
  {
    if (Bulk_mesh_pt->boundary_coordinate_exists(ibound))
    {
      Bulk_mesh_pt->Mesh::template 
	doc_boundary_coordinates<ELEMENT>(ibound,some_file);
    }
  }
  some_file.close();


  // Output boundaries
  //------------------
  sprintf(filename,"%s/boundaries%i.dat",Doc_info.directory().c_str(),
	  Doc_info.number());
  some_file.open(filename);
  Bulk_mesh_pt->output_boundaries(some_file);
  some_file.close();


  // Output volumes and areas
  //-------------------------
  std::ofstream volumes_and_areas_file;
  sprintf(filename,"%s/volumes%i.dat",Doc_info.directory().c_str(),
	  Doc_info.number());
  volumes_and_areas_file.open(filename);


  // Output bulk elements in torus region
  //-------------------------------------
  double volume_in_torus_region = 0.0;
  sprintf(filename,"%s/soln_in_torus_region%i.dat",Doc_info.directory().c_str(),
	  Doc_info.number());
  some_file.open(filename);
  unsigned region_id = Torus_region_id;
  unsigned n_el = Bulk_mesh_pt->nregion_element(region_id);
  for (unsigned e=0; e<n_el; e++)
  {
    if (do_bulk_output) 
    {
      Bulk_mesh_pt->region_element_pt(region_id,e)->output(some_file,nplot);
    }
    volume_in_torus_region += Bulk_mesh_pt->
      region_element_pt(region_id,e)->size();
  }
  some_file.close();

  // Output bulk elements in region 0
  //--------------------------------- 
  double volume_in_region0 = 0.0;
  sprintf(filename,"%s/soln_in_zero_region%i.dat",Doc_info.directory().c_str(),
	  Doc_info.number());
  some_file.open(filename);
  region_id = 0;
  n_el = Bulk_mesh_pt->nregion_element(region_id);
  for (unsigned e=0; e<n_el; e++)
  {
    if (do_bulk_output) 
    {
      Bulk_mesh_pt->region_element_pt(region_id,e)->output(some_file,nplot);
    }
    volume_in_region0 += Bulk_mesh_pt->region_element_pt(region_id,e)->size();
  }
  some_file.close();
 

  // Get total mesh volume
  double total_mesh_volume = 0.0;
  n_el = Bulk_mesh_pt->nelement();
  for (unsigned e=0; e<n_el; e++)
  {
    total_mesh_volume += Bulk_mesh_pt->finite_element_pt(e)->size();
  }

  // QUEHACERES come back to this 
  // // Check volumes:
  // oomph_info << "Error in total region volume balance: " <<
  //   abs(total_mesh_volume-(volume_in_above_disk_region+
  // 			   volume_in_below_disk_region+
  // 			   volume_in_torus_region+
  // 			   volume_in_region0))/total_mesh_volume*100.0 
  // 	     << " % " << std::endl;

  // oomph_info << "Error in above/below disk region volume balance: " <<
  //   abs(volume_in_above_disk_region-volume_in_below_disk_region)/
  //   volume_in_above_disk_region*100.0 << " % " << std::endl;
 

  sprintf(filename,"%s/face_elements_on_outer_boundary%i.dat",
	  Doc_info.directory().c_str(),
	  Doc_info.number());
  some_file.open(filename);
  for(unsigned b = First_boundary_id_for_outer_boundary;
      b < First_boundary_id_for_outer_boundary + 6; b++)
  {
    unsigned n_el = Bulk_mesh_pt->nboundary_element(b);
    for (unsigned e=0; e<n_el; e++)
    {
      FiniteElement* el_pt = Bulk_mesh_pt->boundary_element_pt(b,e);
     
      // What is the index of the face of the bulk element at the boundary
      int face_index = Bulk_mesh_pt->face_index_at_boundary(b,e);

      // Build the corresponding face element
      NavierStokesFaceElement<ELEMENT>* surface_element_pt =
      	new NavierStokesFaceElement<ELEMENT>(el_pt,face_index);

      surface_element_pt->output(some_file, nplot);

      delete surface_element_pt;
    }    
  }
  some_file.close();
  
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
      NavierStokesFaceElement<ELEMENT>* surface_element_pt =
      	new NavierStokesFaceElement<ELEMENT>(el_pt,face_index);
     
      // Get surface area
      torus_surface_area += surface_element_pt->size();
     
      // Output
      surface_element_pt->output(some_file, nplot);
     
      // ...and we're done!
      delete surface_element_pt;
    }
  }
  some_file.close();
  oomph_info << "Torus surface area: " <<  torus_surface_area << std::endl;
 
  // Attach face elements to part of disk inside torus
  //--------------------------------------------------
  sprintf(filename,"%s/face_elements_on_disk_in_torus%i.dat",
	  Doc_info.directory().c_str(),
	  Doc_info.number());
  some_file.open(filename);
  double disk_in_torus_surface_area = 0.0;
  region_id = Torus_region_id;
  unsigned nb = One_based_boundary_id_for_disk_within_torus.size();
  for (unsigned i=0; i<nb; i++)
  {
    unsigned b = One_based_boundary_id_for_disk_within_torus[i]-1;
    unsigned nel = Bulk_mesh_pt->nboundary_element_in_region(b,region_id);
    for (unsigned e=0; e<nel; e++)
    {
      FiniteElement* el_pt = 
	Bulk_mesh_pt->boundary_element_in_region_pt(b,region_id,e);
     
      // What is the index of the face of the bulk element at the boundary
      int face_index = Bulk_mesh_pt->
	face_index_at_boundary_in_region(b,region_id,e);

      // Build the corresponding surface power element
      NavierStokesFaceElement<ELEMENT>* surface_element_pt =
      	new NavierStokesFaceElement<ELEMENT>(el_pt, face_index);
     
      // Get surface area
      disk_in_torus_surface_area += surface_element_pt->size();
     
      // Output
      surface_element_pt->output(some_file, nplot);
     
      // ...and we're done!
      delete surface_element_pt;
    }
  }
  some_file.close();
  oomph_info << "Disk in torus surface area: "
	     <<  disk_in_torus_surface_area << std::endl;

 
  // Attach face elements to part of disk outside torus
  //--------------------------------------------------
  sprintf(filename,"%s/face_elements_on_disk_outside_torus%i.dat",
	  Doc_info.directory().c_str(),
	  Doc_info.number());
  some_file.open(filename);
  double disk_outside_torus_surface_area = 0.0;
  region_id = 0; 
  nb = One_based_boundary_id_for_disk_outside_torus.size();
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
      
      // Build the corresponding flux jump element
      NavierStokesFaceElement<ELEMENT>* surface_element_pt =
      	new NavierStokesFaceElement<ELEMENT>(el_pt, face_index);
     
      // Get surface area
      disk_outside_torus_surface_area += surface_element_pt->size();
     
      // Output      
      surface_element_pt->output(some_file, nplot);
     
      // ...and we're done!
      delete surface_element_pt;
    }
  }
  some_file.close();
  oomph_info << "Disk outside torus surface area: "
	     <<  disk_outside_torus_surface_area << std::endl;

  oomph_info << "Total surface area of disk with torus: "
	     <<  disk_in_torus_surface_area+disk_outside_torus_surface_area 
	     << std::endl;

  // QUEHACERES come back to this
  // // Doc volumes and areas
  // volumes_and_areas_file << volume_in_above_disk_region << " " 
  // 			 << volume_in_below_disk_region << " " 
  // 			 << volume_in_torus_region << " " 
  // 			 << total_mesh_volume << " " 
  // 			 << volume_in_region0 << " "
  // 			 << torus_surface_area << " " 
  // 			 << disk_in_torus_surface_area << " " 
  // 			 << disk_outside_torus_surface_area << " " 
  // 			 << disk_upper_layer_surface_area << " " 
  // 			 << disk_lower_layer_surface_area << " " 
  // 			 << free_standing_disk_surface_area << " " 
  // 			 << std::endl;
  // volumes_and_areas_file.close();

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
      NavierStokesFaceElement<ELEMENT>* surface_element_pt =
      	new NavierStokesFaceElement<ELEMENT>(el_pt, face_index);
     
      // Output
      surface_element_pt->output(some_file, nplot);
     
      // ...and we're done!
      delete surface_element_pt;
    }
  }
  some_file.close();

    // Attach face elements to upper disk boundaries
  //--------------------------------------------------
  sprintf(filename,"%s/face_elements_on_upper_disk%i.dat",
	  Doc_info.directory().c_str(),
	  Doc_info.number());
  some_file.open(filename);
 
  for (unsigned b=First_upper_disk_boundary_id; b<=Last_upper_disk_boundary_id; b++)
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
      NavierStokesFaceElement<ELEMENT>* surface_element_pt =
      	new NavierStokesFaceElement<ELEMENT>(el_pt, face_index);
     
      // Output
      surface_element_pt->output(some_file, nplot);
     
      // ...and we're done!
      delete surface_element_pt;
    }
  }
  some_file.close();
  
  // @@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
  // QUEHACERES debug
  sprintf(filename,"%s/nodal_values_on_disk_from_face_elements%i.dat",
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

      // Build the corresponding surface power element
      NavierStokesFaceElement<ELEMENT>* surface_element_pt =
      	new NavierStokesFaceElement<ELEMENT>(el_pt, face_index);
     
      // Output
      for(unsigned j=0; j<surface_element_pt->nnode(); j++)
      {
	Node* node_pt = surface_element_pt->node_pt(j);

	for(unsigned i=0; i<3; i++)
	{
	  some_file << node_pt->x(i) << " ";
	}
	for(unsigned i=0; i<node_pt->nvalue(); i++)
	{
	  some_file << node_pt->value(i) << " ";
	}
	some_file << std::endl;
      }
     
      // ...and we're done!
      delete surface_element_pt;
    }
  }
  some_file.close();

  // -------------------------------------------------------------
  sprintf(filename,
	  "%s/radial_solution_phi0%d.dat",
	  Doc_info.directory().c_str(), Doc_info.number() );
  some_file.open(filename);

  // global coordinates
  Vector<double> x(3);

  // loop over sample points
  for (Vector<std::pair<GeomObject*,Vector<double> > >::iterator it =
	 Radial_sample_point_pt.begin(); it != Radial_sample_point_pt.end(); it++)
  {
    // get element containing current sample point
    NavierStokesFaceElement<ELEMENT>* el_pt =
      dynamic_cast<NavierStokesFaceElement<ELEMENT>*>(it->first);

    // get local coordinate for this sample point
    Vector<double> s = it->second;

    // get global coordinate corresponding to this elements local coord
    el_pt->interpolated_x(s, x);

    // get interpolated solution at this local coordinate
    Vector<double> u(3);

    el_pt->interpolated_u_nst(s, u);
    
    for(unsigned i=0; i<3; i++)
    {      
      some_file << x[0] << " ";
    }

    // add the pressure
    u.push_back(el_pt->interpolated_p_nst(s));

    // output
    for(unsigned i=0; i<4; i++)
    {
      some_file << u[i] << " ";
    }
    some_file << std::endl;
  }
  some_file.close();
  
  // @@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
  
  // -------------------------------------------------------------
  // output the main map to check we've assembled it correctly
  // -------------------------------------------------------------

  sprintf(filename, "%s/upper_disk_elements_from_nodal_map.dat", Doc_info.directory().c_str());
  some_file.open(filename);

  for (unsigned ibound = First_lower_disk_boundary_id;
       ibound <= Last_lower_disk_boundary_id; ibound++)
  {     
    unsigned nboundary_node = Bulk_mesh_pt->nboundary_node(ibound);
    for(unsigned n=0; n<nboundary_node; n++)
    {
      // get a pointer to this boundary node
      Node* node_pt = Bulk_mesh_pt->boundary_node_pt(ibound, n);

      // get the set of all upper disk elements associated with this node
      std::set<std::pair<ELEMENT*, unsigned> > elem_set =
	Disk_node_to_upper_disk_element_and_index_map[node_pt];
	
      typename std::set<std::pair<ELEMENT*, unsigned> >::iterator it;
  
      for(it = elem_set.begin(); it != elem_set.end(); it++)
      {
	// get the element and the nodal index
	std::pair<ELEMENT*, unsigned> element_index_pair = *it;

	// output the element
	element_index_pair.first->output(some_file, nplot);
      }
    }
  }

  some_file.close();

  // Output solution
  //----------------
  sprintf(filename,"%s/soln%i.dat",Doc_info.directory().c_str(),
	  Doc_info.number());
  some_file.open(filename);
  if (do_bulk_output) Bulk_mesh_pt->output(some_file,nplot);
  some_file.close();

  // Output solution showing element outlines
  //-----------------------------------------
  sprintf(filename,"%s/coarse_soln%i.vtu",Doc_info.directory().c_str(),
	  Doc_info.number());
  some_file.open(filename);
  if (do_bulk_output) Bulk_mesh_pt->output_paraview(some_file,2);
  some_file.close();


  // Get norm of solution
  //---------------------
  sprintf(filename,"%s/norm%i.dat",Doc_info.directory().c_str(),
	  Doc_info.number());
  some_file.open(filename);
  double norm_soln = 0.0;
  Bulk_mesh_pt->compute_norm(norm_soln);  
  some_file << sqrt(norm_soln) << std::endl;
  oomph_info << "Norm of computed solution: "   << sqrt(norm_soln)  << endl;
  some_file.close();

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
      
    jac.sparse_indexed_output(some_file);

    some_file.close();

    oomph_info << "Output sparse jacobian matrix to " << filename << "\n\n";
  }

  if (CommandLineArgs::command_line_flag_has_been_set("--describe_dofs"))
  {
    sprintf(filename,"%s/describe_dofs.dat", Doc_info.directory().c_str());
    some_file.open(filename);
    describe_dofs(some_file);
    
    some_file.close();

    oomph_info << "Output description of dofs\n";
  }
  if (CommandLineArgs::command_line_flag_has_been_set("--describe_nodes"))
  {
    sprintf(filename,"%s/describe_nodes.dat", Doc_info.directory().c_str());
    some_file.open(filename);

    for(unsigned j=0; j<Bulk_mesh_pt->nnode(); j++)
    {
      // grab the node
      Node* node_pt = Bulk_mesh_pt->node_pt(j);

      // get it's coordinates
      double x = node_pt->x(0);
      double y = node_pt->x(1);
      double z = node_pt->x(2);

      some_file << j << " " << x << " " << y << " " << z << " " << node_pt << std::endl;
    }
    
    some_file.close();
    oomph_info << "Output description of nodes\n";
  }
  
  //Increment counter for solutions 
  Doc_info.number()++;
  
} // end of doc


//========================================================================
/// Driver
//========================================================================
int main(int argc, char* argv[])
{  
  // set up the multi-processor interface
  MPI_Helpers::init(argc,argv);
  
  // Store command line arguments
  CommandLineArgs::setup(argc,argv);
  
  // length of downstream region occupied by impedance elements
  CommandLineArgs::specify_command_line_flag("--suppress_bulk_output");

  CommandLineArgs::specify_command_line_flag("--describe_dofs");
  CommandLineArgs::specify_command_line_flag("--describe_nodes");
  
  // rigid body velocity of the plate
  CommandLineArgs::specify_command_line_flag("--velocity_x",
					     &Global_Parameters::disk_velocity[0]);
  CommandLineArgs::specify_command_line_flag("--velocity_y",
					     &Global_Parameters::disk_velocity[1]);
  CommandLineArgs::specify_command_line_flag("--velocity_z",
					     &Global_Parameters::disk_velocity[2]);

  // amplitude of the warping of the disk
  CommandLineArgs::specify_command_line_flag("--epsilon",
					     &Global_Parameters::epsilon);

  // wavenumber of the warping of the disk
  CommandLineArgs::specify_command_line_flag("--n",
					     &Global_Parameters::n);
  
  // get output directory
  CommandLineArgs::specify_command_line_flag("--dir", &Global_Parameters::output_directory);
 
#ifndef DO_TETGEN

  // Gmsh command line invocation
  CommandLineArgs::specify_command_line_flag
    ("--gmsh_command_line",
     &Global_Parameters::Gmsh_command_line_invocation);

#endif

  // Parse command line
  CommandLineArgs::parse_and_assign(); 
 
  // Doc what has actually been specified on the command line
  CommandLineArgs::doc_specified_flags();

#ifndef DO_TETGEN

  // Are you suicidal?
  if (!CommandLineArgs::command_line_flag_has_been_set("--gmsh_command_line"))
  {
    std::string error_msg
      ("You haven't specified how gmsh is invoked on the command line\n");
    error_msg += "on your computer, so I'll use the default\n\n" + 
      Global_Parameters::Gmsh_command_line_invocation
      + "\n\nwhich, unless you're mheil, is unlikely to work and I will "
      + "now die...\n";
    throw OomphLibError(error_msg, 
			OOMPH_CURRENT_FUNCTION,
			OOMPH_EXCEPTION_LOCATION);
  }

#endif


  
  // Note that this can make tetgen die!
  //feenableexcept(FE_INVALID | FE_DIVBYZERO | FE_OVERFLOW | FE_UNDERFLOW);

  // Shut up prefix
  oomph_info.output_modifier_pt() = &default_output_modifier;

  // Number of output points per edge
  // QUEHACERES 2 for mesh debug
  unsigned nplot = 2; // 5;

  // Build problem
  FlowAroundDiskProblem <ProjectableTaylorHoodElement<TTaylorHoodElement<3> > > problem;

  // // QUEHACERES for debug
  // problem.newton_solver_tolerance() = 5e-8;
  
  //Output initial guess
  problem.doc_solution(nplot);
  
  unsigned max_adapt = 0; 
  for (unsigned i=0; i<=max_adapt; i++)
  {
    // Solve the bastard!
    problem.newton_solve();

    //Output solution
    problem.doc_solution(nplot);

    if (i != max_adapt)
    {
      problem.adapt();
    }
  }

  // Shut down oomph-lib's MPI
  MPI_Helpers::finalize();

  return 0;
}


