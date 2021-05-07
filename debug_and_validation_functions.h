// main Problem class definition
#include "flow_around_disk_problem.h"

#include "external_src/oomph_superlu_4.3/slu_ddefs.h"

template <class ELEMENT>
void FlowAroundDiskProblem<ELEMENT>::read_and_recompute_singular_drag(
  const std::string& filename) const
{
  oomph_info << "Reading singular amplitudes..." << std::endl;
  std::ifstream some_file;
  some_file.open(filename);

  std::string line;
  
  unsigned line_counter = 0;
  
  while(std::getline(some_file, line))
  {
    // ignore zone info
    if(line_counter % 6 == 0)
    {
      line_counter++;
      continue;
    }
    // skip every other line (assumes nplot=5)
    if(line_counter % 2 == 0)
    {
      line_counter++;
      continue;
    }
    
    std::istringstream iss(line);

    // read the data in
    Vector<double> amplitude_data(12, 0.0);
    for(unsigned i=0; i<12; i++)
    {
      iss >> amplitude_data[i];
    }

    unsigned elem_number = line_counter / 6;
    FiniteElement* elem_pt = dynamic_cast<FiniteElement*>(
      Singular_fct_element_mesh_pt->element_pt(elem_number));

    unsigned node_number = (line_counter % 6) / 2;
    Node* node_pt = elem_pt->node_pt(node_number);

    /* oomph_info << "Element number: " << elem_number << ", " */
    /* 	       << "Node number: " << node_number << std::endl; */
    for(unsigned i=0; i<4; i++)
      node_pt->set_value(i, amplitude_data[4+i]);

    line_counter++;
  }
  
  some_file.close();

  oomph_info << "Re-computing singular drag..." << std::endl;
  
  // compute the contribution of the singular functions to the
  // total drag
  Vector<double> sing_drag_upper(Dim, 0.0), sing_drag_lower(Dim, 0.0),
    sing_torque_upper(Dim, 0.0), sing_torque_lower(Dim, 0.0);

  Vector<double> centre_of_mass(Dim, 0.0);
  
  Singular_drag_elements[0]->compute_total_singular_drag_and_torque(centre_of_mass,
								    sing_drag_upper,
								    sing_torque_upper);
  
  Singular_drag_elements[1]->compute_total_singular_drag_and_torque(centre_of_mass,
								    sing_drag_lower,
								    sing_torque_lower);

  oomph_info << "Total singular drag: ";
  for(unsigned i=0; i<3; i++)
  {    
    oomph_info << sing_drag_upper[i] + sing_drag_lower[i] << " ";
  }
  oomph_info << std::endl;
}


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

//========================================================================
//== start of set_values_to_singular_solution ============================
/// Function to assign the singular solution to all nodes of the mesh
//========================================================================
template <class ELEMENT>
void FlowAroundDiskProblem<ELEMENT>::set_values_to_singular_solution() const
{  
  unsigned nel = Torus_region_mesh_pt->nelement();
  for(unsigned e=0; e<nel; e++)
  {
    ELEMENT* elem_pt = dynamic_cast<ELEMENT*>(Torus_region_mesh_pt->element_pt(e));

    // loop over the nodes in this element
    for(unsigned j=0; j<elem_pt->nnode(); j++)
    {
      // get a pointer to this node
      Node* node_pt = elem_pt->node_pt(j);
      
      // get Cartesian coordinates of this node
      Vector<double> x(3, 0.0);

      for(unsigned i=0; i<3; i++)
	x[i] = node_pt->x(i);

      LagrangianCoordinates lagr_coords;
      std::pair<GeomObject*, Vector<double> > line_element_and_local_coord;
      
      // now get the element which computes the scaled singular functions
      get_lagr_coordinates_and_singular_element(elem_pt,
						x,
						lagr_coords,
						line_element_and_local_coord);

      // cast the GeomObject to a singular line element      
      auto sing_el_pt = dynamic_cast<
	ScalableSingularityForNavierStokesLineElement<SINGULAR_ELEMENT_NNODE_1D>*>
	(line_element_and_local_coord.first);

      // local coordinate in the singular element for the zeta of this plot point
      Vector<double> s_singular_el = line_element_and_local_coord.second;

      // get the total singular contribution at this point
      Vector<double> u_sing_total =
	sing_el_pt->total_singular_contribution(lagr_coords,
						s_singular_el);

      for(unsigned i=0; i<3; i++)
	node_pt->set_value(i, u_sing_total[i]);

      // are we at a vertex where pressure is stored?
      if(j<4)
	node_pt->set_value(3, u_sing_total[3]);
      
    } // end loop over nodes
  }
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
      for(unsigned i=0; i<3; i++)
	x[i] = node_pt->x(i);
      
      // get the total solution at this point (in the bulk, so
      // the FE solution is the total solution)
      Vector<double> u(4, 0.0);
      FlatDiskExactSolutions::total_exact_solution(x, Global_Parameters::u_disk_rigid_translation, 
						   Global_Parameters::omega_disk,
						   Global_Parameters::rotation_matrix,
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
      
      LagrangianCoordinates lagr_coords;
      std::pair<GeomObject*, Vector<double> > line_element_and_local_coord;
      
      // now get the element which computes the scaled singular functions
      get_lagr_coordinates_and_singular_element(elem_pt,
						x,
						lagr_coords,
						line_element_and_local_coord);

      // cast the GeomObject to a singular line element      
      auto sing_el_pt = 
	dynamic_cast<SingularLineElement*>(line_element_and_local_coord.first);

      // local coordinate in the singular element for the zeta of this plot point
      Vector<double> s_singular_el = line_element_and_local_coord.second;

      // now account for any slight effects of the non-infinite cuvature
      // which might cause the sign of z to flip from the lower to the upper disk
      // by catching the case where the Lagrangian xi3 coordinate is negative
      // but z has sneaked over 0
      if(x[2] > 0 && lagr_coords.sign_of_xi3() < 0)
	x[2] = -1e-9;
      
      // get the total solution at this point (in the bulk, so
      // the FE solution is the total solution)
      Vector<double> u_total(4, 0.0);
      FlatDiskExactSolutions::total_exact_solution(x, Global_Parameters::u_disk_rigid_translation, 
						   Global_Parameters::omega_disk,
						   Global_Parameters::rotation_matrix,
						   u_total);

      
      // get the total singular contribution at this point
      Vector<double> u_sing_total =
	sing_el_pt->total_singular_contribution(lagr_coords,
						s_singular_el);

      // QUEHACERES leave velocities at zero for now, i.e. assuming we've
      // got the correct singular amplitude
      /* // set the velocities to the non-singular bit */
      for(unsigned i=0; i<3; i++)
      	node_pt->set_value(i, u_total[i] - u_sing_total[i]);

      // set the pressure if it's there (if there are only velocities and
      // corresponding Lagrange multipliers there will be an odd number of values)

      // vertices are enumerated first, so nodes j <= Dim have pressures
      if(j <= Dim)
      {
	// catch the infinity case
	if( (abs(lagr_coords.xi1 - 1) < 1e-4) && (abs(lagr_coords.xi3) < 1e-6))	  
	{
	  double infinity = 10;
	  node_pt->set_value(3, 0.0); // infinity); //  * cos(lagr_coords.zeta));
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
void FlowAroundDiskProblem<ELEMENT>::validate_singular_stress() const
{  
  oomph_info << "\nValidating singular stress...\n"
	     << "-----------------------------\n" << std::endl;

  double t_start = TimingHelpers::timer();

  // set the singular amplitudes to their analytic values
  impose_fake_singular_amplitude();

  // QUEHACERES delete?
  set_values_to_singular_solution();
  
  char filename[500];

  sprintf(filename, "%s/error_in_singular_stress.dat",
	  Doc_info.directory().c_str());
  
  // open the output file to record the error (tecplot format)
  ofstream stress_error_output(filename);
  
  // open the output file to record the error (plain format)
 
  sprintf(filename, "%s/error_in_singular_stress_plain.dat",
	  Doc_info.directory().c_str());
  
  ofstream stress_error_output_plain(filename);

  // column headers
  stress_error_output << "x,y,z,err_xx,err_xy,err_xz,err_yx,err_yy,err_yz,"
		      << "err_zx,err_zy,err_zz,sing_xx,sing_xy,sing_xz,"
		      << "sing_yx,sing_yy,sing_yz,sing_zx,sing_zy,sing_zz,"
		      <<"fd_xx,fd_xy,fd_xz,fd_yx,fd_yy,fd_yz,fd_zx,fd_zy,fd_zz"
		      << "p_sing";
  
  // number of plot points per side
  unsigned nplot = Global_Parameters::Nplot_for_bulk;

  oomph_info << "Computing singular analytic and finite-diff stress...\n";

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

      std::pair<GeomObject*, Vector<double> > line_elem_and_local_coord = 
	elem_pt->line_element_and_local_coordinate_at_plot_point(iplot);

      auto sing_el_pt = dynamic_cast<
	ScalableSingularityForNavierStokesLineElement<SINGULAR_ELEMENT_NNODE_1D>*>
	(line_elem_and_local_coord.first);

      // get the Lagrangian (xi1,xi2,xi3) coordinates at this knot
      LagrangianCoordinates lagr_coords = elem_pt->lagr_coordinate_at_plot_point(iplot);
      
      DenseMatrix<double> stress_sing(3,3, 0.0);

      // get the singular stress from the line element
      sing_el_pt->total_singular_stress(lagr_coords, line_elem_and_local_coord.second, stress_sing);
      
      Vector<double> u_exact(4, 0.0);
      Analytic_Functions::exact_solution_flat_disk(x, u_exact);

      DenseMatrix<double> du_dx_fd(3,3,0.0);

      const double fd_dx = 1e-8;
      
      // now compute the solution at a small increment in each Cartesian direction
      for(unsigned j=0; j<3; j++)
      {
	// coordinate vectors with increments applied to their components
	Vector<double> x_plus_dx = x;

	// add the increments in each direction
	x_plus_dx[j] += fd_dx;

	Vector<double> u_exact_plus_dx(3, 0.0);
	Analytic_Functions::exact_solution_flat_disk(x_plus_dx, u_exact_plus_dx);

	for(unsigned i=0; i<3; i++)
	  du_dx_fd(i,j) = (u_exact_plus_dx[i] - u_exact[i]) / fd_dx;
      }

      DenseMatrix<double> stress_exact(3,3,0.0);
      for(unsigned i=0; i<3; i++)
      {
	for(unsigned j=0; j<3; j++)
	{
	  if(i==j)
	    stress_exact(i,j) -= u_exact[3];

	  stress_exact(i,j) += du_dx_fd(i,j) + du_dx_fd(j,i);
	}
      }

      // -----------------------------------------
      // Error
      // -----------------------------------------
	
      // compute the error
      for(unsigned i=0; i<dim; i++)
      {
	for(unsigned j=0; j<dim; j++)
	{
	  // compute the error between the interpolated "FE" stress and the exact singular stress
	  double error = /* stress_fe(i,j) */ stress_exact(i,j) - stress_sing(i,j);

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
      	  stress_error_output       << stress_exact(i,j) /* stress_fe(i,j) */ << " ";
      	  stress_error_output_plain << stress_exact(i,j) /* stress_fe(i,j) */ << " ";
      	}
      }

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
Vector<std::pair<unsigned, double>> FlowAroundDiskProblem<ELEMENT>::
  compute_singular_amplitudes_from_disk_velocity(const double& zeta) const
{
  using namespace Global_Parameters;
    
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

  double c_broadside_rotation_no_az = omega_out_of_plane * cos(zeta);
  double c_broadside_rotation_az = omega_out_of_plane * sin(zeta);

  // ### QUEHACERES 
  /* // in-plane amplitude is the in-plane speed modulated by */
  /* // an in-phase function of the angle of this point relative to the */
  /* // translation vector */
  /* double c_in_plane = u_in_plane * sqrt(2) * cos(zeta) * cos(zeta_rotation); */
    
  // in-plane rotation amplitude is the in-plane speed modulated by
  // a function pi/2 out of phase with the angle of this point relative to the
  // translation vector
  double c_in_plane_rotation = u_in_plane * (-sin(zeta) +
					     cos(zeta)) * cos(zeta_translation) +
    omega_disk[2];

  double c_in_plane_az    = sin(zeta) * u_in_plane;
  double c_in_plane_no_az = cos(zeta) * u_in_plane;
    
  // now package them up for return
  Vector<std::pair<unsigned,double>> amplitudes(4);

  // QUEHACERES 
  /* amplitudes[0] = std::make_pair(Sing_fct_id_broadside, c_broadside); */
  amplitudes[0] = std::make_pair(Sing_fct_id_broadside_rotation_no_az,
				 c_broadside_rotation_no_az);
  amplitudes[1] = std::make_pair(Sing_fct_id_broadside_rotation_az_only,
				 c_broadside_rotation_az);
  amplitudes[2] = std::make_pair(Sing_fct_id_in_plane_no_az,
				 c_in_plane_no_az); // QUEHACERES c_in_plane;
  amplitudes[3] = std::make_pair(Sing_fct_id_in_plane_az_only,
				 c_in_plane_az); // QUEHACERES c_in_plane_rotation;
  
  return amplitudes;
}


//==start_of_impose_fake_singular_amplitude===============================
/// Set the singular amplitude to a prescribed value and bypass the proper calculation
//========================================================================
template <class ELEMENT>
void FlowAroundDiskProblem<ELEMENT>::
impose_fake_singular_amplitude(const bool& impose_zero_amplitude) const
{
  // tell all the elements in the singular element mesh about the fake
  // amplitude to impose
  for(unsigned e=0; e<Singular_fct_element_mesh_pt->nelement(); e++)
  {    
    // get a pointer to this singular line element in the upper mesh
    auto sing_el_pt =
      dynamic_cast<SingularLineElement*>(Singular_fct_element_mesh_pt->element_pt(e));

    // if we just want to zero out all the amplitudes, do it and don't mess around
    // computing anything else
    if(impose_zero_amplitude)
    {
      sing_el_pt->impose_zero_singular_amplitude();
      continue;
    }

    // map from singular function ID to a vector of the corresponding nodal amplitudes
    std::map<unsigned, Vector<double>> sing_id_to_nodal_amplitudes_map;
    
    // now loop over the nodes and modulate the in-plane amplitude
    // by the azimuthal angle of each
    for(unsigned j=0; j<sing_el_pt->nnode(); j++)
    {
      // get the zeta for this node
      double zeta = sing_el_pt->zeta_nodal(j,0,0);

      // compute the amplitudes, initially set to zero
      Vector<std::pair<unsigned, double>> ids_and_amplitudes;

      ids_and_amplitudes = compute_singular_amplitudes_from_disk_velocity(zeta);
      
      // add to the list
      for(std::pair<unsigned,double>& id_amplitude_pair : ids_and_amplitudes)
      {
	sing_id_to_nodal_amplitudes_map[id_amplitude_pair.first].push_back(id_amplitude_pair.second);
      }
    }

    // now tell the singular elements about them
    for(auto const& id_amplitudes_pair : sing_id_to_nodal_amplitudes_map)
    {
      sing_el_pt->impose_singular_fct_amplitude(id_amplitudes_pair.first,
						id_amplitudes_pair.second);
    }
  }  
}

//==start_of_pin_singular_function========================================
/// Set (and pin) a particular singular amplitude to a prescribed value
//========================================================================
template <class ELEMENT>
void FlowAroundDiskProblem<ELEMENT>::pin_singular_function(const unsigned& sing_fct_id,
							   const double& val)
{
  unsigned nel = Singular_fct_element_mesh_pt->nelement();

  for(unsigned e=0; e<nel; e++)
  {
    // get a pointer to this singular line element in the upper mesh
    auto sing_el_pt =
      dynamic_cast<SingularLineElement*>(Singular_fct_element_mesh_pt->element_pt(e));

    Vector<double> amplitude(sing_el_pt->nnode(), val);
    
    sing_el_pt->impose_singular_fct_amplitude(sing_fct_id, amplitude);
  }
}

template <class ELEMENT>
void FlowAroundDiskProblem<ELEMENT>::output_augmented_elements_eulerian_and_lagr_coords() const
{
  ofstream some_file;
  std::ostringstream filename;
  
  filename << Doc_info.directory() << "/torus_region_eulerian_and_lagr_coords.dat";
  
  some_file.open(filename.str().c_str());
    
  // do Torus region bulk elements
  // Bulk elements in torus region
  unsigned region_id = Torus_region_id;
  unsigned n_el = Bulk_mesh_pt->nregion_element(region_id);
  for (unsigned e=0; e<n_el; e++)
  {    
    ELEMENT* torus_region_el_pt = dynamic_cast<ELEMENT*>(
      Bulk_mesh_pt->region_element_pt(region_id, e));

    const unsigned nnode = torus_region_el_pt->nnode();
    
    Shape psi(nnode);

    Vector<double> s(3, 0.0);
      
    // number of integration points in this element
    const unsigned n_intpt = torus_region_el_pt->integral_pt()->nweight();
    
    for(unsigned ipt=0; ipt<n_intpt; ipt++)
    {          
      // get the shape functions at this knot
      torus_region_el_pt->shape_at_knot(ipt, psi);
      
      // get the coords at this integration point
      LagrangianCoordinates lagr_coords_at_knot =
	torus_region_el_pt->lagr_coordinate_at_knot(ipt);

      // Eulerian coordinates of this knot
      Vector<double> x(Dim, 0.0);
      
      // loop over the nodes and compute the coords at this knot
      for(unsigned l=0; l<nnode; l++) 
      {
	// grab a pointer to the current node
	Node* node_pt = torus_region_el_pt->node_pt(l);

	// interpolate the Eulerian coordinates
	for(unsigned i=0; i<3; i++)	  
	  x[i] += node_pt->x(i) * psi[l];
      }
      
      // output
      some_file << x[0] << " " << x[1] << " " << x[2] << " "
		<< lagr_coords_at_knot.xi1 << " "
		<< lagr_coords_at_knot.xi2 << " "
		<< lagr_coords_at_knot.xi3 << std::endl;
    }
      
  }

  some_file.close();

  filename.str("");
  
  filename << Doc_info.directory() << "/bc_element_eulerian_and_lagr_coords.dat";
  
  some_file.open(filename.str().c_str());
    
  n_el = Face_mesh_for_bc_pt->nelement();
  for (unsigned e=0; e<n_el; e++)
  {    
    NavierStokesWithSingularityBCFaceElement<ELEMENT>* bc_el_pt =
      dynamic_cast<NavierStokesWithSingularityBCFaceElement<ELEMENT>*>(
      Face_mesh_for_bc_pt->element_pt(e));

    const unsigned nnode = bc_el_pt->nnode();
    
    Shape psi(nnode);

    Vector<double> s(3, 0.0);
      
    // number of integration points in this element
    const unsigned n_intpt = bc_el_pt->integral_pt()->nweight();
    
    for(unsigned ipt=0; ipt<n_intpt; ipt++)
    {          
      // get the shape functions at this knot
      bc_el_pt->shape_at_knot(ipt, psi);
      
      // get the coords at this integration point
      LagrangianCoordinates lagr_coords_at_knot =
	bc_el_pt->lagr_coordinate_at_knot(ipt);

      // Eulerian coordinates of this knot
      Vector<double> x(Dim, 0.0);
      
      // loop over the nodes and compute the coords at this knot
      for(unsigned l=0; l<nnode; l++) 
      {
	// grab a pointer to the current node
	Node* node_pt = bc_el_pt->node_pt(l);

	// interpolate the Eulerian coordinates
	for(unsigned i=0; i<3; i++)	  
	  x[i] += node_pt->x(i) * psi[l];
      }
      
      // output
      some_file << x[0] << " " << x[1] << " " << x[2] << " "
		<< lagr_coords_at_knot.xi1 << " "
		<< lagr_coords_at_knot.xi2 << " "
		<< lagr_coords_at_knot.xi3 << std::endl;
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
							    Global_Parameters::rotation_matrix,
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
    get_total_disk_force_and_torque(centre_of_rotation_dummy,
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
