#ifndef OOMPH_RIGID_DISK_ELEMENTS_HEADER
#define OOMPH_RIGID_DISK_ELEMENTS_HEADER

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
      unsigned num_plot_points = this->nplot_points(nplot);
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

  void interpolated_x(const Vector<double>& s, Vector<double>& x) const
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

  void add_contribution_to_drag_and_torque(const Vector<double>& centre_of_mass,
					   Vector<double>& total_drag,
					   Vector<double>& total_torque) const
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

    //Loop over the integration points
    for(unsigned ipt=0; ipt<n_intpt; ipt++)
    {       
      //Assign values of s
      for(unsigned i=0; i<(Dim-1); i++)
      {
	s[i] = this->integral_pt()->knot(ipt,i);
      }

      // get the Eulerian coordinates at this knot for the torque calculation
      Vector<double> x(Dim, 0.0);
      this->interpolated_x(s, x);
      
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

      // compute the traction t_i = \tau_{ij} n_j
      Vector<double> traction(Dim, 0.0);
      for(unsigned i=0; i<Dim; i++)
      {
	for(unsigned j=0; j<Dim; j++)
	{
	  traction[i] += stress(i,j) * unit_normal[j];
	}
      }

      // now add the weighted traction to the total drag
      for(unsigned i=0; i<Dim; i++)
      {
	total_drag[i] += traction[i] * W;
      }
      
      // compute the radial vector of this knot from the CoM
      Vector<double> r(Dim, 0.0);
      for(unsigned i=0; i<Dim; i++)
	r[i] = x[i] - centre_of_mass[i];

      // now get the torque G = r^F
      Vector<double> torque(Dim, 0.0);
      VectorHelpers::cross(r, traction, torque);

      // add the contribution of the torque at this knot to the integral
      for(unsigned i=0; i<Dim; i++)
      {
	total_torque[i] += torque[i] * W;
      }      
    } // end loop over integration points
  }

  void contribution_to_centre_of_mass(Vector<double>& com) const
  {
    // Find out how many nodes there are
    const unsigned n_node = this->nnode();
    
    // shorthand
    const unsigned Dim = this->Dim;
    
    // Set up memory for the shape and test functions
    Shape psi(n_node), test(n_node);
     
    // Set the value of Nintpt
    const unsigned n_intpt = this->integral_pt()->nweight();
     
    // Set the Vector to hold local coordinates
    Vector<double> s(Dim-1, 0.0);

    // Set the Vector to hold local coordinates in the bulk element
    Vector<double> s_bulk(Dim, 0.0);

    // Eulerian coordinates
    Vector<double> x(Dim, 0.0);
	
    // Saves result of integration
    Vector<double> total_force(Dim, 0.0);

    // Get pointer to assocated bulk element
    ELEMENT* bulk_el_pt = dynamic_cast<ELEMENT*>(bulk_element_pt());
      
    //Loop over the integration points
    for(unsigned ipt=0; ipt<n_intpt; ipt++)
    {       
      //Assign values of s
      for(unsigned i=0; i<(Dim-1); i++)
      {
	s[i] = this->integral_pt()->knot(ipt,i);
      }

      this->get_local_coordinate_in_bulk(s, s_bulk);
     
      //Get x position from bulk
      bulk_el_pt->interpolated_x(s_bulk, x);
	
      //Get the integral weight
      double w = this->integral_pt()->weight(ipt);
       
      //Find the shape and test functions and return the Jacobian
      //of the mapping
      double J = this->J_eulerian(s); 
       
      //Premultiply the weights and the Jacobian
      double W = w*J;

      // add the integrated position vector contribution
      for(unsigned i=0; i<3; i++)
      {
	com[i] += x[i] * W;
      }
    }
  }
      
private:
  
  unsigned Dim;
};

#endif
