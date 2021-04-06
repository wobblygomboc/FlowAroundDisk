//LIC// ====================================================================
//LIC// This file forms part of oomph-lib, the object-oriented, 
//LIC// multi-physics finite-element library, available 
//LIC// at http://www.oomph-lib.org.
//LIC// 
//LIC//    Version 1.0; svn revision $LastChangedRevision: 1097 $
//LIC//
//LIC// $LastChangedDate: 2015-12-17 11:53:17 +0000 (Thu, 17 Dec 2015) $
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
// Header file for elements that are used to ...hierher
#ifndef OOMPH_STOKES_SING_FACE_ELEMENTS_HEADER
#define OOMPH_STOKES_SING_FACE_ELEMENTS_HEADER

#include "additional_maths.h"
#include "navier_stokes.h"
#include "poisson.h"
#include "chebyshev_gauss_integration.h"

namespace oomph
{  
  //----------------------ONE DIMENSIONAL CIRCULAR (LINE) MESH-----------------

  //============================================================================
  /// A simple one dimensional circular mesh: uniformly spaced nodes in a
  /// circle of radius 1
  //============================================================================
  template<class ELEMENT>
    class CircularLineMesh : public Mesh
  {

  public:

    /// Mesh Constructor. The argument is the desired number of elements
    CircularLineMesh(const unsigned& n_element)
    {
      // Resize the vector of pointers to elements: there are n_element elements
      Element_pt.resize(n_element); 

      // Construct the first element (Note the use of the template parameter)
      Element_pt[0] = new ELEMENT;

      // Find the number of nodes per element (N.B. all elements are identical
      // so we can determine this value once and for all). 
      unsigned n_node = finite_element_pt(0)->nnode();
      //Loop over all the nodes of the first element
      for(unsigned n=0; n<n_node; n++)
      {
	// Construct the next node and add it to the Mesh::Node_pt vector
	// Note that these interior nodes need not (and should not)
	// be boundary nodes, so they are created using the construct_node
	// function, which has the same interface as
	// construct_boundary_node()
	Node_pt.push_back(finite_element_pt(0)->construct_boundary_node(n));
	/* Node_pt.push_back(finite_element_pt(0)->construct_node(n)); */
      }
  
      // Loop over the remaining elements apart from the last
      for(unsigned e=1; e<(n_element-1); e++)
      {
	// Construct the e-th element
	Element_pt[e] = new ELEMENT;

	// The first local node of the e-th element is the last local node
	// of the (e-1)-th element. We MUST NOT construct the node twice.
	// Instead, we set the pointer in the e-th element to point to the
	// previously created node in the (e-1)-th element.
	finite_element_pt(e)->node_pt(0) = 
	  finite_element_pt(e-1)->node_pt(n_node-1);

	// Loop over the remaining nodes of the e-th element
	for(unsigned n=1; n<n_node; n++)
	{
	  // Construct the next node and add it to the Mesh::Node_pt vector
	  // Note that these interior nodes need not (and should not)
	  // be boundary nodes, so they are created using the construct_node
	  // function, which has the same interface as
	  // construct_boundary_node()
	  Node_pt.push_back(finite_element_pt(e)->construct_boundary_node(n));
	  /* Node_pt.push_back(finite_element_pt(e)->construct_node(n)); */
	}
      } // End of loop over elements
  
  
      // Construct the final element
      Element_pt[n_element-1] = new ELEMENT;
  
      // The first local node of the final element is the last local node
      // of the penultimate element. We MUST NOT construct the node twice.
      // Instead, we set the pointer in the final element to point to the
      // previously created node in the penultimate element.
      finite_element_pt(n_element-1)->node_pt(0) = 
	finite_element_pt(n_element-2)->node_pt(n_node-1);

      // Loop over the remaining central nodes of the final element
      for(unsigned n=1; n<(n_node-1); n++)
      {
	// Construct the next node and add it to the Mesh::Node_pt vector
	// Note that these interior nodes need not (and should not)
	// be boundary nodes, so they are created using the construct_node
	// function()
	Node_pt.push_back(finite_element_pt(n_element-1)->construct_boundary_node(n));
	/* Node_pt.push_back(finite_element_pt(n_element-1)->construct_node(n)); */
      }

      // Now make the mesh circular - don't want to construct a new final node,
      // but want to set the pointer to the last node in the last element to
      // be the pointer to the first node in the first element
      finite_element_pt(n_element - 1)->node_pt(n_node-1) =
	finite_element_pt(0)->node_pt(0);
      
      // We've now created all the nodes -- let's set their positions:
      // -------------------------------------------------------------
      // Find the total number of nodes
      unsigned n_global_node = nnode();

      // Loop over all nodes
      for(unsigned n=0; n<n_global_node; n++)
      {
	// equal spacing in azimuthal angle
	double phi = double(n) * 2.0 * MathematicalConstants::Pi / (n_global_node);

	// and convert to 2D Cartesian coordinates
	Node_pt[n]->x(0) = cos(phi);
	Node_pt[n]->x(1) = sin(phi);
      }
    } // End of constructor

  private:
    /* double random_var_to_alter_mem_map; */
  }; // End of OneDimMesh class.
  
  ////////////////////////////////////////////////////////////////////////
  ////////////////////////////////////////////////////////////////////////
  ////////////////////////////////////////////////////////////////////////

  //======================================================================
  /// \short A class for 1D line elements which manage the singular
  /// function pointers, store the singular amplitudes, and provide various
  /// functionality such as providing scaled and unscaled singular
  /// velocities / pressures, interpolating amplitudes and their
  /// derivatives, etc.
  /// N.B. this has nothing to do with the Poisson equation! this class
  /// merely inherits from QPoissonElement to hijack useful functionality
  //======================================================================
  
  template <unsigned NNODE_1D>
  class ScalableSingularityForNavierStokesLineElement :
    public virtual QPoissonElement<1, NNODE_1D>    
  {
  public:

    // ------------------------------------------
    // statics and typedefs
    // ------------------------------------------
    
    // expose the template argument so the mesh class can figure this out
    static const unsigned _NNODE_1D_ = NNODE_1D;
    
    // typedefs for the functions which provide the singular solutions and
    // their derivatives
    typedef Vector<double>(*UnscaledSingSolnFctPt) (const LagrangianCoordinates&);

    typedef DenseMatrix<double>(*GradientOfUnscaledSingSolnFctPt)
      (const LagrangianCoordinates&);

    // function which provides dzeta/dx
    typedef Vector<double>(*DzetaDxFctPt)(const LagrangianCoordinates&);

    // function pointer for helper function which computes the
    // constitutive relationship
    typedef DenseMatrix<double> (*StressFctPt)(const DenseMatrix<double>& strain_rate,
					       const double& p);
    
    // ------------------------------------------
    
    /// \short Constructor
    ScalableSingularityForNavierStokesLineElement() :						  
    Dim(1), Bulk_dim(3), Dzeta_dx_fct_pt(nullptr), Stress_fct_pt(nullptr),
      Zeta_has_been_setup(false)
    {
      // QUEHACERES get this programmatically
      this->set_nodal_dimension(3);      
    }
	    
    /// Broken copy constructor
    ScalableSingularityForNavierStokesLineElement(
      const ScalableSingularityForNavierStokesLineElement&)
    {
      BrokenCopy::broken_copy("NavierStokesWithSingularityLineElement");
    }

    /// Broken assignment operator
    void operator=(const ScalableSingularityForNavierStokesLineElement&)
    {
      BrokenCopy::broken_assign("ScalableSingularityForNavierStokesLineElement");
    }

    // Function to get the pointer to the function which computes dzeta/dx
    DzetaDxFctPt& dzeta_dx_fct_pt()
    {
      return Dzeta_dx_fct_pt;
    }

    // get the function pointer to the function which computes the
    // stress from the velocity gradient and the pressure
    StressFctPt& stress_fct_pt()
    {
      return Stress_fct_pt;
    }
    
    void add_unscaled_singular_fct_and_gradient_pt(
      const UnscaledSingSolnFctPt& sing_fct_pt,
      const GradientOfUnscaledSingSolnFctPt& grad_of_sing_fct_pt,
      const unsigned& sing_fct_id)
    {      
      // check if we've already got this ID
      if(Singular_fct_index.find(sing_fct_id) != Singular_fct_index.end())
      {
	throw OomphLibError(
	"Can't add duplicate singular function, functions must have a unique ID",
	OOMPH_CURRENT_FUNCTION,
	OOMPH_EXCEPTION_LOCATION);
      }

#ifdef PARANOID

      if(sing_fct_pt == nullptr)
      {
	throw OomphLibError(
	"Unscaled singular function pointer is null!",
	OOMPH_CURRENT_FUNCTION,
	OOMPH_EXCEPTION_LOCATION);
      }

      if(grad_of_sing_fct_pt == nullptr)
      {
	throw OomphLibError(
	"Gradient of unscaled singular function pointer is null!",
	OOMPH_CURRENT_FUNCTION,
	OOMPH_EXCEPTION_LOCATION);
      }
	    
#endif
      
      // add the ID to the map, index is the current size of the vectors
      // of function pointers
      unsigned nsing = Unscaled_singular_fct_pt.size();
      Singular_fct_index[sing_fct_id] = nsing; // since we're zero-indexed

      // add the ID to the list
      Singular_fct_id.push_back(sing_fct_id);
      
      // now add the function pointers
      Unscaled_singular_fct_pt.push_back(sing_fct_pt);
      Gradient_of_unscaled_singular_fct_pt.push_back(grad_of_sing_fct_pt);

      // and finally, resize the nodes in this element to accommodate
      // the new singular function
      for(unsigned j=0; j<this->nnode(); j++)
      {
	Node* node_pt = this->node_pt(j);

	// check it hasn't already been resized by a neighbouring element
	if(node_pt->nvalue() < (nsing + 1))
	  node_pt->resize(nsing+1);
      }
    }
    
    ///Function to compute unscaled version of the requested unscaled version
    Vector<double> unscaled_singular_fct(const LagrangianCoordinates& lagr_coords,
					 const unsigned& sing_fct_id) const
    {      
      // get the index for this ID from the map
      // (std::map::at() throws an exception if this index doesn't exist)
      unsigned sing_fct_index = Singular_fct_index.at(sing_fct_id);
      
      return Unscaled_singular_fct_pt[sing_fct_index](lagr_coords); 
    }

    ///Compute unscaled version of gradient of the requested singular function
    DenseMatrix<double> gradient_of_unscaled_singular_fct(const LagrangianCoordinates& lagr_coords,
							  const unsigned& sing_fct_id) const
    {
      // get the index for this ID from the map
      // (std::map::at() throws an exception if this index doesn't exist)
      unsigned sing_fct_index = Singular_fct_index.at(sing_fct_id);
      
      return Gradient_of_unscaled_singular_fct_pt[sing_fct_index](lagr_coords); 
    }

    ///Compute scaled version of the requested singular function
    Vector<double> singular_fct(const LagrangianCoordinates& lagr_coords,
				const Vector<double>& s,
				const unsigned& sing_fct_id) const
    {
      // get number of unknowns in the problem; plus one because we want pressure as well
      // as the velocity components
      const unsigned nvalue = lagr_coords.ncoord + 1;

      // storage for the scaled basis functions
      Vector<double> scaled_singular_fct(nvalue, 0.0);

      // get the unscaled functions
      Vector<double> u_sing_unscaled = unscaled_singular_fct(lagr_coords, sing_fct_id);

      // get the singular amplitude associated with this ID and these local coords
      double c = interpolated_amplitude(s, sing_fct_id);
      
      // scale 'em
      for(unsigned i=0; i<nvalue; i++)
      {
	scaled_singular_fct[i] = c * u_sing_unscaled[i];
      }
      
      return scaled_singular_fct;
    }
            
    ///Compute scaled version of gradient of the requested singular function
    DenseMatrix<double> gradient_of_singular_fct(const LagrangianCoordinates& lagr_coords,
						 const Vector<double>& s,
						 const unsigned& sing_fct_id) const
    {
      // get the interpolated singular amplitude
      const double c = interpolated_amplitude(s, sing_fct_id);

      // get the unscaled singular velocity
      Vector<double> u_sing_unscaled = unscaled_singular_fct(lagr_coords, sing_fct_id);
      
      // get the unscaled singular velocity gradient
      DenseMatrix<double> dudx_sing_unscaled =
	gradient_of_unscaled_singular_fct(lagr_coords, sing_fct_id);

      // get the Lagrangian derivative of the singular amplitude, dc/dzeta
      double dc_dzeta = interpolated_dc_dzeta(s, sing_fct_id);

#ifdef PARANOID
      if (Dzeta_dx_fct_pt == nullptr)
      {
	throw OomphLibError("Error: function pointer for dzeta/dx hasn't been set\n",
			    OOMPH_CURRENT_FUNCTION,
			    OOMPH_EXCEPTION_LOCATION);
      }
#endif

      // get the Eulerian derivatives dzeta/dx from the function pointer
      Vector<double> dzeta_dx = Dzeta_dx_fct_pt(lagr_coords);
      
      // number of rows and columns
      const unsigned N = dudx_sing_unscaled.nrow();
      const unsigned M = dudx_sing_unscaled.ncol();

      // the scaled Eulerian velocity gradient
      DenseMatrix<double> dudx_sing(N,M, 0.0);
      
      // now compute it: d(c*u_sing)/dx = dc/dzeta * dzeta/dx * u_sing + c*du_sing/dx      
      for(unsigned i=0; i<N; i++)
      {
	for(unsigned j=0; j<M; j++)
	{
	  dudx_sing(i,j) = u_sing_unscaled[i] * dc_dzeta * dzeta_dx[j] +
	    c * dudx_sing_unscaled(i,j);
	}
      }
      return dudx_sing;
    }

    // wrapper to get the total contribution of the singular part of the solution
    Vector<double> total_singular_contribution(const LagrangianCoordinates& lagr_coords,
					       const Vector<double>& s) const
    {
      // total singular part of the solution
      Vector<double> u_sing_total(lagr_coords.ncoord + 1, 0.0);

      // loop over each singular function and add the contribution to the total
      for(std::map<unsigned,unsigned>::const_iterator it = Singular_fct_index.begin();
	  it != Singular_fct_index.end(); it++)
      {
	// get the singular function ID
	unsigned sing_fct_id = it->first;
	
	// get the contribution of the ith singular function
	Vector<double> u_sing = singular_fct(lagr_coords, s, sing_fct_id);

	// add it to the total
	for(unsigned j=0; j<u_sing.size(); j++)
	  u_sing_total[j] += u_sing[j];
      }
      return u_sing_total;
    }

    // wrapper to get the total contribution of the derivatives of the
    // singular functions
    DenseMatrix<double>
      total_singular_gradient_contribution(const LagrangianCoordinates& rzp_coords,
					   const Vector<double>& s) const
    {
      const unsigned nval = rzp_coords.ncoord;
      
      // total of the singular derivatives      
      DenseMatrix<double> dudx_total(nval, nval, 0.0);
            

      // loop over each singular function and add the contribution to the total
      for(std::map<unsigned,unsigned>::const_iterator it = Singular_fct_index.begin();
	  it != Singular_fct_index.end(); it++)
      {
	// get the singular function ID
	unsigned sing_fct_id = it->first;
	
	// get the contribution of the ith singular function
	DenseMatrix<double> dudx = gradient_of_singular_fct(rzp_coords, s, sing_fct_id);

	// add it to the total
	for(unsigned m=0; m<nval; m++)
	{
	  for(unsigned n=0; n<nval; n++)
	  {
	    dudx_total(m,n) += dudx(m,n);
	  }
	}
      }

      return dudx_total;
    }

    void total_singular_strain_rate(const LagrangianCoordinates& lagr_coords,
				    const Vector<double>& s,
				    DenseMatrix<double>& strain_rate_sing_total) const
    {
      // the sum of all scaled singular functions
      Vector<double> u_sing_total = total_singular_contribution(lagr_coords, s);

      // total singular contributions to velocity gradient tensor and strain rate
      DenseMatrix<double> dudx_sing_total =
	total_singular_gradient_contribution(lagr_coords, s);
      
      // now compute the singular strain rate
      for(unsigned i=0; i<Dim; i++)
      {
	for(unsigned j=0; j<Dim; j++)
	{
	  strain_rate_sing_total(i,j) = 0.5 * (dudx_sing_total(i,j) +
					       dudx_sing_total(j,i));
	}
      }
    }

    void total_singular_stress(const LagrangianCoordinates& lagr_coords,
			       const Vector<double>& s,
			       DenseMatrix<double>& stress_sing_total) const
    {
      // get the total singular strain rate
      DenseMatrix<double> strain_rate_sing_total(3, 3, 0.0);
      total_singular_strain_rate(lagr_coords, s, strain_rate_sing_total);

      // get the total singular pressure      
      Vector<double> u_sing_total = total_singular_contribution(lagr_coords,s);      
      double p_sing_total = u_sing_total[3];

      // now use the function pointer to get the constitutive relation
      if(Stress_fct_pt == nullptr)
      {
	throw OomphLibError(
	  "Stress function pointer hasn't been set for this singular line element",
	  OOMPH_CURRENT_FUNCTION,
	  OOMPH_EXCEPTION_LOCATION);
      }
      else
      {
	stress_sing_total = Stress_fct_pt(strain_rate_sing_total, p_sing_total);
      }
    }
    
    /// \Short function to compute the interpolated amplitude of the
    /// nth singular function at the local coordinate s
    double interpolated_amplitude(const Vector<double>& s, const unsigned& sing_fct_id) const
    {
      // get the index for this ID from the map
      // (std::map::at() throws an exception if this index doesn't exist)
      const unsigned sing_fct_index = Singular_fct_index.at(sing_fct_id);
      
      // make space for the shape functions
      Shape psi(this->nnode());
      
      // get 'em
      this->shape(s, psi);

      // amplitude of this singular function
      double interpolated_c = 0;

      // loop over each node in this element and add its contribution to the
      // interpolated amplitude
      for(unsigned j=0; j<this->nnode(); j++)
      {
	interpolated_c += this->nodal_value(j, sing_fct_index) * psi[j];
      }

      return interpolated_c;
    }

    /// \short Function to compute the interpolated value of dc/dzeta, i.e. the
    /// variation in the singular amplitude w.r.t. zeta
    double interpolated_dc_dzeta(const Vector<double>& s, const unsigned& sing_fct_id) const
    {
      // get the index for this ID from the map
      // (std::map::at() throws an exception if this index doesn't exist)
      const unsigned sing_fct_index = Singular_fct_index.at(sing_fct_id);

      // shorthand
      unsigned nnode = this->nnode();
      
      // make space for the shape functions and their derivatives
      Shape psi(nnode);      
      DShape dpsi_ds(nnode, Dim);

      // Find values of shape functions and their derivatives
      this->dshape_local(s, psi, dpsi_ds);

      // compute interpolated local derivatives dc/ds and dzeta/ds
      Vector<double> interpolated_dc_ds(Dim, 0.0);
      Vector<double> interpolated_dzeta_ds(Dim, 0.0);
      
      for(unsigned l=0; l<nnode; l++) 
      {
	// loop for generality, in reality this is 1D so only 1 local coordinate
	for(unsigned j=0; j<Dim; j++)
	{
	  interpolated_dc_ds[j]    += this->nodal_value(l, sing_fct_index) * dpsi_ds(l,j);
	  interpolated_dzeta_ds[j] += zeta_nodal(l,0,j) * dpsi_ds(l,j);
	}
      }

      // now compute dc/dzeta = dc/ds * 1/(dzeta/ds)
      double interp_dc_dzeta = 0.0;
      for(unsigned j=0; j<Dim; j++)
      {
	interp_dc_dzeta += interpolated_dc_ds[j] / interpolated_dzeta_ds[j];
      }
            
      return interp_dc_dzeta;
    }
    
    // Override to provide the global boundary zeta for the nth node on the
    // edge of the disk. 
    double zeta_nodal(const unsigned& n, const unsigned& k, const unsigned& i) const
    {
      if(!Zeta_has_been_setup)
      {
	throw OomphLibError(
	  "Boundary zeta coordinate hasn't been setup for this line element",
	  OOMPH_CURRENT_FUNCTION,
	  OOMPH_EXCEPTION_LOCATION);
      }
      
      return Zeta[n];
    }

    // 1D wrapper
    double zeta_nodal(const unsigned& n) const
    {
      return zeta_nodal(n, 0, 0);
    }

    // how many singular functions are we subtracting?
    unsigned nsingular_fct() const
    {
      return Unscaled_singular_fct_pt.size();
    }

    /// Output with various contributions
    void output(std::ostream& outfile, const unsigned& nplot)
    {
      // dimensionality of the nodes (not the element)
      unsigned node_dim = this->node_pt(0)->ndim();
      
      //Vector of local coordinates
      Vector<double> s(this->Dim);

      // Tecplot header info
      outfile << this->tecplot_zone_string(nplot);
   
      // Loop over plot points
      unsigned num_plot_points = this->nplot_points(nplot);
      for (unsigned iplot=0; iplot<num_plot_points; iplot++)
      {
	// Get local coordinates of plot point
	this->get_s_plot(iplot, nplot, s);
	
	Vector<double> x(node_dim);
	for(unsigned i=0; i<node_dim; i++) 
	{
	  x[i] = this->interpolated_x(s,i);
	  outfile << x[i] << " ";
	}

	// output the interpolated value of zeta at this plot point so that
	// we can plot stuff easily in azimuthal space
	double zeta = interpolated_zeta(s);

	outfile << zeta << " ";
	
	// now loop over the singular functions and get their
	// interpolated amplitudes at this plot point
	for(std::map<unsigned,unsigned>::iterator it = Singular_fct_index.begin();
	  it != Singular_fct_index.end(); it++)
	{
	  // get the singular function ID
	  unsigned sing_fct_id = it->first;
	  
	  double amplitude = interpolated_amplitude(s, sing_fct_id);
	  outfile << amplitude << " ";
	}
	
	// and the gradients
	for(std::map<unsigned,unsigned>::iterator it = Singular_fct_index.begin();
	  it != Singular_fct_index.end(); it++)
	{
	  // get the singular function ID
	  unsigned sing_fct_id = it->first;
	  
	  double dc_dzeta = interpolated_dc_dzeta(s, sing_fct_id);
	  outfile << dc_dzeta << " ";
	}
	
	outfile << std::endl;
      }
      
      // Write tecplot footer (e.g. FE connectivity lists)
      this->write_tecplot_zone_footer(outfile, nplot);   
    }

    /// Call this to bypass the correct computation of the
    /// residual for the nth singular function and replace it by r_c = C-ampl
    void impose_singular_fct_amplitude(const unsigned& sing_fct_id,
				       const Vector<double>& ampl) const
    {
      unsigned nnode = this->nnode();
      
      // check the user has provided enough values for the whole element
      if(ampl.size() != nnode)
      {
	ostringstream error_message;

	error_message << "Error: need an amplitude for each node in this element;\n"
		      << ampl.size() << " provided, but this element has " << nnode
		      << " nodes.\n";
	
	throw OomphLibError(error_message.str(),
			    OOMPH_CURRENT_FUNCTION,
			    OOMPH_EXCEPTION_LOCATION);
      }

      // get the index for this ID from the map
      // (std::map::at() throws an exception if this index doesn't exist)
      unsigned sing_fct_index = Singular_fct_index.at(sing_fct_id);
      
      for(unsigned j=0; j<nnode; j++)
      {	
	// pin it so that the amplitude is no longer a dof
	this->node_pt(j)->pin(sing_fct_index);

	// set its value to the imposed amplitude
	this->node_pt(j)->set_value(sing_fct_index, ampl[j]);	
      }      
    } 

    /// Reset all singular amplitudes to compute r_c properly via integral
    void dont_impose_singular_fct_amplitude() const
    {
      for(unsigned j=0; j<this->nnode(); j++)
      {
	for(unsigned i=0; i<this->node_pt(i)->nvalue(); i++)
	{
	  this->node_pt(j)->unpin_value(i);
	}
      }
    } 

    // unpin a particular singular amplitude
    void dont_impose_singular_fct_amplitude(const unsigned& sing_fct_id) const
    {
      // get the index for this ID from the map
      // (std::map::at() throws an exception if this index doesn't exist)
      unsigned sing_fct_index = Singular_fct_index.at(sing_fct_id);
      
      for(unsigned j=0; j<this->nnode(); j++)
      {
	this->node_pt(j)->unpin_value(sing_fct_index);
      }
    }

    /// Add the element's contribution to its residual vector
    inline void fill_in_contribution_to_residuals(Vector<double>& residuals)
    {  	      
      // No local contributions - all done by external elements!	    
    }

    // overload the Jacobian contributions (don't want Poisson!)
    inline void fill_in_contribution_to_jacobian(Vector<double>& residuals,
    						 DenseMatrix<double>& jacobian)
    {
      // No local contributions - all done by external elements!
    }
    
    // function to return whether we're imposing the amplitude of the nth singular fct
    bool is_singular_fct_amplitude_imposed(const unsigned& sing_fct_id)
    {
      // get the index for this ID from the map
      // (std::map::at() throws an exception if this index doesn't exist)
      unsigned sing_fct_index = Singular_fct_index.at(sing_fct_id);
      
      // we're not allowing some nodes to be pinned and others unpinned at the moment,
      // so it's enough to check the first node
      return this->node_pt(0)->is_pinned(sing_fct_index);
    }

    // returns a list of the IDs, i.e. the keys to the ID->index map 
    Vector<unsigned> singular_fct_ids() const
    {
      return Singular_fct_id;
    }
     

    /// \short Function to compute the shape and test functions and to return 
    /// the Jacobian of mapping between local and global (Eulerian)
    /// coordinates. N.B. this is public, as the contributions to this elements
    /// residuals are all made via external elements which need access to these
    /// shape functions
    inline double shape_and_test(const Vector<double>& s, Shape& psi, Shape& test)
      const
    {
      //Find number of nodes
      unsigned n_node = this->nnode();

      //Get the shape functions
      this->shape(s, psi);

      //Set the test functions to be the same as the shape functions
      for(unsigned i=0; i<n_node; i++)
      {
	test[i] = psi[i];
      }

      //Return the value of the jacobian
      return this->J_eulerian(s);
    }

    /// \short Function to compute the Eulerian derivatives of the shape
    /// functions. N.B. this is public, as the contributions to this elements
    /// residuals are all made via external elements which need access to these
    /// shape functions
    inline void dshape_eulerian(const LagrangianCoordinates& lagr_coords,
				const Vector<double>& s,
				DShape& dpsi_dx) const
    {      
      //Find number of nodes
      unsigned nnode = this->nnode();

      Shape psi_dummy(nnode);
      DShape dpsi_ds(nnode, Dim);

      // Get the derivatives of the shape functions
      // w.r.t. the local (1D) coordinate
      this->dshape_local(s, psi_dummy, dpsi_ds);

      // get the derivative of zeta (which parameterises the line mesh of
      // these elements) w.r.t. this elements local coordinate
      double interpolated_dzeta_ds = 0.0;

      for(unsigned l=0; l<nnode; l++) 
      {
	interpolated_dzeta_ds += zeta_nodal(l) * dpsi_ds(l,0);
      }
      
#ifdef PARANOID
      if (Dzeta_dx_fct_pt == nullptr)
      {
	throw OomphLibError("Error: function pointer for dzeta/dx hasn't been set\n",
			    OOMPH_CURRENT_FUNCTION,
			    OOMPH_EXCEPTION_LOCATION);
      }
#endif

      // get the Eulerian derivatives dzeta/dx from the function pointer
      Vector<double> dzeta_dx = Dzeta_dx_fct_pt(lagr_coords);

      // now compute the Eulerian derivatives via chain rule:
      // dshape/dx = dshape/ds * 1/(dzeta/ds) * dzeta/dx
      // (this works because we're in 1D so don't have to compute the inverse
      // of a matrix that would be needed for inverting functions of
      // multiple variables)
      for(unsigned n=0; n<nnode; n++)
      {
	// loop over the problem dimensions (not this elements dimensions)
	for(unsigned i=0; i<this->node_pt(0)->ndim(); i++)
	{
	  dpsi_dx(n,i) = dpsi_ds(n,0) * dzeta_dx[i] / interpolated_dzeta_ds;
	}
      }
    }

    double interpolated_zeta(const Vector<double>& s) const
    {
      // setup memory for shape functions
      Shape psi(this->nnode());
      
      // Get 'em
      this->shape(s, psi);

      // interpolated zeta coordinate
      double zeta = 0;

      // loop over the nodes and add their contributions
      for(unsigned j=0; j<this->nnode(); j++)
      {
	zeta += zeta_nodal(j) * psi[j];
      }

      return zeta;
    }
    
    // compute the boundary zeta, which is the initial azimuthal angle
    void setup_zeta_nodal(const bool& use_zeta_2pi_instead_of_0 = false)
    {
      // make sure we've got enough space
      Zeta.resize(this->nnode(), 0.0);

      // tolerance on floating-point 
      const double tol = 1e-10;
      for(unsigned j=0; j<this->nnode(); j++)
      {
	// get the azimthal angle of this node
	double zeta = atan2pi(this->node_pt(j)->x(1), this->node_pt(j)->x(0));

	// are we approaching from the negative or positive half plane?
	if((abs(zeta) < tol) && use_zeta_2pi_instead_of_0)
	  zeta = 2.0 * MathematicalConstants::Pi;
	else if((abs(zeta - 2.0 * MathematicalConstants::Pi) < tol) &&
		 !use_zeta_2pi_instead_of_0)
	  zeta = 0;

	if(zeta < 0)
	{
	  throw OomphLibError(
	    "something's gone wrong, zeta shouldn't be negative",
	    OOMPH_CURRENT_FUNCTION,
	    OOMPH_EXCEPTION_LOCATION);
	}
	
	Zeta[j] = zeta;
      }

      // set the flag to say we've done it
      Zeta_has_been_setup = true;
    }
    
  private:

    /// Dimensionality of this element
    unsigned Dim;

    /// Dimensionality of the whole problem
    unsigned Bulk_dim;
    
    /// Pointers to the singular functions
    Vector<UnscaledSingSolnFctPt> Unscaled_singular_fct_pt;

    /// Pointers to gradients of the singular funcions
    Vector<GradientOfUnscaledSingSolnFctPt> Gradient_of_unscaled_singular_fct_pt;

    /// Function pointer to a function which computes dzeta/dx_j
    DzetaDxFctPt Dzeta_dx_fct_pt;

    /// Function pointer to a pointer which computes the constitutive relationship,
    /// i.e. stress as a function of strain rate
    StressFctPt Stress_fct_pt;
    
    /// Object which maps an ID for a given singular function to it's
    // corresponding index in the vectors of singular function pointers
    std::map<unsigned, unsigned> Singular_fct_index;

    /// The IDs of the singular functions, i.e. the keys for the above map
    Vector<unsigned> Singular_fct_id;
    
    /// Flag to check that the nodal zeta values have been setup appropriately
    // so that we don't need to recompute
    bool Zeta_has_been_setup;
    
    /// The boundary zeta of each node in this element
    Vector<double> Zeta;
  };

  //////////////////////////////////////////////////////////////////////////////
  //////////////////////////////////////////////////////////////////////////////
  //////////////////////////////////////////////////////////////////////////////
  //////////////////////////////////////////////////////////////////////////////

  //===========================================================================
  // Generic base class for augmented face elements
  //===========================================================================
  template <class ELEMENT>
    class NavierStokesWithSingularityFaceElement :
    public virtual FaceGeometry<ELEMENT>, public virtual FaceElement
  {
  public:

    /* // Extract the number of nodes per side from the templated element */
    // QUEHACERES sort this out, don't get this from a global #define
    static const unsigned NNODE_1D = SINGULAR_ELEMENT_NNODE_1D; // ELEMENT::_NNODE_1D_;

    // shorthand
    typedef ScalableSingularityForNavierStokesLineElement<NNODE_1D> SingularLineElement;

    typedef void (*ExactTractionFctPt)(const double& t,
				       const Vector<double>& x,
				       const Vector<double>& unit_normal,
				       Vector<double>& traction);

    // Function which computes the Eulerian derivatives of the edge coordinates
    typedef double (*DxiDxFctPt)(const LagrangianCoordinates&);
    
    /// \short Constructor, takes the pointer to the "bulk" element and the 
    /// index of the face to which the element is attached. Optional final
    /// arg is the identifier for the additional unknowns (Lagrange multipliers)
    NavierStokesWithSingularityFaceElement(ELEMENT* const& bulk_el_pt, 
					     const int& face_index, 
					     const unsigned& id=0);

    NavierStokesWithSingularityFaceElement()
    {
      	throw OomphLibError(
	  "Don't call empty constructor for NavierStokesWithSingularityFaceElement",
	  OOMPH_CURRENT_FUNCTION,
	  OOMPH_EXCEPTION_LOCATION);
    }

     /// Broken copy constructor
    NavierStokesWithSingularityFaceElement(
      const NavierStokesWithSingularityFaceElement& dummy) 
    { 
	BrokenCopy::broken_copy("NavierStokesWithSingularityFaceElement");
    } 
    
    /// Broken assignment operator
    void operator=(const NavierStokesWithSingularityFaceElement&) 
      {
	BrokenCopy::broken_assign("NavierStokesWithSingularityFaceElement");
      }
    
    // set the edge coordinates at each knot
    void set_lagr_coordinates_at_knot(const Vector<LagrangianCoordinates>& coords)
    {
      // number of knot points in this element
      unsigned nknot = this->integral_pt()->nweight();

      // number of coordinates supplied
      unsigned ncoords = coords.size();

      // check the right number of coordinates have been passed in
      if(ncoords != nknot)
      {
	ostringstream error_message;
	error_message << "Number of sets of coordinates provided is not consistent "
		      << "with the number of knots.\n"
		      << "Number of supplied coordinate sets: " << ncoords << "\n"
		      << "Number of knots:                    " << nknot;
	
	throw OomphLibError(
	  error_message.str(),
	  OOMPH_CURRENT_FUNCTION,
	  OOMPH_EXCEPTION_LOCATION);
      }

      // make enough space
      Lagrangian_coordinates_at_knot.resize(nknot);
	
      // set 'em
      for(unsigned i=0; i<nknot; i++)
      	Lagrangian_coordinates_at_knot[i] = coords[i];          
    }

    // ========================================================================
    // get the edge coordinates at the ith knot
    LagrangianCoordinates lagr_coordinate_at_knot(const unsigned& i) const
    {
      return Lagrangian_coordinates_at_knot[i];
    }

    // ========================================================================
    // set the line element and local coordinate for the ith knot
    void set_line_element_and_local_coordinate_at_knot(
      const Vector<std::pair<GeomObject*, Vector<double> > >&
      line_element_and_local_coordinate_at_knot)
    {
      // number of knot points in this element
      unsigned nknot = this->integral_pt()->nweight();

      // number of coordinates supplied
      unsigned ncoords = line_element_and_local_coordinate_at_knot.size();

      // check the right number of coordinates have been passed in
      if(ncoords != nknot)
      {
	ostringstream error_message;
	error_message << "Number of element-coordinate pairs provided is not consistent "
		      << "with the number of knots.\n"
		      << "Number of supplied coordinate sets: " << ncoords << "\n"
		      << "Number of knots:                    " << nknot;
	
	throw OomphLibError(
	  error_message.str(),
	  OOMPH_CURRENT_FUNCTION,
	  OOMPH_EXCEPTION_LOCATION);
      }

      // make enough space
      Line_element_and_local_coordinate_at_knot.resize(nknot);
      
      // set 'em
      for(unsigned i=0; i<nknot; i++)
      {
	Line_element_and_local_coordinate_at_knot[i] =
	  line_element_and_local_coordinate_at_knot[i];
      }

      // now we have the singular element pointers, set their dofs
      // as external data to this element
      set_singular_amplitudes_as_external_data();
    }

    // ========================================================================
    // get the line element and local coordinate for the ith knot
    std::pair<GeomObject*, Vector<double> >
      line_element_and_local_coordinate_at_knot(const unsigned& i) const
    {
      if(Line_element_and_local_coordinate_at_knot.size() <= i)
      {
	GeomObject* dummy_gom_obj_pt = nullptr;
	std::pair<GeomObject*, Vector<double> > dummy_pair =
	  std::make_pair(dummy_gom_obj_pt, Vector<double>(0));
	  
	return dummy_pair;
      }
      
      return Line_element_and_local_coordinate_at_knot[i];
    }

    // ========================================================================
    // get the edge coordinates and (a pointer to) the element which stores
    // the singular amplitudes
    void lagr_coords_and_singular_element_at_knot(
      const unsigned& ipt,
      LagrangianCoordinates& lagr_coords,
      SingularLineElement*& sing_el_pt,
      Vector<double>& s_singular_el) const
    {
      // get the line element and local coordinate which corresponds to the
      // singular amplitude for this knot
      std::pair<GeomObject*, Vector<double> > line_elem_and_local_coord = 
	this->line_element_and_local_coordinate_at_knot(ipt);

      // cast the GeomObject to a singular line element      
      sing_el_pt = dynamic_cast<SingularLineElement*>
	(line_elem_and_local_coord.first);

      // check we've actually got one, or we're in trouble
      if(sing_el_pt == nullptr)
      {
	if(this->No_throw_if_no_singular_elem_pt)
	{
	  return;
	}
	else
	{ 
	  ostringstream error_message;

	  error_message << "Error: this singular face element has no "
			<< "singular line element pointer\n";
	    
	  throw OomphLibError(error_message.str().c_str(),
			      OOMPH_CURRENT_FUNCTION,
			      OOMPH_EXCEPTION_LOCATION);
	}
      }

      // local coordinate in the singular element for the zeta of this knot
      s_singular_el = line_elem_and_local_coord.second;

      // get the \rho,\zeta,\phi coordinates at this knot
      lagr_coords = this->lagr_coordinate_at_knot(ipt);
    }

    // ========================================================================
    // Assign the singular line elements associated with this face elements
    // integration points as external data for this element, since it will
    // make a contribution to the singular amplitudes
    void set_singular_amplitudes_as_external_data()
    { 
      // number of integration points in this element
      const unsigned n_intpt = this->integral_pt()->nweight();

      // Make space in our list of external equation indices
      C_external_data_index_at_knot.resize(n_intpt);
      
      // set of unique singular line element nodes to who's dofs (singular amplitudes)
      // this face element is responsible for making a contribution to.
      std::set<Node*> external_singular_node_set;
      
      // need to loop over the integration points, because different integration points
      // in this same element might correspond to different singular line elements; 
      // find the singular line element associated with each, and add it to our set
      for(unsigned ipt=0; ipt<n_intpt; ipt++)
      {  
	SingularLineElement* sing_el_pt = nullptr;
	LagrangianCoordinates lagr_coords_dummy;
	Vector<double> s_singular_el_dummy;
	lagr_coords_and_singular_element_at_knot(ipt, lagr_coords_dummy,
						 sing_el_pt, s_singular_el_dummy);

	// for each knot in this face element, we want to store the
	// index to each external node
	C_external_data_index_at_knot[ipt].resize(sing_el_pt->nnode());
	
	// now loop over this singular elements nodes
	for(unsigned j=0; j<sing_el_pt->nnode(); j++)
	{
	  C_external_data_index_at_knot[ipt][j] = this->add_external_data(sing_el_pt->node_pt(j));
	}	
      }
    }
    
    // ========================================================================
    /// \short Specify the value of nodal zeta from the face geometry
    /// The "global" intrinsic coordinate of the element when
    /// viewed as part of a geometric object should be given by
    /// the FaceElement representation, by default (needed to break
    /// indeterminacy if bulk element is SolidElement)
    double zeta_nodal(const unsigned& n, const unsigned& k,
		      const unsigned& i) const 
    {
      return FaceElement::zeta_nodal(n,k,i);
    }

    /// Pin Lagrange multiplier associated with ith coordinate at specified local node
    void pin_lagrange_multiplier_at_specified_local_node(const unsigned& j,
							 const unsigned& i,
							 const int& id = -1) const
    {
      // get the face IDs map for this node
      map<unsigned, unsigned> map_l = *(
	dynamic_cast<BoundaryNodeBase*>(this->node_pt(j))->
	index_of_first_value_assigned_by_face_element_pt() );

      unsigned lambda_index;

      // if no id specified, just take the index for the first (and probably only)
      // boundary in the map
      if(id == -1)
      {
	lambda_index = map_l.begin()->second + i;
      }	
      else
      {
	// otherwise, get the nodal index for the specified boundary ID
	lambda_index = map_l.at(id) + i;
      }

      // pin and set to zero
      this->node_pt(j)->pin(lambda_index);
      this->node_pt(j)->set_value(lambda_index, 0.0);
    }

    /// Unpin ith component of FE part of the solution at specified local node
    void unpin_u_fe_at_specified_local_node(const unsigned& j, const unsigned& i) const
    {   
      this->node_pt(j)->unpin(i);	  	
    }
      
    /// C-style output function -- forward to broken version in FiniteElement
    /// until somebody decides what exactly they want to plot here...
    void output(FILE* file_pt)
    {
      FiniteElement::output(file_pt);
    }

    /// \short C-style output function -- forward to broken version in 
    /// FiniteElement until somebody decides what exactly they want to plot 
    /// here...
    void output(FILE* file_pt, const unsigned& n_plot)
    {
      FiniteElement::output(file_pt, n_plot);
    }

    // function to get the pressure at local coordinates s in the
    // face element via the bulk element
    double interpolated_p_nst(const Vector<double>& s) const
    {
      // number of dimensions in the bulk
      unsigned dim_bulk = Dim + 1;
      
      // local coordinates in bulk element
      Vector<double> s_bulk(dim_bulk, 0.0);
      
      // Get pointer to assocated bulk element
      ELEMENT* bulk_el_pt = dynamic_cast<ELEMENT*>(bulk_element_pt());

      this->get_local_coordinate_in_bulk(s, s_bulk);
      
      //Get pressure from bulk
      return bulk_el_pt->interpolated_p_nst(s_bulk);
    }
    
    // function to get the strain rate at local coordinates s in the face element
    // via the bulk element
    void strain_rate(const Vector<double>& s, DenseMatrix<double>& _strain_rate) const
    {
      // number of dimensions in the bulk
      unsigned dim_bulk = Dim + 1;
      
      // local coordinates in bulk element
      Vector<double> s_bulk(dim_bulk, 0.0);
      
      // Get pointer to assocated bulk element
      ELEMENT* bulk_el_pt = dynamic_cast<ELEMENT*>(bulk_element_pt());

      this->get_local_coordinate_in_bulk(s, s_bulk);
      
      //Get pressure from bulk
      return bulk_el_pt->strain_rate(s_bulk, _strain_rate);
    }

    void integrated_force_and_moment(const Vector<double>& centre_of_rotation,
				     Vector<double>& integrated_force,
				     Vector<double>& integrated_moment) const
    {
      //Find out how many nodes there are
      const unsigned n_node = this->nnode();

      // shorthand
      const unsigned Dim = this->Dim;
    
      //Set up memory for the shape and test functions
      Shape psi(n_node), test(n_node);
     
      // The number of integration points (knots)
      const unsigned n_intpt = this->integral_pt()->nweight();
     
      //Set the Vector to hold local coordinates
      Vector<double> s(Dim-1);
	
      // shorthand
      ELEMENT* bulk_el_pt = dynamic_cast<ELEMENT*>(this->bulk_element_pt());

      // Make sure we've got enough space and intialised everything to zero
      integrated_force.resize(Dim, 0.0);
      integrated_force.initialise(0.0);

      integrated_moment.resize(3, 0.0);
      integrated_moment.initialise(0.0);
  
      //Loop over the integration points
      for(unsigned ipt=0; ipt<n_intpt; ipt++)
      {       
	//Assign values of s
	for(unsigned i=0; i<(Dim-1); i++)
	{
	  s[i] = this->integral_pt()->knot(ipt,i);
	}

	// Eulerian position of this integration point
	Vector<double> x(Dim, 0.0);
	
	for(unsigned i=0; i<Dim; i++)
	  x[i] = this->interpolated_x(s,i);

	// get the body force 
	Vector<double> body_force(3, 0.0);
	bulk_el_pt->get_body_force(x, body_force);
      
	// radial distance from this Eulerian position to the specified
	// centre of rotation
	Vector<double> r(Dim, 0.0);
	
	for(unsigned i=0; i<Dim; i++)
	  r[i] = x[i] - centre_of_rotation[i];
	  
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

	// FE strain rate tensor
	DenseMatrix<double> strain_rate_fe(Dim, Dim, 0.0);

	this->strain_rate(s, strain_rate_fe);

	// get the FE pressure
	double p_fe = this->interpolated_p_nst(s);
	  
	// get contribution to the FE stress on this element
	DenseMatrix<double> stress_fe(Dim, Dim, 0.0);

	stress_fe = (*bulk_el_pt->stress_fct_pt())(strain_rate_fe, p_fe);

	// Now get the total singular contributions
	// --------------------------------------------------------

	// total singular stress
	DenseMatrix<double> stress_sing_total(Dim, Dim, 0.0);
	  
	// get the line element and local coordinate which corresponds to the
	// singular amplitude for this knot
	std::pair<GeomObject*, Vector<double> > line_elem_and_local_coord = 
	  this->line_element_and_local_coordinate_at_knot(ipt);

	// cast the GeomObject to a singular line element      
	SingularLineElement* sing_el_pt =
	  dynamic_cast <SingularLineElement*>( line_elem_and_local_coord.first );

	if(sing_el_pt != nullptr)
	{
	  // local coordinate in the singular element for the zeta of this knot
	  Vector<double> s_singular_el = line_elem_and_local_coord.second;
      
	  // get the \rho,\zeta,\phi coordinates at this knot
	  LagrangianCoordinates lagr_coords_at_knot = this->lagr_coordinate_at_knot(ipt);

	  // the sum of all scaled singular functions
	  Vector<double> u_sing_total(Dim+1, 0.0);

	  // total singular contributions to strain rate tensor
	  DenseMatrix<double> strain_rate_sing_total(Dim, Dim, 0.0);
      
	  u_sing_total = sing_el_pt->total_singular_contribution(lagr_coords_at_knot,
								 s_singular_el);

	  sing_el_pt->total_singular_strain_rate(lagr_coords_at_knot,
						 s_singular_el,
						 strain_rate_sing_total);

	  // total singular pressure
	  double p_sing_total = u_sing_total[this->P_index_nst];

	  // total scaled singular stress
	  stress_sing_total = (*bulk_el_pt->stress_fct_pt())(strain_rate_sing_total,
							     p_sing_total);	  
	}

	Vector<double> total_traction(Dim, 0.0);
	
	// now add the weighted contribution of the traction at this knot to the total
	for(unsigned i=0; i<Dim; i++)
	{
	  for(unsigned j=0; j<Dim; j++)
	  {
	    total_traction[i] += (stress_fe(i,j) + stress_sing_total(i,j)) * unit_normal[j];
	  }	  
	}

	for(unsigned i=0; i<Dim; i++)
	{
	  // QUEHACERES experimental - half the body force, since the
	  // top and bottom surfaces are integrated over separately so don't
	  // want double counting //  + 0.5*body_force[i]
	  integrated_force[i] += (total_traction[i] ) * W;
	}

	integrated_moment[2] += (r[0]*total_traction[1] - r[1]*total_traction[0]) * W;
	
	if (Dim == 3)
	{
	  integrated_moment[0] += (r[1]*total_traction[2] - r[2]*total_traction[1]) * W;
	  integrated_moment[1] += (r[2]*total_traction[0] - r[0]*total_traction[2]) * W;	 
	}
	
      } // end loop over integration points
    }
    
    ExactTractionFctPt& exact_traction_fct_pt()
    {
      return Exact_traction_fct_pt;
    }

    DxiDxFctPt& dxi_dx_fct_pt()
    {
      return Dxi_dx_fct_pt;
    }

    // allows these elements to be used in non-augmented regions without dying
    // trying to look for associated singular line elements
    void no_throw_if_no_singular_elem_pt()
    {
      No_throw_if_no_singular_elem_pt = true;
    }

    void throw_if_no_singular_elem_pt()
    {
      No_throw_if_no_singular_elem_pt = false;
    }

    // wrapper to allow the protected FiniteElement::fill_in... to be
    // called from the outside (for debug)
    void fill_in_jacobian_by_fd(Vector<double>& residuals,
				DenseMatrix<double>& jacobian)
    {
      FiniteElement::fill_in_contribution_to_jacobian(residuals, jacobian);
    }
      
  protected:

    // derivative of Eulerian position w.r.t.
    // local coordinates in this face element
    void interpolated_dx_ds(const Vector<double>& s,
			    DenseMatrix<double>& interpolated_dxds) const
    {
      //Find the number of nodes in the bulk element
      const unsigned n_node_bulk = Bulk_element_pt->nnode();
   
      //Find the number of position types in the bulk element
      const unsigned n_position_type_bulk = 
	Bulk_element_pt->nnodal_position_type();

      //Find the spatial dimension of the FaceElement
      const unsigned element_dim = dim();

      //Find the overall dimension of the problem 
      //(assume that it's the same for all nodes)
      const unsigned spatial_dim = nodal_dimension();
 
      //Construct the local coordinate in the bulk element
      Vector<double> s_bulk(2, 0.0);
      
      //Get the local coordinates in the bulk element
      this->get_local_coordinate_in_bulk(s, s_bulk);
    
      //Allocate storage for the shape functions and their derivatives wrt
      //local coordinates
      Shape psi_dummy(n_node_bulk, n_position_type_bulk);
      DShape dpsids(n_node_bulk, n_position_type_bulk, 2);
      
      //Get the value of the shape functions at the given local coordinate
      this->Bulk_element_pt->dshape_local(s_bulk, psi_dummy, dpsids);
 
      //Calculate all derivatives of the spatial coordinates 
      //wrt local coordinates
      interpolated_dxds.resize(2, spatial_dim, 0.0);

      // zero out to be sure since we're adding
      interpolated_dxds.initialise(0.0);
      
      //Loop over all parent nodes
      for(unsigned l=0; l<n_node_bulk; l++)
      {
	//Loop over all position types in the bulk
	for(unsigned k=0; k<n_position_type_bulk; k++)
	{
	  //Loop over derivative direction
	  for(unsigned j=0; j<2; j++)
	  {
	    //Loop over coordinate directions
	    for(unsigned i=0; i<spatial_dim; i++)
	    {
	      //Compute the spatial derivative
	      interpolated_dxds(i,j) += 
		this->Bulk_element_pt->nodal_position_gen(l,k,i) * dpsids(l,k,j);
	    }
	  }
	}
      }
    }
    
    /// \short Function to compute the shape and test functions and to return 
    /// the Jacobian of mapping between local and global (Eulerian)
    /// coordinates
    inline double shape_and_test(const Vector<double>& s,
				 Shape& psi,
				 Shape& test) const
    {
      //Find number of nodes
      unsigned n_node = nnode();

      //Get the shape functions
      shape(s, psi);

      //Set the test functions to be the same as the shape functions
      for(unsigned i=0; i<n_node; i++)
      {
	test[i] = psi[i];
      }

      //Return the value of the jacobian
      return J_eulerian(s);
    }


    /// \short Function to compute the shape and test functions and to return 
    /// the Jacobian of mapping between local and global (Eulerian)
    /// coordinates
    inline double shape_and_test_at_knot(const unsigned& ipt,
					 Shape& psi,
					 Shape& test) const
    {
      //Find number of nodes
      unsigned n_node = nnode();

      //Get the shape functions
      shape_at_knot(ipt,psi);

      //Set the test functions to be the same as the shape functions
      for(unsigned i=0; i<n_node; i++)
      {
	test[i] = psi[i];
      }

      //Return the value of the jacobian
      return J_eulerian_at_knot(ipt);
    }

  protected:
    
    /// Number of spatial dimensions in the problem
    unsigned Dim;

    /// The index at which the unknown is stored at the nodes
    unsigned P_index_nst;

    /// ID of the boundary this face element sits on
    unsigned Boundary_id;
    
    /// \short Indices of external Data that store the values of the amplitudes of
    /// the singular functions - the outer index is the integration point number,
    /// and the inner index is the node number of the singular line element
    /// which coresponds to the integration point
    Vector<Vector<unsigned> > C_external_data_index_at_knot;

    /// Don't throw an error if there is no pointer to a singular line element
    bool No_throw_if_no_singular_elem_pt;
    
  private:
    
    // QUEHACERES exact solution fct pointer, for debug
    ExactTractionFctPt Exact_traction_fct_pt;

    /// \short Pointer to a function which computes the Eulerian derivatives of the
    /// edge coordinates, i.e. dxi/dx
    DxiDxFctPt Dxi_dx_fct_pt;
    
    /// \short Edge coordinates (\rho, \zeta, \phi) of each of this element's
    /// knot points
    Vector<LagrangianCoordinates> Lagrangian_coordinates_at_knot;

    /// \short Edge coordinates (\rho, \zeta, \phi) of each of this element's nodes
    Vector<LagrangianCoordinates> Nodal_lagr_coordinates;
    
    /// \short The line element and its local coordinate(s) which correspond to the zeta
    /// values at each knot point in this bulk element. N.B. the element
    /// is 1D so only has 1 local coordinate, but we keep this as a vector
    /// for consistency with the general case
    Vector<std::pair<GeomObject*, Vector<double> > >
      Line_element_and_local_coordinate_at_knot;

    /// \short The line element and its local coordinate(s) which correspond to the zeta
    /// values at eachnode in this bulk element. N.B. the element
    /// is 1D so only has 1 local coordinate, but we keep this as a vector
    /// for consistency with the general case
    Vector<std::pair<GeomObject*, Vector<double> > >
      Line_element_and_local_coordinate_at_node;
  };

  //===========================================================================
  /// Constructor, takes the pointer to the "bulk" element, the 
  /// index of the fixed local coordinate and its value represented
  /// by an integer indicating which face we're on.
  /// Optional final arg is the identifier for the new values created
  /// by this face element
  //===========================================================================
  template <class ELEMENT>
    NavierStokesWithSingularityFaceElement<ELEMENT>::
    NavierStokesWithSingularityFaceElement(ELEMENT* const& bulk_el_pt, 
					     const int& face_index, 
					     const unsigned& id) : 
  FaceGeometry<ELEMENT>(), FaceElement(), Boundary_id(id),
    /* C_external_data_index_at_knot(0), */ No_throw_if_no_singular_elem_pt(false),
    Exact_traction_fct_pt(nullptr), Dxi_dx_fct_pt(nullptr)
  {
    // Let the bulk element build the FaceElement, i.e. setup the pointers 
    // to its nodes (by referring to the appropriate nodes in the bulk
    // element), etc.
    bulk_el_pt->build_face_element(face_index, this);
 
#ifdef PARANOID
    {
      //Check that the element is not a refineable 3d element
      ELEMENT* elem_pt = dynamic_cast<ELEMENT*>(bulk_el_pt);
      //If it's three-d
      if(elem_pt->dim() == 3)
      {
	//Is it refineable
	RefineableElement* ref_el_pt = dynamic_cast<RefineableElement*>(elem_pt);
	if(ref_el_pt != nullptr)
	{
	  if (this->has_hanging_nodes())
	  {
	    throw OomphLibError(
	      "This face element will not work correctly if nodes are hanging\n",
	      OOMPH_CURRENT_FUNCTION,
	      OOMPH_EXCEPTION_LOCATION);
	  }
	}
      }
    }
#endif   

    // Extract the dimension of the problem from the dimension of 
    // the first node
    // QUEHACERES this fixes the weird link error when building with -O3
    Dim = this->node_pt(0)->ndim();

    // Set up P_index_nst. Initialise to Dim, (since we have Dim velocity components indexed
    // from zero, followed by the pressure) which probably won't change
    // in most cases, oh well, the price we pay for generality
    P_index_nst = Dim;

    // Cast to the appropriate NavierStokesEquation so that we can
    // find the index at which the variable is stored
    // We assume that the dimension of the full problem is the same
    // as the dimension of the node, if this is not the case you will have
    // to write custom elements, sorry
    switch(Dim)
    {
      //One dimensional problem
      case 1:
      {
	throw OomphLibError("1D fluid dynamics makes no sense!\n",
			    OOMPH_CURRENT_FUNCTION,
			    OOMPH_EXCEPTION_LOCATION);
      }
      break;
    
      //Two dimensional problem
      case 2:
      {
	NavierStokesEquations<2>* eqn_pt = 
	  dynamic_cast<NavierStokesEquations<2>*>(bulk_el_pt);
	//If the cast has failed die
	if(eqn_pt == nullptr)
	{
	  std::string error_string =
	    "Bulk element must inherit from NavierStokesEquations.";
	  error_string += 
	    "Nodes are two dimensional, but cannot cast the bulk element to\n";
	  error_string += "NavierStokesEquations<2>\n.";
	  error_string += 
	    "If you desire this functionality, you must implement it yourself\n";
       
	  throw OomphLibError(error_string,
			      OOMPH_CURRENT_FUNCTION,
			      OOMPH_EXCEPTION_LOCATION);
	}
	else
	{
	  //Read the index from the (cast) bulk element.
	  P_index_nst = eqn_pt->p_nodal_index_nst();
	}
      }
      break;
    
      //Three dimensional problem
      case 3:
      {
	NavierStokesEquations<3>* eqn_pt = 
	  dynamic_cast<NavierStokesEquations<3>*>(bulk_el_pt);
	
	//If the cast has failed die
	if(eqn_pt == nullptr)
	{
	  std::string error_string =
	    "Bulk element must inherit from NavierStokesEquations.";
	  error_string += 
	    "Nodes are three dimensional, but cannot cast the bulk element to\n";
	  error_string += "NavierStokesEquations<3>\n.";
	  error_string += 
	    "If you desire this functionality, you must implement it yourself\n";
       
	  throw OomphLibError(error_string,
			      OOMPH_CURRENT_FUNCTION,
			      OOMPH_EXCEPTION_LOCATION);       
	}
	else
	{
	  //Read the index from the (cast) bulk element.
	  P_index_nst = eqn_pt->p_nodal_index_nst();
	}
      }
      break;

      //Any other case is an error
      default:
	std::ostringstream error_stream; 
	error_stream <<  "Dimension of node is " << Dim 
		     << ". It should be 2, or 3!" << std::endl;
     
	throw OomphLibError(error_stream.str(),
			    OOMPH_CURRENT_FUNCTION,
			    OOMPH_EXCEPTION_LOCATION);
	break;
    }
  }


  
  ////////////////////////////////////////////////////////////////////////
  ////////////////////////////////////////////////////////////////////////
  ////////////////////////////////////////////////////////////////////////

  template<unsigned DIM, unsigned NNODE_1D>
    class QSingularDragIntegralDiskElement : public QPoissonElement<DIM, NNODE_1D>
  {
  public:

    // shorthand for a singular line element which provides all the
    // singular function goodies
    typedef ScalableSingularityForNavierStokesLineElement<NNODE_1D> SingularLineElement;

    // Function pointer to a function which provides the outer unit normal
    // at a given set of Lagrangian coordinates
    typedef void (*OuterUnitNormalFctPt)(const LagrangianCoordinates&, Vector<double>&);
    
    // Default constructor
  QSingularDragIntegralDiskElement() :
    Sign_of_outer_unit_normal(1), Outer_unit_normal_fct_pt(nullptr)
    {
      // set the integration scheme of this element to be our custom
      // mixed Gauss/Cebyshev-Gauss scheme
      this->set_integration_scheme(&Mixed_integration_scheme);

      // set the nodal dimension, since this is a surface embedded in a higher
      // dimensional domain
      this->set_nodal_dimension(DIM+1);
	
      // QUEHACERES do some zeta shit
    }
    
    QSingularDragIntegralDiskElement(const QSingularDragIntegralDiskElement&)
    {
      BrokenCopy::broken_copy("QSingularDragIntegralDiskElement");
    }

    void operator=(const QSingularDragIntegralDiskElement&)
    {
      BrokenCopy::broken_assign("QSingularDragIntegralDiskElement");
    }
    
    /// Override this with a broken version to prevent these elements
    /// being used accidentally for the problems residuals
    inline void fill_in_contribution_to_residuals(Vector<double>& residuals)
    {
      std::stringstream error;
      error << "These elements should be used only for computing drag integrals "
	    << "and should not be contributing directly to the problem's residuals\n";
	  
      throw OomphLibError(error.str(),
			  OOMPH_CURRENT_FUNCTION,
			  OOMPH_EXCEPTION_LOCATION);
    }

    // ========================================================================
    // set the line element and local coordinate for the ith knot
    void set_line_element_and_local_coordinate_at_knot(
      const Vector<std::pair<GeomObject*, Vector<double> > >&
      line_element_and_local_coordinate_at_knot)
    {
      // number of knot points in this element
      unsigned nknot = this->integral_pt()->nweight();

      // number of coordinates supplied
      unsigned ncoords = line_element_and_local_coordinate_at_knot.size();

      // check the right number of coordinates have been passed in
      if(ncoords != nknot)
      {
	ostringstream error_message;
	error_message << "Number of element-coordinate pairs provided is not consistent "
		      << "with the number of knots.\n"
		      << "Number of supplied coordinate sets: " << ncoords << "\n"
		      << "Number of knots:                    " << nknot;
	
	throw OomphLibError(
	  error_message.str(),
	  OOMPH_CURRENT_FUNCTION,
	  OOMPH_EXCEPTION_LOCATION);
      }

      // make enough space
      Line_element_and_local_coordinate_at_knot.resize(nknot);
      
      // set 'em
      for(unsigned i=0; i<nknot; i++)
      {
	Line_element_and_local_coordinate_at_knot[i] =
	  line_element_and_local_coordinate_at_knot[i];
      }
    }

    // set the edge coordinates of each knot point
    void set_lagr_coordinates_at_knot(const Vector<LagrangianCoordinates>& coords)
    {
      // Number of integration points
      const unsigned n_intpt = this->integral_pt()->nweight();
      
      if(coords.size() != n_intpt)
      {
	ostringstream error_message;
	error_message << "number of sets of coordinates provided is not consistent "
		      << "with number of integration points in this element\n"
		      << "Number of supplied coordinate sets: " << coords.size() << "\n"
		      << "Number of integration points:       " << n_intpt;
	throw OomphLibError(
	  error_message.str(),
	  OOMPH_CURRENT_FUNCTION,
	  OOMPH_EXCEPTION_LOCATION);
      }

      // make sure we've got enough space
      Lagrangian_coordinates_at_knot.resize(n_intpt);

      // set 'em
      for(unsigned i=0; i<n_intpt; i++)
      {
	Lagrangian_coordinates_at_knot[i]  = coords[i];
      }
    }

    // get the edge coordinates at the ith knot point
    LagrangianCoordinates lagr_coordinate_at_knot(const unsigned& i) const
    {
      return Lagrangian_coordinates_at_knot[i];
    }

    OuterUnitNormalFctPt& outer_unit_normal_fct_pt()
    {
      return Outer_unit_normal_fct_pt;
    }
      
    void outer_unit_normal(const LagrangianCoordinates& lagr_coords,
			   Vector<double>& unit_normal) const
    {
      // do we actually have a function pointer?
      if(Outer_unit_normal_fct_pt == nullptr)
      {
	throw OomphLibError(
	  "Outer unit normal function pointer hasn't been set for this drag integration element",
	  OOMPH_CURRENT_FUNCTION,
	  OOMPH_EXCEPTION_LOCATION);
      }
      else
      {
	// get the unit normal from the function pointer
	Outer_unit_normal_fct_pt(lagr_coords, unit_normal);
	
	// apply the correct sign (to account for upper/lower side of disk)
	for(auto& ni : unit_normal)
	  ni *= (double)(Sign_of_outer_unit_normal);
      }
    }

    /// \short main point of this element - function to compute this element's
    /// contribution to the singular hydrodynamic drag
    Vector<double> compute_total_singular_drag() const
    {
      Vector<double> drag(DIM+1, 0.0);
      
      // local coordinates of knot
      Vector<double> s(DIM, 0.0);

      for(unsigned ipt=0; ipt<this->integral_pt()->nweight(); ipt++)
      {
	//Assign values of s
	for(unsigned i=0; i<DIM; i++)
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

	// get the Lagrangian coordinates of this point
	LagrangianCoordinates lagr_coords = Lagrangian_coordinates_at_knot[ipt];

	// cast the GeomObject to a singular line element      
	SingularLineElement* sing_el_pt = dynamic_cast<SingularLineElement*>
	  (Line_element_and_local_coordinate_at_knot[ipt].first);

	// local coordinate in the singular element for the zeta of this knot
	Vector<double> s_singular_el =
	  Line_element_and_local_coordinate_at_knot[ipt].second;

	// the sum of all scaled singular functions
	Vector<double> u_sing_total(DIM+2, 0.0);

	DenseMatrix<double> stress_sing_total(DIM+1, DIM+1, 0.0);
	sing_el_pt->total_singular_stress(lagr_coords, s_singular_el, stress_sing_total);

	// get the outer unit normal
	Vector<double> unit_normal(DIM+1, 0.0);
	outer_unit_normal(lagr_coords, unit_normal);

	// now compute traction
	Vector<double> traction(DIM+1, 0.0);
	for(unsigned i=0; i<DIM+1; i++)
	{
	  for(unsigned j=0; j<DIM+1; j++)
	  {
	    traction[i] = stress_sing_total(i,j) * unit_normal[j];
	  }
	}

	// add the contribution of the traction at this knot to the integral
	for(unsigned i=0; i<DIM+1; i++)
	  drag[i] += traction[i] * W;
      }

      return drag;
    }
    
    // if we're on the lower side of the disk, flip the sign of the unit normal
    void is_lower_disk_element()
    {
      Sign_of_outer_unit_normal = -1;
    }
    
  private:

    // mixed Gauss/Chebyshev-Gauss integration scheme
    static TwoDGaussTensorProductChebyshevGauss<NNODE_1D, NNODE_1D> Mixed_integration_scheme;

    /// \short Edge coordinates (\rho, \zeta, \phi) of each of this element's
    /// knot points
    Vector<LagrangianCoordinates> Lagrangian_coordinates_at_knot;

    /// \short The line element and its local coordinate(s) which correspond to the zeta
    /// values at each knot point in this bulk element. N.B. the element
    /// is 1D so only has 1 local coordinate, but we keep this as a vector
    /// for consistency with the general case
    Vector<std::pair<GeomObject*, Vector<double> > >
      Line_element_and_local_coordinate_at_knot;

    /// \short these elements sit on the surface of the disk, so cannot determine
    /// the sign of the unit normal from the Lagrangian coordinates
    int Sign_of_outer_unit_normal;

    /// \short Function pointer to a function which computes the outer unit normal
    /// at given Lagrangian coordinates
    OuterUnitNormalFctPt Outer_unit_normal_fct_pt;
      
    // ********
    // QUEHACERES start here - need same machinery as with the face elements,
    // i.e. need to store sing elem and local coords for each knot point
    
  };

  ////////////////////////////////////////////////////////////////////////
  ////////////////////////////////////////////////////////////////////////
  ////////////////////////////////////////////////////////////////////////
  
  //======================================================================
  /// \short A class for elements that imposes Dirichlet boundary 
  /// conditions on complete solution (such that u_fe + C u_sing = u_bc) using a
  /// Lagrange multiplier. Thus the element introduce an additional
  /// unknown at the nodes it's attached to. C and u_sing are specified
  /// via a ScalableSingularityForNavierStokesElement.
  //======================================================================
  template <class ELEMENT>
    class NavierStokesWithSingularityBCFaceElement : 
    public virtual NavierStokesWithSingularityFaceElement<ELEMENT>
  {
    
  public:
    
    // shorthand typedef for the baseclass typedef 
      typedef typename
	NavierStokesWithSingularityFaceElement<ELEMENT>::SingularLineElement SingularLineElement;
      
      /// \short Constructor, takes the pointer to the "bulk" element and the 
      /// index of the face to which the element is attached. Optional final
      /// arg is the identifier for the additional unknowns multiplier
      NavierStokesWithSingularityBCFaceElement(ELEMENT* const& bulk_el_pt, 
					       const int& face_index,					       
					       const unsigned& id = 0,
					       const bool& is_on_upper_disk_surface = false); 
  
      ///\short  Broken empty constructor
      NavierStokesWithSingularityBCFaceElement()
      {
	throw OomphLibError(
	  "Don't call empty constructor for NavierStokesWithSingularityBCFaceElement",
	  OOMPH_CURRENT_FUNCTION,
	  OOMPH_EXCEPTION_LOCATION);
      }
  
      /// Broken copy constructor
      NavierStokesWithSingularityBCFaceElement(
	const NavierStokesWithSingularityBCFaceElement& dummy) 
      { 
	BrokenCopy::broken_copy("NavierStokesWithSingularityBCFaceElement");
      } 
  
      /// Broken assignment operator
      void operator=(const NavierStokesWithSingularityBCFaceElement&) 
	{
	  BrokenCopy::broken_assign("NavierStokesWithSingularityBCFaceElement");
	}

      /// Add the element's contribution to its residual vector
      inline void fill_in_contribution_to_residuals(Vector<double>& residuals)
      {
	//Call the generic residuals function with flag set to 0
	//using a dummy matrix argument
	fill_in_generic_residual_contribution_navier_stokes_bc(
	  residuals, GeneralisedElement::Dummy_matrix, 0);
      }
      
#ifndef USE_FD_JACOBIAN
      /// \short Add the element's contribution to its residual vector and its
      /// Jacobian matrix
      inline void fill_in_contribution_to_jacobian(Vector<double>& residuals,
      						   DenseMatrix<double>& jacobian)
      {
      	//Call the generic routine with the flag set to 1
      	fill_in_generic_residual_contribution_navier_stokes_bc(residuals, jacobian, 1);
      }
#endif
      
      /// Output function
      void output(std::ostream& outfile) const
      {
	const unsigned n_plot = 5;
	output(outfile, n_plot);
      }

      /// \short Output function
      void output(std::ostream& outfile, const unsigned& nplot) const
      {
	//Vector of local coordinates
	Vector<double> s(this->Dim-1);
   
	// Tecplot header info
	outfile << this->tecplot_zone_string(nplot);
   
	// Loop over plot points
	unsigned num_plot_points = this->nplot_points(nplot);
	
	for (unsigned iplot=0; iplot<num_plot_points; iplot++)
	{
	  // Get local coordinates of plot point
	  this->get_s_plot(iplot, nplot, s);
     
	  Vector<double> x(this->Dim);
	  for(unsigned i=0; i<this->Dim; i++) 
	  {
	    x[i] = this->interpolated_x(s,i);
	    outfile << x[i] << " ";
	  }

	  // get the local bulk coordinates    
	  Vector<double> s_bulk = this->local_coordinate_in_bulk(s);

	  // pointer to the bulk element this face element is attached to
	  ELEMENT* bulk_el_pt = dynamic_cast<ELEMENT*>(this->bulk_element_pt());
	  
	  Vector<double> u(3, 0.0);
	  bulk_el_pt->interpolated_u_nst(s_bulk, u);
	  
	  for(unsigned i=0; i<this->Dim; i++)
	  {
	    outfile << u[i] << " ";
	  }
	  
	  outfile << this->interpolated_p_nst(s) << std::endl;
	}

	// Write tecplot footer (e.g. FE connectivity lists)
	this->write_tecplot_zone_footer(outfile, nplot);
      }

      /// \short Provide nodal values of desired boundary values.
      /// They're imposed by Lagrange multipliers.
      void set_nodal_boundary_values(const DenseMatrix<double>& nodal_boundary_value)
      {
#ifdef PARANOID
	if (nodal_boundary_value.nrow() != this->nnode())
	{
	  std::stringstream error;
	  error << "nodel_boundary_value is a matrix with " 
		<< nodal_boundary_value.nrow() 
		<< " rows, but should have the same number of rows as the number of nodes, "
		<< this->nnode();
	  throw OomphLibError(error.str(),
			      OOMPH_CURRENT_FUNCTION,
			      OOMPH_EXCEPTION_LOCATION);
	}
#endif
	Nodal_boundary_value = nodal_boundary_value;
      }

      // QUEHACERES for debug, output the value of the Lagrange multipliers on the Dirichlet
      // boundaries and the associated stress
      void output_lagrange_multiplier_and_stress(std::ostream& outfile,
						 const unsigned& nplot) const
      {
	ostringstream error_message;
	  error_message << "Don't use this at the moment, broken since we switched over to "
			<< "line singularity - this element contains a vector of singular "
			<< "elements and local coordinates at *knots*, but for this function"
			<< "we need those elements and coordinates at the plot points";
	  
	throw OomphLibError(error_message.str(),
	  OOMPH_CURRENT_FUNCTION,
	      OOMPH_EXCEPTION_LOCATION);
      }

      // function to get the total normal stress acting on this element
      // (useful for convergence studies of total force, etc.)
      Vector<double> get_contribution_to_normal_stress() const;

      bool is_on_upper_disk_surface() const
      {
	return Is_on_upper_disk_surface;
      }

      void set_upper_disk_surface()
      {
	Is_on_upper_disk_surface = true;
      }
	
    private:

      /// \short Add the element's contribution to its residual vector.
      /// flag=1(or 0): do (or don't) compute the contribution to the
      /// Jacobian as well. 
      void fill_in_generic_residual_contribution_navier_stokes_bc(
	Vector<double>& residuals, DenseMatrix<double>& jacobian, 
	const unsigned& flag);

      /// Desired boundary values at nodes
      DenseMatrix<double> Nodal_boundary_value;

      /// are we attached to the upper or lower surface?
      bool Is_on_upper_disk_surface;
      
  }; // end of NavierStokesWithSingularityBCFaceElement class 

  //////////////////////////////////////////////////////////////////////
  //////////////////////////////////////////////////////////////////////
  //////////////////////////////////////////////////////////////////////
  
  //===========================================================================
  /// Constructor, takes the pointer to the "bulk" element, the 
  /// index of the fixed local coordinate and its value represented
  /// by an integer indicating which face we're on.
  /// Optional final arg is the identifier for the new values created
  /// by this face element
  //===========================================================================
  template <class ELEMENT>
    NavierStokesWithSingularityBCFaceElement<ELEMENT>::
    NavierStokesWithSingularityBCFaceElement(ELEMENT* const& bulk_el_pt, 
					     const int& face_index,					     
					     const unsigned& id,
					     const bool& is_on_upper_disk_surface) : 
  NavierStokesWithSingularityFaceElement<ELEMENT>(bulk_el_pt, face_index, id),
    Is_on_upper_disk_surface(is_on_upper_disk_surface)
  { 
    unsigned n_node = this->nnode();

    // Make space for Dim Lagrange multipliers
    Vector<unsigned> n_additional_values(n_node, this->Dim);

    // add them (this takes care of duplicate IDs so no additional
    // checks needed here to avoid redundant nodal values)
    this->add_additional_values(n_additional_values, this->Boundary_id);

    // flip the sign of the normal, since outer_unit_normal() gives the normal
    // which points out of the bulk element we're attached to, but since
    // we're attached to the disk we want to point into the bulk element,
    // i.e. if we're on the top surface of the disk the bulk element is above
    // and so the default is for the normal to point downwards, but the upper
    // disk normal points upwards
    this->normal_sign() = -this->normal_sign();
    
  } // end NavierStokesWithSingularityBCFaceElement constructor


  template <class ELEMENT>
  Vector<double> NavierStokesWithSingularityBCFaceElement<ELEMENT>::
    get_contribution_to_normal_stress() const
  {
    //Find out how many nodes there are
    const unsigned n_node = this->nnode();

    // shorthand
    const unsigned Dim = this->Dim;
    
    //Set up memory for the shape and test functions
    Shape psi(n_node), test(n_node);
     
    // The number of integration points (knots)
    const unsigned n_intpt = this->integral_pt()->nweight();
     
    //Set the Vector to hold local coordinates
    Vector<double> s(Dim-1);
	
    // shorthand
    ELEMENT* bulk_el_pt = dynamic_cast<ELEMENT*>(this->bulk_element_pt());

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

      this->strain_rate(s, strain_rate_fe);

      // get the FE pressure
      double p_fe = this->interpolated_p_nst(s);
	  
      // get contribution to the FE stress on this element
      DenseMatrix<double> stress_fe(Dim, Dim, 0.0);

      stress_fe = (*bulk_el_pt->stress_fct_pt())(strain_rate_fe, p_fe);

      // Now get the total singular contributions
      // --------------------------------------------------------

      // get the line element and local coordinate which corresponds to the
      // singular amplitude for this knot
      std::pair<GeomObject*, Vector<double> > line_elem_and_local_coord = 
	this->line_element_and_local_coordinate_at_knot(ipt);

      // cast the GeomObject to a singular line element      
      SingularLineElement* sing_el_pt =
	dynamic_cast <SingularLineElement*>( line_elem_and_local_coord.first );

      // local coordinate in the singular element for the zeta of this knot
      Vector<double> s_singular_el = line_elem_and_local_coord.second;
      
      // get the \rho,\zeta,\phi coordinates at this knot
      LagrangianCoordinates lagr_coords_at_knot = this->lagr_coordinate_at_knot(ipt);

      // the sum of all scaled singular functions
      Vector<double> u_sing_total(Dim+1, 0.0);

      // total singular contributions to velocity gradient tensor and strain rate      
      DenseMatrix<double> strain_rate_sing_total(Dim, Dim, 0.0);
      
      if(sing_el_pt != nullptr)
      {
	u_sing_total = sing_el_pt->total_singular_contribution(lagr_coords_at_knot,
							       s_singular_el);

	sing_el_pt->total_singular_strain_rate(lagr_coords_at_knot,
					       s_singular_el,
					       strain_rate_sing_total);
      }
      else
      {
	oomph_info << "Weird shit happened, we've got a BC face element with no "
		   << "associated singular element\n" << std::endl;
	abort();
      }
      
      // total singular pressure
      double p_sing_total = u_sing_total[this->P_index_nst];

      // total singular stress
      DenseMatrix<double> stress_sing_total(Dim, Dim, 0.0);
      stress_sing_total = (*bulk_el_pt->stress_fct_pt())(strain_rate_sing_total,
							 p_sing_total);

      // add the FE and singular bits together to get the total stress
      DenseMatrix<double> stress_total(Dim, Dim, 0.0);
      
      for(unsigned i=0; i<Dim; i++)
      {
	for(unsigned j=0; j<Dim; j++)
	{
	  stress_total(i,j) = stress_fe(i,j) + stress_sing_total(i,j);
	}
      }
      
      // now add the weighted contribution to the normal stress to the
      // integral over this element
      for(unsigned i=0; i<Dim; i++)
      {
	for(unsigned j=0; j<Dim; j++)
	{	      
	  total_force[i] += stress_total(i,j) * unit_normal[j] * W;
	}
      }
    } // end loop over integration points

    return total_force;
  }

  
  //===========================================================================
  /// \short Compute the element's residual vector and the Jacobian matrix.
  /// Adds this boundary face element's contribution to the equation which
  /// determines the Lagrange multipliers, and adds the Lagrange multiplier
  /// contribution to the bulk equations
  //===========================================================================
  template <class ELEMENT>
    void NavierStokesWithSingularityBCFaceElement<ELEMENT>::
    fill_in_generic_residual_contribution_navier_stokes_bc(
      Vector<double>& residuals, DenseMatrix<double>& jacobian, 
      const unsigned& flag) 
  {
    // shorthands
    const unsigned Dim = this->Dim;
    const unsigned Boundary_id = this->Boundary_id;

    // pointer to the bulk element this face element is attached to
    ELEMENT* bulk_el_pt = dynamic_cast<ELEMENT*>(this->bulk_element_pt());
    
    //Find out how many nodes there are
    const unsigned n_node = this->nnode();
    const unsigned n_node_bulk = bulk_el_pt->nnode();
    
    //Set up memory for the shape and test functions
    Shape psi(n_node), test(n_node);

    Shape psi_bulk(n_node_bulk);
	  
    DShape dpsidx(n_node_bulk, Dim);
    
    //Set the value of Nintpt
    const unsigned n_intpt = this->integral_pt()->nweight();
     
    //Set the Vector to hold local coordinates
    Vector<double> s(Dim-1);
  
    // loop over the integration points
    //--------------------------------
    for(unsigned ipt=0; ipt<n_intpt; ipt++)
    {       
      //Assign values of s
      for(unsigned i=0; i<(Dim-1); i++) 
      {
	s[i] = this->integral_pt()->knot(ipt, i);
      }

      // get the local bulk coordinates    
      Vector<double> s_bulk = this->local_coordinate_in_bulk(s);
      
      // get the integral weight
      double w = this->integral_pt()->weight(ipt);
       
      // find the shape and test functions and return the Jacobian
      //of the mapping from Eulerian -> local coords
      double J = this->shape_and_test(s, psi, test);
	
      // get the derivatives of the shape functions from the bulk element      
      bulk_el_pt->dshape_eulerian(s_bulk, psi_bulk, dpsidx);
      
      // premultiply the weights and the Jacobian
      double W = w * J;
      
      // calculate stuff at integration point
      Vector<double> u_fe(Dim, 0.0);
      Vector<double> u_bc(Dim, 0.0);	

      // get the interpolated BC-enforcing Lagrange mutliplier field at
      // this Gauss point
      Vector<double> lambda_bc(Dim, 0.0);

      // get the momentum-enforcing Lagrange multiplier field from the bulk
      Vector<double> lambda_momentum = bulk_el_pt->interpolated_lambda(s_bulk);
      
      for(unsigned l=0; l<n_node; l++)
      {
	// grab a pointer to the current node
	Node* node_pt = this->node_pt(l);

	// get the map which gives the starting nodal index for
	// the Lagrange multipliers associated with each boundary ID
	std::map<unsigned, unsigned> first_index = *(
	  dynamic_cast<BoundaryNodeBase*>(node_pt)->
	  index_of_first_value_assigned_by_face_element_pt() );
		  
	for(unsigned i=0; i<Dim; i++)
	{
	  // get the nodal values of the FE solution and the boundary conditions
	  u_fe[i] += this->nodal_value(l,i)    * psi[l];
	  u_bc[i] += Nodal_boundary_value(l,i) * psi[l];

	  unsigned lambda_bc_index = first_index.at(this->Boundary_id) + i;
	  lambda_bc[i] += this->nodal_value(l, lambda_bc_index) * psi[l];
	}	  
      }

      // cast the GeomObject to a singular line element      
      SingularLineElement* sing_el_pt = nullptr;
      
      LagrangianCoordinates lagr_coords_at_knot;
      Vector<double> s_singular_el(1, 0.0);
      
      // get the edge coordinates at this knot, the corresponding singular
      // line element and it's local coordinates
      this->lagr_coords_and_singular_element_at_knot(ipt,						       
						     lagr_coords_at_knot,
						     sing_el_pt,
						     s_singular_el);

       // Stuff related to singular fct
      Vector<double> u_sing(Dim+1, 0.0);

      // the total contribution of the singular velocity gradients
      DenseMatrix<double> dudx_sing_total(Dim, Dim, 0.0);

      // get contribution of total singular pressure and
      // total singular strain-rate to total stress tensor
      DenseMatrix<double> stress_sing_total(Dim, Dim, 0.0);
      
      // check we've actually got a singular line element pointer	
      if (sing_el_pt == nullptr)
      {
	if(!this->No_throw_if_no_singular_elem_pt)
	{
	  throw OomphLibError(
	    "This NavierStokesWithSingularityBCFaceElement doesn't have a singular line element pointer",
	    OOMPH_CURRENT_FUNCTION,
	    OOMPH_EXCEPTION_LOCATION);
	}
      }
      else
      {
	// total singular contribution
	u_sing = sing_el_pt->total_singular_contribution(lagr_coords_at_knot,
							 s_singular_el);
	
	// the total contribution of the singular velocity gradients
	dudx_sing_total = sing_el_pt->
	  total_singular_gradient_contribution(lagr_coords_at_knot,
					       s_singular_el);

	// total singular pressure
	double p_sing_total = u_sing[this->P_index_nst];

	// total singular contribution to the strain-rate
	DenseMatrix<double> strain_rate_sing_total(Dim, Dim, 0.0);

	// compute the strain-rate from the velocity gradients,
	// \epsion_{ij} = 1/2 (du_i/dx_j + du_j/dx_i)
	for (unsigned i=0; i<Dim; i++)
	{
	  for(unsigned j=0; j<Dim; j++)
	  {
	    strain_rate_sing_total(i,j) = 0.5*(dudx_sing_total(i,j) + dudx_sing_total(j,i));
	  }
	}

#ifdef PARANOID
	if(bulk_el_pt->stress_fct_pt() == 0)
	{
	  throw OomphLibError(
	    "Error: the stress function pointer has not been set for the augmented elements\n",
	    OOMPH_CURRENT_FUNCTION,
	    OOMPH_EXCEPTION_LOCATION);
	}
#endif 
	stress_sing_total =
	  (*bulk_el_pt->stress_fct_pt())(strain_rate_sing_total, p_sing_total);
      }

      // ======================================================================
      //Now add to the appropriate equations

      
      // contribution to the singular amplitude equations
      // --------------------------------------------------
      
      // get a list of the singular function IDs this element has
      Vector<unsigned> sing_ids = sing_el_pt->singular_fct_ids();

      const unsigned nnode_sing = sing_el_pt->nnode();

      // get the shape functions from the singular element at the local coordinate
      // which corresponds to this knot in this face element
      Shape psi_sing(nnode_sing);
      Shape test_dummy(nnode_sing);
      sing_el_pt->shape_and_test(s_singular_el, psi_sing, test_dummy);
      
      // loop over the nodes in the singular element associated with this integration point
      for(unsigned ising_node=0; ising_node<nnode_sing; ising_node++)
      {	
	// external data index for this singular node
	unsigned ext_index = this->C_external_data_index_at_knot[ipt][ising_node];

	// loop over the singular functions
	for(unsigned ising=0; ising < sing_el_pt->nsingular_fct(); ising++)
	{
	  // get the ID
	  unsigned sing_fct_id = sing_ids[ising];

	  // and get the unscaled singular function associated with this ID
	  Vector<double> u_sing_unscaled =
	    sing_el_pt->unscaled_singular_fct(lagr_coords_at_knot, sing_fct_id);

	  // external equation number which determines the singular amplitude
	  int external_eqn_c = this->external_local_eqn(ext_index, ising);

	  // if this singular amplitude isn't pinned
	  if(external_eqn_c >= 0)	      
	  {	  
	    for(unsigned i=0; i<Dim; i++)
	    {		
	      residuals[external_eqn_c] +=
		lambda_bc[i] * u_sing_unscaled[i] * psi_sing[ising_node] * W;

	      // Jacobian?
	      if (flag == 1)
	      {		  
		//Loop over the test functions in this face element
		for(unsigned l=0; l<n_node; l++)
		{
		  // grab a pointer to the current node
		  Node* node_pt = this->node_pt(l);
	  
		  // get the map which gives the starting nodal index for
		  // the Lagrange multipliers associated with each boundary ID
		  std::map<unsigned, unsigned> first_index = *(
		    dynamic_cast<BoundaryNodeBase*>(node_pt)->
		    index_of_first_value_assigned_by_face_element_pt() );
	    
		  // get the nodal index of the Lagrange multiplier for this
		  // coordinate direction and boundary ID
		  unsigned lambda_bc_index = first_index.at(Boundary_id) + i;

		  // get the local Lagrange multiplier equation number 
		  int local_eqn_lagr = this->nodal_local_eqn(l, lambda_bc_index);

		  if(local_eqn_lagr >= 0)
		  {
		    // dC/dlambda-tilde
		    jacobian(external_eqn_c, local_eqn_lagr) +=
		      u_sing_unscaled[i] * psi[l] * psi_sing[ising_node] * W;

		    // symmetric jacobian entry, so chuck it in here too
		    jacobian(local_eqn_lagr, external_eqn_c) +=
		      u_sing_unscaled[i] * psi[l] * psi_sing[ising_node] * W;
		  }		    
		}
	      }
	    } // end loop over dimensions
	  } // end pinned c check
	} // end loop over singular functions
	
      }

      // Now do BC-enforcing Lagrange multiplier
      // ---------------------------------------------------------------
      
      // get the outer unit normal on this boundary point
      Vector<double> unit_normal(Dim, 0.0);
      this->outer_unit_normal(s, unit_normal);

      // FE traction in augmented region along Dirchlet boundary
      Vector<double> traction_fe(Dim, 0.0);
      
      // get the FE traction along this Dirchlet boundary
      bulk_el_pt->get_traction(s_bulk, unit_normal, traction_fe);
	
      //Loop over the test functions
      for(unsigned l=0; l<n_node; l++)
      {
	// grab a pointer to the current node
	Node* node_pt = this->node_pt(l);
	  
	// get the map which gives the starting nodal index for
	// the Lagrange multipliers associated with each boundary ID
	std::map<unsigned, unsigned> first_index = *(
	  dynamic_cast<BoundaryNodeBase*>(node_pt)->
	  index_of_first_value_assigned_by_face_element_pt() );

	// loop over the coordinate directions
	for(unsigned d=0; d<Dim; d++)
	{
	  // get the nodal index of the Lagrange multiplier for this
	  // coordinate direction and boundary ID
	  unsigned lambda_bc_index = first_index[Boundary_id] + d;

	  // get the local Lagrange multiplier equation number 
	  int local_eqn_lagr = this->nodal_local_eqn(l, lambda_bc_index);
	      
	  // equation for FE velocity
	  int local_eqn_u_fe = this->nodal_local_eqn(l, d);
	    
#ifdef PARANOID
	  // Lagrange multiplier active but u_fe pinned won't work!
	  if ( (local_eqn_lagr >= 0) && (local_eqn_u_fe < 0) )
	  {
	    throw OomphLibError(
	      "Lagrange multiplier active but u_fe pinned won't work!",
	      OOMPH_CURRENT_FUNCTION,
	      OOMPH_EXCEPTION_LOCATION);
	  }
#endif	    
	    	      
	  // Lagrange multiplier for BC residual: It's determined by enforcing
	  // that u_fe + C u_sing = u_bc
	  if(local_eqn_lagr >= 0)
	  {	   
	    residuals[local_eqn_lagr] += (u_fe[d] + u_sing[d] - u_bc[d]) * psi[l]*W;
	    
	    // Jacobian?
	    if (flag == 1)
	    {
	      for(unsigned l2=0; l2<n_node; l2++)
	      {
		// QUEHACERES again, get this index more systematically
		int local_unknown_u_fe = this->nodal_local_eqn(l2, d);
		if (local_unknown_u_fe >= 0)
		{
		  // dlambda_tilde/du
		  jacobian(local_eqn_lagr, local_unknown_u_fe) += psi[l2] * psi[l]*W;
		}
	      }
	    }
	  }

	  // ========================================================
	  // Contribution of Lagrange multipliers to augmented velocity eqn:
	  // ========================================================
	  if (local_eqn_u_fe >= 0)
	  { 
	    residuals[local_eqn_u_fe] += lambda_bc[d] * psi[l] * W;

	    if (flag == 1)
	    {
	      for(unsigned l2=0; l2<n_node; l2++)
	      {
		// grab a pointer to the second node
		Node* node2_pt = this->node_pt(l2);

		// get the map which gives the starting nodal index for
		// the Lagrange multipliers associated with each boundary ID
		std::map<unsigned, unsigned> first_index2 = *(
		  dynamic_cast<BoundaryNodeBase*>(node2_pt)->
		  index_of_first_value_assigned_by_face_element_pt() );
		
		// get the index of the Lagrange multiplier of the second node
		// associated with this face ID and direction 
		unsigned lambda_bc_index2 = first_index2[Boundary_id] + d;
		int local_unknown_lambda_bc = this->nodal_local_eqn(l2, lambda_bc_index2);
		    
		if (local_unknown_lambda_bc >= 0)
		{
		  // QUEHACERES write-up transpose of Eqn. 1.16
		  jacobian(local_eqn_u_fe, local_unknown_lambda_bc) += psi[l2] * psi[l] * W;
		}   
	      }
	    }	// end Jacobian check    
	  }

	  // no contributions from momentum-enforcing LMs, since they are zero
	  // on Dirchlet boundaries
	  
	} // end loop over directions	
      } // end loop over nodes

      // Boundary contributions of momentum LMs to the pressure residual
      // --------------------------------------------------------------------

      // number of pressure nodes in this face element
      const unsigned n_pres = Dim;
      
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
	unsigned k_in_bulk = this->bulk_node_number(k);

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
	      lambda_momentum[j] * unit_normal[j] * psip[k_in_bulk] * W;
	  }

	  // no Jacobian contributions to momentum-enforcing LMs,
	  // since they are zero on Dirchlet boundaries
	  
	}
      }
    } // end loop over integration points
  } // end of fill_in_generic_residual_contribution_navier_stokes_bc()

////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////

//======================================================================
/// \short A class for FaceElements that impose the jump in flux
/// at the interface between regions where a "singular" solution is
/// added to the FE solution and where it's not. The element attaches
/// itself to the faces of the outermost bulk elements in which
/// the singular solution is added. To facilitate the imposition 
/// in the jump in the FE solution, the FaceElement creates a new
/// set of nodes and makes the bulk element point to these. 
///
/// NOTE: Also need to update node pointers of any other bulk elements
///       that are not formally part of this boundary but have isolated
///       nodes on it! Use the existing_duplicate_node_pt map to figure
///       which nodes must be replaced (only for elements in region
///       that contains singularity)
/// 
/// Element retains pointers to the original nodes (which are now only
/// pointed to by the adjacent element in which the singular solution 
/// is not added. Continutity of the solution such that 
/// u_fe + u_sing on the "left" = u_fe on the "right" is enforced
/// by a Lagrange multiplier (stored as additional nodal Data in the
/// FaceElement; the flux is imposed by adding he appropriate contribution
/// to the bulk residual in the element on the "left". 
/// The final complication arises because the amplitude of the
/// "singular" solution is in general an unknown which is handled
/// by a ScalableSingularityForPoissonElement, so we store a pointer to that.
/// The data in that element acts as external Data for the current
/// element.
//======================================================================
  template <class ELEMENT>
    class NavierStokesWithSingularityStressJumpFaceElement :
    public NavierStokesWithSingularityFaceElement<ELEMENT>
  {
 
  public:

    // shorthand typedef for the baseclass typedef 
    typedef typename
      NavierStokesWithSingularityFaceElement<ELEMENT>::SingularLineElement SingularLineElement;
    
    /// \short Constructor, takes the pointer to the "bulk" element and the 
    /// index of the face to which the element is attached. 
    /// Map keeps a running count of duplicate nodes already created;
    /// existing_duplicate_node_pt[orig_node_pt]=new_node_pt.
    /// Optional final arg is the identifier for the lagrange multiplier
    NavierStokesWithSingularityStressJumpFaceElement(
      ELEMENT* const& augmented_bulk_el_pt, 
      const int& face_index,
      ELEMENT* const& non_augmented_bulk_el_pt,
      std::map<Node*,Node*>& existing_duplicate_node_pt,
      const unsigned& lambda_hat_id,
      const unsigned& lambda_hat_hat_id); 

    ///\short  Broken empty constructor
    NavierStokesWithSingularityStressJumpFaceElement()
    {
      std::string error_string = "Don't call empty constructor for ";
      error_string+="NavierStokesWithSingularityStressJumpFaceElement";
      throw OomphLibError(error_string,
			  OOMPH_CURRENT_FUNCTION,
			  OOMPH_EXCEPTION_LOCATION);
    }

    /// Broken copy constructor
    NavierStokesWithSingularityStressJumpFaceElement(
      const NavierStokesWithSingularityStressJumpFaceElement& dummy) 
    { 
      BrokenCopy::broken_copy("NavierStokesWithSingularityStressJumpFaceElement");
    } 
 
    /// Broken assignment operator
    void operator=(const NavierStokesWithSingularityStressJumpFaceElement&) 
      {
	BrokenCopy::broken_assign("NavierStokesWithSingularityStressJumpFaceElement");
      }

    /// Add the element's contribution to its residual vector
    inline void fill_in_contribution_to_residuals(Vector<double>& residuals)
    {
      //Call the generic residuals function with flag set to 0
      //using a dummy matrix argument
      fill_in_generic_residual_contribution_nst_sing_jump(
	residuals, GeneralisedElement::Dummy_matrix, 0);
    }
    
#ifndef USE_FD_JACOBIAN
      
    /// \short Add the element's contribution to its residual vector and its
    /// Jacobian matrix
    inline void fill_in_contribution_to_jacobian(Vector<double>& residuals,
						 DenseMatrix<double>& jacobian)
    {
      //Call the generic routine with the flag set to 1
      fill_in_generic_residual_contribution_nst_sing_jump
	(residuals, jacobian, 1);
    }
    
#endif
      
    /// Output function
    void output(std::ostream& outfile) const
    {
      const unsigned n_plot=5;
      output(outfile,n_plot);
    }

    /// \short Output function
    void output(std::ostream& outfile, const unsigned& nplot) const
    {
      // shorthand
      const unsigned Dim = this->Dim;
	  
      // Dimension of element 
      const unsigned el_dim = this->dim();
    
      //Vector of local coordinates
      Vector<double> s(el_dim);

      // Number of nodes
      unsigned n_node = this->nnode();
      Shape psi(n_node);

      // Tecplot header info
      outfile << this->tecplot_zone_string(nplot);
    
      // Loop over plot points
      unsigned num_plot_points = this->nplot_points(nplot);
      for (unsigned iplot=0; iplot<num_plot_points; iplot++)
      {      
	// Get local coordinates/shape fcts at plot point
	this->get_s_plot(iplot, nplot, s);
      
	this->shape(s, psi);
      
	//Calculate stuff at integration point
	Vector<double> x(Dim, 0.0);
	Vector<double> u_augmented(Dim, 0.0);
	Vector<double> u_bulk(Dim, 0.0);

	// get the interpolated Lagrange mutliplier field at this point,
	// taking into account any other LM fields which contribute to this one
	Vector<double> lambda(Dim, 0.0);

	for(unsigned l=0; l<n_node; l++) 
	{
	  // grab a pointer to the current node
	  Node* node_pt = this->node_pt(l);

	  // get the map which gives the starting nodal index for
	  // the Lagrange multipliers associated with each boundary ID
	  std::map<unsigned, unsigned> first_index = *(
	    dynamic_cast<BoundaryNodeBase*>(node_pt)->
	    index_of_first_value_assigned_by_face_element_pt() );
	
	  for(unsigned i=0; i<Dim; i++)
	  {
	    u_augmented[i]  += this->nodal_value(l,i) * psi[l];
	    u_bulk[i] += Orig_node_pt[l]->value(i) * psi[l];

	    x[i] += this->nodal_position(l,i) * psi[l];

	    // get the nodal index, accounting for the dimension offset
	    unsigned lambda_index = first_index[this->Boundary_id] + i;

	    lambda[i] += this->nodal_value(l, lambda_index) * psi[l];
	  }
	}
	  
	for(unsigned i=0; i<Dim; i++) 
	{
	  outfile << x[i] << " ";
	}
	for(unsigned i=0; i<Dim; i++) 
	{
	  outfile << u_augmented[i] << " ";
	}
	for(unsigned i=0; i<Dim; i++) 
	{
	  outfile << u_bulk[i] << " ";
	}
	for(unsigned i=0; i<Dim; i++) 
	{
	  outfile << u_bulk[i] - u_augmented[i] << " ";
	}
	for(unsigned i=0; i<Dim; i++) 
	{
	  outfile << lambda[i] << " ";
	}
	
	outfile << std::endl;
	  
      } // end loop over plot points
	
    
      // Write tecplot footer (e.g. FE connectivity lists)
      this->write_tecplot_zone_footer(outfile, nplot);    
    }

    void output_traction_and_lagrange_multipliers(std::ofstream& outfile) const
    {
      // shorthand
      const unsigned Dim = this->Dim;
	  
      // Dimension of element 
      const unsigned el_dim = this->dim();
    
      //Vector of local coordinates
      Vector<double> s(el_dim, 0.0);

      // Number of nodes
      const unsigned n_node = this->nnode();

      // make space for the shape and test functions
      Shape psi(n_node), test(n_node);;
      
      // number of integration points
      const unsigned n_intpt = this->integral_pt()->nweight();
      
      //Loop over the integration points
      //--------------------------------
      for(unsigned ipt=0; ipt<n_intpt; ipt++)
      {
	// get local coordinates of this integration point
	for(unsigned i=0; i<(this->Dim-1); i++)
	{
	  s[i] = this->integral_pt()->knot(ipt,i);
	}
   
	//Calculate stuff at integration point
	Vector<double> x(Dim, 0.0);
	Vector<double> u_augmented(Dim, 0.0);
	Vector<double> u_bulk(Dim, 0.0);
	
	// get the interpolated Lagrange mutliplier field at this Gauss point
	Vector<double> lambda(Dim, 0.0);

	Vector<double> traction_augmented(Dim, 0.0);
	Vector<double> traction_bulk(Dim, 0.0);

	// get the outer unit normal on this boundary point
	Vector<double> unit_normal(Dim, 0.0);
	this->outer_unit_normal(s, unit_normal);

	// pointer to the bulk element this face element is attached to
	ELEMENT* bulk_el_pt = dynamic_cast<ELEMENT*>(this->bulk_element_pt());

	// Get the local bulk coordinates    
	Vector<double> s_bulk = this->local_coordinate_in_bulk(s);
        
	// now compute the Lagrange multipliers
	for(unsigned l=0; l<n_node; l++) 
	{
	  // grab a pointer to the current node
	  Node* node_pt = this->node_pt(l);

	  // get the map which gives the starting nodal index for
	  // the Lagrange multipliers associated with each boundary ID
	  std::map<unsigned, unsigned> first_index = *(
	    dynamic_cast<BoundaryNodeBase*>(node_pt)->
	    index_of_first_value_assigned_by_face_element_pt() );
	
	  for(unsigned i=0; i<Dim; i++)
	  {	    
	    x[i] += this->nodal_position(l,i) * psi[l];

	    unsigned lambda_index = first_index.at(this->Boundary_id) + i;
	    lambda[i] += this->nodal_value(l, lambda_index) * psi[l];
	  }
	}

	// get the traction onto the augmented region	
	bulk_el_pt->get_traction(s_bulk, unit_normal, traction_augmented);

	// now get the traction in the bulk
	// ----------------------------------

	// first, find the local coordinates of this point in the
	// bulk element
	Vector<double> s_bulk_region(Dim, 0.0);
	
	GeomObject* geom_object_pt = nullptr;

	this->non_augmented_bulk_elem_pt()->locate_zeta(x, geom_object_pt, s_bulk_region);

	// now get the traction (passing in the same unit normal as before,
	// so will need to flip the sign of the traction)
	this->non_augmented_bulk_element_pt()->get_traction(s_bulk_region,
							    unit_normal,
							    traction_bulk);
	
	for(unsigned i=0; i<Dim; i++) 
	{
	  outfile << x[i] << ",";
	}
	for(unsigned i=0; i<Dim; i++) 
	{
	  outfile << traction_augmented[i] << ",";
	}
	for(unsigned i=0; i<Dim; i++)
	{
	  outfile << traction_bulk[i] << ",";
	}
	for(unsigned i=0; i<Dim; i++) 
	{
	  outfile << lambda[i] << ",";
	}
	
	outfile << std::endl;
	  
      } // end loop over knot points
    }
 
    /// Pin Lagrange multipliers and set to zero
    void pin_lagrange_multipliers_and_set_to_zero() const
    {
      unsigned nnod = this->nnode();
      for (unsigned j=0; j<nnod; j++)
      {
	// grab a pointer to the current node
	Node* node_pt = this->node_pt(j);
	  
	// get the map which gives the starting nodal index for
	// the Lagrange multipliers associated with each boundary ID
	std::map<unsigned, unsigned> first_index = *(
	  dynamic_cast<BoundaryNodeBase*>(node_pt)->
	  index_of_first_value_assigned_by_face_element_pt() );
	
	for(unsigned i=0; i<this->Dim; i++)
	{
	  // get the nodal index, accounting for the dimension offset
	  unsigned lambda_index = first_index[this->Boundary_id] + i;

	  node_pt->pin(lambda_index);
	  node_pt->set_value(lambda_index, 0.0);
	}	  
      }
    }

    // function to return the pointer to the non-augmented bulk element
    // (not the bulk element in the augmented region to which this element is attached)
    ELEMENT*& non_augmented_bulk_element_pt()
    {
      return Non_augmented_bulk_element_pt;
    }

    // switch on the warning about reducing locate zeta tolerance for finding
    // this face element's knot points in the (non-augmented) bulk element
    void warn_about_reducing_tolerance_for_bulk_locate_zeta()
    {
      Shut_up_about_reducing_tolerance_for_bulk_locate_zeta = false;
    }
      
  private:   

    /// \short Add the element's contribution to its residual vector.
    /// flag=1(or 0): do (or don't) compute the contribution to the
    /// Jacobian as well. 
    void fill_in_generic_residual_contribution_nst_sing_jump(
      Vector<double>& residuals, DenseMatrix<double>& jacobian, 
      const unsigned& flag);
  
    /// Vector of pointers to orig nodes
    Vector<Node*> Orig_node_pt;

    /// \short Pointer to the non-augmented element which also shares this face
    /// on the augmented region boundary
    ELEMENT* Non_augmented_bulk_element_pt;

    // This is the external index for the non-augmented nodes which
    // correspond directly to the augmented nodes of this face element.
    // This vector is indexed by the node number in this face element.
    Vector<unsigned> External_data_index_for_non_aug_node;

    // This is the external data index for all the nodes in the non-augmented
    // bulk element corresponding to this face element. This vecto is
    // indexed by the node number in the non-augmented bulk element.
    Vector<unsigned> External_data_index_all_non_aug_nodes;

    /// \short Map which takes the index of a node in the bulk element
    /// to which this face element is attached, and returns the external data
    /// index for this face element. These bulk nodes make contributions to
    /// this elements Jacobian because it depends on derivatives of bulk
    /// quantities (Lagrange multipliers)
    std::map<unsigned, unsigned> External_data_index_bulk_aug_node_map;

    /// \short ID given to the additional nodal values added to store
    /// the Lagrange multipliers which enforce continuity of the
    /// Lagrange multipliers which enforce the governing momentum PDEs
    unsigned Lambda_hat_hat_id;

    /// Suppress warning message about having to reduce the locate zeta
    /// tolerance after a failed attempt to locate one of this face element's
    /// knot points in the corresponding (non-augmented) bulk element
    bool Shut_up_about_reducing_tolerance_for_bulk_locate_zeta;
  }; 

  //===========================================================================
  /// Constructor, takes the pointer to the "bulk" element, the 
  /// index of the fixed local coordinate and its value represented
  /// by an integer indicating where the face is located.
  /// Map keeps a running count of duplicate nodes already created;
  /// existing_duplicate_node_pt[orig_node_pt]=new_node_pt.
  /// Optional final arg is the identifier for the lagrange multiplier
  //===========================================================================
  template <class ELEMENT>
    NavierStokesWithSingularityStressJumpFaceElement<ELEMENT>::
    NavierStokesWithSingularityStressJumpFaceElement(
      ELEMENT* const& augmented_bulk_el_pt, 
      const int& face_index,
      ELEMENT* const& non_augmented_bulk_el_pt,
      std::map<Node*,Node*>& existing_duplicate_node_pt,
      const unsigned& lambda_hat_id,
      const unsigned& lambda_hat_hat_id) : 
  NavierStokesWithSingularityFaceElement<ELEMENT>(augmented_bulk_el_pt,
						  face_index,
						  lambda_hat_id),
    Non_augmented_bulk_element_pt(non_augmented_bulk_el_pt),
    Lambda_hat_hat_id(lambda_hat_hat_id), 
    Shut_up_about_reducing_tolerance_for_bulk_locate_zeta(true)
  {   
    // Back up original nodes and make new ones
    unsigned nnod = this->nnode();
    Orig_node_pt.resize(nnod);
    External_data_index_for_non_aug_node.resize(nnod);
      
    // have we definitely got a non-augmented bulk element pointer?
    if(non_augmented_bulk_el_pt == nullptr)
    {
      throw OomphLibError("Can't pass empty non-augmented bulk element pointer to this "
			  "stress-jump face element\n",
			    OOMPH_CURRENT_FUNCTION,
			    OOMPH_EXCEPTION_LOCATION);
    }

    // make enough space 
    unsigned nnode_non_aug = Non_augmented_bulk_element_pt->nnode();
    External_data_index_all_non_aug_nodes.resize(nnode_non_aug, 0);

    for (unsigned j=0; j<nnod; j++)
    {
      // Here's the node that used to be shared between the
      // adjacent bulk elements. Note: this may be updated
      // below if it turns out that this bulk node
      // was already replaced.
      Node* nod_pt = this->node_pt(j);
      Orig_node_pt[j] = nod_pt;
     
      // Find this original node in the map; if we find it
      // it's already been duplicated earlier when creating another
      // FaceElement; in that case use that duplicate rather than
      // creating another one.
      std::map<Node*,Node*>::iterator it = existing_duplicate_node_pt.find(nod_pt);

      bool is_replacement = false;
      Node* orig_for_replaced_node_pt = nullptr;
      if (it != existing_duplicate_node_pt.end())
      {
	// Use the existing duplicate node
	this->node_pt(j) = (*it).second;
      }
      // Make a new one (as boundary node) hierher update comment when done
      else // @@
      {       
	// See if it is a replacement (can happen if a very 
	// deformed bulk element has multiple faces on this boundary)
	for (std::map<Node*,Node*>::iterator it =
	       existing_duplicate_node_pt.begin();
	     it != existing_duplicate_node_pt.end(); it++)
        {
	  if ((*it).second == nod_pt)
          {
	    is_replacement = true;
	    orig_for_replaced_node_pt = (*it).first;
	    Orig_node_pt[j] = orig_for_replaced_node_pt;
	    break;
          }
        }
       
	// Stick with the existing replacement
	if (is_replacement)
        {
	  this->node_pt(j) = nod_pt;
        }
	// Make new node
	else //--
        {
	  unsigned n_dim = nod_pt->ndim();
	  unsigned n_position_type = nod_pt->nposition_type();
	  unsigned n_value = nod_pt->nvalue();
	  this->node_pt(j) = new BoundaryNode<Node>(n_dim, n_position_type, n_value);
         
	  // It has the same coordinate; hierher add history values too
	  // when imlementing as Navier Stokes
	  for (unsigned i=0; i<n_dim; i++)
          {
	    this->node_pt(j)->x(i) = nod_pt->x(i);
          }
         
	  // ...and the same values
	  for (unsigned i=0; i<n_value; i++)
          {
	    this->node_pt(j)->set_value(i, nod_pt->value(i));
          }
         
	  // It is on the same boundaries
	  std::set<unsigned>* boundaries_pt;
	  nod_pt->get_boundaries_pt(boundaries_pt);
	  for (std::set<unsigned>::iterator it = boundaries_pt->begin();
	       it != boundaries_pt->end(); it++)
          {
	    // Get/set boundary ID
	    unsigned new_boundary_id = (*it);
	    this->node_pt(j)->add_to_boundary(new_boundary_id);
           
	    // Get/set boundary coordinates
	    if (nod_pt->boundary_coordinates_have_been_set_up())
            {
	      // QUEHACERES come back to this!
	      /* unsigned n = nod_pt->ncoordinates_on_boundary(new_boundary_id); */
	      /* Vector<double> boundary_zeta(n); */
	      /* nod_pt->get_coordinates_on_boundary(new_boundary_id, */
	      /* 					  boundary_zeta); */
	      /* this->node_pt(j)->set_coordinates_on_boundary(new_boundary_id, */
	      /* 					      boundary_zeta); */
            }
	    else
            {
	      // hierher throw? (Doesn't happen at the moment, i.e. 
	      // when this diagnostic was finally commented out)
             
	      /*   oomph_info << "No boundary coordinates have been set up" */
	      /*              << " for new local node " << j  */
	      /*              << " at : " */
	      /*              << node_pt(j)->x(0) << " "  */
	      /*              << node_pt(j)->x(1) << " "  */
	      /*              << node_pt(j)->x(2) << " "  */
	      /*              << std::endl; */
            }
          }
	  
	  // Copy across index map for additional values
	  std::map<unsigned, unsigned>* index_pt =
	    dynamic_cast<BoundaryNodeBase*>(nod_pt)->
	    index_of_first_value_assigned_by_face_element_pt();
	  
	  if (index_pt != nullptr)
          {
	    std::map<unsigned, unsigned>* new_index_pt =
	      new std::map<unsigned, unsigned>;
	    
	    dynamic_cast<BoundaryNodeBase*>(this->node_pt(j))->
	      index_of_first_value_assigned_by_face_element_pt() =
	      new_index_pt;
	    
	    for (std::map<unsigned, unsigned>::iterator it =
		   (*index_pt).begin(); it != (*index_pt).end(); it++)
            {
	      (*new_index_pt)[(*it).first] = (*it).second;
            }           
          }
         
	  // Keep track 
	  existing_duplicate_node_pt[nod_pt] = this->node_pt(j);
         
	  // Switch over node for bulk element that we attached ourselves
	  // to
	  unsigned j_in_bulk = this->bulk_node_number(j);
	  augmented_bulk_el_pt->node_pt(j_in_bulk) = this->node_pt(j);
         
        } // end existing node is already replacement vs make new one
      } 

      // The original node now acts as external data for this element
      // (we still need it to enforce continuity)
      External_data_index_for_non_aug_node[j] = this->add_external_data(Orig_node_pt[j]);
    }

    // Now loop over the nodes in the corresponding non-augmented bulk element
    // and add its nodes as external data, since they all contribute to the
    // non-augmented traction which this elements residuals depends on        
    for(unsigned j=0; j<nnode_non_aug; j++)
    {
      External_data_index_all_non_aug_nodes[j] =
	this->add_external_data(Non_augmented_bulk_element_pt->node_pt(j));
    }
    
    Vector<Node*> non_boundary_aug_nodes;
    Vector<unsigned> node_nums_in_bulk(this->nnode());
    
    for(unsigned j=0; j<this->nnode(); j++)
      node_nums_in_bulk[j] = this->bulk_node_number(j);

    for(unsigned j=0; j<augmented_bulk_el_pt->nnode(); j++)
    {
      // skip if this node is on this face
      if(std::find(node_nums_in_bulk.begin(), node_nums_in_bulk.end(), j)
	 != node_nums_in_bulk.end() )
	continue;
    
      // if it isn't on this face, then add it as external data
      // and add the external index to our list
      External_data_index_bulk_aug_node_map[j] = 
	this->add_external_data(augmented_bulk_el_pt->node_pt(j));
    }

    // double check we've got the right number
#ifdef PARANOID
    
    // number of nodes in the bulk element we're attached to which aren't
    // shared with this face face element
    unsigned nnode_non_boundary_aug = augmented_bulk_el_pt->nnode() - this->nnode();

    if(External_data_index_bulk_aug_node_map.size() != nnode_non_boundary_aug)
    {
      throw OomphLibError("Wrong number of augmented non-face bulk nodes added "
			  "as external data",
			  OOMPH_CURRENT_FUNCTION,
			  OOMPH_EXCEPTION_LOCATION);
    }
#endif
    
    // Make space for Dim Lagrange multipliers which enforce continuity of
    // the solution across the augmented boundary
    Vector<unsigned> n_additional_values(nnod, this->Dim);
    this->add_additional_values(n_additional_values, lambda_hat_id);

    // ...and add another Dim LMs to handle the continuity of the
    // momentum-enforcing LMs
    this->add_additional_values(n_additional_values, lambda_hat_hat_id);
    
  } // end NavierStokesWithSingularityStressJumpFaceElement constructor

  //===========================================================================
  /// Compute the element's residual vector and the Jacobian matrix.
  //===========================================================================
  template <class ELEMENT>
    void NavierStokesWithSingularityStressJumpFaceElement<ELEMENT>::
    fill_in_generic_residual_contribution_nst_sing_jump(
      Vector<double>& residuals, DenseMatrix<double>& jacobian, 
      const unsigned& flag)
  {
    // shorthand
    ELEMENT* bulk_el_pt = dynamic_cast<ELEMENT*>(this->bulk_element_pt());
    
    //Find out how many nodes there are
    const unsigned n_node = this->nnode();
    const unsigned n_node_bulk = bulk_el_pt->nnode();
    
    // shorthand
    const unsigned Dim = this->Dim;

    //Set up memory for the shape and test functions and their derivatives
    Shape psi(n_node), test(n_node);

    Shape psi_bulk(n_node_bulk);
	  
    DShape dpsidx(n_node_bulk, Dim);

    // pressure basis functions from the bulk element
    Shape psip(bulk_el_pt->npres_nst());    
      
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

      // Get the shape and test functions and return
      // the Jacobian of the mapping
      double J = this->shape_and_test(s, psi, test);

      // get the derivatives of the shape functions from the bulk element      
      bulk_el_pt->dshape_eulerian(s_bulk, psi_bulk, dpsidx);

      // get the pressure basis functions from the bulk element
      bulk_el_pt->pshape_nst(s_bulk, psip);
      
      //Premultiply the weights and the Jacobian
      double W = w * J;

      // the velocity fields in the augmented and non-augmented (bulk)
      // regions, and the Lagrange multipliers which enforce continuity
      // of the solution across the boundary of the augmented region
      Vector<double> u_augmented(Dim, 0.0);
      Vector<double> u_non_aug(Dim, 0.0);

      // interpolated continuity-enforcing Lagrange mutliplier field at this point
      Vector<double> lambda_cont(Dim, 0.0);
      
      // Eulerian coordinates
      Vector<double> x(Dim, 0.0);

      // LM which enforces continuity of the momentum-enforcing LMs
      Vector<double> lambda_hat_hat(Dim, 0.0);
      
      // loop over the nodes and compute the above at this integration point
      for(unsigned l=0; l<n_node; l++) 
      {
	// grab a pointer to the current node
	Node* node_pt = this->node_pt(l);

	// get the map which gives the starting nodal index for
	// the Lagrange multipliers associated with each boundary ID
	std::map<unsigned, unsigned> first_index = *(
	  dynamic_cast<BoundaryNodeBase*>(node_pt)->
	  index_of_first_value_assigned_by_face_element_pt() );
	
	for(unsigned i=0; i<Dim; i++)
	{
	  u_augmented[i] += this->nodal_value(l,i) * psi[l];
	  u_non_aug[i]   += Orig_node_pt[l]->value(i) * psi[l];
	  
	  // interpolate the Eulerian coordinates
	  x[i] += node_pt->x(i) * psi[l];
	  
	  unsigned lambda_cont_index = first_index.at(this->Boundary_id) + i;

	  lambda_cont[i] += this->nodal_value(l, lambda_cont_index) * psi[l];
	  
	  // get the nodal index of the Lagrange multiplier for this
	  // coordinate direction and boundary ID
	  unsigned lambda_hat_hat_index = first_index.at(Lambda_hat_hat_id) + i;

	  // interpolated lambda_hat_hat
	  lambda_hat_hat[i] += this->nodal_value(l, lambda_hat_hat_index) * psi[l];
	}
      }

      // cast the GeomObject to a singular line element      
      SingularLineElement* sing_el_pt = nullptr;

      LagrangianCoordinates lagr_coords_at_knot;
      Vector<double> s_singular_el(1, 0.0);

      // get the edge coordinates at this knot, the corresponding singular
      // line element and it's local coordinates
      this->lagr_coords_and_singular_element_at_knot(ipt,
						     lagr_coords_at_knot,
						     sing_el_pt,
						     s_singular_el);
      
      // check we've actually got one, or we're in trouble
      if(sing_el_pt == nullptr)
      {
	ostringstream error_message;

	error_message << "Error: this stress jump element has no "
		      << "singular line element pointer, so it "
		      << "cannot compute the correct jump in stress\n";
	    
	throw OomphLibError(error_message.str().c_str(),
			    OOMPH_CURRENT_FUNCTION,
			    OOMPH_EXCEPTION_LOCATION);
      }
      
      // get the list of singular function IDs
      const Vector<unsigned> sing_ids = sing_el_pt->singular_fct_ids();

      // unscaled stuff. These are stored in an array so that they can be
      // looped over when implementing analytic jacobian
      Vector<Vector<double> > u_sing_unscaled(sing_el_pt->nsingular_fct());
      Vector<DenseMatrix<double> > dudx_sing_unscaled(sing_el_pt->nsingular_fct());

      // unscaled singular pressures
      Vector<double> p_sing_unscaled(sing_el_pt->nsingular_fct());
      
      // the sum of all scaled singular functions
      Vector<double> u_sing_total =
	sing_el_pt->total_singular_contribution(lagr_coords_at_knot,
                                        	s_singular_el);

      // the total contribution of the singular velocity gradients
      DenseMatrix<double> dudx_sing_total = sing_el_pt->
	total_singular_gradient_contribution(lagr_coords_at_knot,
					     s_singular_el);
      
      // total singular pressure
      double p_sing_total = u_sing_total[this->P_index_nst];
      
      // loop over all the singular functions this element knows about and
      // compute the sum of their contributions to the total solution      
      for(unsigned ising=0; ising < sing_el_pt->nsingular_fct(); ising++)
      {
	// get the ID
	unsigned sing_fct_id = sing_ids[ising];
	
	// unscaled stuff
	// ----------------
	// unscaled singular function 
	u_sing_unscaled[ising] =	    
	  sing_el_pt->unscaled_singular_fct(lagr_coords_at_knot, sing_fct_id);	  
	 
	// unscaled singular gradient 
	dudx_sing_unscaled[ising] =
	  sing_el_pt->gradient_of_unscaled_singular_fct(lagr_coords_at_knot, sing_fct_id);
	 
	// unscaled singular pressure 
	p_sing_unscaled[ising] = u_sing_unscaled[ising][this->P_index_nst];
	
      }      

      // Compute outer unit normal at the specified local coordinate
      // to compute scaled and unscaled flux of singular solution
      Vector<double> unit_normal(Dim, 0.0);
      this->outer_unit_normal(s, unit_normal);

      // total singular contribution to the strain-rate
      DenseMatrix<double> strain_rate_sing_total(Dim, Dim, 0.0);

      // array of strain rate tensors associated with each singular function
      Vector<DenseMatrix<double> > strain_rate_sing_unscaled
	(sing_el_pt->nsingular_fct(), DenseMatrix<double>(Dim, Dim, 0.0));

      // compute the strain-rate from the velocity gradients,
      // \epsion_{ij} = 1/2 (du_i/dx_j + du_j/dx_i)
      for (unsigned i=0; i<Dim; i++)
      {
	for(unsigned j=0; j<Dim; j++)
	{
	  strain_rate_sing_total(i,j) = 0.5*(dudx_sing_total(i,j) + dudx_sing_total(j,i));

	  // and compute the individual contributions of each singular function
	  for(unsigned ising=0; ising<sing_el_pt->nsingular_fct(); ising++)
	  {
	    strain_rate_sing_unscaled[ising](i,j) +=
	      0.5*(dudx_sing_unscaled[ising](i,j) + dudx_sing_unscaled[ising](j,i));
	  }
	}
      }
	
      // get contribution of total singular pressure and
      // total singular strain-rate to total stress tensor
      DenseMatrix<double> stress_sing_total(Dim, Dim);

#ifdef PARANOID
      if(bulk_el_pt->stress_fct_pt() == nullptr)
      {
	throw OomphLibError(
	  "Error: the stress function pointer has not been set for the augmented elements\n",
	  OOMPH_CURRENT_FUNCTION,
	  OOMPH_EXCEPTION_LOCATION);
      }
#endif      
      stress_sing_total =
	(*bulk_el_pt->stress_fct_pt())(strain_rate_sing_total, p_sing_total);

      // get stress associated with each singular function
      Vector<DenseMatrix<double> > stress_sing_unscaled(sing_el_pt->nsingular_fct());
      
      for(unsigned ising=0; ising < sing_el_pt->nsingular_fct(); ising++)
      {
	stress_sing_unscaled[ising] =
	  (*bulk_el_pt->stress_fct_pt())(strain_rate_sing_unscaled[ising],
					 p_sing_unscaled[ising]);
      }

      // get the momentum-enforcing Lagrange multiplier field from the bulk
      Vector<double> lambda_momentum = bulk_el_pt->interpolated_lambda(s_bulk);
            
      // now get the traction from the non-augmented side
      // -------------------------------------------------

      // first, find the local coordinates of this point in the
      // non-augmented bulk element
      Vector<double> s_bulk_non_aug(Dim, 0.0);
	
      GeomObject* geom_object_pt = nullptr;
	
      Non_augmented_bulk_element_pt->locate_zeta(x, geom_object_pt, s_bulk_non_aug);

      if (geom_object_pt == nullptr)
      {
	if(!Shut_up_about_reducing_tolerance_for_bulk_locate_zeta)
	{
	  oomph_info << "Warning: didn't find point: ";
	
	  for(double xi : x)
	    oomph_info << xi << " ";

	  oomph_info << "\nin the non-augmented bulk element corresponding to "
		     << "this StressJumpFaceElement with "
		     << "Locate_zeta_helpers::Newton_tolerance = "
		     << Locate_zeta_helpers::Newton_tolerance << ".\n"
		     << "Trying again with tolerance: 1e-6..." << std::endl;
	}
	
	// backup the current tolerance
	double locate_zeta_tol_backup = Locate_zeta_helpers::Newton_tolerance;

	// try again with a slightly looser tolerance
	Locate_zeta_helpers::Newton_tolerance = 1e-6;

	Non_augmented_bulk_element_pt->locate_zeta(x, geom_object_pt, s_bulk_non_aug);

	// if it still doesn't converge, then throw
	if (geom_object_pt == nullptr)
	{
	  ostringstream error_message;

	  error_message << "Couldn't find corresponding local coordinates "
			<< "in the non-augmented bulk element corresponding "
			<< "to this StressJumpFaceElement\n for Eulerian coordinates: ";

	  for(double xi : x)
	    error_message << xi << " ";
	
	  error_message << "\n\n Nodes in augmented element: \n";

	  for(unsigned n=0; n<bulk_el_pt->nnode(); n++)
	  {
	    Node* node_pt = bulk_el_pt->node_pt(n);

	    for(unsigned i=0; i<Dim; i++)
	    {
	      error_message << node_pt->x(i) << " ";
	    }
	    error_message << "\n";
	  }
	  error_message << "\nNodes in non-augmented element: \n";

	  for(unsigned n=0; n<Non_augmented_bulk_element_pt->nnode(); n++)
	  {
	    Node* node_pt = Non_augmented_bulk_element_pt->node_pt(n);

	    for(unsigned i=0; i<Dim; i++)
	    {
	      error_message << node_pt->x(i) << " ";
	    }
	    error_message << "\n";
	  }
	
	  throw OomphLibError(error_message.str().c_str(),
			      OOMPH_CURRENT_FUNCTION,
			      OOMPH_EXCEPTION_LOCATION);
	}

	// reset the tolerance
	Locate_zeta_helpers::Newton_tolerance = locate_zeta_tol_backup;
      }
      
      // now get the traction from the non-augmented bulk element
      Vector<double> traction_non_aug(Dim, 0.0);
      Non_augmented_bulk_element_pt->get_traction(s_bulk_non_aug,
					       unit_normal,
					       traction_non_aug);

      Vector<double> lambda_momentum_non_aug = 
	Non_augmented_bulk_element_pt->interpolated_lambda(s_bulk_non_aug);
      
      Vector<double> traction_sing_total(Dim, 0.0);

      for(unsigned i=0; i<Dim; i++)
      {
	for(unsigned j=0; j<Dim; j++)
	{	  
	  traction_sing_total[i] += stress_sing_total(i,j) * unit_normal[j];
	}
      }

      // get non-augmented shape functions and derivatives 
      Shape psi_non_aug(Non_augmented_bulk_element_pt->nnode());
      Shape psip_non_aug(Non_augmented_bulk_element_pt->npres_nst());
      DShape dpsidx_non_aug(Non_augmented_bulk_element_pt->nnode(), Dim);
      
      Non_augmented_bulk_element_pt->
	dshape_eulerian(s_bulk_non_aug, psi_non_aug, dpsidx_non_aug);

      Non_augmented_bulk_element_pt->pshape_nst(s_bulk_non_aug, psip_non_aug);
      
      // ======================================================================
      // Now add to the appropriate equations
      // ======================================================================

      
      // =================================================
      // contribution to the singular amplitude equations
      // =================================================
      
      const unsigned nnode_sing = sing_el_pt->nnode();

      // get the shape functions and their Eulerian derivatives from the
      // singular element at the local coordinate which corresponds to this
      // knot in this face element
      Shape psi_sing(nnode_sing);
      Shape test_dummy(nnode_sing);
      DShape dpsi_sing_dx(nnode_sing, Dim);

      sing_el_pt->shape_and_test(s_singular_el, psi_sing, test_dummy);
      sing_el_pt->dshape_eulerian(lagr_coords_at_knot, s_singular_el, dpsi_sing_dx);

      // loop over the nodes in the singular element associated with this integration point
      for(unsigned ising_node=0; ising_node<nnode_sing; ising_node++)
      {
	// external data index for this singular node
	unsigned ext_index = this->C_external_data_index_at_knot[ipt][ising_node];

	// loop over the singular functions
	for(unsigned ising=0; ising < sing_el_pt->nsingular_fct(); ising++)
	{
	  // external equation number which determines the singular amplitude
	  int external_eqn_c = this->external_local_eqn(ext_index, ising);
	  
	  // if this singular amplitude isn't pinned
	  if(external_eqn_c >= 0)
	  {
	    for(unsigned i=0; i<Dim; i++)
	    {
	      // contribution of the continuity-enforcing LM to the C equations
	      residuals[external_eqn_c] +=
		lambda_cont[i] * u_sing_unscaled[ising][i] * psi_sing[ising_node] * W;

	      // Jacobian?
	      if (flag == 1)
	      {
		//Loop over the test functions
		for(unsigned l=0; l<n_node; l++)
		{
		  Node* node_pt = this->node_pt(l);

		  // get the map which gives the starting nodal index for
		  // the Lagrange multipliers associated with each boundary ID
		  std::map<unsigned, unsigned> first_index = *(
		    dynamic_cast<BoundaryNodeBase*>(node_pt)->
		    index_of_first_value_assigned_by_face_element_pt() );
	
		  // get the nodal index of the Lagrange multiplier for this
		  // coordinate direction and boundary ID
		  unsigned lambda_cont_index = first_index.at(this->Boundary_id) + i;
		  
		  int local_eqn_lagr = this->nodal_local_eqn(l, lambda_cont_index);
	  
		  if (local_eqn_lagr >= 0)
		  {		      
		    // dC/dlambda-hat
		    jacobian(external_eqn_c, local_eqn_lagr) +=
		      u_sing_unscaled[ising][i] * psi[l] * psi_sing[ising_node] * W;

		    // symmetric jacobian entry, so chuck it in here too
		    jacobian(local_eqn_lagr, external_eqn_c) +=
		      u_sing_unscaled[ising][i] * psi[l] * psi_sing[ising_node] * W;
		  }
		}
	      }
	    } // end loop over dimensions
	      
	  } // end external_eqn_c >= check
	} // end ext_index >= 0 check
	
      } //  end loop over singular nodes
      
      //Loop over the test functions
      for(unsigned l=0; l<n_node; l++)
      {
	Node* node_pt = this->node_pt(l);

	Node* nonaug_node_pt = Orig_node_pt[l];

	// get the map which gives the starting nodal index for
	// the Lagrange multipliers associated with each boundary ID
	std::map<unsigned, unsigned> first_index = *(
	  dynamic_cast<BoundaryNodeBase*>(node_pt)->
	  index_of_first_value_assigned_by_face_element_pt() );

	Vector<double> x_l(Dim, 0.0);
	for(unsigned i=0; i<Dim; i++)
	  x_l[i] = node_pt->x(i);
	
	for(unsigned i=0; i<Dim; i++)
	{
	  // =======================================================================
	  // Contributions to the bulk Lagrange multiplier equations which
	  // enforce the Stokes momentum PDEs
	  // =======================================================================
	  
	  // Index of the bulk Lagrange multiplier which enforces momentum PDEs 
	  unsigned lambda_momentum_index =
	    bulk_el_pt->index_of_lagrange_multiplier(node_pt, i);

	  int local_eqn_lagr_mom = this->nodal_local_eqn(l, lambda_momentum_index);
	  
	  if(local_eqn_lagr_mom >= 0)
	  {
	    residuals[local_eqn_lagr_mom] += lambda_hat_hat[i] * psi[l] * W;
	    
	    if (flag == 1)
	    {
	      // Now add contributions from lambda-hat-hat to lambda
	      for(unsigned l2=0; l2<n_node; l2++)
	      {
		Node* node2_pt = this->node_pt(l2);

		// get the map which gives the starting nodal index for
		// the Lagrange multipliers associated with each boundary ID
		std::map<unsigned, unsigned> first_index2 = *(
		  dynamic_cast<BoundaryNodeBase*>(node2_pt)->
		  index_of_first_value_assigned_by_face_element_pt() );
		
		unsigned lambda_hat_hat_index = first_index2.at(Lambda_hat_hat_id) + i;

		// get lambda-hat-hat equation number at this node
	      	int local_unknown = this->nodal_local_eqn(l2, lambda_hat_hat_index);
		
	      	if (local_unknown >= 0)
	      	{
	      	  jacobian(local_eqn_lagr_mom, local_unknown) +=
	      	    psi[l2] * psi[l] * W;
	      	}
	      }
	    
	    } // end Jacobian flag check	    
	  } // end momentum LM equations

	  // =======================================================================
	  // Continuity of momentum-enforcing Lagrange multiplier equations:
	  // Determined from continuity of the momentum-enforcing LMs from the
	  // augmented to non-augmented region
	  // =======================================================================
	  
	  // Index of the bulk Lagrange multiplier which enforces momentum PDEs	    
	  unsigned lambda_mom_non_aug_index =
	    Non_augmented_bulk_element_pt->index_of_lagrange_multiplier(nonaug_node_pt, i);

	  unsigned ext_index = External_data_index_for_non_aug_node[l];
	  
	  int local_eqn_lagr_mom_non_aug =
	    this->external_local_eqn(ext_index, lambda_mom_non_aug_index);
	  
	  if(local_eqn_lagr_mom_non_aug >= 0)
	  {
	    residuals[local_eqn_lagr_mom_non_aug] -= lambda_hat_hat[i] * psi[l] * W;

	    if(flag)
	    {	      
	      for(unsigned l2=0; l2<n_node; l2++)
	      {
		Node* node2_pt = this->node_pt(l2);

		// get the map which gives the starting nodal index for
		// the Lagrange multipliers associated with each boundary ID
		std::map<unsigned, unsigned> first_index2 = *(
		  dynamic_cast<BoundaryNodeBase*>(node2_pt)->
		  index_of_first_value_assigned_by_face_element_pt() );
		
		unsigned lambda_hat_hat_index = first_index2.at(Lambda_hat_hat_id) + i;

		// get lambda-hat-hat equation number at this node
	      	int local_unknown = this->nodal_local_eqn(l2, lambda_hat_hat_index);
		
	      	if (local_unknown >= 0)
	      	{
	      	  jacobian(local_eqn_lagr_mom_non_aug, local_unknown) -=
	      	    psi[l2] * psi[l] * W;
	      	}
	      }
	    }
	  }
	  
	  // =======================================================================
	  // Continuity Lagrange multiplier equations: Determined from continuity of
	  // solution with (scaled!) singular solution in the augmented region.
	  // =======================================================================
	  
	  // get the nodal index of the Lagrange multiplier for this
	  // coordinate direction and boundary ID
	  unsigned lambda_cont_index = first_index.at(this->Boundary_id) + i;

	  if(lambda_cont_index < Dim)
	  {
	    ostringstream error_message;

	    error_message << "wrong lambda_cont index! Apparently index for "
			  << "Boundary_id: " << this->Boundary_id << " is: "
			  << first_index[this->Boundary_id] << "\n";
	    
	    throw OomphLibError(error_message.str().c_str(),
				OOMPH_CURRENT_FUNCTION,
				OOMPH_EXCEPTION_LOCATION);
	  }
	  
	  int local_eqn_lagr = this->nodal_local_eqn(l, lambda_cont_index);
	  
	  if (local_eqn_lagr >= 0)
	  {
	    // Equation enforcing continuity of soln across augmented boundary
	    residuals[local_eqn_lagr] +=
	      ((u_augmented[i] + u_sing_total[i]) - u_non_aug[i]) * psi[l] * W;

	    // compute Jacobian
	    if (flag == 1)
	    {
	      for(unsigned l2=0; l2<n_node; l2++)
	      {               
		int local_unknown_augmented = this->nodal_local_eqn(l2, i);
		if (local_unknown_augmented >= 0)
		{
		  jacobian(local_eqn_lagr, local_unknown_augmented) +=
		    psi[l2] * psi[l] * W;
		}

		int local_unknown_non_aug = this->external_local_eqn(
		  External_data_index_for_non_aug_node[l2], i);
		
		if (local_unknown_non_aug >= 0)
		{
		  // QUEHACERES write-up Eqn. 1.18
		  jacobian(local_eqn_lagr, local_unknown_non_aug) -=
		    psi[l2] * psi[l] * W;
		}
	      }	      
	    } // end Jacobian flag check
	  }

	  // =======================================================================
	  // Add the equation which determines lambda-hat-hat, the LMs which
	  // enforce the continuity of the momentum-enforcing LMs
	  // =======================================================================
	  {
	    // get the nodal index from the map of additional values with the
	    // dimension offset
	    unsigned lambda_hat_hat_index = first_index.at(Lambda_hat_hat_id) + i;

	    // get the equation number for this nodal dof
	    int local_eqn_lambda_hat_hat = this->nodal_local_eqn(l, lambda_hat_hat_index);

	    // is it pinned?
	    if(local_eqn_lambda_hat_hat >= 0)
	    {
	      // equation is simply the difference between the momentum-enforcing LMs
	      // evaluated either side of the augmented boundary
	      residuals[local_eqn_lambda_hat_hat] +=
		(lambda_momentum[i] - lambda_momentum_non_aug[i]) * psi[l] * W;

	      if(flag == 1)
	      {
		for(unsigned l2=0; l2<n_node; l2++)
		{
		  Node* node2_pt = this->node_pt(l2);

		  Node* nonaug_node2_pt = Orig_node_pt[l2];
		    
		  // Index of the bulk Lagrange multiplier which enforces momentum PDEs	    
		  unsigned lambda_mom_non_aug2_index =
		    Non_augmented_bulk_element_pt->index_of_lagrange_multiplier(nonaug_node2_pt, i);

		  unsigned ext_index = External_data_index_for_non_aug_node[l2];
	  
		  int local_eqn_lagr_mom_non_aug =
		    this->external_local_eqn(ext_index, lambda_mom_non_aug2_index);

		  if(local_eqn_lagr_mom_non_aug >= 0)
		  {
		    jacobian(local_eqn_lambda_hat_hat, local_eqn_lagr_mom_non_aug) -=
		      psi[l] * psi[l2] * W;
		  }

		  // Index of the bulk Lagrange multiplier which enforces momentum PDEs	    
		  unsigned lambda_mom2_index =
		    bulk_el_pt->index_of_lagrange_multiplier(node2_pt, i);

		  int local_eqn_lagr_mom = this->nodal_local_eqn(l2, lambda_mom2_index);

		  if(local_eqn_lagr_mom >= 0)
		  {
		    jacobian(local_eqn_lambda_hat_hat, local_eqn_lagr_mom) +=
		      psi[l] * psi[l2] * W;
		  }
		}
	      } // end Jacobian check
	    }
	  }
	  
	  // ==================================================================
	  // Contribution of momentum- and continuity-enforcing Lagrange
	  // multipliers to augmented velocity equations
	  // ==================================================================
	  
	  int local_eqn_augmented = this->nodal_local_eqn(l, i);
	  if (local_eqn_augmented >= 0)
	  {
	    residuals[local_eqn_augmented] += lambda_cont[i] * psi[l]*W;

	    // compute Jacobian
	    if (flag == 1)
	    {
	      // loop over this face element's nodes again
	      for(unsigned l2=0; l2<n_node; l2++)
	      {
		Node* node2_pt = this->node_pt(l2);

		// get the map which gives the starting nodal index for
		// the Lagrange multipliers associated with each boundary ID
		std::map<unsigned, unsigned> first_index2 = *(
		  dynamic_cast<BoundaryNodeBase*>(node2_pt)->
		  index_of_first_value_assigned_by_face_element_pt() );

		// get the nodal index of the Lagrange multiplier for this
		// coordinate direction and boundary ID
		unsigned lambda_cont_index = first_index2.at(this->Boundary_id) + i;
		
		int local_unknown_lambda_cont = this->nodal_local_eqn(l2, lambda_cont_index);

		// if lambda-hat isn't pinned
		if (local_unknown_lambda_cont >= 0)
		{
		  // du/dlambda_p
		  jacobian(local_eqn_augmented, local_unknown_lambda_cont) +=
		    psi[l2] * psi[l] * W;
		}		
	      }

	    } // end Jacobian flag check
          } // end of augmented velocity equations

	  // ==================================================================
	  // Contribution of Lagrange multiplier to velocity eqn in
	  // non-augmented region
	  // ==================================================================

	  // get the external equation number for the non-augmented nodal dof
	  int local_eqn_non_aug =
	    this->external_local_eqn(External_data_index_for_non_aug_node[l], i);

// QUEHACERES delete? 
/* #ifdef PARANOID	   */
/* 	  if(l_in_non_aug_bulk < 0) */
/* 	  { */
/* 	    oomph_info << "couldn't find the node number of this node in the non-aug" */
/* 		       << " bulk element\n" << std::endl; */
/* 	    abort(); */
/* 	  } */

/* #endif */
	  if (local_eqn_non_aug >= 0)
	  { 
	    // QUEHACERES write-up Eq. 1.4 term 2
	    residuals[local_eqn_non_aug] -= lambda_cont[i] * psi[l] * W;

	    // compute Jacobian
	    if (flag == 1)
	    {
	      for(unsigned l2=0; l2<n_node; l2++)
	      {
		Node* node2_pt = this->node_pt(l2);

		// get the map which gives the starting nodal index for
		// the Lagrange multipliers associated with each boundary ID
		std::map<unsigned, unsigned> first_index2 = *(
		  dynamic_cast<BoundaryNodeBase*>(node2_pt)->
		  index_of_first_value_assigned_by_face_element_pt() );

		// get the nodal index of the Lagrange multiplier for this
		// coordinate direction and boundary ID
		unsigned lambda_cont_index = first_index2.at(this->Boundary_id) + i;
		  
		int local_unknown_lambda_cont = this->nodal_local_eqn(l2, lambda_cont_index);
		if (local_unknown_lambda_cont >= 0) 
		{
		  jacobian(local_eqn_non_aug, local_unknown_lambda_cont) -=
		    psi[l2] * psi[l]*W;
		}

	      }

	    } // end Jacobian check
	  } // end check for non-aug velocity pinned

	} // end loop over dimensions
      } // end loop over test functions      
    } // end loop over integration points
    
  } // end of fill_in_generic_residual_contribution_nst_sing_jump

  //====================================================================
  /// Class for TractionElements in a PDE-constrained optimisation
  // framework - really just a wrapper to override the residual contributions
  //====================================================================  
  template <class ELEMENT>
    class NavierStokesPdeConstrainedOptimisationTractionElement :
    public virtual NavierStokesTractionElement <ELEMENT>
  {
  public:

    ///Constructor, which takes a "bulk" element and the value of the index
    ///and its limit - forward to parent constructor
    NavierStokesPdeConstrainedOptimisationTractionElement(
      ELEMENT* const& element_pt, 
      const int& face_index,
      const bool& called_from_refineable_constructor=false) :
    NavierStokesTractionElement<ELEMENT>(element_pt, face_index,
					 called_from_refineable_constructor)
    { }

    NavierStokesPdeConstrainedOptimisationTractionElement(
      const NavierStokesPdeConstrainedOptimisationTractionElement&)
    {
      BrokenCopy::broken_copy("NavierStokesPdeConstrainedOptimisationTractionElement");
    }

    void operator=(const NavierStokesPdeConstrainedOptimisationTractionElement&)
    {
      BrokenCopy::broken_assign("NavierStokesPdeConstrainedOptimisationTractionElement");
    }
    
    /// This function returns just the residuals
    inline void fill_in_contribution_to_residuals(Vector<double>& residuals)
    {
      // Call the generic residuals function with flag set to 0
      // using a dummy matrix argument
      fill_in_generic_residual_contribution_fluid_traction_constrained_opt(
	residuals, GeneralisedElement::Dummy_matrix, 0);
    }

    /// This function returns the residuals and the jacobian
    inline void fill_in_contribution_to_jacobian(Vector<double>& residuals,
						 DenseMatrix<double>& jacobian)
    {
      //Call the generic routine with the flag set to 1
      fill_in_generic_residual_contribution_fluid_traction_constrained_opt(residuals, jacobian, 1);
    }

    void fill_in_generic_residual_contribution_fluid_traction_constrained_opt(
      Vector<double>& residuals,
      DenseMatrix<double>& jacobian,
      const bool& flag = false)
    {
      // Get continuous time from timestepper of first node
      const double time = this->node_pt(0)->time_stepper_pt()->time_pt()->time();
      
      //Find out how many nodes there are
      const unsigned n_node = this->nnode();
 
      //Set up memory for the shape and test functions
      Shape psif(n_node), testf(n_node);
 
      //Set the value of n_intpt
      const unsigned n_intpt = this->integral_pt()->nweight();
       
      // shorthand
      const unsigned Dim = this->Dim;

      // shorthand for bulk element this face element is attached to
      ELEMENT* bulk_el_pt = dynamic_cast<ELEMENT*>(this->bulk_element_pt());
      
      //Loop over the integration points
      for(unsigned ipt=0; ipt<n_intpt; ipt++)
      {
	//Get the integral weight
	double w = this->integral_pt()->weight(ipt);
   
	//Find the shape and test functions and return the Jacobian
	//of the mapping
	double J = this->shape_and_test_at_knot(ipt, psif, testf);
   
	//Premultiply the weights and the Jacobian
	double W = w*J;
   
	//Need to find position to feed into Traction function
	Vector<double> interpolated_x(Dim, 0.0);
   
	//Calculate velocities and derivatives
	for(unsigned l=0; l<n_node; l++) 
	{
	  //Loop over velocity components
	  for(unsigned i=0; i<this->Dim; i++)
	  {
	    interpolated_x[i] += this->nodal_position(l,i) * psif[l];
	  }
	}
   
	// Get the outer unit normal
	Vector<double> interpolated_n(Dim, 0.0);
	this->outer_unit_normal(ipt, interpolated_n);

	//Get the user-defined traction terms
	Vector<double> traction(Dim, 0.0);
	this->get_traction(time, interpolated_x, interpolated_n, traction);

	// ====================================================================
	//Now add to the appropriate equations
	// ====================================================================
	
	//Loop over the test functions
	for(unsigned l=0; l<n_node; l++)
	{
	  //Loop over the velocity components
	  for(unsigned i=0; i<Dim; i++)
	  {
	    // get the bulk node number corresponding to this face element node
	    const unsigned node_number_in_bulk = this->bulk_node_number(l);
	    
	    // Index of the bulk Lagrange multiplier which enforces momentum PDEs 
	    unsigned lambda_momentum_index =
	      bulk_el_pt->index_of_lagrange_multiplier(node_number_in_bulk, i);

	    int local_eqn_lagr_mom = this->nodal_local_eqn(l, lambda_momentum_index);
	  
	    // if this dof isn't pinned
	    if(local_eqn_lagr_mom >= 0)
	    {
	      // contribution of the prescribed traction to the bulk LM field
	      residuals[local_eqn_lagr_mom] += traction[i] * psif[l] * W;

	      // no contribution to Jacobian since the traction is prescribed
	      // and so doesn't depend on velocities, pressures or other LMs
	    }
	  }
	}
      }
    }
  };


  
  
  ////////////////////////////////////////////////////////////////////////////
  ////////////////////////////////////////////////////////////////////////////
  ////////////////////////////////////////////////////////////////////////////

  // hierher really need to tidy this up! Should only need one class 
  // for T and Q
  //
  //====================================================================
  /// New class. Mainly overloads output-related functions to add
  /// "singular function" (which is assumed to satisfy the Laplace
  /// equation; therefore no change to the governing (bulk) equations) 
  /// to the FE solution. 
  //====================================================================
  template <unsigned DIM>
    class TNavierStokesWithSingularityPdeConstrainedMinElement :
    public virtual TTaylorHoodElement<DIM>
  {
    
  public:
    
    typedef void (*ExactNonSingularFctPt)
      (const Vector<double>& x, Vector<double>& u, DenseMatrix<double>& grad_u);

    // function pointer for helper function which computes the
    // constitutive relationship
    typedef DenseMatrix<double> (*StressFctPt)(const DenseMatrix<double>& strain_rate,
					       const double& p);

    /// \short Function pointer to a functional which is to be minimised
    /// (subject to constraints imposed by the governing Stokes PDEs)
    typedef double (*FunctionalToMinimiseFctPt)(const Vector<double>& u);

    /// \short Function pointer to the derivative of the functional 
    /// w.r.t. the solution
    typedef Vector<double> (*DfunctionalDuFctPt)(const Vector<double>& u);

    /// Function pointer to a function which computes the body force
    typedef void (*BodyForceFctPt) (const Vector<double>& x, Vector<double>& body_force);
      
    // expose the template argument so this can be referenced externally
    static const unsigned _NNODE_1D_ = TTaylorHoodElement<DIM>::_NNODE_1D_;

    // expose the dimension template argument
    static const unsigned _DIM_ = DIM;
    
    // shorthand for the singular line element
    // QUEHACERES hacky! get rid of the #define at some point
    typedef ScalableSingularityForNavierStokesLineElement<SINGULAR_ELEMENT_NNODE_1D>
      SingularLineElement;

    unsigned nnode_on_face() const
    {
      switch (DIM)
      {
	case 2:
	  return _NNODE_1D_;
	  
	case 3:
	  return 3 * _NNODE_1D_ - 3;
      }
    }
    
    /// Constructor                            
  TNavierStokesWithSingularityPdeConstrainedMinElement() :
    Nplot(0), Is_lower_disk_element(false), Is_augmented_element(false),     
      Exact_non_singular_fct_pt(nullptr), Dfunctional_du_fct_pt(nullptr),
      Stress_fct_pt(nullptr), Body_force_fct_pt(nullptr), Nsingular_fct(0),
      Pressure_regularisation_factor(0.0), Velocity_regularisation_factor(0.0)
      , random_var_to_alter_mem_map(0.0)
    { }

    TNavierStokesWithSingularityPdeConstrainedMinElement(
      const TNavierStokesWithSingularityPdeConstrainedMinElement&)
    {
      BrokenCopy::broken_copy("TNavierStokesWithSingularityPdeConstrainedMinElement");
    }

    void operator=(const TNavierStokesWithSingularityPdeConstrainedMinElement&)
    {
      BrokenCopy::broken_assign("TNavierStokesWithSingularityPdeConstrainedMinElement");
    }

    // add extra values at each node and record the index of the first one.
    // These are the Lagrange multipliers which weakly enforce the Stokes momentum eqs
    void add_lagrange_multiplier_dofs(std::map<Node*,unsigned>& node_to_first_lm_index_map)
    {      
      // make enough space in the vector of maps for each node
      Lambda_index.resize(this->nnode());

      unsigned nnod = this->nnode();
      for(unsigned n=0; n<nnod; n++)
      {
	// number of Lagrange multipliers in the augmented nodes -
	// one for each momentum equation
	unsigned nlagrange_multiplier = DIM;

	// if we're on a vertex node (these are enumerated first, so it's safe
	// to check if we're at an index less than DIM+1) then we also want
	// to add space for a continuity-enforcing Lagrange multiplier
	if(n < DIM + 1)
	  nlagrange_multiplier++;
	
	Node* node_pt = this->node_pt(n);
	unsigned nvalue = node_pt->nvalue();

	// has this node already had the dofs added?
	if(node_to_first_lm_index_map.find(node_pt) == node_to_first_lm_index_map.end())
	{	
	  // resize the node to accomodate new LM dofs
	  node_pt->resize(nvalue + nlagrange_multiplier);

	  // set the index of the first LM
	  Lambda_index[n] = nvalue;

	  // and add this index to the complete node to index map
	  node_to_first_lm_index_map[node_pt] = nvalue;
	}
	else
	{
	  // we've already done it for this node (from a neighbouring element)
	  // so just need to set the index
	  Lambda_index[n] = node_to_first_lm_index_map[node_pt];
	}
      }
    }
    
    // get the nodal index of the first LM at node n
    unsigned index_of_lagrange_multiplier(const unsigned& n, const unsigned& i) const
    {      
      return Lambda_index[n] + i;
    }

    // index of continuity-enforcing LM
    unsigned index_of_lagrange_multiplier_p(const unsigned& n) const
    {
#ifdef PARANOID      
      if(n > DIM)
      {
	// lambda_p is only stored at the vertices
	throw OomphLibError("Error: No continuity Lagrange multiplier for this node\n",
			    OOMPH_CURRENT_FUNCTION,
			    OOMPH_EXCEPTION_LOCATION);
      }
#endif
      return Lambda_index[n] + DIM;
    }

    // get the nodal index of the first LM at a given node pointer
    // (useful for face elements attached to this bulk element where the nodal
    // enumeration is different in general)
    unsigned index_of_lagrange_multiplier(Node* const& node_pt, const unsigned& i) const
    {
      for(unsigned n=0; n<this->nnode(); n++)
      {	
	if(this->node_pt(n) == node_pt)	  
	  return Lambda_index[n] + i;
      }

      ostringstream error_message;
      error_message << "Error, couldn't get bulk Lagrange multiplier index "
		    << "for the requested node " << node_pt << " (";
      for(unsigned i=0; i<DIM; i++)
	error_message << node_pt->x(i) << " ";
      
      error_message << ") as this "
		    << "node isn't in this element!\n\n"
		    << "Nodes in this element: \n";
      
      for(unsigned n=0; n<this->nnode(); n++)
      {
	error_message << this->node_pt(n) << " ("
		      << this->node_pt(n)->x(0) << " "
		      << this->node_pt(n)->x(1) << " "
		      << this->node_pt(n)->x(2) << ")"
		      << std::endl;
      }
      
      // if we didn't return in the above loop then something's gone wrong
      throw OomphLibError(error_message.str().c_str(),
			  OOMPH_CURRENT_FUNCTION,
			  OOMPH_EXCEPTION_LOCATION);
    }

    // get the edge coordinates and (a pointer to) the element which stores
    // the singular amplitudes
    void lagr_coords_and_singular_element_at_knot(
      const unsigned& ipt,
      LagrangianCoordinates& lagr_coords,
      SingularLineElement*& sing_el_pt,
      Vector<double>& s_singular_el) const
    {
      if(!Is_augmented_element)
      {
	throw OomphLibError("Error: don't call lagr_coords_and_singular_element_at_knot() "
			    "from a non-augmented element, this makes no sense!",
			    OOMPH_CURRENT_FUNCTION,
			    OOMPH_EXCEPTION_LOCATION);
      }
      
      // get the line element and local coordinate which corresponds to the
      // singular amplitude for this knot
      std::pair<GeomObject*, Vector<double> > line_elem_and_local_coord = 
	this->line_element_and_local_coordinate_at_knot(ipt);

      // cast the GeomObject to a singular line element      
      sing_el_pt = dynamic_cast<SingularLineElement*>
	(line_elem_and_local_coord.first);

      // check we've actually got one if we're an augmented element, otherwise
      // we're in trouble
      if(sing_el_pt == nullptr)
      {
	throw OomphLibError("Error: this is an augmented element "
			    "but lagr_coords_and_singular_element_at_knot() "
			    "has been called before the singular line element pointer "
			    "has been set\n",
			    OOMPH_CURRENT_FUNCTION,
			    OOMPH_EXCEPTION_LOCATION);
      }

      // local coordinate in the singular element for the zeta of this knot
      s_singular_el = line_elem_and_local_coord.second;

      // get the \rho,\zeta,\phi coordinates at this knot
      lagr_coords = this->lagr_coordinate_at_knot(ipt);
    }
    
    // ========================================================================
    // Assign the singular line elements associated with this face elements
    // integration points as external data for this element, since it will
    // make a contribution to the singular amplitudes
    void set_singular_amplitudes_as_external_data()
    { 
      // number of integration points in this element
      const unsigned n_intpt = this->integral_pt()->nweight();

      // Make space in our list of external equation indices
      C_external_data_index_at_knot.resize(n_intpt);
      
      // set of unique singular line element nodes to who's dofs (singular amplitudes)
      // this face element is responsible for making a contribution to.
      std::set<Node*> external_singular_node_set;
      
      // need to loop over the integration points, because different integration points
      // in this same element might correspond to different singular line elements; 
      // find the singular line element associated with each, and add it to our set
      for(unsigned ipt=0; ipt<n_intpt; ipt++)
      {  
	SingularLineElement* sing_el_pt = nullptr;
	LagrangianCoordinates lagr_coords_dummy;
	Vector<double> s_singular_el_dummy;
	lagr_coords_and_singular_element_at_knot(ipt, lagr_coords_dummy,
						 sing_el_pt, s_singular_el_dummy);

	// for each knot in this face element, we want to store the
	// index to each external node
	C_external_data_index_at_knot[ipt].resize(sing_el_pt->nnode());
	
	// now loop over this singular elements nodes
	for(unsigned j=0; j<sing_el_pt->nnode(); j++)
	{
	  C_external_data_index_at_knot[ipt][j] =
	    this->add_external_data(sing_el_pt->node_pt(j));
	}	
      }
    }

    // compute the body force at this Eulerian position
    void get_body_force(const Vector<double>& x, Vector<double>& body_force) const
    {
      if (Body_force_fct_pt == nullptr)
      {
	body_force.resize(3);
	body_force.initialise(0.0);
      }
      else
      {
	Body_force_fct_pt(x, body_force);
      }
    }

    BodyForceFctPt& body_force_fct_pt()
    {
      return Body_force_fct_pt;
    }
    
    void compute_error(std::ofstream& outfile,
		       FiniteElement::SteadyExactSolutionFctPt exact_soln_fct_pt,
		       double& v_error, double& p_error, double& norm) const
    {      
      v_error = 0.0;
      p_error = 0.0;
      norm    = 0.0;

      //Vector of local coordinates
      Vector<double> s(DIM);

      // Vector for coordintes
      Vector<double> x(DIM);

      //Set the value of n_intpt
      unsigned n_intpt = this->integral_pt()->nweight();
   
      outfile << "ZONE" << std::endl;
 
      // Exact solution Vector (u,v,[w],p)
      Vector<double> exact_soln(DIM+1);
   
      //Loop over the integration points
      for(unsigned ipt=0; ipt<n_intpt; ipt++)
      {
      	//Assign values of s
      	for(unsigned i=0; i<DIM; i++)
      	{
      	  s[i] = this->integral_pt()->knot(ipt,i);
      	}

      	//Get the integral weight
      	double w = this->integral_pt()->weight(ipt);

      	// Get jacobian of mapping
      	double J = this->J_eulerian(s);

      	//Premultiply the weights and the Jacobian
      	double W = w*J;

      	// Get x position as Vector
      	this->interpolated_x(s,x);

	// get the side of the disk right
	const double tol = 1e-8;
	if((abs(x[2]) < tol) && Is_lower_disk_element)
	  x[2] = -tol;
	
      	// Get exact solution at this point
      	(*exact_soln_fct_pt)(x, exact_soln);

      	// get the total solution u = u_{fe} + c \hat u
      	Vector<double> u_fe_plus_sing = interpolated_u_total_navier_stokes(s, ipt);
	
      	// Velocity error
      	for(unsigned i=0; i<DIM; i++)
      	{
      	  norm    += pow(exact_soln[i], 2) * W;
      	  v_error += pow(exact_soln[i] - u_fe_plus_sing[i], 2) * W;
      	}

	p_error = pow(exact_soln[DIM] - u_fe_plus_sing[DIM], 2) * W;
	
      	//Output x,y,...,u_exact
      	for(unsigned i=0; i<DIM; i++)
      	{
      	  outfile << x[i] << " ";
      	}

      	//Output x,y,[z],u_error,v_error,[w_error]
      	for(unsigned i=0; i<DIM+1; i++)
      	{
      	  outfile << exact_soln[i] - u_fe_plus_sing[i] << " ";
      	}
      	outfile << std::endl;
      }
    }

    // get the function pointer to the function which computes the
    // exact non-singular solution as a function of the Eulerian position
    ExactNonSingularFctPt& exact_non_singular_fct_pt()
    {
      return Exact_non_singular_fct_pt;
    }  

    // get the function pointer to the function which computes the
    // stress from the velocity gradient and the pressure
    StressFctPt& stress_fct_pt()
    {
      return Stress_fct_pt;
    }

    // get the function pointer to the function which computes the
    // derivative of the functional w.r.t. the solution
    DfunctionalDuFctPt& dfunctional_du_fct_pt()
    {
      return Dfunctional_du_fct_pt;
    }
    
    // tell this element that it's on the lower side of the disk
    void set_lower_disk_element()
    {
      Is_lower_disk_element = true;
    }

    void total_singular_contribution_and_gradient_at_knot(const unsigned& ipt,
							  Vector<double>& u_sing_total,
							  DenseMatrix<double>& dudx_sing_total) const
    {
      // clear eveything out
      u_sing_total.resize(DIM+1, 0.0);
      dudx_sing_total.resize(DIM, DIM, 0.0);

      u_sing_total.initialise(0.0);
      dudx_sing_total.initialise(0.0);
      
      // get the line element and local coordinate which corresponds to the
      // singular amplitude for this point
      std::pair<GeomObject*, Vector<double> > line_elem_and_local_coord;

      line_elem_and_local_coord = 
	this->line_element_and_local_coordinate_at_knot(ipt);

      // cast the GeomObject to a singular line element      
      SingularLineElement* sing_el_pt =
	dynamic_cast <SingularLineElement*>( line_elem_and_local_coord.first );

      // check if we're subtracting the singularity or not
      if (sing_el_pt == nullptr)
      {
#ifdef PARANOID
	if(Is_augmented_element)
	{	  
	  throw OomphLibError("Error, this is apparently an augmented element but has "
			      "no singular line element pointer, so can't get total singular "
			      "solution and gradient",
			      OOMPH_CURRENT_FUNCTION,
			      OOMPH_EXCEPTION_LOCATION);
	}
#endif
	// if it's a non-augmented element then nothing to do, just return zeros
	return;
      }
      
      // local coordinate in the singular element for the zeta of this plot point
      Vector<double> s_singular_el = line_elem_and_local_coord.second;
      
      // get the \rho,\zeta,\phi coordinates at this knot
      LagrangianCoordinates lagr_coords_at_point;
      
      lagr_coords_at_point = this->lagr_coordinate_at_knot(ipt);

      // get 'em
      u_sing_total =
	sing_el_pt->total_singular_contribution(lagr_coords_at_point,
						s_singular_el);

      dudx_sing_total =
	sing_el_pt->total_singular_gradient_contribution(lagr_coords_at_point,
							 s_singular_el);
    }
    
    // QUEHACERES make these time dependent
    /// \short Return FE representation of function value u_fe
    inline Vector<double> interpolated_u_fe_navier_stokes(const Vector<double>& s) const
    {
      // FE solution vector
      Vector<double> u_fe(DIM);

      // get FE velocity
      this->interpolated_u_nst(s, u_fe);

      // get FE pressure
      double p_fe = this->interpolated_p_nst(s);
     
      // add pressure to the solution vector
      u_fe.push_back(p_fe);
      
      return u_fe;
    } 

    // ------------------------------------------------------------------------
    
    /// \short Return FE representation of function value u_navier_stokes(s) 
    /// plus scaled singular fct (if provided) at local coordinate s
    inline Vector<double> interpolated_u_total_navier_stokes(
      const Vector<double>& s, const unsigned& ipt,
      const bool& use_plot_points = false) const
    {
      // get the interpolated FE bit
      Vector<double> u_fe = interpolated_u_fe_navier_stokes(s);

      // get the line element and local coordinate which corresponds to the
      // singular amplitude for this point
      std::pair<GeomObject*, Vector<double> > line_elem_and_local_coord;
	
      if(use_plot_points)
      {
	line_elem_and_local_coord = 
	  this->line_element_and_local_coordinate_at_plot_point(ipt);
      }
      else
      {
	line_elem_and_local_coord = 
	  this->line_element_and_local_coordinate_at_knot(ipt);
      }
      
      // cast the GeomObject to a singular line element      
      SingularLineElement* sing_el_pt =
	dynamic_cast <SingularLineElement*>( line_elem_and_local_coord.first );

      // check if we're subtracting the singularity or not
      if (sing_el_pt != nullptr)
      {
	// local coordinate in the singular element for the zeta of this plot point
	Vector<double> s_singular_el = line_elem_and_local_coord.second;
      
	// get the Lagrangian (xi1,xi2,xi3) coordinates at this knot
	LagrangianCoordinates lagr_coords_at_point;

	if(use_plot_points)
	  lagr_coords_at_point = this->lagr_coordinate_at_plot_point(ipt);
	else
	  lagr_coords_at_point = this->lagr_coordinate_at_knot(ipt);

	Vector<double> u_sing =
	  sing_el_pt->total_singular_contribution(lagr_coords_at_point,
						  s_singular_el);

	// add singular part of the solution to the FE part to give the total
	// computed solution
	for(unsigned i=0; i<DIM+1; i++)
	{
	  u_fe[i] += u_sing[i];
	}
      }
      return u_fe;
    } 

    void output_error_at_plot_points(std::ostream& outfile, 
				     const Vector<double>& s,
				     const unsigned& iplot,
				     FiniteElement::SteadyExactSolutionFctPt exact_soln_pt) const
    {
      // by definition of this function!
      const bool use_plot_points = true;
      
      // get the total computed solution at these local coords
      Vector<double> u_total = 
	this->interpolated_u_total_navier_stokes(s, iplot, use_plot_points);

      // exact solution at these coords
      Vector<double> exact_soln(4, 0.0);

      // Vector for coordinates
      Vector<double> x(DIM);
      
      // Get x position as Vector
      this->interpolated_x(s, x);
	
      // if we have an exact solution pointer, then get the exact solution
      if (exact_soln_pt != nullptr)
	(*exact_soln_pt)(x, exact_soln);

      // output the difference between the exact and the computed solutions
      // ------------------------------------------------------------------
      
      // coordinates
      for(unsigned i=0; i<DIM; i++) 
      {
	outfile << x[i] << " ";
      }
      
      // output the total solution
      for(unsigned i=0; i<DIM+1; i++)
      {
	outfile << exact_soln[i] - u_total[i] << " ";
      }
      outfile << std::endl;
    }

    /// Output with various contributions
    void output_error_at_plot_points(std::ostream& outfile, 
				     const unsigned& nplot,
				     FiniteElement::SteadyExactSolutionFctPt exact_soln_pt) const
    {
      //Vector of local coordinates
      Vector<double> s(DIM);
   
      // Tecplot header info
      outfile << this->tecplot_zone_string(nplot);
   
      // Loop over plot points
      unsigned num_plot_points = this->nplot_points(nplot);
      for (unsigned iplot=0; iplot < num_plot_points; iplot++)
      {
	// Get local coordinates of plot point
	this->get_s_plot(iplot, nplot, s);

	// do the output
	output_error_at_plot_points(outfile, s, iplot, exact_soln_pt);
      }
      outfile << std::endl;
      
      // Write tecplot footer (e.g. FE connectivity lists)
      this->write_tecplot_zone_footer(outfile, nplot);   
    }

    void output_divergence(std::ostream& outfile,
			   const unsigned& nplot) const
    {
      //Vector of local coordinates
      Vector<double> s(DIM);
   
      // Tecplot header info
      outfile << this->tecplot_zone_string(nplot);
   
      // Loop over plot points
      unsigned num_plot_points = this->nplot_points(nplot);
      for (unsigned iplot=0; iplot < num_plot_points; iplot++)
      {
	// Get local coordinates of plot point
	this->get_s_plot(iplot, nplot, s);

	// get the interpolated position
	Vector<double> x(DIM);
	for(unsigned i=0; i<DIM; i++) 
	{
	  outfile << this->interpolated_x(s,i) << " ";	
	}
	
	// get the line element and local coordinate which corresponds to the
	// singular amplitude for this knot
	std::pair<GeomObject*, Vector<double> > line_elem_and_local_coord = 
	  this->line_element_and_local_coordinate_at_plot_point(iplot);
      
	// cast the GeomObject to a singular line element      
	SingularLineElement* sing_el_pt =
	  dynamic_cast <SingularLineElement*>( line_elem_and_local_coord.first );

	// local coordinate in the singular element for the zeta of this plot point
	Vector<double> s_singular_el = line_elem_and_local_coord.second;

	// do the divergence, div u = du_i/dx_i
	double div_u_fe    = 0;
	double div_u_sing  = 0.0;
	double div_u_total = 0.0;
			
	for(unsigned i=0; i<DIM; i++)
	{
	  // get derivative du_i/dx_i and add to divergence
	  double div_u_fe_i = this->interpolated_dudx_nst(s, i, i);
	  
	  div_u_fe    += div_u_fe_i;
	  div_u_total += div_u_fe_i;
	  
	  if(sing_el_pt != nullptr)
	  {
	    LagrangianCoordinates lagr_coords = lagr_coordinate_at_plot_point(iplot);
	    
	    DenseMatrix<double> dudx_sing =
	      sing_el_pt->total_singular_gradient_contribution(lagr_coords,
							       s_singular_el);

	    div_u_sing += dudx_sing(i,i);

	    div_u_total += dudx_sing(i,i);
	  }
	}

	outfile << div_u_total << " "
		<< div_u_fe << " "
		<< div_u_sing << std::endl;
      }
      
      outfile << std::endl;
      
      // Write tecplot footer (e.g. FE connectivity lists)
      this->write_tecplot_zone_footer(outfile, nplot);   
    }
    
    void output_with_various_contributions(std::ostream& outfile, 
					   const Vector<double>& s,
					   const unsigned& iplot,
					   FiniteElement::SteadyExactSolutionFctPt exact_soln_pt) const
    {
      // since this is just the solution output (not integration)
      // we want to use plot points not knot points
      const bool use_plot_points = true;
      
      Vector<double> x(DIM);
      for(unsigned i=0; i<DIM; i++) 
      {
	x[i] = this->interpolated_x(s,i);	
      }

      // regular part of the solution
      Vector<double> u_exact_non_sing(DIM+1, 0.0);

      // get the regular FE solution, and the full computed solution u = u_FE + u_sing
      Vector<double> u_fe(DIM+1, 0.0);
      Vector<double> u_fe_plus_sing(DIM+1, 0.0);

      u_fe           = this->interpolated_u_fe_navier_stokes(s);
      u_fe_plus_sing = this->interpolated_u_total_navier_stokes(s, iplot, use_plot_points);


      // ==========================================
      // output total solution and FE bits
      // ==========================================
      // coordinates
      for(unsigned i=0; i<DIM; i++) 
      { 
	outfile << x[i] << " ";
      }
      
      // output the total solution
      for(unsigned i=0; i<DIM+1; i++)
      {
	outfile << u_fe_plus_sing[i] << " ";
      }	
      // output the FE bit
      for(unsigned i=0; i<DIM+1; i++)
      {
	outfile << u_fe[i] << " ";
      }

      // ==========================================
      // output the singular bits
      // ==========================================

      // get the line element and local coordinate which corresponds to the
      // singular amplitude for this knot
      std::pair<GeomObject*, Vector<double> > line_elem_and_local_coord = 
	this->line_element_and_local_coordinate_at_plot_point(iplot);
      
      // cast the GeomObject to a singular line element      
      SingularLineElement* sing_el_pt =
	dynamic_cast <SingularLineElement*>( line_elem_and_local_coord.first );

      // local coordinate in the singular element for the zeta of this plot point
      Vector<double> s_singular_el = line_elem_and_local_coord.second;
      
      // do we actually have a pointer to a singular element?
      if(sing_el_pt == nullptr)
      {
	// output the regular part of the exact solution, which is the same as
	// the total exact solution if there's no singular part (i.e. we're in the bulk)
	if(exact_soln_pt != nullptr)
	{
	  Vector<double> exact_soln(4, 0.0);
	  
	  // Get exact solution at this point
	  (*exact_soln_pt)(x, exact_soln);

	  for(unsigned i=0; i<DIM+1; i++)
	  {
	    outfile << exact_soln[i] << " ";
	  }	  
	}
	else
	{
	  // if no exact solution function pointer has been provided,
	  // just output zeros to keep oomph-convert happy
	  for(unsigned i=0; i<DIM+1; i++)
	    outfile << "0 ";
	}

	// -----------------
	
	// now for the singular contributions, just output zeros if there are
	// no singular functions so the oomph-convert script doesn't die
	// when the number of columns isn't the same for all elements
	for(unsigned ising=0; ising < Nsingular_fct; ising++)
	{	
	  for(unsigned i=0; i<DIM+1; i++)
	    outfile << "0 ";
	}		
      }
      else
      {
	// get the Lagrangian coordinates at this plot point
	LagrangianCoordinates lagr_coords_at_plot =
	  this->lagr_coordinate_at_plot_point(iplot);

	// get the total singular contribution at this point
	Vector<double> u_sing_total = 
	  sing_el_pt->total_singular_contribution(lagr_coords_at_plot,
						  s_singular_el);

	// Get exact solution at this point
	Vector<double> exact_soln(4, 0.0);
	if(exact_soln_pt != nullptr)
	{
	  (*exact_soln_pt)(x, exact_soln);
	}

	// now subtract off the singular contributions to get the regular part of
	// the exact solution
	for(unsigned i=0; i<DIM+1; i++)
	{
	  outfile << exact_soln[i] - u_sing_total[i] << " ";
	}
      
	// get a list of the singular function IDs this element has
	Vector<unsigned> sing_ids = sing_el_pt->singular_fct_ids();

	// now loop over the singular functions and output each contribution
	for(unsigned sing_fct_id : sing_ids) 
	{	
	  // singular part of the solution
	  Vector<double> u_sing = sing_el_pt->singular_fct(lagr_coords_at_plot,
							   s_singular_el,
							   sing_fct_id);
	  for(unsigned i=0; i<DIM+1; i++)
	  {
	    outfile << u_sing[i] << " ";
	  }
	}
      }
      
      {
	Vector<double> lambda = interpolated_lambda(s);
	double lambda_p = interpolated_lambda_p(s);

	for(Vector<double>::iterator lambda_it = lambda.begin();
	    lambda_it != lambda.end(); lambda_it++)
	  outfile << *lambda_it << " ";

	outfile << lambda_p << " ";
      }
      
      outfile << std::endl;
    }
    
    /// Output with various contributions
    void output_with_various_contributions(
      std::ostream& outfile, 
      const unsigned& nplot,
      FiniteElement::SteadyExactSolutionFctPt exact_soln_pt = nullptr) const
    {
      //Vector of local coordinates
      Vector<double> s(DIM);
   
      // Tecplot header info
      outfile << this->tecplot_zone_string(nplot);
   
      // Loop over plot points
      unsigned num_plot_points = this->nplot_points(nplot);
      for (unsigned iplot=0; iplot < num_plot_points; iplot++)
      {
	// Get local coordinates of plot point
	this->get_s_plot(iplot, nplot, s);

	// do the output
	output_with_various_contributions(outfile, s, iplot, exact_soln_pt);
      }
      outfile << std::endl;
      
      // Write tecplot footer (e.g. FE connectivity lists)
      this->write_tecplot_zone_footer(outfile, nplot);   
    }
    
    // ========================================================================
    // ========================================================================
    
    // set the edge coordinates of each plot point
    void set_lagr_coordinates_at_plot_point(const Vector<LagrangianCoordinates>& coords)
    {
      unsigned ncoordinates = coords.size();
      
      // check if this has already been set
      if(Nplot == 0)
      {
	// set the number of plot points, used during output
	Nplot = ncoordinates;
      }
      else if(ncoordinates != Nplot)
      {
	ostringstream error_message;
	error_message << "number of sets of coordinates provided is not consistent "
		      << "with previously set Nplot\n"
		      << "Number of supplied coordinate sets: " << ncoordinates << "\n"
		      << "Nplot:                              " << Nplot;
	throw OomphLibError(
	  error_message.str(),
	  OOMPH_CURRENT_FUNCTION,
	  OOMPH_EXCEPTION_LOCATION);
      }

      // make sure we've got enough space
      Lagrangian_coordinates_at_plot_point.resize(Nplot);

      // set 'em
      for(unsigned i=0; i<Nplot; i++)
      {
	Lagrangian_coordinates_at_plot_point[i]  = coords[i];
      }
    }

    // set the edge coordinates of each knot point
    void set_lagr_coordinates_at_knot(const Vector<LagrangianCoordinates>& coords)
    {
      // Number of integration points
      const unsigned n_intpt = this->integral_pt()->nweight();
      
      if(coords.size() != n_intpt)
      {
	ostringstream error_message;
	error_message << "number of sets of coordinates provided is not consistent "
		      << "with number of integration points in this element\n"
		      << "Number of supplied coordinate sets: " << coords.size() << "\n"
		      << "Number of integration points:       " << n_intpt;
	throw OomphLibError(
	  error_message.str(),
	  OOMPH_CURRENT_FUNCTION,
	  OOMPH_EXCEPTION_LOCATION);
      }

      // make sure we've got enough space
      Lagrangian_coordinates_at_knot_point.resize(n_intpt);

      // set 'em
      for(unsigned i=0; i<n_intpt; i++)
      {
	Lagrangian_coordinates_at_knot_point[i]  = coords[i];
      }
    }

    // set the edge coordinates at each knot
    void set_nodal_lagr_coordinates(const Vector<LagrangianCoordinates>& lagr_coords)
    {
      // number of coordinates supplied
      unsigned ncoords = lagr_coords.size();

      const unsigned n_node = this->nnode();
      
      // check the right number of coordinates have been passed in
      if(ncoords != n_node)
      {
	ostringstream error_message;
	error_message << "Number of sets of coordinates provided is not consistent "
		      << "with the number of nodes.\n"
		      << "Number of supplied coordinate sets: " << ncoords << "\n"
		      << "Number of nodes:                    " << n_node;
	
	throw OomphLibError(
	  error_message.str(),
	  OOMPH_CURRENT_FUNCTION,
	  OOMPH_EXCEPTION_LOCATION);
      }

      // make enough space
      Nodal_lagr_coordinates.resize(n_node);
      
      // set 'em
      for(unsigned i=0; i<n_node; i++)
	Nodal_lagr_coordinates[i] = lagr_coords[i];
    }

    // ========================================================================
    // ========================================================================
    
    // get the edge coordinates at the ith plot point
    LagrangianCoordinates lagr_coordinate_at_plot_point(const unsigned& i) const
    {
      return Lagrangian_coordinates_at_plot_point[i];
    }

    // get the edge coordinates at the ith knot point
    LagrangianCoordinates lagr_coordinate_at_knot(const unsigned& i) const
    {
      return Lagrangian_coordinates_at_knot_point[i];
    }

    // get the edge coordinates at the ith node
    LagrangianCoordinates nodal_lagr_coordinates(const unsigned& i) const
    {
      return Nodal_lagr_coordinates[i];
    }

    LagrangianCoordinates interpolated_lagr_coordinates(const Vector<double>& s) const
    {
      LagrangianCoordinates interpolated_lagr_coords(0.0, 0.0, 0.0);

      const unsigned n_node = this->nnode();
      
      //Set up memory for the shape functions
      Shape psi(n_node);

      // get 'em
      this->shape(s, psi);
      
      // The number of integration points (knots)
      const unsigned n_intpt = this->integral_pt()->nweight();

      // loop over each node in this element and add its contribution to the
      // interpolated coordinates
      for(unsigned j=0; j<n_node; j++)
      {
	// get the edge coordinates of this node
	LagrangianCoordinates nodal_lagr_coords = Nodal_lagr_coordinates[j];

	// interpolated with the shape functions
	interpolated_lagr_coords.xi1 += nodal_lagr_coords.xi1 * psi[j];
	interpolated_lagr_coords.xi2 += nodal_lagr_coords.xi2 * psi[j];
	interpolated_lagr_coords.xi3 += nodal_lagr_coords.xi3 * psi[j];
      }     
      
      return interpolated_lagr_coords;
    }

    // ========================================================================
    // ========================================================================

    void lagr_coords_and_singular_element_at_node(const unsigned& inode,
						  LagrangianCoordinates& lagr_coords,
						  SingularLineElement*& sing_el_pt,
						  Vector<double>& s_singular_el) const
    {
      // get the line element and local coordinate which corresponds to the
      // singular amplitude for this knot
      std::pair<GeomObject*, Vector<double> > line_elem_and_local_coord = 
	this->line_element_and_local_coordinate_at_node(inode);

      // cast the GeomObject to a singular line element      
      sing_el_pt = dynamic_cast<SingularLineElement*>(line_elem_and_local_coord.first);

      // check we've actually got one, or we're in trouble
      if(sing_el_pt == nullptr)
      {
	ostringstream error_message;

	error_message << "Error: this singular face element has no "
		      << "singular line element pointer\n";
	    
	throw OomphLibError(error_message.str().c_str(),
			    OOMPH_CURRENT_FUNCTION,
			    OOMPH_EXCEPTION_LOCATION);
      }

      // local coordinate in the singular element for the zeta of this knot
      s_singular_el = line_elem_and_local_coord.second;

      // get the \rho,\zeta,\phi coordinates at this knot
      lagr_coords = this->nodal_lagr_coordinates(inode);
    }

    // ========================================================================
    // ========================================================================
    
    void set_nsingular_fct(const unsigned& nsing)
    {
      Nsingular_fct = nsing;
    }
    
    // set the line element and local coordinate at each plot point
    void set_line_element_and_local_coordinate_at_plot_point(
      const Vector<std::pair<GeomObject*, Vector<double> > >&
      line_element_and_local_coordinate_at_plot_point)
    {
      unsigned ncoordinates = line_element_and_local_coordinate_at_plot_point.size();
      
      // check if this has already been set
      if(Nplot == 0)
      {
	// set the number of plot points, used during output
	Nplot = ncoordinates;
      }
      else if(ncoordinates != Nplot)
      {
	ostringstream error_message;
	error_message << "number of sets of coordinates provided is not consistent "
		      << "with previously set Nplot\n"
		      << "Number of supplied coordinate sets: " << ncoordinates << "\n"
		      << "Nplot:                              " << Nplot;
	throw OomphLibError(
	  error_message.str(),
	  OOMPH_CURRENT_FUNCTION,
	  OOMPH_EXCEPTION_LOCATION);
      }
      
      // make sure we've got enough space
      Line_element_and_local_coordinate_at_plot_point.resize(Nplot);

      // set 'em
      for(unsigned i=0; i<Nplot; i++)
      {
	Line_element_and_local_coordinate_at_plot_point[i] =
	  line_element_and_local_coordinate_at_plot_point[i];
      }
    }

    // set the line element and local coordinate at each knot point
    void set_line_element_and_local_coordinate_at_knot(
      const Vector<std::pair<GeomObject*, Vector<double> > >&
      line_element_and_local_coordinate_at_knot_point)
    {
      // Number of integration points
      const unsigned n_intpt = this->integral_pt()->nweight();
      
      unsigned ncoordinates = line_element_and_local_coordinate_at_knot_point.size();
      
      // check if this has already been set
      if(ncoordinates != n_intpt)
      {
	ostringstream error_message;
	error_message << "number of coordinates provided is not consistent "
		      << "with the number of integration points in this element\n"
		      << "Number of supplied coordinate sets: " << ncoordinates << "\n"
		      << "Number of integration points:       " << n_intpt;
	throw OomphLibError(
	  error_message.str(),
	  OOMPH_CURRENT_FUNCTION,
	  OOMPH_EXCEPTION_LOCATION);
      }
      
      // make sure we've got enough space
      Line_element_and_local_coordinate_at_knot.resize(n_intpt);

      // set 'em
      for(unsigned i=0; i<n_intpt; i++)
      {
	Line_element_and_local_coordinate_at_knot[i] =
	  line_element_and_local_coordinate_at_knot_point[i];
      }

      // now set the nodal dofs of the singular line elements as external data
      // for this bulk element
      if(Is_augmented_element)
      {
	set_singular_amplitudes_as_external_data();
      }
    }

    // set the line element and local coordinate for the ith knot
    void set_line_element_and_local_coordinate_at_node(
      const Vector<std::pair<GeomObject*, Vector<double> > >&
      line_element_and_local_coordinate_at_node)
    {
      const unsigned n_node = this->nnode();
      
      // number of coordinates supplied
      unsigned ncoords = line_element_and_local_coordinate_at_node.size();

      // check the right number of coordinates have been passed in
      if(ncoords != n_node)
      {
	ostringstream error_message;
	error_message << "Number of element-coordinate pairs provided is not consistent "
		      << "with the number of nodes.\n"
		      << "Number of supplied coordinate sets: " << ncoords << "\n"
		      << "Number of nodes:                    " << n_node;
	
	throw OomphLibError(
	  error_message.str(),
	  OOMPH_CURRENT_FUNCTION,
	  OOMPH_EXCEPTION_LOCATION);
      }

      // make enough space
      Line_element_and_local_coordinate_at_node.resize(n_node);
      
      // set 'em
      for(unsigned i=0; i<n_node; i++)
      {
	Line_element_and_local_coordinate_at_node[i] =
	  line_element_and_local_coordinate_at_node[i];
      }
    }

    // ========================================================================
    // ========================================================================
    
    // get the line element and local coordinate for the ith plot point
    std::pair<GeomObject*, Vector<double> >
      line_element_and_local_coordinate_at_plot_point(const unsigned& i) const
    {
      std::pair<GeomObject*, Vector<double> > return_pair;

      // check if we have an element-coordinate pair for this plot point
      // (since this bulk element is used everywhere in the bulk not just
      // in the augmented region, there may be no singular functions
      // associated with this element, and so we just want to return an
      // empty pointer and vector)
      if(i < Line_element_and_local_coordinate_at_plot_point.size())
      {
	return_pair = Line_element_and_local_coordinate_at_plot_point[i];
      }
      else
      {
	GeomObject* dummy_pt = nullptr;
	Vector<double> dummy_vec(1,0);
	return_pair = std::make_pair(dummy_pt, dummy_vec);
      }
      return return_pair;
    }

    // get the line element and local coordinate for the ith knot point
    std::pair<GeomObject*, Vector<double> >
      line_element_and_local_coordinate_at_knot(const unsigned& i) const
    {
      std::pair<GeomObject*, Vector<double> > return_pair;

      // check if we have an element-coordinate pair for this plot point
      // (since this bulk element is used everywhere in the bulk not just
      // in the augmented region, there may be no singular functions
      // associated with this element, and so we just want to return an
      // empty pointer and vector)
      if(i < Line_element_and_local_coordinate_at_knot.size())
      {
	return_pair = Line_element_and_local_coordinate_at_knot[i];
      }
      else
      {
	GeomObject* dummy_pt = nullptr;
	Vector<double> dummy_vec(1,0);
	return_pair = std::make_pair(dummy_pt, dummy_vec);
      }
      return return_pair;
    }

    // get the interpolated momentum Lagrange multiplier components at given local coordinates
    Vector<double> interpolated_lambda(const Vector<double>& s) const
    {	
      unsigned n_node = this->nnode();

      // make space for the basis functions
      Shape psif(n_node);

      // get 'em
      this->shape(s, psif);

      Vector<double> lambda_interp(DIM, 0.0);
      
      for(unsigned n=0; n<n_node; n++)
      {
	for(unsigned i=0; i<DIM; i++)
	{
	  // get the index of the ith LM which enforces the momentum eqs
	  unsigned lambda_index = index_of_lagrange_multiplier(n, i);
	  
	  lambda_interp[i] += this->raw_nodal_value(n, lambda_index) * psif[n];
	}
      }

      return lambda_interp;
    }

    double interpolated_lambda_p(const Vector<double>& s) const
    {
      const unsigned npressure_basis_fn = this->npres_nst();
	
      // make space for the pressure basis functions
      Shape psip(npressure_basis_fn);
      
      // get 'em
      this->pshape_nst(s, psip);

      double lambda_p_interp = 0.0;
       
      for(unsigned n=0; n<npressure_basis_fn; n++)
      {
	// get the index of the LM which enforces the continuity eqn
	unsigned lambda_p_index = index_of_lagrange_multiplier_p(n);
	  
	lambda_p_interp += this->raw_nodal_value(n, lambda_p_index) * psip[n];

      }

      return lambda_p_interp;
    }
    
    // tell this element it's augmented, i.e. use the pde-constrained residuals
    // rather than standard Navier-Stokes
    void set_augmented_element()
    {
      Is_augmented_element = true;
    }

    // Helper function to interrogate whether this element is augmented or not
    bool is_augmented_element() const
    {
      return Is_augmented_element;
    }

    void set_pressure_regularisation_factor(const double& p_reg)
    {
      Pressure_regularisation_factor = p_reg;
    }

    void set_velocity_regularisation_factor(const double& u_reg)
    {
      Velocity_regularisation_factor = u_reg;
    }
    
    /// Add the element's contribution to its residual vector
    inline void fill_in_contribution_to_residuals(Vector<double>& residuals)
    {
      // if this is an augmented element, it's residuals are determined
      // by PDE-constrained minimisation
      fill_in_generic_residual_contribution_pde_constrained_min(residuals,
								GeneralisedElement::Dummy_matrix,
								0);
    }

    /// \short Add the element's contribution to its residual vector and its
    /// Jacobian matrix
    inline void fill_in_contribution_to_jacobian(Vector<double>& residuals,
    						 DenseMatrix<double>& jacobian)
    {
      // N.B. this #ifdef is inside the function not outside, since
      // this class inherits from TaylorHood but if we want to do an FD Jacobian
      // we can't just not overload fill_in...jacobian() because then we'll
      // drop down into the analytic TaylorHood Jacobian
#ifdef USE_FD_JACOBIAN

      FiniteElement::fill_in_contribution_to_jacobian(residuals, jacobian);
      
#else
      fill_in_generic_residual_contribution_pde_constrained_min(residuals,
								jacobian,
								1);
#endif
    }

    // pin the ith LM which enforces momentum at the nth node of this element
    void pin_lagrange_multiplier(const unsigned& n, const unsigned& i,
				 double val = 0.0) const
    {
      unsigned nodal_index = index_of_lagrange_multiplier(n, i);

      this->node_pt(n)->pin(nodal_index);
      this->node_pt(n)->set_value(nodal_index, val);
    }

    // pin the LM which enforces continuity at the nth node of this element
    void pin_lagrange_multiplier_p(const unsigned& n, double val = 0.0) const
    {
      unsigned nodal_index = index_of_lagrange_multiplier_p(n);

      this->node_pt(n)->pin(nodal_index);
      this->node_pt(n)->set_value(nodal_index, val);
    }
    
    /// \short Wrapper function to pin all DIM + 1 PDE-constraint LMs at the
    /// nth node in this element
    void pin_pde_lagrange_multipliers(const unsigned& n, double val = 0.0) const
    {
      // pin the DIM momentum-enforcing LMs
      for(unsigned i=0; i<DIM; i++)
      {
	pin_lagrange_multiplier(n, i, val);
      }

      // and pin the continuity-enforcing LM if there is one here
      if(n <= DIM)
	pin_lagrange_multiplier_p(n, val);     
    }

    /// \short Wrapper function to pin the DIM + 1 momentum-enforcign LMs at the
    /// nth node in this element
    void pin_momentum_lagrange_multipliers(const unsigned& n, double val = 0.0) const
    {
      // pin the DIM momentum-enforcing LMs
      for(unsigned i=0; i<DIM; i++)
      {
	pin_lagrange_multiplier(n, i, val);
      }   
    }

    // compute the volume integral of the functional
    // \Pi = 1/2 |u|^2 over this element
    double integral_of_functional() const
    {      
      //Number of integration points
      unsigned n_intpt = this->integral_pt()->nweight();
   
      //Set the Vector to hold local coordinates
      Vector<double> s(DIM, 0.0);

      // the integral of the functional over this element
      double integral = 0.0;
      
      //Loop over the integration points
      for(unsigned ipt=0; ipt<n_intpt; ipt++)
      {
	//Assign values of s
	for(unsigned i=0; i<DIM; i++)
	{
	  s[i] = this->integral_pt()->knot(ipt, i);
	}
	
	//Get the integral weight
	double w = this->integral_pt()->weight(ipt);
   
	//Call the derivatives of the shape and test functions
	double J = this->J_eulerian_at_knot(ipt);

	double W = w * J;
	
	double functional = 0.0;

	// compute the functional \Pi = 1/2 |u|^2 at this integration point
	for(unsigned i=0; i<DIM; i++)
	  functional += 0.5 * pow(this->interpolated_u_nst(s, i), 2);

	// add the weighted value to the integral over this element
	integral += functional * W;
      }

      return integral;
    }
    
    // override the NavierStokesEquations residual function, as we
    // want to take advantage of the other goodies (source/body functions, du/dx etc.)
    // but don't actually want to solve the standard Stokes equations here
    // override the NavierStokesEquations residual function, as we
    // want to take advantage of the other goodies (source/body functions, du/dx etc.)
    // but don't actually want to solve the standard Stokes equations here
    void fill_in_generic_residual_contribution_pde_constrained_min(
      Vector<double>& residuals,
      DenseMatrix<double>& jacobian,
      const bool& flag = false)
    {
      // Return immediately if there are no dofs
      if (this->ndof()==0) return;

      //Find out how many nodes there are
      unsigned n_node = this->nnode();
       
      // Find out how many pressure dofs there are
      unsigned n_pres = this->npres_nst();

      //Find the indices at which the local velocities are stored
      unsigned u_nodal_index[DIM];
      for(unsigned i=0; i<DIM; i++)
      {
	u_nodal_index[i] = this->u_index_nst(i);
      }

      //Set up memory for the shape and test functions
      Shape psif(n_node), testf(n_node);
      DShape dpsifdx(n_node, DIM), dtestfdx(n_node, DIM);
      DShape dpsipdx(n_node, DIM), dtestpdx(n_node, DIM);
      
      //Set up memory for pressure shape and test functions
      Shape psip(n_pres), testp(n_pres);

      //Number of integration points
      unsigned n_intpt = this->integral_pt()->nweight();
   
      //Set the Vector to hold local coordinates
      Vector<double> s(DIM, 0.0);
 
      //Integers to store the local equations and unknowns
      int local_eqn = 0, local_unknown = 0;

      //Loop over the integration points
      for(unsigned ipt=0; ipt<n_intpt; ipt++)
      {
	//Assign values of s
	for(unsigned i=0; i<DIM; i++)
	{
	  s[i] = this->integral_pt()->knot(ipt, i);
	}
	
	//Get the integral weight
	double w = this->integral_pt()->weight(ipt);
   
	//Call the derivatives of the shape and test functions
	double J = 
	  this->dshape_and_dtest_eulerian_at_knot_nst(ipt, psif, dpsifdx, testf, dtestfdx);

	// Get the the pressure shape and test functions and their derivatives
	this->dpshape_and_dptest_eulerian_nst(s, psip, dpsipdx, testp, dtestpdx);
	  
	//Premultiply the weights and the Jacobian
	double W = w*J;
   
	// Calculate local values of the pressure and velocity components,
	// Lagrange multipliers and the functional to minimise
	// Allocate
	double interpolated_p = 0.0;
	Vector<double> interpolated_dp_dx(DIM, 0.0);
	
	Vector<double> interpolated_u(DIM, 0.0);
	DenseMatrix<double> interpolated_dudx(DIM, DIM, 0.0);
	
	Vector<double> interpolated_x(DIM, 0.0);
	
	Vector<double> interpolated_lambda(DIM, 0.0);
	double interpolated_lambda_p = 0.0;
	DenseMatrix<double> interpolated_dlambda_dx(DIM, DIM, 0.0);

	DenseMatrix<double> stress(DIM, DIM, 0.0);
	
	// Calculate pressure and pressure gradient
	for(unsigned l=0; l<n_pres; l++)
	{
	  interpolated_p += this->p_nst(l) * psip[l];

	  for(unsigned i=0; i<DIM; i++)
	  {
	    interpolated_dp_dx[i] += this->p_nst(l) * dpsipdx(l,i);
	  }
	}
   
	//Calculate velocities, derivatives and Lagrange multipliers:
  
	// Loop over nodes
	for(unsigned l=0; l<n_node; l++) 
	{  
	  //Loop over directions
	  for(unsigned i=0; i<DIM; i++)
	  {
	    // get the index of the ith LM which enforces the momentum eqs
	    unsigned lambda_index =
	      index_of_lagrange_multiplier(l, i);
	    
	    //Get the nodal value of the velocity component
	    double u_value = this->raw_nodal_value(l, u_nodal_index[i]);

	    interpolated_u[i] += u_value * psif[l];
	    interpolated_x[i] += this->raw_nodal_position(l,i) * psif[l];

	    // get the nodal LM component
	    double lambda_value = this->raw_nodal_value(l, lambda_index);
	    interpolated_lambda[i] += lambda_value * psif[l];
	    
	    //Loop over derivative directions and compute the derivatives
	    for(unsigned j=0; j<DIM; j++)
	    {                               
	      interpolated_dudx(i,j) += u_value * dpsifdx(l, j);
	      interpolated_dlambda_dx(i,j) += lambda_value * dpsifdx(l, j);
	    }
	  }
	}

	// Compute lambda_p, which is interpolated with the pressure shape functions
	unsigned Pconv_size = DIM + 1;
	for (unsigned k=0; k<Pconv_size; k++)
	{
	  // get the actual node number from the pressure node numbering
	  unsigned vertex_index = this->Pconv[k];

	  // get the nodal index of \lambda_p
	  unsigned lambda_p_index = index_of_lagrange_multiplier_p(vertex_index);

	  interpolated_lambda_p += this->raw_nodal_value(k, lambda_p_index) * psip[k];
	}

	// get the strain rate at this knot
	DenseMatrix<double> strain_rate(DIM, DIM, 0.0);
	this->strain_rate(s, strain_rate);
#ifdef PARANOID
	if(Stress_fct_pt == nullptr)
	{
	  throw OomphLibError(
	    "This TNavierStokesWithSingularityPdeConstrainedMinElement has no stress function pointer.",
	    OOMPH_CURRENT_FUNCTION,
	    OOMPH_EXCEPTION_LOCATION);
	}
#endif
	// compute the stress tensor
	stress = Stress_fct_pt(strain_rate, interpolated_p);
	
	// compute the derivatives of the functional we're trying to minimise
	Vector<double> dfunctional_du(DIM, 0.0);

	// have we actually got a function pointer for dPi/du? */	
	if(Dfunctional_du_fct_pt != nullptr)
	{
	  dfunctional_du = Dfunctional_du_fct_pt(interpolated_u);
	}

	// total scaled contribution of the singular functions and their
	// Eulerian derivatives
	Vector<double> u_sing_total(DIM+1, 0.0);
	DenseMatrix<double> dudx_sing_total(DIM, DIM, 0.0);
	DenseMatrix<double> stress_sing_total(DIM, DIM, 0.0);
	DenseMatrix<double> strain_rate_sing_total(DIM, DIM, 0.0);
	
	// get 'em
	total_singular_contribution_and_gradient_at_knot(ipt, u_sing_total,
							 dudx_sing_total);
	
	for(unsigned i=0; i<DIM; i++)
	{
	  for(unsigned j=0; j<DIM; j++)
	  {
	    strain_rate_sing_total(i,j) = 0.5 * (dudx_sing_total(i,j) +
						 dudx_sing_total(j,i));
	  }
	}

	// and finally compute the total singular stress
	double p_sing_total = u_sing_total[this->p_index_nst()];
	stress_sing_total = Stress_fct_pt(strain_rate_sing_total, p_sing_total);

	// Get the Eulerian position
	Vector<double> x(DIM, 0.0);
	for(unsigned i=0; i<DIM; i++) 
	{
	  x[i] = this->interpolated_x(s, i);
	}

	// Get the body force 
	Vector<double> body_force(3, 0.0);
	get_body_force(x, body_force);
	
	// ====================================================================
	// Now add to the relevant equations
	// ====================================================================
	
	// ================================================
	// LAGRANGE MULTIPLIER EQUATIONS
	// ================================================

	// momentum-enforcing LMs \lambda_i
	// ---------------------------------
	for(unsigned k=0; k<n_node; k++)
	{
	  for(unsigned j=0; j<DIM; j++)
	  {
	    unsigned lambda_index = index_of_lagrange_multiplier(k, j);
	    local_eqn = this->nodal_local_eqn(k, lambda_index);

	    // if this dof isn't pinned
	    if(local_eqn >= 0)
	    {
	      for(unsigned i=0; i<DIM; i++)
	      {
		// the Stokes momentum equations div\tau = 0 are the eqns
		// for the momentum-enforcing LMs	
		residuals[local_eqn] -=
		  (stress(i,j) + stress_sing_total(i,j)) * dpsifdx(k,i) * W;

		// add contribution from the body force
		// (not integrated by parts)
		residuals[local_eqn] += body_force[j] * psif[k] * W;
	      }

	      if(flag)
	      {
		//Loop over the velocity shape functions again
		for(unsigned l2=0; l2<n_node; l2++)
		{ 
		  //Loop over the velocity components again
		  for(unsigned i2=0; i2<DIM; i2++)
		  {		    
		    local_unknown = this->nodal_local_eqn(l2, u_nodal_index[i2]);
		    
		    // If at a non-zero degree of freedom add in the entry
		    if(local_unknown >= 0)
		    {
		      // dlambda/du
		      jacobian(local_eqn, local_unknown) 
			-= dpsifdx(l2,j) * dpsifdx(k,i2) * W;
                 
		      // Extra component if i2 = j
		      if(i2 == j)
		      {      
			/* Loop over dimensions again */
			for(unsigned n=0; n<DIM; n++)
			{ 
			  jacobian(local_eqn, local_unknown)
			    -= dpsifdx(l2, n) * dpsifdx(k,n) * W;
			}
		      }
		    }
		  } // end loop over LM components
		}

		// add contribution from p
		// -------------------------------
		
		/* Loop over pressure basis functions */
		for(unsigned l2=0; l2<n_pres; l2++)
		{
		  local_unknown = this->p_local_eqn(l2);
		    
		  /* If we are at a non-zero degree of freedom in the entry */
		  if(local_unknown >= 0)
		  {
		    // dlambda/dp
		    jacobian(local_eqn, local_unknown)
		      += psip[l2] * dpsifdx(k,j) * W;
		  }
		}

		// contributions from singular amplitudes (for augmented elements)
		// added in the singular amplitude residuals section below
	      }
	    }	    
	  }
	}

	// contiuity-enforcing LM \lambda_p
	// --------------------------------
	  
	// Loop over all vertex nodes	
	for (unsigned k=0; k<Pconv_size; k++)
	{
	  // get the actual node number from the pressure node numbering
	  unsigned vertex_index = this->Pconv[k];

	  // get the nodal index of \lambda_p
	  unsigned lambda_p_index = index_of_lagrange_multiplier_p(vertex_index);

	  // and it's corresponding local eqn number
	  local_eqn = this->nodal_local_eqn(k, lambda_p_index);

	  // if this dof isn't pinned
	  if(local_eqn >= 0)
	  {
	    for(unsigned j=0; j<DIM; j++)
	    {
	      residuals[local_eqn] +=
		(interpolated_dudx(j,j) + dudx_sing_total(j,j)) * psip[k] * W;
	    }

	    if(flag)
	    {	      
	      /*Loop over the velocity shape functions*/
	      for(unsigned l2=0; l2<n_node; l2++)
	      { 
		/*Loop over momentum LM components*/
		for(unsigned i2=0; i2<DIM; i2++)
		{		  
		  local_unknown = this->nodal_local_eqn(l2, u_nodal_index[i2]);

		  // If at a non-zero degree of freedom add in the entry
		  if(local_unknown >= 0)
		  {
		    // dlambda_p/du
		    jacobian(local_eqn, local_unknown)
		      += dpsifdx(l2,i2)* psip[k] * W;
		  }
		} /*End of loop over i2*/
	      } /*End of loop over l2*/
	    }
	  }
	  
	}

	// ================================================
	// Velocity Equations
	// ================================================
   
	// Loop over the test functions
	for(unsigned k=0; k<n_node; k++)
	{
	  // Loop over the velocity components
	  for(unsigned i=0; i<DIM; i++)
	  {	    
	    local_eqn = this->nodal_local_eqn(k, u_nodal_index[i]);

	    // if this dof isn't pinned
	    if(local_eqn >= 0)
	    {
	      //Add the functional minimisation term
	      residuals[local_eqn] +=
		Velocity_regularisation_factor * dfunctional_du[i] * psif[k] * W;

	      //Add the contribution from the momentum-enforcing LM
	      for(unsigned j=0; j<DIM; j++)
	      {
		residuals[local_eqn] -=
		  (interpolated_dlambda_dx(i,j) + interpolated_dlambda_dx(j,i)) * dpsifdx(k,j) * W;
	      }

	      // contribution from \lambda_p
	      residuals[local_eqn] += interpolated_lambda_p * dpsifdx(k,i) * W;
	      	      
	      //CALCULATE THE JACOBIAN
	      if(flag)
	      {
		//Loop over the velocity shape functions again
		for(unsigned l2=0; l2<n_node; l2++)
		{
		  if(Is_augmented_element)
		  {
		    local_unknown = this->nodal_local_eqn(l2, i);

		    if(local_unknown >= 0)
		    {
		      // du_k/du_l2
		      jacobian(local_eqn, local_unknown) +=
			Velocity_regularisation_factor * psif[l2] * psif[k] * W;
		    }
		  }

		  //Loop over the velocity components again
		  for(unsigned i2=0; i2<DIM; i2++)
		  {
		    unsigned lm_index = index_of_lagrange_multiplier(l2, i2);
		    local_unknown = this->nodal_local_eqn(l2, lm_index);
		    
		    // If at a non-zero degree of freedom add in the entry
		    if(local_unknown >= 0)
		    {
		      // du/dlambda
		      jacobian(local_eqn, local_unknown) 
			-= dpsifdx(l2,i) * dpsifdx(k,i2) * W;
                 
		      // Extra component if i2 = i
		      if(i2 == i)
		      {      
			/* Loop over dimensions again */
			for(unsigned n=0; n<DIM; n++)
			{
			  // du/dlambda
			  jacobian(local_eqn, local_unknown)
			    -= dpsifdx(l2, n) * dpsifdx(k,n) * W;
			}
		      }
		    }
		  } // end loop over LM components
		}

		// add contribution from Lambda_p
		// -------------------------------
		
		/* Loop over pressure basis functions */
		for(unsigned l2=0; l2<n_pres; l2++)
		{
		  unsigned lm_p_index = index_of_lagrange_multiplier_p(l2);
		  local_unknown = this->nodal_local_eqn(l2, lm_p_index);
		    
		  /* If we are at a non-zero degree of freedom in the entry */
		  if(local_unknown >= 0)
		  {
		    // du/dlambda_p
		    jacobian(local_eqn, local_unknown)
		      += psip[l2] * dpsifdx(k,i) * W;
		  }
		}
	      } /*End of Jacobian calculation*/
         
	    } //End of if not boundary condition statement
       
	  } //End of loop over dimensions
	} //End of loop over shape functions

	
	// ================================================
	// Pressure Equations
	// ================================================
	
	//Loop over the shape functions
	for(unsigned k=0; k<n_pres; k++)
	{
	  local_eqn = this->p_local_eqn(k);
	  
	  //If not a boundary conditions
	  if(local_eqn >= 0)
	  {
	    if(Is_augmented_element)
	    {
	      // add the pressure penalty term from the functional to be minimised
	      residuals[local_eqn] +=
		Pressure_regularisation_factor * interpolated_p * psip[k] * W;
	    }
	    
	    //Loop over velocity components
	    for(unsigned j=0; j<DIM; j++)
	    {
	      // contribution from momentum-enforcing LM
	      residuals[local_eqn] += interpolated_dlambda_dx(j,j) * psip[k] * W;
	    }
	    
	    // Jacobian?
	    if(flag)
	    {
	      if(Is_augmented_element)
	      {
		for(unsigned k2=0; k2<n_pres; k2++)
		{
		  local_unknown = this->p_local_eqn(k2);

		  if(local_unknown >= 0)
		  {
		    // dp/dp
		    jacobian(local_eqn, local_unknown) +=
		      Pressure_regularisation_factor * psip[k2] * psip[k] * W;
		  }
		}
	      }
	      
	      /*Loop over the velocity shape functions*/
	      for(unsigned l2=0; l2<n_node; l2++)
	      { 
		/*Loop over momentum LM components*/
		for(unsigned i2=0; i2<DIM; i2++)
		{
		  unsigned lm_index = index_of_lagrange_multiplier(l2, i2);
		  local_unknown = this->nodal_local_eqn(l2, lm_index);

		  // If at a non-zero degree of freedom add in the entry
		  if(local_unknown >= 0)
		  {
		    // dp/dlambda
		    jacobian(local_eqn, local_unknown)
		      += dpsifdx(l2,i2) * psip[k] * W;
		  }
		} /*End of loop over i2*/
	      } /*End of loop over l2*/
	      
	    } /*End of Jacobian calculation*/
       
	  } //End of if not boundary condition
	  
	} //End of loop over p shape functions

	// =================================================
	// contribution to the singular amplitude equations
	// =================================================

	// If we're in the non-augmented bulk so there's no singular line element
	// by definition, so we're done on this knot's residuals
	if(!Is_augmented_element)
	  continue;
	
	// the (pointer to the) singular line element associated with this knot
	SingularLineElement* sing_el_pt = nullptr;

	LagrangianCoordinates lagr_coords_at_knot;
	Vector<double> s_singular_el(1, 0.0);

	// get the edge coordinates at this knot, the corresponding singular
	// line element and it's local coordinates
	this->lagr_coords_and_singular_element_at_knot(ipt,
						       lagr_coords_at_knot,
						       sing_el_pt,
						       s_singular_el);
	// have we actually got a pointer?
	if(sing_el_pt == nullptr)
	{	  
	  // something's gone wrong here, we're in an augmented element
	  // which has no singular line element associated with one of its knots
	  ostringstream error_message;

	  error_message << "Error: this augmented bulk element has no "
			<< "singular line element pointer for knot # "
			<< ipt << "\n" << std::endl;
	    
	  throw OomphLibError(error_message.str().c_str(),
			      OOMPH_CURRENT_FUNCTION,
			      OOMPH_EXCEPTION_LOCATION);	  
	}

	// number of singular functions we're subtracting
	const unsigned nsing_fct = sing_el_pt->nsingular_fct();

	// get the list of singular function IDs
	const Vector<unsigned> sing_ids = sing_el_pt->singular_fct_ids();

	// unsclaed singular functions and their gradients
	Vector<Vector<double> > u_sing_unscaled(nsing_fct);
	Vector<DenseMatrix<double> > dudx_sing_unscaled(nsing_fct);

	// unslcaed pressure
	Vector<double> p_sing_unscaled(nsing_fct);
	  
	// loop over all the singular functions this element knows about and
	// compute the sum of their contributions to the total solution      
	for(unsigned ising=0; ising < nsing_fct; ising++)
	{
	  // get the ID
	  unsigned sing_fct_id = sing_ids[ising];
	
	  // unscaled stuff
	  // ----------------
	  // unscaled singular function 
	  u_sing_unscaled[ising] =	    
	    sing_el_pt->unscaled_singular_fct(lagr_coords_at_knot, sing_fct_id);	  
	 
	  // unscaled singular gradient 
	  dudx_sing_unscaled[ising] =
	    sing_el_pt->gradient_of_unscaled_singular_fct(lagr_coords_at_knot, sing_fct_id);
	 
	  // unscaled singular pressure 
	  p_sing_unscaled[ising] = u_sing_unscaled[ising][this->p_index_nst()];
	}

	// number of nodes in this singular line element
	const unsigned nnode_sing = sing_el_pt->nnode();

	// get the shape functions and their Eulerian derivatives from the
	// singular element at the local coordinate which corresponds to this
	// knot in this face element
	Shape psi_sing(nnode_sing);
	Shape test_dummy(nnode_sing);
	DShape dpsi_sing_dx(nnode_sing, DIM);

	sing_el_pt->shape_and_test(s_singular_el, psi_sing, test_dummy);
	sing_el_pt->dshape_eulerian(lagr_coords_at_knot, s_singular_el, dpsi_sing_dx);

	// loop over the nodes in the singular element associated with this integration point
	for(unsigned ising_node=0; ising_node<nnode_sing; ising_node++)
	{
	  // external data index for this singular node
	  unsigned ext_index = this->C_external_data_index_at_knot[ipt][ising_node];

	  // loop over the singular functions
	  for(unsigned ising=0; ising < sing_el_pt->nsingular_fct(); ising++)
	  {
	    // external equation number which determines the singular amplitude
	    int external_eqn_c = this->external_local_eqn(ext_index, ising);
	  
	    // if this singular amplitude isn't pinned
	    if(external_eqn_c >= 0)	      
	    {
	      for(unsigned i=0; i<DIM; i++)
	      {
		residuals[external_eqn_c] += interpolated_lambda_p *
		  ( dudx_sing_unscaled[ising](i,i) * psi_sing[ising_node]
		    + u_sing_unscaled[ising][i] * dpsi_sing_dx(ising_node, i) ) * W;
		  
		  
		for(unsigned j=0; j<DIM; j++)
		{
		  // contribution of the momentum-enforcing LMs to the C equations
		  residuals[external_eqn_c] += interpolated_dlambda_dx(i,j) *
		    ( p_sing_unscaled[ising] *
		      psi_sing[ising_node] * delta(i,j)
		      - ( psi_sing[ising_node] * (dudx_sing_unscaled[ising](i,j) +
						  dudx_sing_unscaled[ising](j,i)) +
			  dpsi_sing_dx(ising_node,j) * u_sing_unscaled[ising][i] +
			  dpsi_sing_dx(ising_node,i) * u_sing_unscaled[ising][j] )
		      ) * W;		  
		}
	      } // end loop over dimensions

	      // Jacobian?
	      if (flag == 1)
	      {		
		for(unsigned k=0; k<n_node; k++)
		{
		  for(unsigned i=0; i<DIM; i++)
		  {
		    for(unsigned j=0; j<DIM; j++)
		    {			
		      unsigned lambda_index = index_of_lagrange_multiplier(k, i);
		      local_eqn = this->nodal_local_eqn(k, lambda_index);

		      if(local_eqn >= 0)
		      {
			double jac_entry = dpsifdx(k,j) *
			  (delta(i,j) * p_sing_unscaled[ising] * psi_sing[ising_node]
			   - (u_sing_unscaled[ising][i] * dpsi_sing_dx(ising_node, j) +
			      u_sing_unscaled[ising][j] * dpsi_sing_dx(ising_node, i) +
			      psi_sing[ising_node] * (dudx_sing_unscaled[ising](i,j) +
						      dudx_sing_unscaled[ising](j,i) )
			     ) ) * W;

			// dC/dlambda-mom
			jacobian(external_eqn_c, local_eqn) += jac_entry;

			// symmetric Jacobian, so add the corresponding term here
			jacobian(local_eqn, external_eqn_c) += jac_entry;
		      }		      
		    }
		  }
		}
		
		// Now the contributions from the continuity-enforcing LM \lambda_p
		for (unsigned k=0; k<Pconv_size; k++)
		{
		  // get the actual node number from the pressure node numbering
		  unsigned vertex_index = this->Pconv[k];

		  // get the nodal index of \lambda_p
		  unsigned lambda_p_index = index_of_lagrange_multiplier_p(vertex_index);

		  // and it's corresponding local eqn number
		  local_eqn = this->nodal_local_eqn(k, lambda_p_index);

		  // if this dof isn't pinned
		  if(local_eqn >= 0)
		  {
		    for(unsigned j=0; j<DIM; j++)
		    {
		      double jac_entry = psip[k] *
			(dpsi_sing_dx(ising_node, j) * u_sing_unscaled[ising][j] +
			 dudx_sing_unscaled[ising](j,j) * psi_sing[ising_node] ) * W;
		    		    
		      // dC/dlambda_p
		      jacobian(external_eqn_c, local_eqn) += jac_entry;

		      // and add the symmetric entry
		      jacobian(local_eqn, external_eqn_c) += jac_entry;
		    }
		  }
		}
	      }
	    } // end external_eqn_c >= check
	  } // end ext_index >= 0 check
	
	} // end loop over singular nodes
      
      } // end loop over integration points

    } // end fill_in_generic_residual_contribution_pde_constrained_min()

    void compute_centroid(Vector<double>& centroid) const
    {
      // resize and zero out
      centroid.resize(DIM, 0.0);
      centroid.initialise(0.0);

      const unsigned nvertex = DIM + 1;
      
      // loop over the vertices
      for(unsigned n=0; n<nvertex; n++)
      {
	for(unsigned i=0; i<DIM; i++)
	  centroid[i] += this->node_pt(n)->x(i);
      }

      // now average
      for(double& xi : centroid)
	xi /= nvertex;
    }
    
    // override to output scalar values for paraview format
    void scalar_value_fct_paraview(std::ofstream& file_out,
				   const unsigned& i,
				   const unsigned& nplot,
				   FiniteElement::SteadyExactSolutionFctPt exact_soln_pt) const
    {
      //Vector of local coordinates
      Vector<double> s(DIM);
   
      // Vector for coordinates
      Vector<double> x(DIM);
    
      // Exact solution Vector
      Vector<double> exact_soln(4,0.0);

      // tolerance on node sitting on the disk surface
      const double tol = 1e-8;
	
      // Loop over plot points
      unsigned num_plot_points = this->nplot_points_paraview(nplot);
      
      for (unsigned iplot=0; iplot<num_plot_points; iplot++)
      {
	// Get local coordinates of plot point
	this->get_s_plot(iplot,nplot,s);
         
	// Get x position as Vector
	this->interpolated_x(s, x);

	// catch the case where the node should be sitting on the lower
	// side of the disk
	if((abs(x[2]) < tol) && Is_lower_disk_element)
	  x[2] = -tol;
	 
	// Get exact solution at this point
	(*exact_soln_pt)(x, exact_soln);

	// output this component
	file_out << exact_soln[i] << std::endl;
      }
    }
    
    // specific override for column names for paraview output
    std::string scalar_name_paraview(const unsigned& i) const
    {
      std::string name;

      switch(i)
      {
	case 0:
	  name = "Total Velocity x";
	  break;
	case 1:
	  name = "Total Velocity y";
	  break;
	case 2:
	  name = "Total Velocity z";
	  break;
	case 3:
	  name = "Total Pressure";
	  break;
	case 4:
	  name = "FE Velocity x";
	  break;
	case 5:
	  name = "FE Velocity y";
	  break;
	case 6:
	  name = "FE Velocity z";
	  break;
	case 7:
	  name = "FE Pressure";
	  break;
	case 8:
	  name = "Singular Velocity (broadside) x";
	  break;
	case 9:
	  name = "Singular Velocity (broadside) y";
	  break;
	case 10:
	  name = "Singular Velocity (broadside) z";
	  break;
	case 11:
	  name = "Singular Pressure (broadside)";
	  break;
	case 12:
	  name = "Singular Velocity (in-plane) x";
	  break;
	case 13:
	  name = "Singular Velocity (in-plane) y";
	  break;
	case 14:
	  name = "Singular Velocity (in-plane) z";
	  break;
	case 15:
	  name = "Singular Pressure (in-plane)";
	  break;
	case 16:
	  name = "Singular Velocity (in-plane rotation) x";
	  break;
	case 17:
	  name = "Singular Velocity (in-plane rotation) y";
	  break;
	case 18:
	  name = "Singular Velocity (in-plane rotation) z";
	  break;
	case 19:
	  name = "Singular Pressure (in-plane rotation)";
	  break;
	default:
	  name = "V"+StringConversion::to_string(i);
	  break;	  
      }

      return name;
    }
    
  private:

    /// \short Number of plot points - this isn't a free choice at point where
    /// we want to output, because we need to have pre-computed the (\rho,\zeta,\phi)
    /// coordinates of each plot point, so this number is set when
    /// set_line_element_and_local_coordinate_at_plot_point() is called
    unsigned Nplot;
    
    /// \short Flag to store whether this is a lower disk element
    /// (useful for output stuff where we're using Eulerian coordinates and need to
    /// decide whether a point on the disk surface has phi=0 or phi=2pi to get
    /// the pressure jump correct)
    bool Is_lower_disk_element;
    
    /// \short Boolean value to determine if this element is augmented with singular 
    // functions (and whose residuals are determined by PDE-constrained minimisation)
    // or whether this is a non-augmented element governed by the standard
    // (Navier-)Stokes equations
    bool Is_augmented_element;
    
    // vector which gives the nodal index of the first momentum-enforcing
    // Lagrange multiplier for each node in this element
    Vector<unsigned> Lambda_index;
 
    /// Pointer to exact non-singular fct (only for post-processing!)
    ExactNonSingularFctPt Exact_non_singular_fct_pt;
    
    /// \short Function pointer to the derivative of the functional 
    /// w.r.t. the solution
    DfunctionalDuFctPt Dfunctional_du_fct_pt;

    /// Pointer to function which computes the stress
    StressFctPt Stress_fct_pt;

    /// Pointer to function which computes the body force acting at a given point
    BodyForceFctPt Body_force_fct_pt;
		    
    /// \short Number of singular functions we're subtracting (only used for
    /// output purposes, for bulk elements which have no singular line element pointer
    /// but need to get the number of output columns to match the elements
    /// with singular functions)
    unsigned Nsingular_fct;
            
    /// \short Edge coordinates (\rho, \zeta, \phi) of each of this element's
    /// plot points
    Vector<LagrangianCoordinates> Lagrangian_coordinates_at_plot_point;

    // same as above but for each knot point (to compute error)
    Vector<LagrangianCoordinates> Lagrangian_coordinates_at_knot_point;

    // same as above but for each node
    Vector<LagrangianCoordinates> Nodal_lagr_coordinates;
    
    /// \short The line element and its local coordinate(s) which correspond to the zeta
    /// values at each plot point in this bulk element. N.B. the element
    /// is 1D so only has 1 local coordinate, but we keep this as a vector
    /// for consistency with the general case
    Vector<std::pair<GeomObject*, Vector<double> > >
      Line_element_and_local_coordinate_at_plot_point;

    // same as above but for knot points (to compute error)
    Vector<std::pair<GeomObject*, Vector<double> > >
      Line_element_and_local_coordinate_at_knot;

     // same as above but for nodes
    Vector<std::pair<GeomObject*, Vector<double> > >
      Line_element_and_local_coordinate_at_node;

    /// \short Indices of external Data that store the values of the amplitudes of
    /// the singular functions - the outer index is the integration point number,
    /// and the inner index is the node number of the singular line element
    /// which coresponds to the integration point
    Vector<Vector<unsigned> > C_external_data_index_at_knot;

    /// \short scaling factor for additional regularisation of the FE pressure
    /// in the augmented region in the functional to be minimised
    double Pressure_regularisation_factor;

    /// \short scaling factor for 
    double Velocity_regularisation_factor;

    double random_var_to_alter_mem_map;
  };


  //////////////////////////////////////////////////////////////////////
  //////////////////////////////////////////////////////////////////////
  //////////////////////////////////////////////////////////////////////


  //=======================================================================
  /// Face geometry for the TNavierStokesWithSingularityPdeConstrainedMinElement elements: The spatial 
  /// dimension of the face elements is one lower than that of the
  /// bulk element but they have the same number of points
  /// along their 1D edges.
  //=======================================================================
  template <unsigned DIM>
    class FaceGeometry<TNavierStokesWithSingularityPdeConstrainedMinElement<DIM> > :
    public virtual TElement<DIM-1,TNavierStokesWithSingularityPdeConstrainedMinElement<DIM>::_NNODE_1D_>
  {

  public:
 
    /// \short Constructor: Call the constructor for the
    /// appropriate lower-dimensional TElement
  FaceGeometry() : TElement<DIM-1,TNavierStokesWithSingularityPdeConstrainedMinElement<DIM>::_NNODE_1D_>() {}

    FaceGeometry(const FaceGeometry&)
    {
      BrokenCopy::broken_copy(
    	"FaceGeometry<TNavierStokesWithSingularityPdeConstrainedMinElement<DIM> >");
    }

    void operator=(const FaceGeometry&)
    {
      BrokenCopy::broken_assign(
	"FaceGeometry<TNavierStokesWithSingularityPdeConstrainedMinElement<DIM> >");
    }
    
    /* FaceGeometry(const FaceGeometry&&) */
    /* { */
    /*   	throw OomphLibError( */
    /* 	"Move assignment operator not implemented for ", */
    /* 	OOMPH_CURRENT_FUNCTION, */
    /* 	OOMPH_EXCEPTION_LOCATION); */
    /* } */
  };

  ////////////////////////////////////////////////////////////////////////
  ////////////////////////////////////////////////////////////////////////
  ////////////////////////////////////////////////////////////////////////


  //=======================================================================
  /// Face geometry for the 1D TNavierStokesWithSingularityPdeConstrainedMinElement elements: Point elements
  //=======================================================================
  template <>
    class FaceGeometry<TNavierStokesWithSingularityPdeConstrainedMinElement<1> >: 
    public virtual PointElement
    {

    public:
 
      /// \short Constructor: Call the constructor for the
      /// appropriate lower-dimensional TElement
    FaceGeometry() : PointElement() {} 

    };
  
}

#endif
