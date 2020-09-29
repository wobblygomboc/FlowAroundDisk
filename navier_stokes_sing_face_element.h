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

#include "navier_stokes.h"

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

namespace oomph
{

  // type to specify we're using the (\rho,\zeta,\phi) edge coordinate system
  // not the global Cartesian system
  struct EdgeCoordinates
  {
    // default the coordinates to zero
    EdgeCoordinates()
      {
	rho  = 0.0;
	zeta = 0.0;
	phi  = 0.0;
	
	// QUEHACERES
	p0 = 0;
      }
    
    double rho;
    double zeta;
    double phi;

    // QUEHACERES
    double p0;
    
    static const unsigned ncoord = 3;
  };
  
  ////////////////////////////////////////////////////////////////////////
  ////////////////////////////////////////////////////////////////////////
  ////////////////////////////////////////////////////////////////////////
  // QUEHACERES put useful comment here once we've hashed out the details
  
  template <unsigned NNODE_1D>
    class ScalableSingularityForNavierStokesLineElement : 
    public virtual FaceGeometry<TElement<2,NNODE_1D> >, public virtual FaceElement
  {
  public:

    // typedefs for the functions which provide the singular solutions and
    // their derivatives
    typedef Vector<double>(*UnscaledSingSolnFctPt) (const EdgeCoordinates& rzp_coords);

    // QUEHACERES probably want to delete this if we don't need it for the C residuals
    typedef DenseMatrix<double>(*GradientOfUnscaledSingSolnFctPt)
      (const EdgeCoordinates& rzp_coords);

    // function which provides dzeta/dx
    typedef Vector<double>(*DzetaDxFctPt)(const EdgeCoordinates& edge_coords);
	
    /// \short Constructor, takes a pointer to the 2D disk element and the
    /// index of the face to which this element is attached. 
    ScalableSingularityForNavierStokesLineElement(FiniteElement* const& disk_el_pt, 
						  const int& face_index,
						  // QUEHACERES delete const unsigned& nsingular_function,
						  std::map<Node*,Node*>& existing_duplicate_node_pt,
						  const bool& use_zeta_2pi_instead_of_0=false);
    
    /// Broken empty constructor
    ScalableSingularityForNavierStokesLineElement()
    {
       throw OomphLibError(
	"Don't call empty constructor for NavierStokesWithSingularityLineElement",
	OOMPH_CURRENT_FUNCTION,
	OOMPH_EXCEPTION_LOCATION);
    }

    // Destructor, need to clear up the nodes since these were created separately
    // with new
    ~ScalableSingularityForNavierStokesLineElement()
    {
      for(unsigned j=0; j<nnode(); j++)
	delete node_pt(j);
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
    
    // QUEHACERES delete, now using add_unscaled_sing...()
    /* ///Function to get pointer to the unscaled version of a singular function */
    /* UnscaledSingSolnFctPt& unscaled_singular_fct_pt(const unsigned& sing_fct_id) */
    /* { */
    /*   // get the index for this ID from the map */
    /*   // (std::map::at() throws an exception if this index doesn't exist) */
    /*   unsigned sing_fct_index = Singular_fct_index.at(sing_fct_id); */
      
    /*   return Unscaled_singular_fct_pt[sing_fct_index]; */
    /* } */

    /* ///Function to get pointer to unscaled version of the gradient of */
    /* // a singular function */
    /* GradientOfUnscaledSingSolnFctPt& */
    /*   gradient_of_unscaled_singular_fct_pt(const unsigned& sing_fct_id) */
    /* { */
    /*   // get the index for this ID from the map */
    /*   // (std::map::at() throws an exception if this index doesn't exist) */
    /*   unsigned sing_fct_index = Singular_fct_index.at(sing_fct_id); */
      
    /*   return Gradient_of_unscaled_singular_fct_pt[sing_fct_index]; */
    /* } */

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

      if(sing_fct_pt == 0)
      {
	throw OomphLibError(
	"Unscaled singular function pointer is null!",
	OOMPH_CURRENT_FUNCTION,
	OOMPH_EXCEPTION_LOCATION);
      }

      if(grad_of_sing_fct_pt == 0)
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
    Vector<double> unscaled_singular_fct(const EdgeCoordinates& edge_coords,
					 const unsigned& sing_fct_id) const
    {
      // QUEHACERES delete
      /* if(Unscaled_singular_fct_pt[ising] == 0) */
      /* { */
      /* 	Vector<double> empty(rzp_coords.ncoord, 0.0); */
      /* 	return empty; */
      /* } */
      
      // get the index for this ID from the map
      // (std::map::at() throws an exception if this index doesn't exist)
      unsigned sing_fct_index = Singular_fct_index.at(sing_fct_id);
      
      return Unscaled_singular_fct_pt[sing_fct_index](edge_coords); 
    }

    ///Compute unscaled version of gradient of the requested singular function
    DenseMatrix<double> gradient_of_unscaled_singular_fct(const EdgeCoordinates& edge_coords,
							  const unsigned& sing_fct_id) const
    {
      // get the index for this ID from the map
      // (std::map::at() throws an exception if this index doesn't exist)
      unsigned sing_fct_index = Singular_fct_index.at(sing_fct_id);
      
      return Gradient_of_unscaled_singular_fct_pt[sing_fct_index](edge_coords); 
    }

    ///Compute scaled version of the requested singular function
    Vector<double> singular_fct(const EdgeCoordinates& edge_coords,
				const Vector<double>& s,
				const unsigned& sing_fct_id) const
    {
      // get number of unknowns in the problem; plus one because we want pressure as well
      // as the velocity components
      const unsigned nvalue = edge_coords.ncoord + 1;

      // storage for the scaled basis functions
      Vector<double> scaled_singular_fct(nvalue, 0.0);

      // get the unscaled functions
      Vector<double> u_sing_unscaled = unscaled_singular_fct(edge_coords, sing_fct_id);

      // get the singular amplitude associated with this ID and these local coords
      double c = interpolated_amplitude(s, sing_fct_id);
      
      // scale 'em
      for(unsigned i=0; i<nvalue; i++)
      {
	scaled_singular_fct[i] = c * u_sing_unscaled[i];
      }

      /* // QUEHACERES for debug */
      /* if(ising == 1) */
      /* 	scaled_singular_fct[0] += 1; */
      
      return scaled_singular_fct;
    }
            
    ///Compute scaled version of gradient of the requested singular function
    DenseMatrix<double> gradient_of_singular_fct(const EdgeCoordinates& edge_coords,
						 const Vector<double>& s,
						 const unsigned& sing_fct_id) const
    {
      // get the interpolated singular amplitude
      const double c = interpolated_amplitude(s, sing_fct_id);

      // get the unscaled singular velocity
      Vector<double> u_sing_unscaled = unscaled_singular_fct(edge_coords, sing_fct_id);
      
      // get the unscaled singular velocity gradient
      DenseMatrix<double> dudx_sing_unscaled =
	gradient_of_unscaled_singular_fct(edge_coords, sing_fct_id);

      // get the Lagrangian derivative of the singular amplitude, dc/dzeta
      double dc_dzeta = interpolated_dc_dzeta(s, sing_fct_id);

#ifdef PARANOID
      if (Dzeta_dx_fct_pt == 0)
      {
	throw OomphLibError("Error: function pointer for dzeta/dx hasn't been set\n",
			    OOMPH_CURRENT_FUNCTION,
			    OOMPH_EXCEPTION_LOCATION);
      }
#endif

      // get the Eulerian derivatives dzeta/dx from the function pointer
      Vector<double> dzeta_dx = Dzeta_dx_fct_pt(edge_coords);
      
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
    Vector<double> total_singular_contribution(const EdgeCoordinates& edge_coords,
					       const Vector<double>& s) const
    {
      // total singular part of the solution
      Vector<double> u_sing_total(edge_coords.ncoord + 1, 0.0);

      // loop over each singular function and add the contribution to the total
      for(std::map<unsigned,unsigned>::const_iterator it = Singular_fct_index.begin();
	  it != Singular_fct_index.end(); it++)
      // QUEHACERES delete
      /* for(unsigned ising=0; ising<Nsingular_fct; ising++) */
      {
	// get the singular function ID
	unsigned sing_fct_id = it->first;
	
	// get the contribution of the ith singular function
	Vector<double> u_sing = singular_fct(edge_coords, s, sing_fct_id);

	// add it to the total
	for(unsigned j=0; j<u_sing.size(); j++)
	  u_sing_total[j] += u_sing[j];
      }
      return u_sing_total;
    }

    // wrapper to get the total contribution of the derivatives of the
    // singular functions
    DenseMatrix<double>
      total_singular_gradient_contribution(const EdgeCoordinates& rzp_coords,
					   const Vector<double>& s) const
    {
      const unsigned nval = rzp_coords.ncoord;
      
      // total of the singular derivatives      
      DenseMatrix<double> dudx_total(nval, nval, 0.0);
            

      // loop over each singular function and add the contribution to the total
      for(std::map<unsigned,unsigned>::const_iterator it = Singular_fct_index.begin();
	  it != Singular_fct_index.end(); it++)
      // QUEHACERES delete
      /* for(unsigned ising=0; ising<Nsingular_fct; ising++) */
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
    
    /// \Short function to compute the interpolated amplitude of the
    /// nth singular function at the local coordinate s
    double interpolated_amplitude(const Vector<double>& s, const unsigned& sing_fct_id) const
    {
      // get the index for this ID from the map
      // (std::map::at() throws an exception if this index doesn't exist)
      const unsigned sing_fct_index = Singular_fct_index.at(sing_fct_id);
      
      // make space for the shape functions
      Shape psi(nnode());
      
      // get 'em
      shape(s, psi);

      // amplitude of this singular function
      double interpolated_c = 0;

      // loop over each node in this element and add its contribution to the
      // interpolated amplitude
      for(unsigned j=0; j<nnode(); j++)
      {
	interpolated_c += nodal_value(j, sing_fct_index) * psi[j];
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

      // make space for the shape functions and their derivatives
      Shape psi(nnode());      
      DShape dpsi_ds(nnode(), Dim);

      // Find values of shape functions and their derivatives
      this->dshape_local(s, psi, dpsi_ds);

      // compute interpolated local derivatives dc/ds and dzeta/ds
      Vector<double> interpolated_dc_ds(Dim, 0.0);
      Vector<double> interpolated_dzeta_ds(Dim, 0.0);
      
      for(unsigned l=0; l<nnode(); l++) 
      {
	for(unsigned j=0; j<Dim; j++)
	{
	  interpolated_dc_ds[j]    += nodal_value(l, sing_fct_index) * dpsi_ds(l,j);
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
    double zeta_nodal(const unsigned &n, const unsigned &k, const unsigned &i) const
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

	// now loop over the singular functions and get their
	// interpolated amplitudes at this plot point
	for(std::map<unsigned,unsigned>::iterator it = Singular_fct_index.begin();
	  it != Singular_fct_index.end(); it++)
	// QUEHACERES delete
	/* for(unsigned ising=0; ising<Nsingular_fct; ising++) */
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
				       const Vector<double>& ampl)
    {
      // check the user has provided enough values for the whole element
      if(ampl.size() != nnode())
      {
	ostringstream error_message;

	error_message << "Error: need an amplitude for each node in this element;\n"
		      << ampl.size() << " provided, but this element has " << nnode()
		      << " nodes.\n";
	
	throw OomphLibError(error_message.str(),
			    OOMPH_CURRENT_FUNCTION,
			    OOMPH_EXCEPTION_LOCATION);
      }

      // QUEHACERES delete
      /* // check we have enough singular functions */
      /* if(n >= nsingular_fct()) */
      /* { */
      /* 	ostringstream error_message; */
      /* 	error_message << "Error: provided amplitudes for the singular function ising =" */
      /* 		      << n << ", but this element only has " << Nsingular_fct */
      /* 		      << " singular functions\n"; */
	
      /* 	throw OomphLibError(error_message.str(), */
      /* 			    OOMPH_CURRENT_FUNCTION, */
      /* 			    OOMPH_EXCEPTION_LOCATION); */
      /* } */
      
      // get the index for this ID from the map
      // (std::map::at() throws an exception if this index doesn't exist)
      unsigned sing_fct_index = Singular_fct_index.at(sing_fct_id);
      
      for(unsigned j=0; j<nnode(); j++)
      {	
	// pin it so that the amplitude is no longer a dof
	node_pt(j)->pin(sing_fct_index);

	// set its value to the imposed amplitude
	node_pt(j)->set_value(sing_fct_index, ampl[j]);	
      }      
    } 

    /// Reset all singular amplitudes to compute r_c properly via integral
    void dont_impose_singular_fct_amplitude()
    {
      for(unsigned j=0; j<nnode(); j++)
      {
	for(unsigned i=0; i<node_pt(i)->nvalue(); i++)
	{
	  node_pt(j)->unpin_value(i);
	}
      }
    } 

    // unpin a particular singular amplitude
    void dont_impose_singular_fct_amplitude(const unsigned& sing_fct_id)
    {
      // get the index for this ID from the map
      // (std::map::at() throws an exception if this index doesn't exist)
      unsigned sing_fct_index = Singular_fct_index.at(sing_fct_id);
      
      for(unsigned j=0; j<nnode(); j++)
      {
	  node_pt(j)->unpin_value(sing_fct_index);
      }
    }
    
    /// Add the element's contribution to its residual vector
    inline void fill_in_contribution_to_residuals(Vector<double>& residuals)
    {
      //Call the generic residuals function with flag set to 0
      //using a dummy matrix argument
      fill_in_generic_residual_contribution_navier_stokes_sing_fct(
	residuals,GeneralisedElement::Dummy_matrix, 0);
    }

    // function to return whether we're imposing the amplitude of the nth singular fct
    bool is_singular_fct_amplitude_imposed(const unsigned& sing_fct_id)
    {
      // get the index for this ID from the map
      // (std::map::at() throws an exception if this index doesn't exist)
      unsigned sing_fct_index = Singular_fct_index.at(sing_fct_id);
      
      // we're not allowing some nodes to be pinned and others unpinned at the moment,
      // so it's enough to check the first node
      return node_pt(0)->is_pinned(sing_fct_index);
    }

    // returns a list of the IDs, i.e. the keys to the ID->index map 
    Vector<unsigned> singular_fct_ids() const
    {
      return Singular_fct_id;
    }
      
      
  private:

    /// Add the element's contribution to its residual vector
    inline void fill_in_generic_residual_contribution_navier_stokes_sing_fct(
      Vector<double>& residuals,
      DenseMatrix<double>& jacobian,
      const unsigned& flag);

    /// \short Function to compute the shape and test functions and to return 
    /// the Jacobian of mapping between local and global (Eulerian)
    /// coordinates
    inline double shape_and_test(const Vector<double> &s, Shape &psi, Shape &test)
      const
    {
      //Find number of nodes
      unsigned n_node = nnode();

      //Get the shape functions
      shape(s,psi);

      //Set the test functions to be the same as the shape functions
      for(unsigned i=0; i<n_node; i++)
      {
	test[i] = psi[i];
      }

      //Return the value of the jacobian
      return J_eulerian(s);
    }

    // compute the boundary zeta, which is the initial azimuthal angle
    void setup_zeta_nodal(const bool& use_zeta_2pi_instead_of_0 = false)
    {
      // make sure we've got enough space
      Zeta.resize(nnode(), 0.0);

      double tol = 1e-10;
      for(unsigned j=0; j<nnode(); j++)
      {
	// get the azimthal angle of this node
	double zeta = atan2pi(node_pt(j)->x(1), node_pt(j)->x(0));

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

    /// Pointers to the singular functions
    Vector<UnscaledSingSolnFctPt> Unscaled_singular_fct_pt;

    /// Pointers to gradients of the singular funcions
    Vector<GradientOfUnscaledSingSolnFctPt> Gradient_of_unscaled_singular_fct_pt;

    /// Function pointer to a function which computes dzeta/dx_j
    DzetaDxFctPt Dzeta_dx_fct_pt;
      
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

    // QUEHACERES delete
    /* /// Number of singular functions to subtract */
    /* unsigned Nsingular_fct; */

    /// Dimensionality of this element
    unsigned Dim;

    /// Dimensionality of the whole problem
    unsigned Bulk_dim;
  };
  
  
  /// \short Constructor, takes a pointer to the 2D disk element and the
  /// index of the face to which this element is attached
  // QUEHACERES don't hard code the dimensions
  template <unsigned NNODE_1D>
  ScalableSingularityForNavierStokesLineElement<NNODE_1D>::
    ScalableSingularityForNavierStokesLineElement(FiniteElement* const& disk_el_pt, 
						  const int& face_index,
						  // QUEHACERES delete: const unsigned& nsingular_fct,
						  std::map<Node*,Node*>& existing_duplicate_node_pt,
						  const bool& use_zeta_2pi_instead_of_0) :
    FaceGeometry<TElement<2,NNODE_1D> >(), FaceElement(), // QUEHACERES delete: Nsingular_fct(nsingular_fct),
    Zeta_has_been_setup(false), Dim(1), Bulk_dim(3), Dzeta_dx_fct_pt(0)
  {
    // QUEHACERES delete
    /* // make enough space for the singular function pointers */
    /* Unscaled_singular_fct_pt.resize(Nsingular_fct, 0); */
    /* Gradient_of_unscaled_singular_fct_pt.resize(Nsingular_fct, 0); */
      
    // ------------------------------------------------------------------------
    // The nodes of this element are constructed by duplicating the nodes of a
    // line element attached to the plate. This duplication process is simplified
    // because we want a stand-alone mesh - these new nodes don't need the
    // Navier-Stokes dofs, or any additional dofs (Lagrange multipliers, etc.),
    // they don't need to be added to appropriate boundaries, and they do not
    // depend on external data. 
    
    // Let the bulk element build the FaceElement, i.e. setup the pointers 
    // to its nodes (by referring to the appropriate nodes in the bulk
    // element), etc.
    disk_el_pt->build_face_element(face_index, this);
    
    for(unsigned j=0; j<nnode(); j++)
    {
      Node* original_node_pt = node_pt(j);

      // boundary node which stores the singular amplitude(s)
      Node* singular_amplitude_node_pt = 0;
	  
      // check if we've already duplicated this one
      if(existing_duplicate_node_pt.find(original_node_pt)
	 == existing_duplicate_node_pt.end())	
      {
	// if we haven't make a new one and setup it's values, coordinates,
	// boundary info etc. and add it to the singular mesh
	  
	// get key attributes from the old node
	unsigned n_dim           = original_node_pt->ndim();
	unsigned n_position_type = original_node_pt->nposition_type();

	// QUEHACERES not time dependent for now

	// create a new node - this node stores as many dofs as there are singular
	// functions, but will be resized later
	singular_amplitude_node_pt = new Node(n_dim, n_position_type, 0); // QUEHACERES delete: Nsingular_fct);

	// It has the same coordinates
	for (unsigned i=0; i<n_dim; i++)
	{
	  singular_amplitude_node_pt->x(i) = original_node_pt->x(i);
	}

	// QUEHACERES delete
	/* // initialise the amplitudes to zero */
	/* for (unsigned i=0; i<Nsingular_fct; i++) */
	/* { */
	/*   singular_amplitude_node_pt->set_value(i, 0.0); */
	/* } */
      }
      else
      {
	// use the one we already have
	singular_amplitude_node_pt = existing_duplicate_node_pt[original_node_pt];
      }
      
      // switch over the node
      node_pt(j) = singular_amplitude_node_pt;
      
      // Keep track 
      existing_duplicate_node_pt[original_node_pt] = node_pt(j);      	  
    }
    
    // ------------------------------------------------------------------------
    
    // Now setup the nodal values of zeta
    setup_zeta_nodal(use_zeta_2pi_instead_of_0);
  }


  //////////////////////////////////////////////////////////////////////////////

  /// Add the element's contribution to its residual vector
  template <unsigned NNODE_1D>
    inline void ScalableSingularityForNavierStokesLineElement<NNODE_1D>::
    fill_in_generic_residual_contribution_navier_stokes_sing_fct(Vector<double>& residuals,
								 DenseMatrix<double>& jacobian,
								 const unsigned& flag)
  {    
      if (ndof() == 0)
      {
	return;
      }

#ifdef PARANOID
      // hierher paranoid check null pointers and zero sized vectors    
#endif

      // Number of integration points
      const unsigned n_intpt = this->integral_pt()->nweight();
     
      //Set the Vector to hold local coordinates
      Vector<double> s(Dim);

      //Set up memory for the shape and test functions
      Shape psi(nnode()), test(nnode());
    
      //Loop over the integration points
      //--------------------------------
      for(unsigned ipt=0; ipt<n_intpt; ipt++)
      {
	 //Assign values of s
	for(unsigned i=0; i<Dim; i++) 
	{
	  s[i] = this->integral_pt()->knot(ipt, i);
	}
 
	//Get the integral weight
	double w = this->integral_pt()->weight(ipt);
       
	//Find the shape and test functions and return the Jacobian
	//of the mapping
	double J = this->shape_and_test(s, psi, test);
       
	//Premultiply the weights and the Jacobian
	double W = w*J;

	for(unsigned l=0; l<nnode(); l++)
	{
	  // grab a pointer to the current node
	  Node* nod_pt = node_pt(l);

	  // loop over the nodal values (the singular amplitudes)
	  for(unsigned i=0; i<nod_pt->nvalue(); i++)
	  {
	    // Get eqn number of residual that determines C
	    int local_eqn_c = nodal_local_eqn(l, i);
      
	    if (local_eqn_c >= 0)
	    {
	      // QUEHACERES come back to this, shouldn't get here if
	      // the dofs are pinned
	      throw OomphLibError("QUEHACERES not implemented yet",
				  OOMPH_CURRENT_FUNCTION,
				  OOMPH_EXCEPTION_LOCATION);
	    }
	  } // end loop over singular amplitudes
	}
	
      } // end loop over knots
    }
  
  //////////////////////////////////////////////////////////////////////////////
  //////////////////////////////////////////////////////////////////////////////
  //////////////////////////////////////////////////////////////////////////////
  //////////////////////////////////////////////////////////////////////////////

  //===========================================================================
  // QUEHACERES generic base class
  //===========================================================================
  template <class ELEMENT>
    class NavierStokesWithSingularityFaceElement :
    public virtual FaceGeometry<ELEMENT>, public virtual FaceElement
  {
  public:

    // Extract the number of nodes per side from the templated element
    static const unsigned NNODE_1D = ELEMENT::_NNODE_1D_;

    // shorthand
    typedef ScalableSingularityForNavierStokesLineElement<NNODE_1D> SingularLineElement;
    
    /// \short Constructor, takes the pointer to the "bulk" element and the 
    /// index of the face to which the element is attached. Optional final
    /// arg is the identifier for the additional unknowns multiplier
    NavierStokesWithSingularityFaceElement(FiniteElement* const& bulk_el_pt, 
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
    void set_edge_coordinates_at_knot(const Vector<EdgeCoordinates>& coords)
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
      Edge_coordinates_at_knot.resize(nknot);
	
      // set 'em
      for(unsigned i=0; i<nknot; i++)
      {
	Edge_coordinates_at_knot[i].rho  = coords[i].rho;
	Edge_coordinates_at_knot[i].zeta = coords[i].zeta;
	Edge_coordinates_at_knot[i].phi  = coords[i].phi;
      }     
    }

    // get the edge coordinates at the ith knot
    EdgeCoordinates edge_coordinate_at_knot(const unsigned& i) const
    {
      return Edge_coordinates_at_knot[i];
    }

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

    // get the line element and local coordinate for the ith knot
    std::pair<GeomObject*, Vector<double> >
      line_element_and_local_coordinate_at_knot(const unsigned& i) const
      {
	return Line_element_and_local_coordinate_at_knot[i];
      }

    void edge_coords_and_singular_element_at_knot(
      const unsigned& ipt,
      EdgeCoordinates& edge_coords,
      SingularLineElement*& sing_el_pt,
      Vector<double>& s_singular_el) const
    {
      // get the line element and local coordinate which corresponds to the
      // singular amplitude for this knot
      std::pair<GeomObject*, Vector<double> > line_elem_and_local_coord = 
	this->line_element_and_local_coordinate_at_knot(ipt);

      // QUEHACERES avoid the hard coded template arg here
      // cast the GeomObject to a singular line element      
      sing_el_pt = dynamic_cast<SingularLineElement*>
	(line_elem_and_local_coord.first);

      // check we've actually got one, or we're in trouble
      if(sing_el_pt == 0)
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
      edge_coords = this->edge_coordinate_at_knot(ipt);
    }

    /// \short Specify the value of nodal zeta from the face geometry
    /// The "global" intrinsic coordinate of the element when
    /// viewed as part of a geometric object should be given by
    /// the FaceElement representation, by default (needed to break
    /// indeterminacy if bulk element is SolidElement)
    double zeta_nodal(const unsigned &n, const unsigned &k,           
		      const unsigned &i) const 
    {
      return FaceElement::zeta_nodal(n,k,i);
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

    const unsigned& boundary_id()
    {
      return Boundary_id;
    }

    // QUEHACERES delete
    /* // how many singular functions does this face element know about? */
    /* unsigned nsingular_fct() const */
    /* { */
    /*   return Nsingular_fct; */
    /* } */

    // function to get the pressure at local coordinates s in the
    // face element via the bulk element
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
    
    // function to get the strain rate at local coordinates s in the face element
    // via the bulk element
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

    // function to return the pointer to the bulk element in the bulk region
    // (not the bulk element in the augmented region to which this element is attached)
    ELEMENT*& bulk_region_bulk_element_pt()
    {
      return Bulk_region_bulk_element_pt;
    }
 
  protected:

    /// \short Function to compute the shape and test functions and to return 
    /// the Jacobian of mapping between local and global (Eulerian)
    /// coordinates
    inline double shape_and_test(const Vector<double> &s, Shape &psi, Shape &test)
      const
    {
      //Find number of nodes
      unsigned n_node = nnode();

      //Get the shape functions
      shape(s,psi);

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
    inline double shape_and_test_at_knot(const unsigned &ipt,
					 Shape &psi, Shape &test)
      const
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

    /// The index at which the unknown is stored at the nodes
    unsigned P_index_nst;

    // QUEHACERES delete
    /* /// Number of singular functions to subtract */
    /* unsigned Nsingular_fct; */

    /// Number of spatial dimensions in the problem
    unsigned Dim;

    /// ID of the boundary this face element sits on
    unsigned Boundary_id;
    
  private:

    /// \short Edge coordinates (\rho, \zeta, \phi) of each of this element's
    /// knot points
    Vector<EdgeCoordinates> Edge_coordinates_at_knot;

    /// \short Edge coordinates (\rho, \zeta, \phi) of each of this element's nodes
    Vector<EdgeCoordinates> Nodal_edge_coordinates;
    
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

    /// \short Pointer to the element in the bulk region which also shares this face
    /// (since this->bulk_element_pt() will return an element in the augmented region)
    ELEMENT* Bulk_region_bulk_element_pt;
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
    NavierStokesWithSingularityFaceElement(FiniteElement* const& bulk_el_pt, 
					     const int& face_index, 
					     const unsigned& id) : 
  FaceGeometry<ELEMENT>(), FaceElement(), /* QUEHACERES delete: Nsingular_fct(0), */ Boundary_id(id)
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
	if(ref_el_pt != 0)
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
	NavierStokesEquations<1>* eqn_pt = 
	  dynamic_cast<NavierStokesEquations<1>*>(bulk_el_pt);
	
	//If the cast has failed die
	if(eqn_pt == 0)
	{
	  std::string error_string =
	    "Bulk element must inherit from NavierStokesEquations.";
	  error_string += 
	    "Nodes are one dimensional, but cannot cast the bulk element to\n";
	  error_string += "NavierStokesEquations<1>\n.";
	  error_string += 
	    "If you desire this functionality, you must implement it yourself\n";
       
	  throw OomphLibError(error_string,
			      OOMPH_CURRENT_FUNCTION,
			      OOMPH_EXCEPTION_LOCATION);
	}
	//Otherwise read out the value
	else
	{
	  //Read the index from the (cast) bulk element
	  P_index_nst = eqn_pt->p_nodal_index_nst();
	}
      }
      break;
    
      //Two dimensional problem
      case 2:
      {
	NavierStokesEquations<2>* eqn_pt = 
	  dynamic_cast<NavierStokesEquations<2>*>(bulk_el_pt);
	//If the cast has failed die
	if(eqn_pt == 0)
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
	if(eqn_pt == 0)
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
		     << ". It should be 1,2, or 3!" << std::endl;
     
	throw OomphLibError(error_stream.str(),
			    OOMPH_CURRENT_FUNCTION,
			    OOMPH_EXCEPTION_LOCATION);
	break;
    }
  }


  
  ////////////////////////////////////////////////////////////////////////
  ////////////////////////////////////////////////////////////////////////
  ////////////////////////////////////////////////////////////////////////
    

  // QUEHACERES probably delete this at some point once Z2 stuff working
  //===========================================================================
  /// NavierStokesWithSingularityBoundaryIntegralFaceElement is a class of face elements 
  ///used to compute the contribution to the residuals from the the singular function
  //===========================================================================
  template <class ELEMENT>
    class NavierStokesWithSingularityBoundaryIntegralFaceElement :
    public NavierStokesWithSingularityFaceElement<ELEMENT>
  {
      
  public:

    /// \short Function pointer to the "exact" non-singular function fct(x,u,grad u)
    typedef void (*ExactNonSingularFctPt)(const Vector<double>& x, Vector<double>& u,
					  DenseMatrix<double>& grad_u);
      
    /// \short Constructor, takes the pointer to the "bulk" element and the 
    /// index of the face to which the element is attached.
    NavierStokesWithSingularityBoundaryIntegralFaceElement(FiniteElement* const& bulk_el_pt, 
							   const int& face_index,
							   const int& boundary_id=0);

    ///\short  Broken empty constructor
    NavierStokesWithSingularityBoundaryIntegralFaceElement()
    {
      throw OomphLibError(
	"Don't call empty constructor for NavierStokesWithSingularityBoundaryIntegralFaceElement",
	OOMPH_CURRENT_FUNCTION,
	OOMPH_EXCEPTION_LOCATION);
    }

    /// Broken copy constructor
    NavierStokesWithSingularityBoundaryIntegralFaceElement(
      const NavierStokesWithSingularityBoundaryIntegralFaceElement& dummy) 
    { 
      BrokenCopy::broken_copy("NavierStokesWithSingularityBoundaryIntegralFaceElement");
    } 
      
    /// Broken assignment operator
    void operator=(const NavierStokesWithSingularityBoundaryIntegralFaceElement&) 
      {
	BrokenCopy::broken_assign("NavierStokesWithSingularityBoundaryIntegralFaceElement");
      }

    /// Add the element's contribution to its residual vector
    inline void fill_in_contribution_to_residuals(Vector<double> &residuals)
    {
      //Call the generic residuals function with flag set to 0
      //using a dummy matrix argument
      fill_in_generic_residual_contribution_to_singular_amplitude(
	residuals,GeneralisedElement::Dummy_matrix,0);
    }

#ifndef USE_FD_JACOBIAN
    /// \short Add the element's contribution to its residual vector and its
    /// Jacobian matrix
    inline void fill_in_contribution_to_jacobian(Vector<double> &residuals,
						 DenseMatrix<double> &jacobian)
    {
      //Call the generic routine with the flag set to 1
      fill_in_generic_residual_contribution_to_singular_amplitude(residuals,jacobian,1);
    }
#endif
    
    /// Output function -- forward to broken version in FiniteElement
    /// until somebody decides what exactly they want to plot here...
    void output(std::ostream& outfile)
    {
      FiniteElement::output(outfile);
    }

    /// Output with various contributions
    void output(std::ostream &outfile, const unsigned &nplot)
    {
      // shorthand
      unsigned Dim = this->Dim;
      
      //Vector of local coordinates
      Vector<double> s(Dim-1);
   
      // Tecplot header info
      outfile << this->tecplot_zone_string(nplot);
   
      // Loop over plot points
      unsigned num_plot_points = this->nplot_points(nplot);
      for (unsigned iplot=0; iplot<num_plot_points; iplot++)
      {
	// Get local coordinates of plot point
	this->get_s_plot(iplot,nplot,s);
     
	Vector<double> x(Dim);
	for(unsigned i=0; i<Dim; i++) 
	{
	  x[i]=this->interpolated_x(s,i);
	  outfile << x[i] << " ";
	}

	// Compute outer unit normal at the specified local coordinate
	Vector<double> unit_normal(Dim);
	this->outer_unit_normal(s, unit_normal);

	outfile << unit_normal[0] << " " << unit_normal[1] << " " << endl;

	// QUEHACERES output something useful here
      }
   
      // Write tecplot footer (e.g. FE connectivity lists)
      this->write_tecplot_zone_footer(outfile,nplot);
   
    }

    /// \short Compute this element's contribution to the integral that determines C
    double get_contribution_integral();
 
    /// Pointer to exact non singular fct (and its gradient) only used
    /// to validate the computation of the integral. Ignored if null 
    /// which is the default
    // hierher make read only and use set/unset fct to enable/disable;
    // currently we'd have to reset this to null!
    ExactNonSingularFctPt& exact_non_singular_fct_pt()
    {
      return Exact_non_singular_fct_pt;
    }

  private:

    /// \short Add the element's contribution to its residual vector.
    /// flag=1(or 0): do (or don't) compute the contribution to the
    /// Jacobian as well. 
    void fill_in_generic_residual_contribution_to_singular_amplitude(
      Vector<double> &residuals, DenseMatrix<double> &jacobian, 
      const unsigned& flag);

    /// Pointer to exact non singular fct (and its gradient) only used
    /// to validate the computation of the integral. Ignored if null 
    /// which is the default
    ExactNonSingularFctPt Exact_non_singular_fct_pt;
  }; 

  ////////////////////////////////////////////////////////////////////////////
  ////////////////////////////////////////////////////////////////////////////
  ////////////////////////////////////////////////////////////////////////////

  //===========================================================================
  /// Constructor, takes the pointer to the "bulk" element, the 
  /// index of the fixed local coordinate and its value represented
  /// by an integer (+/- 1), indicating that the face is located
  /// at the max. or min. value of the "fixed" local coordinate
  /// in the bulk element.
  //===========================================================================
  template <class ELEMENT>
    NavierStokesWithSingularityBoundaryIntegralFaceElement<ELEMENT>::
    NavierStokesWithSingularityBoundaryIntegralFaceElement(FiniteElement* const& bulk_el_pt, 
							   const int& face_index,
							   const int& boundary_id) : 
  NavierStokesWithSingularityFaceElement<ELEMENT>(bulk_el_pt, face_index, boundary_id),
    Exact_non_singular_fct_pt(0) { }

  ////////////////////////////////////////////////////////////////////////////
  ////////////////////////////////////////////////////////////////////////////
  ////////////////////////////////////////////////////////////////////////////

  //===========================================================================
  /// Compute the element's residual vector and the (zero) Jacobian matrix.
  //===========================================================================
  template <class ELEMENT>
    void NavierStokesWithSingularityBoundaryIntegralFaceElement<ELEMENT>::
    fill_in_generic_residual_contribution_to_singular_amplitude(
      Vector<double> &residuals, DenseMatrix<double> &jacobian, 
      const unsigned& flag)
  {

    // hierher populate when contribution is split up
    oomph_info << "This shouldn't be called at the moment\n";
    abort();
  }

  ////////////////////////////////////////////////////////////////////////////
  ////////////////////////////////////////////////////////////////////////////
  ////////////////////////////////////////////////////////////////////////////

  
  //===========================================================================
  /// Calculate the contribution of the face element to the integral that
  /// determines the amplitude, via the Lorentz reciprocity theorem. 
  //===========================================================================
  template <class ELEMENT>
    double NavierStokesWithSingularityBoundaryIntegralFaceElement<ELEMENT>::
    get_contribution_integral()
  {
    // shorthand
    const unsigned Dim = this->Dim;
    
    //Find out how many nodes there are
    const unsigned n_node = this->nnode();
  
    //Set up memory for the shape and test functions
    Shape psif(n_node), testf(n_node);
 
    //Set the value of Nintpt
    const unsigned n_intpt = this->integral_pt()->nweight();
 
    // Set the Vector to hold local coordinates (in this face element, not the
    // bulk element this is attached to)
    Vector<double> s(Dim-1);
 
    // Saves result of integration
    double integral_result = 0.0;

    //Loop over the integration points
    //--------------------------------
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
      double J = this->shape_and_test(s, psif, testf);
   
      //Premultiply the weights and the Jacobian
      double W = w*J;

      // compute outer normal unit vector
      Vector<double> unit_normal(Dim);
      this->outer_unit_normal(s, unit_normal);

      // local coordinates in bulk element this face element is attached to
      Vector<double> s_bulk(Dim);

      // global coordinates
      Vector<double> x(Dim); 

      // get global coordinates
      for(unsigned i=0; i<Dim; i++)
      { 
	x[i] = this->interpolated_x(s,i); 
      } 
   
      // Get gradient of scaled/unscaled singular velocity functions
      DenseMatrix<double> dudx_sing(Dim, Dim);
      DenseMatrix<double> dudx_sing_unscaled(Dim, Dim);

      // ### UPDATE! don't delete entirely
      /* dudx_sing          = Navier_stokes_sing_el_pt->gradient_of_singular_fct(x); */
      /* dudx_sing_unscaled = Navier_stokes_sing_el_pt->gradient_of_unscaled_singular_fct(x); */
      
      // Get the values of the singular functions at our current location
      Vector<double> u_sing(Dim+1);
      Vector<double> u_sing_unscaled(Dim+1);

      // ### UPDATE! don't delete entirely
      /* u_sing          = Navier_stokes_sing_el_pt->singular_fct(x); */
      /* u_sing_unscaled = Navier_stokes_sing_el_pt->unscaled_singular_fct(x); */

      // get singular pressure
      double p_sing          = u_sing[this->P_index_nst];
      double p_sing_unscaled = u_sing_unscaled[this->P_index_nst];
      
      // shorthand
      ELEMENT* bulk_el_pt = dynamic_cast<ELEMENT*>(this->bulk_element_pt());

      // compute the singular contribution to the strain-rate
      DenseMatrix<double>strain_rate_sing(Dim, Dim, 0.0);
      DenseMatrix<double>strain_rate_sing_unscaled(Dim, Dim, 0.0);
      
      for (unsigned i=0; i<Dim; i++)
      {
	for(unsigned j=0; j<Dim; j++)
	{
	  strain_rate_sing(i,j) = 0.5*(dudx_sing(i,j) + dudx_sing(j,i));
	  
	  strain_rate_sing_unscaled(i,j) =
	    0.5*(dudx_sing_unscaled(i,j) + dudx_sing_unscaled(j,i));
	}
      }
	
      // get contribution of singular pressure and singular velocity gradients to stress tensor
      DenseMatrix<double> stress_sing(Dim, Dim);
      DenseMatrix<double> stress_sing_unscaled(Dim, Dim);
      stress_sing = (*bulk_el_pt->stress_fct_pt())(strain_rate_sing, p_sing);
      stress_sing_unscaled = (*bulk_el_pt->stress_fct_pt())(strain_rate_sing_unscaled, p_sing_unscaled);
      
      // Get the local bulk coordinates    
      s_bulk = this->local_coordinate_in_bulk(s);
      
      Vector<double> u_fe(Dim);
      
      // Get FE part of the velocity
      for(unsigned i=0; i<Dim; i++)
      {
	u_fe[i] = bulk_el_pt->interpolated_u_nst(s_bulk, i);
      }
      
      // get FE part of pressure
      double p_fe = bulk_el_pt->interpolated_p_nst(s_bulk);

      // get FE part of velocity gradient tensor
      DenseMatrix<double> dudx_fe(Dim, Dim);
      
      for(unsigned i=0; i<Dim; i++)
      {
	for(unsigned j=0; j<Dim; j++)
	{
	  // get derivative du_i/dx_j
	  dudx_fe(i,j) = bulk_el_pt->interpolated_dudx_nst(s_bulk, i, j);
	}
      }
      
      // get the FE strain rate 1/2(du_i/dx_j + du_j/dx_i)
      DenseMatrix<double> strain_rate_fe(Dim, Dim, 0.0);
      
      bulk_el_pt->strain_rate(s_bulk, strain_rate_fe);
      
      // FE part of the stress
      DenseMatrix<double> stress_fe(Dim, Dim, 0.0);

      // compute it from consitutive equation
      stress_fe = (*bulk_el_pt->stress_fct_pt())(strain_rate_fe, p_fe);
            
      // get FE part of the traction
      Vector<double> traction_fe(Dim, 0.0);

      // get FE traction from bulk element
      bulk_el_pt->get_traction(s_bulk, unit_normal, traction_fe);

      // ================================================
      // Now we compute the contibution of the integral

      for(unsigned i=0; i<Dim; i++)
      {
	for(unsigned j=0; j<Dim; j++)
	{
	  // Lorentz reciprocity theorem
	  integral_result += W * unit_normal[j] * (stress_fe(i,j) * u_sing_unscaled[i]
						   - stress_sing_unscaled(i,j) * u_fe[i] );
	}
      }      
    } // end loop over integration points
        
    return integral_result;
  }

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
      NavierStokesWithSingularityBCFaceElement(FiniteElement* const &bulk_el_pt, 
					       const int& face_index,
					       const unsigned &id = 0); 
  
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

      // ### <- but keep around until we've sorted out the external data stuff
      
      /* /// \short Set pointer to element that stores singular fct. Data that stores */
      /* /// the amplitude of the singular fct and its index is retrieved from */
      /* /// that element so the Data can be used as external Data in this */
      /* /// element. */
      /* void set_navier_stokes_sing_el_pt( */
      /* 	Vector<ScalableSingularityForNavierStokesElement<ELEMENT>*> navier_stokes_sing_el_pt)  */
      /* { */
      /* 	// shorthand, get number of singlar functions */
      /* 	unsigned nsing = navier_stokes_sing_el_pt.size(); */

      /* 	// set the number of singular functions */
      /* 	this->Nsingular_fct = nsing; */
	
      /* 	// make sure we've got enough space */
      /* 	Navier_stokes_sing_el_pt.resize(nsing); */
      /* 	C_external_data_index.resize(nsing); */
      /* 	C_external_data_value_index.resize(nsing); */

      /* 	// loop over the singular functions and add their amplitudes as external data */
      /* 	for(unsigned ising=0; ising<nsing; ising++) */
      /* 	{ */
      /* 	  Navier_stokes_sing_el_pt[ising] = navier_stokes_sing_el_pt[ising]; */
      /* 	  C_external_data_index[ising] = this->add_external_data( */
      /* 	    navier_stokes_sing_el_pt[ising]->data_that_stores_amplitude_of_singular_fct()); */
      /* 	  C_external_data_value_index[ising] = */
      /* 	    navier_stokes_sing_el_pt[ising]->index_of_value_that_stores_amplitude_of_singular_fct(); */
      /* 	} */
      /* } */

      /// Add the element's contribution to its residual vector
      inline void fill_in_contribution_to_residuals(Vector<double> &residuals)
      {
	//Call the generic residuals function with flag set to 0
	//using a dummy matrix argument
	fill_in_generic_residual_contribution_navier_stokes_bc(
	  residuals, GeneralisedElement::Dummy_matrix, 0);
      }

/* #ifndef USE_FD_JACOBIAN */
      /// \short Add the element's contribution to its residual vector and its
      /// Jacobian matrix
      inline void fill_in_contribution_to_jacobian(Vector<double> &residuals,
      						   DenseMatrix<double> &jacobian)
      {
      	//Call the generic routine with the flag set to 1
      	fill_in_generic_residual_contribution_navier_stokes_bc(residuals, jacobian, 1);
      }
/* #endif   */

      /// Output function
      void output(std::ostream &outfile)
      {
	const unsigned n_plot = 5;
	output(outfile, n_plot);
      }

      /// \short Output function
      void output(std::ostream &outfile, const unsigned &nplot)
      {
	//oomph_info << "hierher need to update output fct" << std::endl;
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
	  outfile << endl;
	}
	return;
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

      /// Pin Lagrange multiplier associated with ith coordinate at specified local node
      void pin_lagrange_multiplier_at_specified_local_node(const unsigned& j,
							   const unsigned& i,
							   const int& id = -1)
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
	  lambda_index = map_l.begin()->second;
	}	
	else
	{
	  // otherwise, get the nodal index for the specified boundary ID
	  lambda_index = map_l[id];
	}
	this->node_pt(j)->pin(lambda_index+i);
      }

      /// Unpin ith component of FE part of the solution at specified local node
      void unpin_u_fe_at_specified_local_node(const unsigned& j, const unsigned& i)
      {   
	this->node_pt(j)->unpin(i);	  	
      }

      // QUEHACERES for debug, output the value of the Lagrange multipliers on the Dirichlet
      // boundaries and the associated stress
      void output_lagrange_multiplier_and_stress(std::ostream& outfile,
						 const unsigned& nplot )
      {
	ostringstream error_message;
	  error_message << "Don't use this at the moment, broken since we switched over to "
			<< "line singularity - this element contains a vector of singular "
			<< "elements and local coordinates at *knots*, but for this function"
			<< "we need those elements and coordinates at the plot points";
	  
	throw OomphLibError(error_message.str(),
	  OOMPH_CURRENT_FUNCTION,
	      OOMPH_EXCEPTION_LOCATION);
	
	/* // shorthand */
	/* unsigned Dim = this->Dim; */
	
	/* // pointer to the bulk element this face element is attached to */
	/* ELEMENT* bulk_el_pt = dynamic_cast<ELEMENT*>(this->bulk_element_pt()); */
      
	/* // get the number of nodes in this face element */
	/* unsigned nnode = this->nnode(); */
    
	/* //Set up memory for the shape and test functions */
	/* Shape psi(nnode), test(nnode); */
    
	/* unsigned num_plot_points = this->nplot_points(nplot); */
	/* for (unsigned iplot=0; iplot < num_plot_points; iplot++) */
	/* { */
	/*   // Get local coordinates of this plot point */
	/*   Vector<double> s(this->dim());       */
	/*   this->get_s_plot(iplot, nplot, s); */
      
	/*   //Find the shape and test functions and return the Jacobian */
	/*   //of the mapping */
	/*   double J = this->shape_and_test(s, psi, test); */
         
	/*   // local coordinates in bulk element this face element is attached to */
	/*   Vector<double> s_bulk(this->Dim); */

	/*   // Get the local bulk coordinates     */
	/*   s_bulk = this->local_coordinate_in_bulk(s); */

	/*   Vector<double> x(Dim); */

	/*   for(unsigned i=0; i<Dim; i++) */
	/*   { */
	/*     x[i] = this->interpolated_x(s,i); */
	/*   } */
	  
	/*   // Lagrange multipliers */
	/*   Vector<double> lambda(Dim); */
      
	/*   unsigned nnode = this->nnode(); */
	/*   for(unsigned j=0; j<nnode; j++) */
	/*   { */
	/*     // get the map which gives the starting nodal index for */
	/*     // the Lagrange multipliers associated with each boundary ID */
	/*     std::map<unsigned, unsigned> first_index = *( */
	/*       dynamic_cast<BoundaryNodeBase*>(this->node_pt(j))-> */
	/*       index_of_first_value_assigned_by_face_element_pt() ); */
	
	/*     for(unsigned i=0; i<Dim; i++) */
	/*     { */
	/*       unsigned lambda_index = first_index[this->Boundary_id] + i; */
	  
	/*       // get the Lagrange multiplier */
	/*       lambda[i] += this->nodal_value(j, lambda_index) * psi[j]; */
	/*     } */
	/*   } */

	/*   // get FE part of pressure */
	/*   double p_fe = bulk_el_pt->interpolated_p_nst(s_bulk); */

	/*   // total pressure */
	/*   double p_total = p_fe; */

	/*   // get the singular parts of the pressure */
	/*   Vector<double> p_sing(this->Nsingular_fct, 0.0); */
	  
	/*   for(unsigned ising=0; ising<this->Nsingular_fct; ising++) */
	/*   { */
	/*     // Get the values of the scaled singular functions at our current location */
	/*     Vector<double> u_sing(Dim+1); */
	/*     u_sing = Navier_stokes_sing_el_pt[ising]->singular_fct(x); */
	    
	/*     // get scaled singular pressure */
	/*     p_sing[ising] = u_sing[this->P_index_nst]; */

	/*     // add to the total pressure */
	/*     p_total += p_sing[ising]; */
	/*   } */
	  
	/*   // get the FE strain rate 1/2(du_i/dx_j + du_j/dx_i) */
	/*   DenseMatrix<double> strain_rate_fe(Dim, Dim);       */
	/*   bulk_el_pt->strain_rate(s_bulk, strain_rate_fe); */

	/*   // get the stress from the constitutive law */
	/*   DenseMatrix<double> stress_fe(Dim, Dim); */
	/*   stress_fe = (*bulk_el_pt->stress_fct_pt())(strain_rate_fe, p_fe); */
	  	  

	/*   // ==================================================== */
	/*   // output stuff */
	/*   // ==================================================== */
	  
	/*   // output the coordinates of this point */
	/*   for(unsigned i=0; i<Dim; i++) */
	/*   {	     */
	/*     outfile << x[i] << " "; */
	/*   } */
	  
	/*   // output the Lagrange multipliers at this plot point */
	/*   for(unsigned i=0; i<Dim; i++) */
	/*   { */
	/*     outfile << lambda[i] << " "; */
	/*   } */

	/*   // output total pressure */
	/*   outfile << p_total << " "; */
	  
	/*   // output the traction at this point */
	/*   for(unsigned i=0; i<Dim; i++)	 */
	/*   { */
	/*     for(unsigned j=0; j<Dim; j++)	 */
	/*     { */
	/*       outfile << stress_fe(i,j) << " "; */
	/*     } */
	/*   } */

	/*   for(unsigned ising=0; ising<this->Nsingular_fct; ising++) */
	/*   { */
	/*     // Get gradient of scaled singular velocity functions */
	/*     DenseMatrix<double> dudx_sing(Dim, Dim); */
	/*     dudx_sing = Navier_stokes_sing_el_pt[ising]->gradient_of_singular_fct(x); */

	/*     // compute the unscaled singular contribution to the strain-rate */
	/*     DenseMatrix<double>strain_rate_sing(Dim, Dim, 0.0); */
      
	/*     for (unsigned i=0; i<Dim; i++) */
	/*     { */
	/*       for(unsigned j=0; j<Dim; j++) */
	/*       { */
	/* 	strain_rate_sing(i,j) = 0.5*(dudx_sing(i,j) + dudx_sing(j,i)); */
	/*       } */
	/*     } */

	/*     // get contribution of singular pressure and singular velocity gradients to stress tensor */
	/*     DenseMatrix<double> stress_sing(Dim, Dim); */
	/*     stress_sing = (*bulk_el_pt->stress_fct_pt())(strain_rate_sing, p_sing[ising]); */
	  
	/*     // output the traction at this point */
	/*     for(unsigned i=0; i<Dim; i++)	 */
	/*     { */
	/*       for(unsigned j=0; j<Dim; j++)	 */
	/*       { */
	/* 	outfile << stress_sing(i,j) << " "; */
	/*       } */
	/*     } */
	/*   } */
	/*     outfile << std::endl;	   */
	/* } */
      }

      // function to get the total normal stress acting on this element
      // (useful for convergence studies of total force, etc.)
      double get_contribution_to_normal_stress() const;

      // add a nodal index to the set of torus boundary nodes
      void add_nodal_index_on_torus(const unsigned& i)
      {
	Index_of_nodes_on_torus.insert(i);
      }
      
    private:

      /// \short Add the element's contribution to its residual vector.
      /// flag=1(or 0): do (or don't) compute the contribution to the
      /// Jacobian as well. 
      void fill_in_generic_residual_contribution_navier_stokes_bc(
	Vector<double> &residuals, DenseMatrix<double> &jacobian, 
	const unsigned& flag);

      /// Desired boundary values at nodes
      DenseMatrix<double> Nodal_boundary_value;
      
      /// \short Indices of external Data that store the values of the amplitudes of
      /// the singular functions
      Vector<unsigned> C_external_data_index;
  
      /// \short Indices of values (within external Data) that store the
      /// values of the amplitudes of the singular functions
      Vector<unsigned> C_external_data_value_index;

      /// \short set of nodal indices for nodes on the torus boundary.
      /// These nodes are shared by the stress-jump elements which impose
      /// their own boundary conditions, but only one set of BCs can be
      /// applied, so these BC elements will not contribute to the residuals
      /// of any nodes whose indices are in this set.
      std::set<unsigned> Index_of_nodes_on_torus;
	
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
    NavierStokesWithSingularityBCFaceElement(FiniteElement* const& bulk_el_pt, 
					     const int& face_index, 
					     const unsigned& id) : 
  NavierStokesWithSingularityFaceElement<ELEMENT>(bulk_el_pt, face_index, id)
  { 
    unsigned n_node = this->nnode();

    // Make space for Dim Lagrange multipliers
    Vector<unsigned> n_additional_values(n_node, this->Dim);

    // add them (this takes care of duplicate IDs so no additional
    // checks needed here to avoid redundant nodal values)
    this->add_additional_values(n_additional_values, this->Boundary_id);
    
  } // end NavierStokesWithSingularityBCFaceElement constructor


  template <class ELEMENT>
  double NavierStokesWithSingularityBCFaceElement<ELEMENT>::
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
    double normal_stress_integral = 0.0;
    
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

      // QUEHACERES avoid the hard coded template arg here
      // cast the GeomObject to a singular line element      
      SingularLineElement* sing_el_pt =
	dynamic_cast <SingularLineElement*>( line_elem_and_local_coord.first );

      // local coordinate in the singular element for the zeta of this knot
      Vector<double> s_singular_el = line_elem_and_local_coord.second;
      
      // get the \rho,\zeta,\phi coordinates at this knot
      EdgeCoordinates edge_coords_at_knot = this->edge_coordinate_at_knot(ipt);

      // the sum of all scaled singular functions
      Vector<double> u_sing_total(Dim+1, 0.0);

      // total singular contributions to velocity gradient tensor and strain rate
      DenseMatrix<double> dudx_sing_total(Dim, Dim, 0.0);
      DenseMatrix<double> strain_rate_sing_total(Dim, Dim, 0.0);
      
      if(sing_el_pt != 0)
      {
	u_sing_total = sing_el_pt->total_singular_contribution(edge_coords_at_knot,
							       s_singular_el);
      
	dudx_sing_total = sing_el_pt->
	  total_singular_gradient_contribution(edge_coords_at_knot,
					       s_singular_el);
      }
      
      // total singular pressure
      double p_sing_total = u_sing_total[this->P_index_nst];

      // now compute the singular strain rate
      for(unsigned i=0; i<Dim; i++)
      {
	for(unsigned j=0; j<Dim; j++)
	{
	  strain_rate_sing_total(i,j) = 0.5 * (dudx_sing_total(i,j) +
					       dudx_sing_total(j,i));
	}
      }

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
	  normal_stress_integral += W * unit_normal[j] * stress_total(i,j);
	}
      }
    } // end loop over integration points

    return normal_stress_integral;    
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
      Vector<double> &residuals, DenseMatrix<double> &jacobian, 
      const unsigned& flag) 
  {
    // shorthands
    const unsigned Dim = this->Dim;
    const unsigned Boundary_id = this->Boundary_id;
    
    //Find out how many nodes there are
    const unsigned n_node = this->nnode();
     
    //Set up memory for the shape and test functions
    Shape psi(n_node), test(n_node);
     
    //Set the value of Nintpt
    const unsigned n_intpt = this->integral_pt()->nweight();
     
    //Set the Vector to hold local coordinates
    Vector<double> s(Dim-1);

    //Loop over the integration points
    //--------------------------------
    for(unsigned ipt=0; ipt<n_intpt; ipt++)
    {       
      //Assign values of s
      for(unsigned i=0; i<(Dim-1); i++) 
      {
	s[i] = this->integral_pt()->knot(ipt, i);
      }
       
      //Get the integral weight
      double w = this->integral_pt()->weight(ipt);
       
      //Find the shape and test functions and return the Jacobian
      //of the mapping
      double J = this->shape_and_test(s, psi, test);
       
      //Premultiply the weights and the Jacobian
      double W = w*J;
       
      //Calculate stuff at integration point
      Vector<double> u_fe(Dim, 0.0);
      Vector<double> u_bc(Dim, 0.0);	
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
	  // get the nodal index, accounting for the dimension offset
	  unsigned lambda_index = first_index[Boundary_id] + i;
	    
	  // get the nodal values of the FE solution and the boundary conditions
	  u_fe[i] += this->nodal_value(l,i)    * psi[l];
	  u_bc[i] += Nodal_boundary_value(l,i) * psi[l];

	  // get the interpolated Lagrange multipliers
	  lambda[i] += this->nodal_value(l, lambda_index) * psi[l];
	}
	  
      }

      // cast the GeomObject to a singular line element      
      SingularLineElement* sing_el_pt = 0;
      
      EdgeCoordinates edge_coords_at_knot;
      Vector<double> s_singular_el(1, 0.0);
      
      // get the edge coordinates at this knot, the corresponding singular
      // line element and it's local coordinates
      this->edge_coords_and_singular_element_at_knot(ipt,						       
						     edge_coords_at_knot,
						     sing_el_pt,
						     s_singular_el);
      
      // QUEHACERES delete
      /* // get the line element and local coordinate which corresponds to the */
      /* // singular amplitude for this knot */
      /* std::pair<GeomObject*, Vector<double> > line_elem_and_local_coord =  */
      /* 	this->line_element_and_local_coordinate_at_knot(ipt); */

      /* // cast the GeomObject to a singular line element */
      /* ScalableSingularityForNavierStokesLineElement<3>* sing_el_pt = */
      /* 	dynamic_cast<ScalableSingularityForNavierStokesLineElement<3>*> */
      /* 	(line_elem_and_local_coord.first); */

      /* // local coordinate in the singular element for the zeta of this knot */
      /* Vector<double> s_singular_el = line_elem_and_local_coord.second; */
      
      /* // get the \rho,\zeta,\phi coordinates at this knot */
      /* EdgeCoordinates edge_coords_at_knot = this->edge_coordinate_at_knot(ipt); */

       // Stuff related to singular fct
      Vector<double> u_sing(Dim+1, 0.0);

      // check we've actually got a singular line element pointer	
      if (sing_el_pt == 0)
      { 
	throw OomphLibError(
	"This NavierStokesWithSingularityBCFaceElement doesn't have a singular line element pointer",
	OOMPH_CURRENT_FUNCTION,
	OOMPH_EXCEPTION_LOCATION);
      }
      
      Vector<Vector<double> > u_sing_unscaled(sing_el_pt->nsingular_fct(), Vector<double>(1,0));
      
      // get a list of the singular function IDs this element has
      Vector<unsigned> sing_ids = sing_el_pt->singular_fct_ids();
      
      u_sing = sing_el_pt->total_singular_contribution(edge_coords_at_knot,
						       s_singular_el);
      // loop over the singular functions
      for(unsigned ising=0; ising < sing_el_pt->nsingular_fct(); ising++)
      {
	// get the ID
	unsigned sing_fct_id = sing_ids[ising];

	// and get the unscaled singular function associated with this ID
	u_sing_unscaled[ising] = sing_el_pt->unscaled_singular_fct(edge_coords_at_knot,
								   sing_fct_id);	 
      }
      

      //Now add to the appropriate equations

      // number of local equation which determines the singular amplitude
      int local_eqn_c = -1;

      // QUEHACERES think about how to loop this properly when we come to implementing the
      // analytic jacobian
      /* if (Navier_stokes_sing_el_pt != 0) */
      /* { */
      /* 	local_eqn_c = external_local_eqn(C_external_data_index, */
      /* 					 C_external_data_value_index); */
      /* } */

      //Loop over the test functions
      for(unsigned l=0; l<n_node; l++)
      {
	// grab a pointer to the current node
	Node* node_pt = this->node_pt(l);
	
	// QUEHACERES experimental:
	// if this node is on a torus boundary, all contributions to it's residuals
	// (velocity, pressure and Lagrange multipliers) are dealt with by the
	// stress jump elements
	if(Index_of_nodes_on_torus.find(l) != Index_of_nodes_on_torus.end())
	{
	  // QUEHACERES debug
	  double r = sqrt(pow(node_pt->x(0),2) + pow(node_pt->x(1),2));
	  
	  continue;
	}
	  
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
	  unsigned lambda_index = first_index[Boundary_id] + d;

	  // get the local Lagrange multiplier equation number 
	  int local_eqn_lagr = this->nodal_local_eqn(l, lambda_index);
	      
	  // QUEHACERES get this nodal index systematically, don't assume it starts at 0
	  int local_eqn_u_fe = this->nodal_local_eqn(l, d);

	  // QUEHACERES debug
	  int global_eqn = node_pt->eqn_number(lambda_index);
	    
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
	    residuals[local_eqn_lagr] += (u_fe[d] + u_sing[d] - u_bc[d]) * test[l]*W;
	      
	    // Jacobian?
	    if (flag == 1)
	    {
	      for(unsigned l2=0; l2<n_node; l2++)
	      {
		// QUEHACERES again, get this index more systematically
		int local_unknown_u_fe = this->nodal_local_eqn(l2, d);
		if (local_unknown_u_fe >= 0)
		{
		  jacobian(local_eqn_lagr, local_unknown_u_fe) += psi[l2] * test[l]*W;
		}
	      }

	      // QUEHACERES come back to this when we have a c equation
	      /* // Deriv. w.r.t. amplitude is simply the unscaled singular fct. */
	      /* for(unsigned ising=0; ising<sing_el_pt->nsingular_fct(); ising++) */
	      /* { */
	      /* 	int local_eqn_c = this->external_local_eqn(C_external_data_index[ising], */
	      /* 						   C_external_data_value_index[ising]); */
	      /* 	if (local_eqn_c >= 0) */
	      /* 	{ */
		
	      /* 	  jacobian(local_eqn_lagr, local_eqn_c) += */
	      /* 	    u_sing_unscaled[ising][d] * test[l]*W; */
	      /* 	}		 */
	      /* } */
	    }
	  }
         
	  // Contribution of Lagrange multiplier to bulk eqn:
	  if (local_eqn_u_fe >= 0)
	  {
	    residuals[local_eqn_u_fe] += lambda[d] * test[l] * W;

	    // QUEHACERES need to review this code, never been tested
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
		unsigned lambda_index2 = first_index2[Boundary_id] + d;
		int local_unknown_lambda = this->nodal_local_eqn(l2, lambda_index2);
		      
		if (local_unknown_lambda >= 0)
		{
		  jacobian(local_eqn_u_fe, local_unknown_lambda) += psi[l2] * test[l] * W;
		}
		    
	      }
	    }
	  }	    
	} // end loop over directions

      } // end loop over nodes
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
      FiniteElement* const& bulk_el_pt, 
      const int& face_index,   
      std::map<Node*,Node*>& existing_duplicate_node_pt,
      const unsigned &id = 0); 

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
    inline void fill_in_contribution_to_residuals(Vector<double> &residuals)
    {
      //Call the generic residuals function with flag set to 0
      //using a dummy matrix argument
      fill_in_generic_residual_contribution_nst_sing_jump(
	residuals,GeneralisedElement::Dummy_matrix,0);
    }

    // QUEHACERES analytic jacobian tested and working without r_c, possibly
    // need to re-enable FD jacobian to test that out.
    // also re-enable once all these different elements have tested analytic
    // jacobians so we can legit have USE_FD_JACOBIAN switched on/off
    /* #ifndef USE_FD_JACOBIAN */
      
    /// \short Add the element's contribution to its residual vector and its
    /// Jacobian matrix
    inline void fill_in_contribution_to_jacobian(Vector<double> &residuals,
						 DenseMatrix<double> &jacobian)
    {
      //Call the generic routine with the flag set to 1
      fill_in_generic_residual_contribution_nst_sing_jump
	(residuals,jacobian,1);
    }

    // QUEHACERES
    /* #endif  */
      
    /// Output function
    void output(std::ostream &outfile)
    {
      const unsigned n_plot=5;
      output(outfile,n_plot);
    }

    /// \short Output function
    void output(std::ostream &outfile, const unsigned &nplot)
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
	    // get the nodal index, accounting for the dimension offset
	    unsigned lambda_index = first_index[this->Boundary_id] + i;
	  
	    u_augmented[i]  += this->nodal_value(l,i) * psi[l];
	    u_bulk[i] += Orig_node_pt[l]->value(i) * psi[l];
	    lambda[i]  += this->nodal_value(l,lambda_index) * psi[l];
	    
	    x[i] += this->nodal_position(l,i) * psi[l];
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

    void output_traction_and_lagrange_multipliers(std::ofstream& outfile)
    {
      // shorthand
      const unsigned Dim = this->Dim;
	  
      // Dimension of element 
      const unsigned el_dim = this->dim();
    
      //Vector of local coordinates
      Vector<double> s(el_dim);

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
   
	// Get the integral weight
	double w = this->integral_pt()->weight(ipt);
   
	//Find the shape and test functions and return the Jacobian
	//of the mapping
	double J = this->shape_and_test(s, psi, test);
   
	//Premultiply the weights and the Jacobian
	double W = w*J;
	
	//Calculate stuff at integration point
	Vector<double> x(Dim, 0.0);
	Vector<double> u_augmented(Dim, 0.0);
	Vector<double> u_bulk(Dim, 0.0);
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
	    // get the nodal index, accounting for the dimension offset
	    unsigned lambda_index = first_index[this->Boundary_id] + i;

	    lambda[i]  += this->nodal_value(l,lambda_index) * psi[l];
	    
	    x[i] += this->nodal_position(l,i) * psi[l];
	  }
	}

	// get the traction onto the augmented region	
	bulk_el_pt->get_traction(s_bulk, unit_normal, traction_augmented);

	// now get the traction in the bulk
	// ----------------------------------

	// first, find the local coordinates of this point in the
	// bulk element
	Vector<double> s_bulk_region(Dim, 0.0);
	
	GeomObject* geom_object_pt = 0x0;
	ELEMENT* bulk_region_bulk_elem_pt = this->bulk_region_bulk_element_pt();
	
	bulk_region_bulk_elem_pt->locate_zeta(x, geom_object_pt, s_bulk_region);

	// now get the traction (passing in the same unit normal as before,
	// so will need to flip the sign of the traction)
	this->bulk_region_bulk_element_pt()->get_traction(s_bulk_region,
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
    // ### <- but keep around until we've sorted out the external data stuff
      
    /* /// \short Set pointer to element that stores singular fcts. Data that stores */
    /* /// the amplitude of the singular fct and its index is retrieved from */
    /* /// that element so the Data can be used as external Data in this */
    /* /// element. */
    /* void set_navier_stokes_sing_el_pt( */
    /* 	Vector<ScalableSingularityForNavierStokesElement<ELEMENT>*> navier_stokes_sing_el_pt)  */
    /* { */
    /* 	// set the number of singular functions */
    /* 	this->Nsingular_fct = navier_stokes_sing_el_pt.size(); */

    /* 	// make sure we've got enough space */
    /* 	Navier_stokes_sing_el_pt.resize(this->Nsingular_fct); */
    /* 	C_external_data_index.resize(this->Nsingular_fct); */
    /* 	C_external_data_value_index.resize(this->Nsingular_fct); */

    /* 	// loop over the singular functions and add their amplitudes as external data */
    /* 	for(unsigned ising=0; ising<this->Nsingular_fct; ising++) */
    /* 	{ */
    /* 	  Navier_stokes_sing_el_pt[ising] = navier_stokes_sing_el_pt[ising]; */
    /* 	  C_external_data_index[ising] = this->add_external_data( */
    /* 	    navier_stokes_sing_el_pt[ising]->data_that_stores_amplitude_of_singular_fct()); */
    /* 	  C_external_data_value_index[ising] = navier_stokes_sing_el_pt[ising]-> */
    /* 	    index_of_value_that_stores_amplitude_of_singular_fct(); */
    /* 	} */
    /* }    */
 
    /// Pin Lagrange multipliers and set to zero
    void pin_lagrange_multipliers_and_set_to_zero()
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

    // compute the contribution of this face element to the mean-squared of the
    // pressure jump across the boundary of the torus
    Vector<double> mean_squared_pressure_jump() const
    {	
      // number of integration points
      const unsigned n_intpt = this->integral_pt()->nweight();

      // number of nodes
      const unsigned n_node = this->nnode();

      // local coordinates of each integration point
      Vector<double> s(this->Dim-1, 0.0);

      // Set up memory for the shape and test functions
      Shape psi(n_node), test(n_node);

      double mean_sqr_dp = 0;
      double mean_sqr_p_bulk = 0;

      //Loop over the integration points
      //--------------------------------
      for(unsigned ipt=0; ipt<n_intpt; ipt++)
      {
	// get local coordinates of this integration point
	for(unsigned i=0; i<(this->Dim-1); i++)
	{
	  s[i] = this->integral_pt()->knot(ipt,i);
	}
   
	// Get the integral weight
	double w = this->integral_pt()->weight(ipt);
   
	//Find the shape and test functions and return the Jacobian
	//of the mapping
	double J = this->shape_and_test(s, psi, test);
   
	//Premultiply the weights and the Jacobian
	double W = w*J;

	// pressure index
	const unsigned p_index = this->Dim;
	  	  
	// interpolated pressure values at this integration point
	double interpolated_p_augmented = 0;
	double interpolated_p_bulk = 0;

	// cast the GeomObject to a singular line element      
	SingularLineElement* sing_el_pt;

	EdgeCoordinates edge_coords_at_knot;
	Vector<double> s_singular_el(1, 0.0);

	// get the edge coordinates at this knot, the corresponding singular
	// line element and it's local coordinates
	this->edge_coords_and_singular_element_at_knot(ipt,
						       edge_coords_at_knot,
						       sing_el_pt,
						       s_singular_el);
	  
	// get the total contribution of the singular functions and their gradients
	// at this integration point

	// the sum of all scaled singular functions
	Vector<double> u_sing_total =
	  sing_el_pt->total_singular_contribution(edge_coords_at_knot,
						  s_singular_el);
  
	// loop over the nodes and compute the interpolated pressure
	// on both sides of the torus boundary
	for(unsigned l=0; l<n_node; l++) 
	{
	  // QUEHACERES for debug
	  unsigned nval = Orig_node_pt[l]->nvalue();
	    
	  // make sure we're not at an mid-edge node which doesn't
	  // contain the pressure
	  if(Orig_node_pt[l]->nvalue() % this->Dim != 0)
	  {
	    // QUEHACERES debug
	    double p_nodal_val_aug = this->nodal_value(l, p_index);
	    double p_nodal_val_bulk = Orig_node_pt[l]->value(p_index);
	    double psi_l = psi[l];
	      
	    interpolated_p_augmented += this->nodal_value(l, p_index)   * psi[l];
	    interpolated_p_bulk      += Orig_node_pt[l]->value(p_index) * psi[l];

	    // QUEHACERES debug
	    if(!isfinite(interpolated_p_augmented) || !isfinite(interpolated_p_bulk))
	      unsigned breakpoint = 1;
	  }
	}

	// add on the singular contributions to the pressure
	interpolated_p_augmented = u_sing_total[p_index];

	double dp = interpolated_p_bulk - interpolated_p_augmented;
	  
	// and add the weighted contribution to the integral
	mean_sqr_p_bulk    += W * pow(interpolated_p_bulk, 2);
	mean_sqr_dp        += W * pow(dp, 2);
      }

      Vector<double> mean_squares(3, 0.0);
      mean_squares[0] = mean_sqr_p_bulk;
      mean_squares[1] = mean_sqr_dp;
	
      return mean_squares;
    }
      
  private:   

    /// \short Add the element's contribution to its residual vector.
    /// flag=1(or 0): do (or don't) compute the contribution to the
    /// Jacobian as well. 
    void fill_in_generic_residual_contribution_nst_sing_jump(
      Vector<double> &residuals, DenseMatrix<double> &jacobian, 
      const unsigned& flag);

    /// Vector of pointers to orig nodes
    Vector<Node*> Orig_node_pt;

    /// \short Index of external Data that stores the value of the amplitude of
    /// the singular function
    Vector<unsigned> C_external_data_index;
      
    /// \short Index of value (within external Data) that stores the
    /// value of the amplitude of the singular function
    Vector<unsigned> C_external_data_value_index;

    // hierher
    Vector<unsigned> External_data_index_for_bulk_node;

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
      FiniteElement* const& bulk_el_pt, 
      const int& face_index, 
      std::map<Node*,Node*>& existing_duplicate_node_pt,
      const unsigned& boundary_id) : 
  NavierStokesWithSingularityFaceElement<ELEMENT>(bulk_el_pt, face_index, boundary_id),
    C_external_data_index(0), C_external_data_value_index(0)
  {   
    // Back up original nodes and make new ones
    unsigned nnod = this->nnode();
    Orig_node_pt.resize(nnod);
    External_data_index_for_bulk_node.resize(nnod);
    
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
      Node* orig_for_replaced_node_pt = 0;
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
	  std::map<unsigned, unsigned>* index_pt=
	    dynamic_cast<BoundaryNodeBase*>(nod_pt)->
	    index_of_first_value_assigned_by_face_element_pt();                 
	  if (index_pt != 0)
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
	  bulk_el_pt->node_pt(j_in_bulk) = this->node_pt(j);
         
        } // end existing node is already replacement vs make new one
      } 

      // The original node now acts as external data for this element
      // (we still need it to enforce continuity)
      External_data_index_for_bulk_node[j] = this->add_external_data(Orig_node_pt[j]);
    }

    // Make space for Dim Lagrange multipliers
    Vector<unsigned> n_additional_values(nnod, this->Dim);
    this->add_additional_values(n_additional_values, boundary_id);

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

    // shorthand
    ELEMENT* bulk_el_pt = dynamic_cast<ELEMENT*>(this->bulk_element_pt());
    
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
      double J = this->shape_and_test(s, psi, test);
       
      //Premultiply the weights and the Jacobian
      double W = w*J;

      // the velocity fields in the augmented and non-augmented (bulk)
      // regions, and the Lagrange multipliers which enforce continuity
      // of the solution across the boundary of the augmented region
      Vector<double> u_augmented(Dim, 0.0);
      Vector<double> u_bulk(Dim, 0.0);
      Vector<double> lambda(Dim, 0.0);

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

	if(first_index.find(this->Boundary_id) == first_index.end())
	{
	  ostringstream error_message;
	  
	  error_message << "Error: Lagrange multiplier index not found for node "
			<< l << "and boundary ID: " << this->Boundary_id << "\n";
	  
	  throw OomphLibError( error_message.str().c_str(),
			       OOMPH_CURRENT_FUNCTION,
			       OOMPH_EXCEPTION_LOCATION);
	}
	
	for(unsigned i=0; i<Dim; i++)
	{
	  // get the nodal index, accounting for the dimension offset
	  unsigned lambda_index = first_index[this->Boundary_id] + i;

	  // QUEHACERES
	  if(lambda_index < 3)
	  {
	    ostringstream error_message;

	    error_message << "wrong lambda index! Apparently index for "
			  << "Boundary_id: " << this->Boundary_id << " is: "
			  << first_index[this->Boundary_id] << "\n";
	    
	    throw OomphLibError(error_message.str().c_str(),
				OOMPH_CURRENT_FUNCTION,
				OOMPH_EXCEPTION_LOCATION);
	  }
	  
	  u_augmented[i] += this->nodal_value(l,i) * psi[l];
	  u_bulk[i]      += Orig_node_pt[l]->value(i) * psi[l];
	  lambda[i]      += this->nodal_value(l, lambda_index) * psi[l];
	}
      }

      // QUEHACERES delete
      /* // get the line element and local coordinate which corresponds to the */
      /* // singular amplitude for this knot */
      /* std::pair<GeomObject*, Vector<double> > line_elem_and_local_coord =  */
      /* 	this->line_element_and_local_coordinate_at_knot(ipt); */

      // QUEHACERES avoid the hard coded template arg here
      // cast the GeomObject to a singular line element      
      SingularLineElement* sing_el_pt;

      EdgeCoordinates edge_coords_at_knot;
      Vector<double> s_singular_el(1, 0.0);

      // get the edge coordinates at this knot, the corresponding singular
      // line element and it's local coordinates
      this->edge_coords_and_singular_element_at_knot(ipt,
						     edge_coords_at_knot,
						     sing_el_pt,
						     s_singular_el);
      
      // QUEHACERES delete;
      /* = */
      /* dynamic_cast<ScalableSingularityForNavierStokesLineElement<3>*> */
      /* (line_elem_and_local_coord.first); */

      // check we've actually got one, or we're in trouble
      if(sing_el_pt == 0)
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
      Vector<unsigned> sing_ids = sing_el_pt->singular_fct_ids();

      // QUEHACERES delete
      /* // local coordinate in the singular element for the zeta of this knot */
      /* Vector<double> s_singular_el = line_elem_and_local_coord.second; */
      
      /* // get the \rho,\zeta,\phi coordinates at this knot */
      /* EdgeCoordinates edge_coords_at_knot = this->edge_coordinate_at_knot(ipt); */

      // unscaled stuff. These are stored in an array so that they can be
      // looped over when implementing analytic jacobian
      Vector<Vector<double> > u_sing_unscaled(sing_el_pt->nsingular_fct());
      Vector<DenseMatrix<double> > dudx_sing_unscaled(sing_el_pt->nsingular_fct());

      // unscaled singular pressures
      Vector<double> p_sing_unscaled(sing_el_pt->nsingular_fct());
      
      // the sum of all scaled singular functions
      Vector<double> u_sing_total =
	sing_el_pt->total_singular_contribution(edge_coords_at_knot,
                                        	s_singular_el);

      // the total contribution of the singular velocity gradients
      DenseMatrix<double> dudx_sing_total = sing_el_pt->
	total_singular_gradient_contribution(edge_coords_at_knot,
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
	  sing_el_pt->unscaled_singular_fct(edge_coords_at_knot, sing_fct_id);	  
	 
	// unscaled singular gradient 
	dudx_sing_unscaled[ising] =
	  sing_el_pt->gradient_of_unscaled_singular_fct(edge_coords_at_knot, sing_fct_id);
	 
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

      if(bulk_el_pt->stress_fct_pt() == 0)
      {
	throw OomphLibError(
	  "Error: the stress function pointer has not been set for the augmented elements\n",
	  OOMPH_CURRENT_FUNCTION,
	  OOMPH_EXCEPTION_LOCATION);
      }
      
      stress_sing_total = (*bulk_el_pt->stress_fct_pt())(strain_rate_sing_total, p_sing_total);

      // get stress associated with each singular function
      Vector<DenseMatrix<double> > stress_sing_unscaled(sing_el_pt->nsingular_fct());
      
      for(unsigned ising=0; ising < sing_el_pt->nsingular_fct(); ising++)
      {
	stress_sing_unscaled[ising] =
	  (*bulk_el_pt->stress_fct_pt())(strain_rate_sing_unscaled[ising],
					 p_sing_unscaled[ising]);
      }
      
      // Now add to the appropriate equations
      // -----------------------------------------------
      
      //Loop over the test functions
      for(unsigned l=0; l<n_node; l++)
      {
	Node* node_pt = this->node_pt(l);

	// get the map which gives the starting nodal index for
	// the Lagrange multipliers associated with each boundary ID
	std::map<unsigned, unsigned> first_index = *(
	  dynamic_cast<BoundaryNodeBase*>(node_pt)->
	  index_of_first_value_assigned_by_face_element_pt() );
	
	for(unsigned d=0; d<Dim; d++)
	{
	  // Lagrange multiplier equations: Determined from continuity of
	  // solution with (scaled!) singular solution in the augmented region.

	  // get the nodal index of the Lagrange multiplier for this
	  // coordinate direction and boundary ID
	  unsigned lambda_index = first_index[this->Boundary_id] + d;

	  // QUEHACERES
	  if(lambda_index < 3)
	  {
	    ostringstream error_message;

	    error_message << "wrong lambda index! Apparently index for "
			  << "Boundary_id: " << this->Boundary_id << " is: "
			  << first_index[this->Boundary_id] << "\n";
	    
	    throw OomphLibError(error_message.str().c_str(),
				OOMPH_CURRENT_FUNCTION,
				OOMPH_EXCEPTION_LOCATION);
	  }
	  
	  int local_eqn_lagr = this->nodal_local_eqn(l, lambda_index);
	  
	  if (local_eqn_lagr >= 0)
	  {
	    residuals[local_eqn_lagr] += ((u_augmented[d] + u_sing_total[d]) - u_bulk[d]) * test[l]*W;

	    // compute Jacobian
	    if (flag == 1)
	    {
	      for(unsigned l2=0; l2<n_node; l2++)
	      {
               
		int local_unknown_augmented = this->nodal_local_eqn(l2, d);
		if (local_unknown_augmented >= 0)
		{
		  jacobian(local_eqn_lagr,local_unknown_augmented) += psi[l2]*test[l]*W;
		}

		int local_unknown_bulk = this->external_local_eqn(
		  External_data_index_for_bulk_node[l2], d);
		
		if (local_unknown_bulk >= 0)
		{
		  jacobian(local_eqn_lagr, local_unknown_bulk) -=
		    psi[l2]*test[l]*W;
		}
	      }

	      // QUEHACERES come back to this when we have a c equation
	      /* for(unsigned ising=0; ising < sing_el_pt->nsingular_fct(); ising++) */
	      /* { */
	      /* 	// Deriv w.r.t. amplitude is simply the unscaled fct */
	      /* 	int local_eqn_c = this->external_local_eqn(C_external_data_index[ising], */
	      /* 						   C_external_data_value_index[ising]); */
	      /* 	if (local_eqn_c >= 0) */
	      /* 	{ */
	      /* 	  jacobian(local_eqn_lagr, local_eqn_c) += */
	      /* 	    u_sing_unscaled[ising][d] * test[l]*W; */
	      /* 	} */
	      /* } */
	    } // end Jacobian flag check
	  }

	  // Contribution of Lagrange multiplier and traction to bulk eqn in the aug region
	  int local_eqn_augmented = this->nodal_local_eqn(l, d);
	  if (local_eqn_augmented >= 0)
	  {	    
	    residuals[local_eqn_augmented] += lambda[d] * test[l]*W;

	    for(unsigned j=0; j<Dim; j++)
	    {
	      residuals[local_eqn_augmented] -= stress_sing_total(d,j)*unit_normal[j] * test[l]*W; 
	    }

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
		unsigned lambda_index = first_index2[this->Boundary_id] + d;
		
		int local_unknown_lambda = this->nodal_local_eqn(l2, lambda_index);
		if (local_unknown_lambda >= 0)
		{
		  jacobian(local_eqn_augmented,local_unknown_lambda) +=
		    psi[l2]*test[l]*W;
		}
	      }

	      // QUEHACERES come back to this when we have a c equation
	      /* for(unsigned ising=0; ising<sing_el_pt->nsingular_fct(); ising++) */
	      /* { */
	      /* 	int local_eqn_c = this->external_local_eqn(C_external_data_index[ising], */
	      /* 						   C_external_data_value_index[ising]); */
	      /* 	if (local_eqn_c >= 0) */
	      /* 	{ */
	      /* 	  // Deriv w.r.t. amplitude is simply the unscaled singular traction */
	      /* 	  // \hat\tau_{ij} n_j */
	      /* 	  for(unsigned j=0; j< Dim; j++) */
	      /* 	  { */
	      /* 	    jacobian(local_eqn_augmented,local_eqn_c) -= */
	      /* 	      stress_sing_unscaled[ising](d,j) * unit_normal[j] *test[l]*W;       */
	      /* 	  } */

	      /* 	} */
	      /* } */
	    } // end Jacobian flag check
         
	    // Contribution of Lagrange multiplier to bulk eqn in non-augmented region
	    int local_eqn_bulk =
	      this->external_local_eqn(External_data_index_for_bulk_node[l], d);
	  
	    if (local_eqn_bulk >= 0)
	    {	      
	      residuals[local_eqn_bulk] -= lambda[d]*test[l]*W;

	      /* // QUEHACERES experimentally adding */
	      /* for(unsigned j=0; j<Dim; j++) */
	      /* {		 */
	      /* 	residuals[local_eqn_bulk] -= stress_sing_total(d,j)*unit_normal[j] * test[l]*W; */
	      /* } */
	      
	      // compute Jacobian
	      if (flag==1)
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
		  unsigned lambda_index = first_index2[this->Boundary_id] + d;
		  
		  int local_unknown_lambda = this->nodal_local_eqn(l2, lambda_index);
		  if (local_unknown_lambda>=0)
		  {
		    jacobian(local_eqn_bulk, local_unknown_lambda) -=
		      psi[l2]*test[l]*W;
		  }
		}
	      }
	    }

	    // QUEHACERES these lagrange multipliers also contribute to r_c if we evaluate
	    // r_c around the augmented region
	  } // end loop over the dimensions
	} // end loop over test functions
      } // end loop over dimensions
    } // end loop over integration points
    
  } // end of fill_in_generic_residual_contribution_nst_sing_jump

  
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
  template <unsigned DIM, unsigned NNODE_1D>
    class TNavierStokesElementWithSingularity : public virtual TTaylorHoodElement<DIM>
  {
    
  public:
    
    typedef void (*ExactNonSingularFctPt)
      (const Vector<double>& x, Vector<double>& u, DenseMatrix<double>& grad_u);

    // function pointer for helper function which computes the stress
    typedef DenseMatrix<double> (*StressFctPt)(const DenseMatrix<double>& du_dx, const double& p);
    
    // shorthand for the singular line element
    typedef ScalableSingularityForNavierStokesLineElement<NNODE_1D> SingularLineElement;
    
    /// Constructor                            
  TNavierStokesElementWithSingularity() :
    Exact_non_singular_fct_pt(0), Nplot(0), Is_lower_disk_element(false)
    { }

    void compute_error(std::ofstream& outfile,
		       FiniteElement::SteadyExactSolutionFctPt exact_soln_fct,
		       double& v_error, double& p_error, double& norm)
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
      	(*exact_soln_fct)(x, exact_soln);

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
    
    ExactNonSingularFctPt& exact_non_singular_fct_pt()
    {
      return Exact_non_singular_fct_pt;
    }  

    StressFctPt& stress_fct_pt()
    {
      return Stress_fct_pt;
    }

    
    // tell this element that it's on the lower side of the disk
    void set_lower_disk_element()
    {
      Is_lower_disk_element = true;
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
      // ### keeping in atm for debug
      // interpolate the position
      Vector<double> x(DIM);
	
      for(unsigned i=0; i<DIM; i++)
      {
      	x[i] = this->interpolated_x(s,i);
      }

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
	  this->line_element_and_local_coordinate_at_knot_point(ipt);
      }
      
      // cast the GeomObject to a singular line element      
      SingularLineElement* sing_el_pt =
	dynamic_cast <SingularLineElement*>( line_elem_and_local_coord.first );

      // check if we're subtracting the singularity or not
      if (sing_el_pt != 0)
      {
	// local coordinate in the singular element for the zeta of this plot point
	Vector<double> s_singular_el = line_elem_and_local_coord.second;
      
	// get the \rho,\zeta,\phi coordinates at this knot
	EdgeCoordinates edge_coords_at_point;

	if(use_plot_points)
	  edge_coords_at_point = this->edge_coordinate_at_plot_point(ipt);
	else
	  edge_coords_at_point = this->edge_coordinate_at_knot_point(ipt);

	Vector<double> u_sing =
	  sing_el_pt->total_singular_contribution(edge_coords_at_point,
						  s_singular_el);

	// add singular part of the solution to the FE part to give the total
	// computed solution
	for(unsigned i=0; i<DIM+1; i++)
	{
	  u_fe[i] += u_sing[i];
	}

	// QUEHACERES delete
	/* for(unsigned ising=0; ising < sing_el_pt->nsingular_fct(); ising++) */
	/* { */
	/*   // singular part of the solution */
	/*   Vector<double> u_sing = sing_el_pt->singular_fct(edge_coords_at_point, */
	/* 						   s_singular_el, */
	/* 						   ising); */
	    
	/*   // add singular part of the solution to the FE part to give the total */
	/*   // computed solution */
	/*   for(unsigned i=0; i<DIM+1; i++) */
	/*   { */
	/*     u_fe[i] += u_sing[i]; */
	/*   } */
	  
	/* } */
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
      Vector<double> exact_soln(4,0.0);

      // Vector for coordinates
      Vector<double> x(DIM);
      
      // Get x position as Vector
      this->interpolated_x(s, x);
	
      // if we have an exact solution pointer, then get the exact solution
      if (exact_soln_pt != 0)
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
	  
	  if(sing_el_pt != 0)
	  {
	    EdgeCoordinates edge_coords = edge_coordinate_at_plot_point(iplot);
	    
	    DenseMatrix<double> dudx_sing =
	      sing_el_pt->total_singular_gradient_contribution(edge_coords,
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

      // QUEHACERES debug @@@@@@
      double tol = 1e-6;
      if(abs(x[0]+0.9) < tol && abs(x[1]-0.0945938) < tol && abs(x[2]) < tol)
      	double breakpoint = 1;
      
      // regular part of the solution
      Vector<double> u_exact_non_sing(DIM+1, 0.0);
            
      // QUEHACERES delete?
      /* DenseMatrix<double> dudx(DIM, DIM); */
      /* // Overwrite with exact version! */
      /* if (Exact_non_singular_fct_pt != 0) */
      /* { */
      /* 	Exact_non_singular_fct_pt(x, u_exact_non_sing, dudx); */
      /* } */
		
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
      if(sing_el_pt == 0)
      {
	// output the regular part of the exact solution, which is the same as
	// the total exact solution if there's no singular part (i.e. we're in the bulk)
	if(exact_soln_pt != 0)
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
	  outfile << "0 0 0 0 ";
	}

	// -----------------
	
	// now for the singular contributions, just output zeros if there are
	// no singular functions so the oomph-convert script doesn't die
	// when the number of columns isn't the same for all elements
	for(unsigned ising=0; ising < Nsingular_fct; ising++)
	{	
	  outfile << "0 0 0 0 ";
	}
	
	outfile << std::endl;
	
	return;
      }

      // get the \rho,\zeta,\phi coordinates at this plot point
      EdgeCoordinates edge_coords_at_plot =
	this->edge_coordinate_at_plot_point(iplot);

      // get the total singular contribution at this point
      Vector<double> u_sing_total = 
	sing_el_pt->total_singular_contribution(edge_coords_at_plot,
						s_singular_el);

      // Get exact solution at this point
      Vector<double> exact_soln(4, 0.0);
      if(exact_soln_pt != 0)
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
      for(unsigned ising=0; ising < sing_ids.size(); ising++)
      {
	// get the ID
	unsigned sing_fct_id = sing_ids[ising];
	
	// singular part of the solution
	Vector<double> u_sing = sing_el_pt->singular_fct(edge_coords_at_plot,
							 s_singular_el,
							 sing_fct_id);

	// QUEHACERES debug
	bool is_finite = isfinite(u_sing[0]) && isfinite(u_sing[1])
	  && isfinite(u_sing[2]) && isfinite(u_sing[3]);

	if(!is_finite)
	{
	  oomph_info << "NOT FINITE!!" << std::endl;
	  abort();
	}
	
	for(unsigned i=0; i<DIM+1; i++)
	{
	  outfile << u_sing[i] << " ";
	}
      }
      // QUEHACERES leave stress out for the time being
      // FE stress and strain-rate tensors
      /* DenseMatrix<double> stress_fe(DIM, DIM); */
      /* DenseMatrix<double> strain_rate_fe(DIM, DIM); */

      /* // singular stress and strain-rate tensors */
      /* DenseMatrix<double> stress_sing(DIM, DIM); */
      /* DenseMatrix<double> stress_total(DIM, DIM); */
	
      /* DenseMatrix<double> strain_rate_sing(DIM, DIM); */
      /* DenseMatrix<double> strain_rate_total(DIM, DIM); */
	
      /* // get the strain rates */
      /* this->strain_rate(s, strain_rate_fe); */
	
      /* // Get gradient of scaled singular velocity functions       */
      /* DenseMatrix<double> dudx_sing(DIM, DIM); */
      /* if(Navier_stokes_sing_el_pt != 0) */
      /* { */
      /* 	dudx_sing = Navier_stokes_sing_el_pt->gradient_of_singular_fct(x); */
      /* } */
      
      /* // compute the unscaled singular contribution to the strain-rate       */
      /* for (unsigned i=0; i<DIM; i++) */
      /* { */
      /* 	for(unsigned j=0; j<DIM; j++) */
      /* 	{ */
      /* 	  strain_rate_sing(i,j)  = 0.5*(dudx_sing(i,j) + dudx_sing(j,i)); */
      /* 	  strain_rate_total(i,j) = strain_rate_fe(i,j) + strain_rate_sing(i,j); */
      /* 	} */
      /* } */
		
      /* // extract pressure from solution */
      /* double p_sing  = u_sing[DIM]; */
      /* double p_fe    = u_fe[DIM];	 */
      /* double p_total = u_fe_plus_sing[DIM]; */
	
      /* // compute stress from constitutive equation */
      /* stress_sing  = (this->stress_fct_pt())(strain_rate_sing, p_sing); */
      /* stress_fe    = (this->stress_fct_pt())(strain_rate_fe, p_fe); */
      /* stress_total = (this->stress_fct_pt())(strain_rate_total, p_total); */
      
      /* for(unsigned i=0; i<DIM; i++) */
      /* { */
      /* 	for(unsigned j=0; j<DIM; j++) */
      /* 	{ */
      /* 	  outfile << stress_total(i,j) << " "; */
      /* 	}  */
      /* } */
      /* for(unsigned i=0; i<DIM; i++) */
      /* { */
      /* 	for(unsigned j=0; j<DIM; j++) */
      /* 	{ */
      /* 	  outfile << stress_fe(i,j) << " "; */
      /* 	}  */
      /* } */
      /* for(unsigned i=0; i<DIM; i++) */
      /* { */
      /* 	for(unsigned j=0; j<DIM; j++) */
      /* 	{ */
      /* 	  outfile << stress_sing(i,j) << " "; */
      /* 	}  */
      /* } */
      outfile << std::endl;
    }
    
    /// Output with various contributions
    void output_with_various_contributions(std::ostream& outfile, 
					   const unsigned& nplot,
					   FiniteElement::SteadyExactSolutionFctPt exact_soln_pt = 0) const
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
    void set_edge_coordinates_at_plot_point(const Vector<EdgeCoordinates>& coords)
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
      Edge_coordinates_at_plot_point.resize(Nplot);

      // set 'em
      for(unsigned i=0; i<Nplot; i++)
      {
	Edge_coordinates_at_plot_point[i]  = coords[i];
      }
    }

    // set the edge coordinates of each knot point
    void set_edge_coordinates_at_knot_point(const Vector<EdgeCoordinates>& coords)
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
      Edge_coordinates_at_knot_point.resize(n_intpt);

      // set 'em
      for(unsigned i=0; i<n_intpt; i++)
      {
	Edge_coordinates_at_knot_point[i]  = coords[i];
      }
    }

    // set the edge coordinates at each knot
    void set_nodal_edge_coordinates(const Vector<EdgeCoordinates>& edge_coords)
    {
      // number of coordinates supplied
      unsigned ncoords = edge_coords.size();

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
      Nodal_edge_coordinates.resize(n_node);
      
      // set 'em
      for(unsigned i=0; i<n_node; i++)
      {
	Nodal_edge_coordinates[i].rho  = edge_coords[i].rho;
	Nodal_edge_coordinates[i].zeta = edge_coords[i].zeta;
	Nodal_edge_coordinates[i].phi  = edge_coords[i].phi;
      }     
    }

    // ========================================================================
    // ========================================================================
    
    // get the edge coordinates at the ith plot point
    EdgeCoordinates edge_coordinate_at_plot_point(const unsigned& i) const
    {
      return Edge_coordinates_at_plot_point[i];
    }

    // get the edge coordinates at the ith knot point
    EdgeCoordinates edge_coordinate_at_knot_point(const unsigned& i) const
    {
      return Edge_coordinates_at_knot_point[i];
    }

    // get the edge coordinates at the ith node
    EdgeCoordinates nodal_edge_coordinates(const unsigned& i) const
    {
      return Nodal_edge_coordinates[i];
    }

    EdgeCoordinates interpolated_edge_coordinates(const Vector<double>& s) const
    {
      EdgeCoordinates interpolated_edge_coords;

      const unsigned n_node = this->nnode();
      
      //Set up memory for the shape functions
      Shape psi(n_node);

      // get 'em
      this->shape(s, psi);
      
      // The number of integration points (knots)
      const unsigned n_intpt = this->integral_pt()->nweight();

      interpolated_edge_coords.rho  = 0.0;
      interpolated_edge_coords.zeta = 0.0;
      interpolated_edge_coords.phi  = 0.0;
	
      // loop over each node in this element and add its contribution to the
      // interpolated coordinates
      for(unsigned j=0; j<n_node; j++)
      {
	// get the edge coordinates of this node
	EdgeCoordinates nodal_edge_coords = Nodal_edge_coordinates[j];

	// interpolated with the shape functions
	interpolated_edge_coords.rho  += nodal_edge_coords.rho  * psi[j];
	interpolated_edge_coords.zeta += nodal_edge_coords.zeta * psi[j];
	interpolated_edge_coords.phi  += nodal_edge_coords.phi  * psi[j];
      }     
      
      return interpolated_edge_coords;
    }

    // ========================================================================
    // ========================================================================

    void edge_coords_and_singular_element_at_node(const unsigned& inode,
						  EdgeCoordinates& edge_coords,
						  SingularLineElement*& sing_el_pt,
						  Vector<double>& s_singular_el) const
    {
      // get the line element and local coordinate which corresponds to the
      // singular amplitude for this knot
      std::pair<GeomObject*, Vector<double> > line_elem_and_local_coord = 
	this->line_element_and_local_coordinate_at_node(inode);

      // QUEHACERES avoid the hard coded template arg here
      // cast the GeomObject to a singular line element      
      sing_el_pt = dynamic_cast<SingularLineElement*> (line_elem_and_local_coord.first);

      // check we've actually got one, or we're in trouble
      if(sing_el_pt == 0)
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
      edge_coords = this->nodal_edge_coordinates(inode);
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
    void set_line_element_and_local_coordinate_at_knot_point(
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
      Line_element_and_local_coordinate_at_knot_point.resize(n_intpt);

      // set 'em
      for(unsigned i=0; i<n_intpt; i++)
      {
	Line_element_and_local_coordinate_at_knot_point[i] =
	  line_element_and_local_coordinate_at_knot_point[i];
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
	GeomObject* dummy_pt = 0;
	Vector<double> dummy_vec(1,0);
	return_pair = std::make_pair(dummy_pt, dummy_vec);
      }
      return return_pair;
    }

    // get the line element and local coordinate for the ith knot point
    std::pair<GeomObject*, Vector<double> >
      line_element_and_local_coordinate_at_knot_point(const unsigned& i) const
    {
      std::pair<GeomObject*, Vector<double> > return_pair;

      // check if we have an element-coordinate pair for this plot point
      // (since this bulk element is used everywhere in the bulk not just
      // in the augmented region, there may be no singular functions
      // associated with this element, and so we just want to return an
      // empty pointer and vector)
      if(i < Line_element_and_local_coordinate_at_knot_point.size())
      {
	return_pair = Line_element_and_local_coordinate_at_knot_point[i];
      }
      else
      {
	GeomObject* dummy_pt = 0;
	Vector<double> dummy_vec(1,0);
	return_pair = std::make_pair(dummy_pt, dummy_vec);
      }
      return return_pair;
    }
    
    /* // override to output scalar values for paraview format */
    /* void scalar_value_fct_paraview(std::ofstream& file_out, */
    /* 				   const unsigned& i, */
    /* 				   const unsigned& nplot) const */
    /* { */
      
    /* } */
    
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

    /// \short Number of singular functions we're subtracting (only used for
    /// output purposes, for bulk elements which have no singular line element pointer
    /// but need to get the number of output columns to match the elements
    /// with singular functions)
    unsigned Nsingular_fct;
    
    /// Pointer to exact non-singular fct (only for post-processing!)
    ExactNonSingularFctPt Exact_non_singular_fct_pt;

    /// Pointer to function which computes the stress
    StressFctPt Stress_fct_pt;
    
    /// \short Edge coordinates (\rho, \zeta, \phi) of each of this element's
    /// plot points
    Vector<EdgeCoordinates> Edge_coordinates_at_plot_point;

    // same as above but for each knot point (to compute error)
    Vector<EdgeCoordinates> Edge_coordinates_at_knot_point;

    // same as above but for each node
    Vector<EdgeCoordinates> Nodal_edge_coordinates;
    
    /// \short The line element and its local coordinate(s) which correspond to the zeta
    /// values at each plot point in this bulk element. N.B. the element
    /// is 1D so only has 1 local coordinate, but we keep this as a vector
    /// for consistency with the general case
    Vector<std::pair<GeomObject*, Vector<double> > >
      Line_element_and_local_coordinate_at_plot_point;

    // same as above but for knot points (to compute error)
    Vector<std::pair<GeomObject*, Vector<double> > >
      Line_element_and_local_coordinate_at_knot_point;

     // same as above but for nodes
    Vector<std::pair<GeomObject*, Vector<double> > >
      Line_element_and_local_coordinate_at_node;
    
    /// \short Flag to store whether this is a lower disk element
    /// (useful for output stuff where we're using Eulerian coordinates and need to
    /// decide whether a point on the disk surface has phi=0 or phi=2pi to get
    /// the pressure jump correct)
    bool Is_lower_disk_element;
  };


  //////////////////////////////////////////////////////////////////////
  //////////////////////////////////////////////////////////////////////
  //////////////////////////////////////////////////////////////////////


  //=======================================================================
  /// Face geometry for the TNavierStokesElementWithSingularity elements: The spatial 
  /// dimension of the face elements is one lower than that of the
  /// bulk element but they have the same number of points
  /// along their 1D edges.
  //=======================================================================
  template <unsigned DIM, unsigned NNODE_1D>
    class FaceGeometry<TNavierStokesElementWithSingularity<DIM,NNODE_1D> > :
    public virtual TElement<DIM-1,NNODE_1D>
  {

  public:
 
    /// \short Constructor: Call the constructor for the
    /// appropriate lower-dimensional TElement
  FaceGeometry() : TElement<DIM-1,NNODE_1D>() {}

  };

  ////////////////////////////////////////////////////////////////////////
  ////////////////////////////////////////////////////////////////////////
  ////////////////////////////////////////////////////////////////////////


  //=======================================================================
  /// Face geometry for the 1D TNavierStokesElementWithSingularity elements: Point elements
  //=======================================================================
  template <unsigned NNODE_1D>
    class FaceGeometry<TNavierStokesElementWithSingularity<1,NNODE_1D> >: 
    public virtual PointElement
    {

    public:
 
      /// \short Constructor: Call the constructor for the
      /// appropriate lower-dimensional TElement
    FaceGeometry() : PointElement() {} 

    };
  
}

#endif
