#ifndef OOMPH_CHEBYSHEV_GAUSS_HEADER
#define OOMPH_CHEBYSHEV_GAUSS_HEADER

namespace oomph
{
  //===================================================================
  /// Class for multidimensional Chebyshev-Gauss integration 
  /// rules
  /// empty - just establishes template parameters
  //===================================================================
  template <unsigned DIM, unsigned NPTS_1D>
    class ChebyshevGauss
  { };

  template<unsigned NPTS_1D>
    class ChebyshevGauss<1, NPTS_1D> : public Integral
  {
  public:

    /// Deafault constructor. Calculates and stores GLL nodes
    ChebyshevGauss()
    {
      //Temporary storage for the integration points
      Vector<double> s(NPTS_1D), w(NPTS_1D);
      
      //call the function that calculates the points
     chebyshev_gauss_nodes(NPTS_1D, s, w);
      
      //Populate the arrays
      for(unsigned i=0; i<NPTS_1D; i++)
      {
	Knot[i][0] = s[i];
	Weight[i]  = w[i];
      }
    }

    /// Number of integration points of the scheme   
    unsigned nweight() const
    {
      return Npts;
    }

    /// Return coordinate s[j] (j=0) of integration point i
    double knot(const unsigned& i, const unsigned& j) const
    {
      return Knot[i][j];
    }

    /// Return weight of integration point i
    double weight(const unsigned& i) const
    {
      return Weight[i];
    }

      
    static void chebyshev_gauss_abscissas(const unsigned& nnode, Vector<double>& x)
    {
      // make space
      x.resize(nnode, 0.0);

      if(nnode < 2)
      {
	throw OomphLibError("Can't do Chebyshev-Gauss integration with < 2 nodes",
			    OOMPH_CURRENT_FUNCTION,
			    OOMPH_EXCEPTION_LOCATION);
      }
      else if (nnode == 2)
      {
	x[0] = -sqrt(2.0)/2.0;
	x[1] = -x[0];
      }
      else
      {      
	unsigned mid = 0;

	// if the number of nodes is odd then the abscissa of the midpoint is 0
	if(nnode % 2 != 0)
	{
	  mid = (nnode-1)/2;
	  x[mid] = 0.0;
	}
	else
	  mid = nnode/2;
      
	for(unsigned i=0; i<mid; i++)
	{
	  x[i] = -std::fabs(
	    cos( (2.0*(double)(i+1) - 1) * MathematicalConstants::Pi / (2.0 * nnode) ) );
	  x[nnode-i-1] = -x[i];
	}
      }
    }

    static void chebyshev_gauss_nodes(const unsigned& nnode,
				      Vector<double>& x, Vector<double>& w)
    {
      // get the knot positions
      chebyshev_gauss_abscissas(nnode, x);

      // now compute weights
      w.resize(nnode, 0.0);

      // nominal weights are pi/n. Here we're also adding in the
      // sqrt(1-s^2) factor into the weighting so that it doesn't
      // need adding separately when evaluating an integral	  
      for(unsigned i=0; i<nnode; i++)
      {
	w[i] = MathematicalConstants::Pi / nnode *
	  sqrt(1.0 - x[i]*x[i]);
      }
    }
    
  private:
    
    /// Number of integration points in scheme
    static const unsigned Npts = NPTS_1D;
    
    /// Array to hold weight and knot points
    //These are not constant, because they are calculated on the fly
    double Knot[NPTS_1D][1], Weight[NPTS_1D];  
  };

  // ==========================================================================
  // quadrature scheme formed by the 2D tensor product of a Chebyshev-Gauss
  // scheme in the 1 direction and a standard Gauss scheme in the 2 direction
  // ==========================================================================
  template<unsigned NPTS_1, unsigned NPTS_2>
    class TwoDChebyshevGaussTensorProductGaussLegendre : public Integral
  {
  public:
    
    TwoDChebyshevGaussTensorProductGaussLegendre()
    {
      // Temporary storage for the Chebyshev-Gauss knot points
      Vector<double> s_cheb_gauss(NPTS_1, 0.0), weight_cheb_gauss(NPTS_1, 0.0);

      // get 'em
      ChebyshevGauss<1, NPTS_1>::chebyshev_gauss_nodes(NPTS_1,
						       s_cheb_gauss,
						       weight_cheb_gauss);

      // now get the standard Gauss knots
      GaussLegendre<1, NPTS_2> gauss_quadrature;

      // set total number of knots and make space
      Npts = NPTS_1 * NPTS_2;
      
      Knot.resize(Npts);
      Weight.resize(Npts, 0.0);

      // keep track of the integration point numbering
      unsigned ipt = 0;
      
      for(unsigned i=0; i<NPTS_1; i++)
      {
	for(unsigned j=0; j<NPTS_2; j++)
	{
	  // make space (this scheme is strictly 2D so this is hard-coded here)
	  Knot[ipt].resize(2, 0.0);
	
	  // we're doing Chebyshev-Gauss spacing in the x1 direction and 
	  // Gauss spacing in the x2 direction
	  Knot[ipt][0] = s_cheb_gauss[i];
	  Knot[ipt][1] = gauss_quadrature.knot(j, 0);

	  // weight is just the tensor product of the weights of the two schemes
	  Weight[ipt] = weight_cheb_gauss[i] * gauss_quadrature.weight(j);

	  ipt++;
	}
      }
    }

     /// Number of integration points of the scheme
    unsigned nweight() const { return Npts; }
 
    /// Return coordinate s[j] of integration point i
    double knot(const unsigned& i, const unsigned& j) const
    {return Knot[i][j];}

    /// Return weight of integration point i
    double weight(const unsigned& i) const { return Weight[i]; }
    
  private:

    // number of integration points
    unsigned Npts;
    
    // dynamically allocated knots and weights
    Vector<Vector<double>> Knot;
    Vector<double> Weight;
  };
}


#endif
