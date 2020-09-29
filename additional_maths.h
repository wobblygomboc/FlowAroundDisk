#ifndef OOMPH_ADDITIONAL_MATHS
#define OOMPH_ADDITIONAL_MATHS

class mVector : public Vector<double>
{
  
public:

mVector() : Vector<double>() {}  
mVector(const unsigned long &n) : Vector<double>(n) {}
mVector(const size_t& n, const double& val) : Vector<double>(n, val) {}
  
  // copy constructor
mVector(const mVector& v) : Vector<double>(v.size(), 0)
  {
    for(unsigned i=0; i<size(); i++)
    {
      (*this)[i] = v[i];
    }	 
  }

  // copy constructor to convert from base class to derived class
mVector(const Vector<double>& v) : Vector<double>(v.size(), 0)
  {
    for(unsigned i=0; i<size(); i++)
    {
      (*this)[i] = v[i];
    }
  }
  
  // dot product
  double operator*(const mVector& v) const
  {      
    // check the dimensions are consistent
    if(size() != v.size())
    {
      ostringstream error_message;
      error_message << "Error, cannot perform vector dot product on vectors of "
		    << "different sizes";
	
      throw OomphLibError(error_message.str(),
			  OOMPH_CURRENT_FUNCTION,
			  OOMPH_EXCEPTION_LOCATION);
    }
      
    double result = 0;

    // dot product a*v = a_i*v_i
    for(unsigned i=0; i<size(); i++)
    {
      result += (*this)[i] * v[i];
    }

    return result;
  }

  // scalar multiplication of a vector
  mVector operator*(const double& a) const
  {      
    mVector result(size(), 0);

    for(unsigned i=0; i<size(); i++)
    {
      result[i] = (*this)[i] * a;
    }
      
    return result;
  }

  // cross product
  mVector operator^(const mVector& b) const
  {      
    // check the dimensions are consistent
    if(size() != b.size())
    {
      ostringstream error_message;
      error_message << "Error, cannot perform vector dot product on vectors of "
		    << "different sizes";
	
      throw OomphLibError(error_message.str(),
			  OOMPH_CURRENT_FUNCTION,
			  OOMPH_EXCEPTION_LOCATION);
    }

    if(size() != 3)
    {
       ostringstream error_message;
       error_message << "Error, can only do a cross product in 3D";
	
      throw OomphLibError(error_message.str(),
			  OOMPH_CURRENT_FUNCTION,
			  OOMPH_EXCEPTION_LOCATION);
    }
    
    // shorthand
    mVector a = (*this);
    
    mVector result(3);
    result[0] = a[1]*b[2] - a[2]*b[1];
    result[1] = a[2]*b[0] - a[0]*b[2];
    result[2] = a[0]*b[1] - a[1]*b[0];
    
    return result;
  }
  
  mVector operator+(const mVector& v) const
  {
    if( size() != v.size() )
    {
      ostringstream error_message;
      error_message << "Error, can't add vectors of different lengths";
	
      throw OomphLibError(error_message.str(),
			  OOMPH_CURRENT_FUNCTION,
			  OOMPH_EXCEPTION_LOCATION);
    }

    mVector result(size());
      
    for(unsigned i=0; i<size(); i++)
    {
      result[i] = (*this)[i] + v[i];	
    }
    return result;
  }

  mVector operator-(const mVector& v) const
  {
    if( size() != v.size() )
    {
      ostringstream error_message;
      error_message << "Error, can't subtract vectors of different lengths";
	
      throw OomphLibError(error_message.str(),
			  OOMPH_CURRENT_FUNCTION,
			  OOMPH_EXCEPTION_LOCATION);
    }

    mVector result(size());
      
    for(unsigned i=0; i<size(); i++)
    {
      result[i] = (*this)[i] - v[i];
    }
    return result;
  }

  mVector operator-() const
  {
    mVector result(size());
      
    for(unsigned i=0; i<size(); i++)
    {
      result[i] = -(*this)[i];	
    }
    return result;
  }

  void operator+=(const mVector& v)
  {
    if( size() != v.size() )
    {
      ostringstream error_message;
      error_message << "Error, can't add vectors of different lengths";
	
      throw OomphLibError(error_message.str(),
			  OOMPH_CURRENT_FUNCTION,
			  OOMPH_EXCEPTION_LOCATION);
    }

    for(unsigned i=0; i<size(); i++)
    {
      (*this)[i] += v[i];	
    }
  }

  void operator-=(const mVector& v)
  {
    if( size() != v.size() )
    {
      ostringstream error_message;
      error_message << "Error, can't add vectors of different lengths";
	
      throw OomphLibError(error_message.str(),
			  OOMPH_CURRENT_FUNCTION,
			  OOMPH_EXCEPTION_LOCATION);
    }

    for(unsigned i=0; i<size(); i++)
    {
      (*this)[i] -= v[i];	
    }
  }
  
  void operator*=(const double& a)
  {
    for(unsigned i=0; i<size(); i++)
    {
      (*this)[i] *= a;	
    }
  }

  double magnitude() const
  {
    double mag = 0;
    for(unsigned i=0; i<size(); i++)
    {
      mag += (*this)[i] * (*this)[i];
    }
    mag = sqrt(mag);
    return mag;
  }

  static mVector e_x()
  {
    mVector x(3,0);
    x[0] = 1;
    return x;
  }

  static mVector e_y()
  {
    mVector x(3,0);
    x[1] = 1;
    return x;
  }
  
  static mVector e_z()
  {
    mVector x(3,0);
    x[2] = 1;
    return x;
  }
      
  friend ostream& operator<<(ostream& output, const mVector& v)
  {
    for(unsigned i=0; i<v.size(); i++)
    {
      output << v[i] << " ";
    }
    return output;
  }
};

///////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////


/// \short Class which implements a mathematical matrix with standard operations
/// such as matrix and vector multiplication
class TransformationMatrix : public DenseDoubleMatrix
{
public:
  
  TransformationMatrix (const unsigned long &n) : DenseDoubleMatrix(n) {}
  
  TransformationMatrix (const unsigned long &n, const unsigned long &m)
    : DenseDoubleMatrix(n,m) {}

  TransformationMatrix (const unsigned long &n, const unsigned long &m,
			const double &initial_val)
    : DenseDoubleMatrix(n,m,initial_val){}

  // copy constructor
TransformationMatrix(const TransformationMatrix& M)
  : DenseDoubleMatrix(M.nrow(), M.ncol())
  {
    for(unsigned i=0; i<nrow(); i++)
    {
      for(unsigned j=0; j<ncol(); j++)
      {
	(*this)(i,j) = M(i,j);
      }
    }	 
  }

  // matrix multiplication
  TransformationMatrix operator*(const TransformationMatrix& B) const
  {
    unsigned nrow_A = nrow();
    unsigned ncol_A = ncol();
      
    unsigned nrow_B = B.nrow();
    unsigned ncol_B = B.ncol();

    // check the dimensions are consistent
    if(ncol_A != nrow_B)
    {
      ostringstream error_message;
      error_message << "Error, cannot perform matrix multiplication A*B if "
		    << "A.ncol() != B.nrow()";
	
      throw OomphLibError(error_message.str(),
			  OOMPH_CURRENT_FUNCTION,
			  OOMPH_EXCEPTION_LOCATION);
    }
      
    TransformationMatrix result(nrow_A, ncol_B);

    // matrix multiplication M = A*B => M_ij = A_ik*B_kj
    for(unsigned i=0; i<nrow_A; i++)
    {
      for(unsigned j=0; j<ncol_B; j++)
      {
	result(i,j) = 0;
	for(unsigned k=0; k<ncol_A; k++)
	{
	  result(i,j) += (*this)(i,k)*B(k,j);
	}
      }
    }

    return result;
  }

  // matrix-vector multiplication
  mVector operator*(const Vector<double>& v) const
  {
    unsigned nrow_A = nrow();
    unsigned ncol_A = ncol();
    unsigned nrow_v = v.size();

    if(ncol_A != nrow_v)
    {
      ostringstream error_message;
      error_message << "Error, cannot perform matrix*vec multiplication A*v if "
		    << "A.ncol() != v.size()";
	
      throw OomphLibError(error_message.str(),
			  OOMPH_CURRENT_FUNCTION,
			  OOMPH_EXCEPTION_LOCATION);
    }
      
    mVector result_vec(nrow_A, 0);

    // matrix-vector multiplication r = A*v => r_i = A_ij*v_j
    for(unsigned i=0; i<nrow_A; i++)
    {
      for(unsigned j=0; j<ncol_A; j++)
      {
	result_vec[i] += (*this)(i,j)*v[j];
      }
    }
      
    return result_vec;
  }

  // matrix addition
  TransformationMatrix operator+(const TransformationMatrix& B) const
  {
    if( (nrow() != B.nrow()) || (ncol() != B.ncol()) )
    {
      ostringstream error_message;
      error_message << "Error adding matrices, dimensions don't match";
	
      throw OomphLibError(error_message.str(),
			  OOMPH_CURRENT_FUNCTION,
			  OOMPH_EXCEPTION_LOCATION);
    }

    TransformationMatrix result(nrow(), ncol() );
      
    for(unsigned i=0; i<nrow(); i++)
    {
      for(unsigned j=0; j<ncol(); j++)
      {
	result(i,j) = (*this)(i,j) + B(i,j);
      }
    }
    return result;
  }

  // normal matrix subtraction
  TransformationMatrix operator-(const TransformationMatrix& B) const
  {
    if( (nrow() != B.nrow()) || (ncol() != B.ncol()) )
    {
      ostringstream error_message;
      error_message << "Error subtracting matrices, dimensions don't match";
	
      throw OomphLibError(error_message.str(),
			  OOMPH_CURRENT_FUNCTION,
			  OOMPH_EXCEPTION_LOCATION);
    }

    TransformationMatrix result(nrow(), ncol() );
      
    for(unsigned i=0; i<nrow(); i++)
    {
      for(unsigned j=0; j<ncol(); j++)
      {
	result(i,j) = (*this)(i,j) - B(i,j);
      }
    }
    return result;
  }

  // operator to negate a matrix without an rvalue
  TransformationMatrix operator-() const
  {
    TransformationMatrix result(nrow(), ncol() );
      
    for(unsigned i=0; i<nrow(); i++)
    {
      for(unsigned j=0; j<ncol(); j++)
      {
	result(i,j) = -(*this)(i,j);
      }
    }
    return result;
  }
  
  friend ostream& operator<<(ostream& output, const TransformationMatrix& M)
  {
    for(unsigned i=0; i<M.nrow(); i++)
    {
      for(unsigned j=0; j<M.ncol(); j++)
      {
	output << M(i,j) << " ";	    
      }
      output << std::endl;
    }
    return output;
  }
};

namespace Analytic_Functions
{ 
  // secant, for mathematica output
  double Sec(const double& theta)
  {
    return 1.0/cos(theta);
  }

  double delta(const unsigned& i, const unsigned& j)
  {
    return static_cast<double>(i == j);
  }

  // function to map any angle (+ve or -ve) to [0,2pi)
  double map_angle_to_range_0_to_2pi(const double& signed_angle)
  {
    // shorthand
    const double two_pi = 2.0 * MathematicalConstants::Pi;

    // effectively implements Python modulo division
    // (C++ would return a negative angle if we did -angle % 2pi)    
    return signed_angle - two_pi * floor(signed_angle/two_pi);    
  }
  
  // reflect an angle the z-axis
  double reflect_angle_wrt_z_axis(const double& theta)
  {
    // flip by subtracting from pi, then make map back to [0,2pi)
    return map_angle_to_range_0_to_2pi(MathematicalConstants::Pi - theta);
  }

  DenseMatrix<double> strain_rate(const DenseMatrix<double>& du_dx)
  {
    DenseMatrix<double> strain_rate(3,3,0);
    
    for(unsigned i=0; i<3; i++)
    {
      for(unsigned j=0; j<3; j++)
      {
	strain_rate(i,j) = 0.5*(du_dx(i,j) + du_dx(j,i));
      }
    }

    return strain_rate;
  }
  
  /// \short Newtonian stress tensor
  DenseMatrix<double> stress(const DenseMatrix<double>& strain_rate,
			     const double& p)
  {
    const unsigned dim = 3;
      
    // stress tensor
    DenseMatrix<double> stress(dim,dim);
    
    for(unsigned i=0; i<dim; i++)
    {
      for(unsigned j=0; j<dim; j++)
      {
	// Newtonian constitutive relation
	stress(i,j) = -p*delta(i,j) + 2.0*strain_rate(i,j);
      }
    }

    return stress;
  }
  
} // end of Analytic_Functions namespace

#endif
