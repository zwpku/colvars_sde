// -*- c++ -*-

// This file is part of the Collective Variables module (Colvars).
// The original version of Colvars and its updates are located at:
// https://github.com/Colvars/colvars
// Please update all Colvars source files before making any changes.
// If you wish to distribute your changes, please submit them to the
// Colvars repository at GitHub.

#ifndef COLVARTYPES_H
#define COLVARTYPES_H

#include <sstream> // TODO specialize templates and replace this with iosfwd
#include <vector>

#include "colvarmodule.h"

#ifndef PI
#define PI 3.14159265358979323846
#endif

// ----------------------------------------------------------------------
/// Linear algebra functions and data types used in the collective
/// variables implemented so far
// ----------------------------------------------------------------------


/// \brief Arbitrary size array (one dimensions) suitable for linear
/// algebra operations (i.e. for floating point numbers it can be used
/// with library functions)
template <class T> class colvarmodule::vector1d
{
protected:

  std::vector<T> data;

public:

  /// Default constructor
  inline vector1d(size_t const n = 0)
  {
    data.resize(n);
    reset();
  }

  /// Constructor from C array
  inline vector1d(size_t const n, T const *t)
  {
    data.resize(n);
    reset();
    size_t i;
    for (i = 0; i < size(); i++) {
      data[i] = t[i];
    }
  }

  /// Explicit Copy constructor
  inline vector1d(const vector1d&) = default;

  /// Explicit Copy assignement
  inline vector1d& operator=(const vector1d&) = default;

  /// Return a pointer to the data location
  inline T * c_array()
  {
    if (data.size() > 0) {
      return &(data[0]);
    } else {
      return NULL;
    }
  }

  /// Return a reference to the data
  inline std::vector<T> &data_array()
  {
    return data;
  }

  /// Return a reference to the data
  inline std::vector<T> const &data_array() const
  {
    return data;
  }

  inline ~vector1d()
  {
    data.clear();
  }

  /// Set all elements to zero
  inline void reset()
  {
    data.assign(data.size(), T(0.0));
  }

  inline size_t size() const
  {
    return data.size();
  }

  inline void resize(size_t const n)
  {
    data.resize(n);
  }

  inline void clear()
  {
    data.clear();
  }

  inline T & operator [] (size_t const i) {
    return data[i];
  }

  inline T const & operator [] (size_t const i) const {
    return data[i];
  }

  inline static void check_sizes(vector1d<T> const &v1, vector1d<T> const &v2)
  {
    if (v1.size() != v2.size()) {
      cvm::error("Error: trying to perform an operation between vectors of different sizes, "+
                 cvm::to_str(v1.size())+" and "+cvm::to_str(v2.size())+".\n");
    }
  }

  inline void operator += (vector1d<T> const &v)
  {
    check_sizes(*this, v);
    size_t i;
    for (i = 0; i < this->size(); i++) {
      (*this)[i] += v[i];
    }
  }

  inline void operator -= (vector1d<T> const &v)
  {
    check_sizes(*this, v);
    size_t i;
    for (i = 0; i < this->size(); i++) {
      (*this)[i] -= v[i];
    }
  }

  inline void operator *= (cvm::real a)
  {
    size_t i;
    for (i = 0; i < this->size(); i++) {
      (*this)[i] *= a;
    }
  }

  inline void operator /= (cvm::real a)
  {
    size_t i;
    for (i = 0; i < this->size(); i++) {
      (*this)[i] /= a;
    }
  }

  inline friend vector1d<T> operator + (vector1d<T> const &v1,
                                        vector1d<T> const &v2)
  {
    check_sizes(v1.size(), v2.size());
    vector1d<T> result(v1.size());
    size_t i;
    for (i = 0; i < v1.size(); i++) {
      result[i] = v1[i] + v2[i];
    }
    return result;
  }

  inline friend vector1d<T> operator - (vector1d<T> const &v1,
                                        vector1d<T> const &v2)
  {
    check_sizes(v1.size(), v2.size());
    vector1d<T> result(v1.size());
    size_t i;
    for (i = 0; i < v1.size(); i++) {
      result[i] = v1[i] - v2[i];
    }
    return result;
  }

  inline friend vector1d<T> operator * (vector1d<T> const &v, cvm::real a)
  {
    vector1d<T> result(v.size());
    size_t i;
    for (i = 0; i < v.size(); i++) {
      result[i] = v[i] * a;
    }
    return result;
  }

  inline friend vector1d<T> operator * (cvm::real a, vector1d<T> const &v)
  {
    return v * a;
  }

  inline friend vector1d<T> operator / (vector1d<T> const &v, cvm::real a)
  {
    vector1d<T> result(v.size());
    size_t i;
    for (i = 0; i < v.size(); i++) {
      result[i] = v[i] / a;
    }
    return result;
  }

  /// Inner product
  inline friend T operator * (vector1d<T> const &v1, vector1d<T> const &v2)
  {
    check_sizes(v1.size(), v2.size());
    T prod(0.0);
    size_t i;
    for (i = 0; i < v1.size(); i++) {
      prod += v1[i] * v2[i];
    }
    return prod;
  }

  /// Squared norm
  inline cvm::real norm2() const
  {
    cvm::real result = 0.0;
    size_t i;
    for (i = 0; i < this->size(); i++) {
      result += (*this)[i] * (*this)[i];
    }
    return result;
  }

  inline cvm::real norm() const
  {
    return cvm::sqrt(this->norm2());
  }

  inline cvm::real sum() const
  {
    cvm::real result = 0.0;
    size_t i;
    for (i = 0; i < this->size(); i++) {
      result += (*this)[i];
    }
    return result;
  }

  /// Slicing
  inline vector1d<T> const slice(size_t const i1, size_t const i2) const
  {
    if ((i2 < i1) || (i2 >= this->size())) {
      cvm::error("Error: trying to slice a vector using incorrect boundaries.\n");
    }
    vector1d<T> result(i2 - i1);
    size_t i;
    for (i = 0; i < (i2 - i1); i++) {
      result[i] = (*this)[i1+i];
    }
    return result;
  }

  /// Assign a vector to a slice of this vector
  inline void sliceassign(size_t const i1, size_t const i2,
                          vector1d<T> const &v)
  {
    if ((i2 < i1) || (i2 >= this->size())) {
      cvm::error("Error: trying to slice a vector using incorrect boundaries.\n");
    }
    size_t i;
    for (i = 0; i < (i2 - i1); i++) {
      (*this)[i1+i] = v[i];
    }
  }

  /// Formatted output

  inline size_t output_width(size_t real_width) const
  {
    return real_width*(this->size()) + 3*(this->size()-1) + 4;
  }

  inline friend std::istream & operator >> (std::istream &is,
                                            cvm::vector1d<T> &v)
  {
    if (v.size() == 0) return is;
    std::streampos const start_pos = is.tellg();
    char sep;
    if ( !(is >> sep) || !(sep == '(') ) {
      is.clear();
      is.seekg(start_pos, std::ios::beg);
      is.setstate(std::ios::failbit);
      return is;
    }
    size_t count = 0;
    while ( (is >> v[count]) &&
            (count < (v.size()-1) ? ((is >> sep) && (sep == ',')) : true) ) {
      if (++count == v.size()) break;
    }
    if (count < v.size()) {
      is.clear();
      is.seekg(start_pos, std::ios::beg);
      is.setstate(std::ios::failbit);
    }
    return is;
  }

  inline friend std::ostream & operator << (std::ostream &os,
                                            cvm::vector1d<T> const &v)
  {
    std::streamsize const w = os.width();
    std::streamsize const p = os.precision();

    os.width(2);
    os << "( ";
    size_t i;
    for (i = 0; i < v.size()-1; i++) {
      os.width(w); os.precision(p);
      os << v[i] << " , ";
    }
    os.width(w); os.precision(p);
    os << v[v.size()-1] << " )";
    return os;
  }

  inline std::string to_simple_string() const
  {
    if (this->size() == 0) return std::string("");
    std::ostringstream os;
    os.setf(std::ios::scientific, std::ios::floatfield);
    os.precision(cvm::cv_prec);
    os << (*this)[0];
    size_t i;
    for (i = 1; i < this->size(); i++) {
      os << " " << (*this)[i];
    }
    return os.str();
  }

  inline int from_simple_string(std::string const &s)
  {
    std::stringstream stream(s);
    size_t i = 0;
    if (this->size()) {
      while ((stream >> (*this)[i]) && (i < this->size())) {
        i++;
      }
      if (i < this->size()) {
        return COLVARS_ERROR;
      }
    } else {
      T input;
      while (stream >> input) {
        if ((i % 100) == 0) {
          data.reserve(data.size()+100);
        }
        data.resize(data.size()+1);
        data[i] = input;
        i++;
      }
    }
    return COLVARS_OK;
  }

};


#endif
