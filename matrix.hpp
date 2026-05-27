#ifndef HAVE_MATRIX_H
#define HAVE_MATRIX_H

#include <iostream>
#include <vector>
#include <cstdlib>
 
class matrix
{
private:
  std::vector<double> data;
  std::vector<int>    p;
  const unsigned int  rows;
  const unsigned int  cols;

  inline unsigned int
  sub2ind(const unsigned int ir, const unsigned int jc) const
  {
    return (ir + jc * rows);
  };

  double &
  value(unsigned int irow, unsigned int jcol)
  {
    return data[sub2ind(irow, jcol)];
  };

  const double &
  value(unsigned int irow, unsigned int jcol) const
  {
    return data[sub2ind(irow, jcol)];
  };

  bool factorized;

public:
  matrix(unsigned int size)
    : rows(size)
    , cols(size)
    , factorized(false)
  {
    data.resize(rows * cols, 0.0);
  };

  matrix(unsigned int rows_, unsigned int cols_)
    : rows(rows_)
    , cols(cols_)
    , factorized(false)
  {
    data.resize(rows * cols, 0.0);
  };

  matrix(matrix const &) = default;

  matrix &
  operator=(matrix const &other)
  {
    if (this == &other)
      return *this;
    assert(rows == other.rows && cols == other.cols);
    data       = other.data;
    p          = other.p;
    factorized = other.factorized;
    return *this;
  }

  void
  randomize()
  {
    for (size_t i = 0; i < data.size(); ++i)
      data[i] = static_cast<double>(std::rand()) / static_cast<double>(RAND_MAX);
    factorized = false;
  }

  unsigned int
  get_rows() const
  {
    return rows;
  }

  unsigned int
  get_cols() const
  {
    return cols;
  }

  double &
  operator()(unsigned int irow, unsigned int jcol)
  {
    return value(irow, jcol);
  };

  const double &
  operator()(unsigned int irow, unsigned int jcol) const
  {
    return value(irow, jcol);
  };

  const double *
  get_data() const
  {
    return &(data[0]);
  };

  double *
  get_data()
  {
    return &(data[0]);
  };

  /// transposed matrix : B = A'
  matrix
  transpose() const;

  void
  solve(matrix &rhs);

  void
  factorize();
};

/// matrix x matrix product : C = A * B
matrix operator*(const matrix &A, const matrix &B);

#endif