/**
 * 	Provides functionality to handle matrix calculations.
 *  A matrix can have real or polynomial coefficients.
 **/
#pragma once

#ifdef _MSC_VER
  using namespace boost;
  #pragma warning( push )
  #pragma warning(disable: 4996)
#else
  #include <unordered_map>
#endif
#include "Polynomial.h"
#include "Complex.h"

#if defined(JDEBUG)
#include <iostream>
#endif

namespace LibHUPF
{

class Matrix
{
public:
  int numRows; //!< number of rows
  int numCols; //!< number of columns
  std::vector<std::vector<double> > m_matrix; //!< real coefficient vector
  std::vector<std::vector<Polynomial> > p_matrix; //!< polynomial coefficient vector
  bool hasPolyElements; //!< true if matrix has polynomial coefficients, false if matrix has real coefficients

  /**
   * Creates empty matrix polynomial coefficient type
   * */
  Matrix()
  {
    hasPolyElements=false;
  };

  /**
   * Copy constructor
   * */
  Matrix(const Matrix &m)
  {
    numCols = m.numCols;
    numRows = m.numRows;
    m_matrix = m.m_matrix;
    p_matrix = m.p_matrix;
    hasPolyElements = m.hasPolyElements;
  };

  /**
  *Square diagonal matrix of width d.size() and diagonal entries d
  **/
  Matrix (std::vector<double> d)
  {
    numRows = int(d.size());
    numCols = int(d.size());
    hasPolyElements=false;

    std::vector<double> tmp;

    tmp.assign(numCols,0.0);
    m_matrix.assign(numRows,tmp);

    for (int i=0; i<numRows; i++)
    {
      m_matrix[i][i]=d[i];
    };
  };


  /**
   * Creates polynomial coefficient or numerical type matrix with n rows and m columns.
   * @param num number of rows
   * @param cols number of columns
   * @param poly indicates polynomial coefficient
   * */
  Matrix(int num, int cols, bool poly = false)
  {
    numRows = num;
    numCols = cols;
    
    if (poly)
    {
      std::vector<Polynomial> tmp;
      tmp.assign(cols,Polynomial());
      p_matrix.assign(num,tmp);
      hasPolyElements=true;
    }
    else
    {
      std::vector<double> tmp;
      tmp.assign(cols,0.0);
      m_matrix.assign(num,tmp);
      hasPolyElements=false;
    }
  };

  /**
   * Get real coefficient at position (n,m)
   * @param n row
   * @param m column
   * @return coefficient
   * */
  double get(int n, int m)
  {
    return m_matrix[n][m];
  };

  /**
   * Set real coefficient at position (n,m)
   * @param n row
   * @param m column
   * @param v value
   * */
  void set(int m, int n, double v)
  {
    m_matrix[m][n] = v;
  };

  /**
  *Transpose real coefficient Matrix (n,m)
  **/

  Matrix Transpose()
  {
    Matrix TM(numCols,numRows);
    for (int i=0; i<numRows; i++)
    {
      for(int j=0; j<numCols; j++)
      {
        TM.m_matrix[j][i]=m_matrix[i][j];
      };
    };
    return TM;
  };

  //inverts an invertible matrix (assuming diagonals are not 0). Achtung: no checking if matrix is invertible!
  Matrix Inverse()
  {
    int cols = 2*numRows;
    Matrix tmp(numRows, cols);
    int i,j,k;

    for(i=0;i<numRows;++i)
    {
      tmp.m_matrix[i][numRows+i] = 1.0;
      for(j=0;j<numRows;++j)
        tmp.m_matrix[i][j] = m_matrix[i][j];
    }

    double ratio = 0.0;
    double a = 1.0;

    for(i=0;i<numRows;++i)
    {
      for(j=0;j<numRows;++j)
      {
        if(i!=j)
        {
          ratio = tmp.m_matrix[j][i]/tmp.m_matrix[i][i];
          for(k=0;k<cols;++k) tmp.m_matrix[j][k] -= ratio * tmp.m_matrix[i][k];
        }
      }
    }

    Matrix out(numRows, numRows);
    for(i=0;i<numRows;++i)
    {
      a = tmp.m_matrix[i][i];
      for(j=0;j<numRows;++j)
      {
        out.m_matrix[i][j] = tmp.m_matrix[i][j+numRows]/a;
      }
    }

    return out;
  }

  static double FrobeniusNorm(const Matrix &m)
  {
    if (m.hasPolyElements) return 0;

    double out = 0;
    for (int r=0; r<m.numRows; r++)
      for (int c=0; c<m.numCols; c++)
        out+=m.m_matrix[r][c]*m.m_matrix[r][c];
    return sqrt(out);
  }

  // function for exchanging two rows of a matrix
  void Swap(int row1, int row2)
  {
      for (int i = 0; i < numCols; i++)
      {
          double temp = m_matrix[row1][i];
          m_matrix[row1][i] = m_matrix[row2][i];
          m_matrix[row2][i] = temp;
      }
  }

  // function for finding rank of a numeric matrix
  static int Rank(const Matrix& M)
  {
    int rank = M.numCols;
    Matrix m(M);

    for (int row = 0; row < rank; row++)
    {
      //rank by Gaussian elimination
      if (fabs(m.m_matrix[row][row])>2.2204460492503131e-16)
      {
        for (int col = 0; col < M.numRows; col++)
        {
          if (col != row)
          {
           double mult = (double)m.m_matrix[col][row]/m.m_matrix[row][row];
           for (int i = 0; i < rank; i++)
             m.m_matrix[col][i] -= mult * m.m_matrix[row][i];
          }
        }
      }
      else
      {
        m.m_matrix[row][row] = 0;
        bool reduce = true;
        // Find the non-zero element in current column
        for (int i = row + 1; i < M.numRows;  i++)
        {
          // Swap the row with non-zero element with this row.
          if (fabs(m.m_matrix[i][row])>2.2204460492503131e-16)
          {
            m.Swap(row, i);
            reduce = false;
            break ;
          }
          else
            m.m_matrix[i][row] = 0;
        }

        // If we did not find any row with non-zero element in current columnm, then all
        // values in this column are 0.
        if (reduce)
        {
          rank--;
          for (int i = 0; i < M.numRows; i++)
            m.m_matrix[i][row] = m.m_matrix[i][rank];
        }
        row--;
      }

    }
    return rank;
  }

  //functon that substitutes the values of a matrix with poly elements
  //Important: matrix should be a polynomial matrix in 2 variables
  Matrix Sub(double va, double vb)
  {
    Matrix m(numRows, numCols);
    if (!hasPolyElements) return m;
    for (int i=0;i<numRows;++i)
      for (int j=0;j<numCols;++j)
      {
        if (p_matrix[i][j].realCoeffType)
          m.m_matrix[i][j] = p_matrix[i][j].eval(vb);
        else m.m_matrix[i][j] = p_matrix[i][j].subs(va).eval(vb);
      }
    return m;
  }

  /**
   * original determinant calculation. Only used if matrix size > 10x10
   * @param m Matrix
   * @return determinant
   * */
  static Polynomial myDeterminant_orig(const Matrix &m)
  {
    using namespace std;

    Polynomial res;
    int n=m.numCols;
    string key;
    string temp;
    //boost::unordered_map<string, Polynomial> minors;
    unordered_map<string, Polynomial> minors;
    //boost::unordered_map<string, Polynomial> newDet;
    unordered_map<string, Polynomial> newDet;
    for(int i=0; i<n; i++)
    {
      temp.clear();
      char tmp1[15];
      sprintf(tmp1,"%i,",i);
      temp.append(tmp1);
      minors.insert(pair<string,Polynomial>(temp,m.p_matrix[i][n-1]));
    }
    for(int c=n-2; c>=0; c--)
    {
      int minorSize = n-c;
      vector<int> minorKey;
      key.clear();
      for(int i=0; i<minorSize; i++)
      {
        minorKey.push_back(i);
        char keytmp[15];
        sprintf(keytmp,"%i,",i);
        key=key+keytmp;
      }
      int fc=0;
      do
      {
        Polynomial det;
        for(int r=0; r<minorSize; r++)
        {
          temp.clear();
          for(int i=0; i<minorSize; i++)
            if(i!=r)
            {
              char tmp1[5];
              sprintf(tmp1,"%i,",minorKey[i]);
              temp=temp+tmp1;
            }


          if(r%2 != 0)
          {
            //det-=m.p_matrix[minorKey[r]][c]*minors[temp];
            det-=minors[temp]*m.p_matrix[minorKey[r]][c];
          }
          else
          {
            //det+=m.p_matrix[minorKey[r]][c]*minors[temp];
            det+=minors[temp]*m.p_matrix[minorKey[r]][c];
          }

        }
        newDet.insert(pair<string,Polynomial>(key,det));
        for(fc=minorSize; fc>0; --fc)
        {
          ++minorKey[fc-1];
          if(minorKey[fc-1] < fc+c)
            break;
        }
        if(fc < minorSize && fc > 0)
        {
          for(int j=fc; j<minorSize; ++j)
            minorKey[j] = minorKey[j-1]+1;
        }
        key.clear();
        for(size_t i=0; i<minorKey.size(); i++)
        {
          char x[15];
          sprintf(x,"%i,",minorKey[i]);
          key=key+x;
        }
      }
      while (fc!=0);
      //minors.clear();
      minors = newDet;
      if(newDet.size() == 1)
      {
        unordered_map<string,Polynomial>::iterator it=newDet.begin();
        res = (*it).second;
      }
      newDet.clear();
    }
    return res;
  };

  /**
   * determinant calculation. Uses laplace expansion to calculate determinant.
   * Only used if matrix size <= 10x10
   * @param m Matrix
   * @return determinant
   * */
  static Polynomial myDeterminant(const Matrix &m)
  {
    using namespace std;
    Polynomial res;
    int POW[] = {1,10,100,1000,10000,100000,1000000,10000000,100000000,1000000000};
    int n=m.numCols;
    if(n>10)
      return myDeterminant_orig(m);
    int key=0;
    unordered_map<int, Polynomial> minors;
    unordered_map<int, Polynomial> newDet;
    for(int i=0; i<n; i++)
    {
      minors[i]=m.p_matrix[i][n-1];
    }
    for(int c=n-2; c>=0; c--)
    {
      int minorSize = n-c;
      vector<int> minorKey;
      minorKey.reserve(minorSize);
      key=0;
      for(int i=0; i<minorSize; i++)
      {
        minorKey.push_back(i);
        key+=POW[i]*i;
      }
      int fc=0;
      do
      {
        Polynomial det;
        for(int r=0; r<minorSize; r++)
        {
          long temp2=0;
          bool iwasr=false;
          for(int i=0; i<minorSize; i++)
          {
            if(i!=r)
            {
              if(iwasr)
              {
                temp2+=POW[i-1]*minorKey[i];
              }
              else
              {
                temp2+=POW[i]*minorKey[i];
              }
            }
            else
            {
              iwasr=true;
            }
          }

          if(r%2 != 0)
          {
            //det-=(m.p_matrix[minorKey[r]][c]*minors[temp2]);
            det-=(minors[temp2]*m.p_matrix[minorKey[r]][c]);
          }
          else
          {
            //det+=(m.p_matrix[minorKey[r]][c]*minors[temp2]);
            det+=(minors[temp2]*m.p_matrix[minorKey[r]][c]);
          }
        }
        newDet[key]=det;
        for(fc=minorSize; fc>0; --fc)
        {
          ++minorKey[fc-1];
          if(minorKey[fc-1] < fc+c)
            break;
        }
        if(fc < minorSize && fc > 0)
        {
          for(int j=fc; j<minorSize; ++j)
            minorKey[j] = minorKey[j-1]+1;
        }
        key=0;
        int count=0;
        for(size_t i=0; i<minorKey.size(); i++)
        {
          key+=POW[count++]*minorKey[i];
        }
      }
      while (fc!=0);

      if(newDet.size() == 1)
      {
        res=newDet.begin()->second;
      }
      minors.swap(newDet);
      newDet.clear();
    }
    return res;
  };

  /**
   * Trys to calculate determinant with upper triangular matrix.
   * Needs Polynomial input and only works if polynomials are
   * divisible. --> Not used in our case, but left for future
   * implementation work.
   **/
  static Polynomial determinantRecursive(Matrix &mat)
  {
    Polynomial subtract_me;
    Polynomial res;
    if(mat.numRows != mat.numCols)
      perror("Only for square matrices.");
    if(mat.numRows==1 && mat.numCols==1)
      return mat.p_matrix[0][0];
    if(mat.numRows==1 && mat.numCols==2)
      return mat.p_matrix[0][0]*mat.p_matrix[1][1] - mat.p_matrix[0][1]*mat.p_matrix[1][0];
    for(int i = 1; i<mat.numRows; i++)
    {
      for(int j = 1; j<mat.numRows; j++)
      {
        subtract_me = mat.p_matrix[i][0] * mat.p_matrix[0][j];
        mat.p_matrix[i][j] = mat.p_matrix[i][j] * mat.p_matrix[0][0] - subtract_me;
      }
    }
    Matrix subMat(mat.numRows-1,mat.numCols-1);
    for(int i=0; i<subMat.numRows; i++)
      for(int j=0; j<subMat.numCols; j++)
        subMat.p_matrix[i][j] = mat.p_matrix[i+1][j+1];
    res = determinantRecursive(subMat);
    for(int i=1; i<=mat.numRows-2; i++)
    {
      res = res / mat.p_matrix[0][0];
    }
    return res;
  };

  /**
   * Calculates determinant of double coefficient matrix.
   * Works with upper triangle matrix. --> Not used in our case
   * (no performance improvement), but left for future implementation work.
   **/
  static double determinantRecursiveDouble(Matrix &mat)
  {
    double subtract_me;
    double res;
    if(mat.numRows != mat.numCols)
      perror("Only for square matrices.");
    if(mat.numRows==1 && mat.numCols==1)
      return mat.m_matrix[0][0];
    if(mat.numRows==1 && mat.numCols==2)
      return mat.m_matrix[0][0]*mat.m_matrix[1][1] - mat.m_matrix[0][1]*mat.m_matrix[1][0];
    for(int i = 1; i<mat.numRows; i++)
    {
      for(int j = 1; j<mat.numRows; j++)
      {
        subtract_me = mat.m_matrix[i][0] * mat.m_matrix[0][j];
        mat.m_matrix[i][j] = mat.m_matrix[i][j] * mat.m_matrix[0][0] - subtract_me;
      }
    }
    Matrix subMat(mat.numRows-1,mat.numCols-1);
    for(int i=0; i<subMat.numRows; i++)
      for(int j=0; j<subMat.numCols; j++)
        subMat.m_matrix[i][j] = mat.m_matrix[i+1][j+1];
    res = determinantRecursiveDouble(subMat);
    for(int i=1; i<=mat.numRows-2; i++)
    {
      res = res / mat.m_matrix[0][0];
    }
    return res;
  };

  static Matrix companionPoly(const Polynomial &a)
  {
    Matrix res(a.degree,a.degree);
    for (int r = 1; r < res.numRows; r++)
      res.m_matrix[r][r-1]=1;
    for (int r = 0; r < res.numRows; r++)
      res.m_matrix[r][res.numCols-1] = -a.coefficient[r] / a.coefficient[a.degree];
    return res;
  };

  static int sizeMatrix(Matrix &m)
  {
    if(m.numCols==m.numRows)
      return m.numCols;
    perror("Size property is valid only for square matrix");
    exit(EXIT_FAILURE);
  };

  /**
   * balance the companion matrix of a polynomial.
   * @param a companion matrix
   * @return balanced companion matrix
   * */
  static Matrix balanceCompanionMatrix(Matrix a)
  {
    int notConverged=1;
    double rowNorm=0;
    double columnNorm=0;
    while (notConverged)
    {
      int i, j;
      double g, f, s;
      int size = a.numCols;
      int radix = 2;
      int radix2 = 4;
      notConverged = 0;
      for (i = 0; i < size; i++)
      {
        // Column Norm, excluding the diagonal
        if (i != size - 1)
          columnNorm = fabs(a.m_matrix[i + 1][i]);
        else
        {
          columnNorm = 0;
          for (j = 0; j < size - 1; j++)
            columnNorm += fabs(a.m_matrix[j][size - 1]);
        }

        // Row Norm, exluding the diagonal
        if (i == 0)
          rowNorm = fabs(a.m_matrix[0][size-1]);
        else if (i == size - 1)
          rowNorm = fabs(a.m_matrix[i][i - 1]);
        else
          rowNorm = fabs(a.m_matrix[i][i - 1]) + fabs(a.m_matrix[i][size - 1]);

        if (columnNorm == 0 || rowNorm == 0)
          continue;

        g = rowNorm / radix;
        f = 1;
        s = rowNorm + columnNorm;

        while (columnNorm < g)
        {
          f *= radix;
          columnNorm *= radix2;
        }
        while (columnNorm > g)
        {
          f /= radix;
          columnNorm /= radix2;
        }
        if ((rowNorm + columnNorm) < 0.95 * s * f)
        {
          notConverged = 1;
          g = 1 / f;
          if (i == 0)
            a.m_matrix[0][size - 1] *= g;
          else
          {
            a.m_matrix[i][i - 1] *= g;
            a.m_matrix[i][size - 1] *= g;
          }
          if (i == size - 1)
            for (j = 0; j < size; j++)
              a.m_matrix[j][i] *= f;
          else
            a.m_matrix[i + 1][i] *= f;

        }
      }
    }
    return a;
  };

  /**
   * this method calculates all eigenvalues of a given matrix.
   * @param a Matrix
   * @return vector<Complex> which contains all eigenvalues
   * */
  static std::vector<Complex> qrCompanion(Matrix a)
  {
    using namespace std;
    vector<Complex> result;
    result.resize(a.numCols);
    /*
    for(int i=0; i<sizeMatrix(a); i++)
    	result.push_back(Complex());*/
    double t = 0.0;
    int iterations = 0, e, i, j, k, m;
    double w, x, y, s, z, a1, a2, a3;
    double epsilon = 2.2204460492503131e-16;
    double p = 0, q = 0, r = 0;

    int notLast;
    int n = sizeMatrix(a);

    int goToNextIteration = 0;
    while (1)
    {
      if (!goToNextIteration)
      {
        // next root
        if (n == 0)
          break;
        iterations = 0;
      }
      goToNextIteration = 0;
      // next iteration
      for (e = n; e >= 2; e--)
      {
        a1 = fabs(a.m_matrix[e - 1][e - 2]);
        a2 = fabs(a.m_matrix[e - 2][e - 2]);
        a3 = fabs(a.m_matrix[e - 1][e - 1]);

        if (a1 <= epsilon * (a2 + a3))
          break;
      }
      x = a.m_matrix[n - 1][n - 1];

      if (e == n)
      {
        result[n - 1].real=x + t;
        n--;
        continue;
      }
      y = a.m_matrix[n - 2][n - 2];
      w = a.m_matrix[n - 2][n - 1] * a.m_matrix[n - 1][n - 2];

      if (e == (n - 1))
      {
        p = (y - x) / 2;
        q = p * p + w;
        y = sqrt(fabs(q));

        x += t;

        if (q > 0)  //two real roots
        {
          if (p < 0)
            y = -y;
          y += p;
          result[n - 1].real=x - w / y;
          result[n - 2].real=x + y;
        }
        else
        {
          result[n - 1].real=x + p;
          result[n-1].imaginary= -y;

          result[n - 2].real=x + p;
          result[n-2].imaginary=y;
        }
        n -= 2;

        continue;
      }

      if (iterations == 60)
      {
        perror("Maximum no. of iterations reached");
        exit(EXIT_FAILURE);
      }

      if (iterations % 10 == 0 && iterations > 0)
      {
        // use an exceptional shift
        t += x;

        for (i = 1; i <= n; i++)
          a.m_matrix[i - 1][i - 1] -= x;
        s = fabs(a.m_matrix[n - 1][n - 2]) + fabs(a.m_matrix[n - 2][n - 3]);
        y = 0.75 * s;
        x = y;
        w = -0.4375 * s * s;
      }

      iterations++;

      for (m = n - 2; m >= e; m--)
      {
        double a1, a2, a3;

        z = a.m_matrix[m - 1][m - 1];
        // TODO: FMAT shift less 1
        r = x - z;
        s = y - z;
        p = a.m_matrix[m - 1][m] + (r * s - w) / a.m_matrix[m][m - 1];
        q = a.m_matrix[m][m] - z - r - s;
        r = a.m_matrix[m + 1][m];
        s = fabs(p) + fabs(q) + fabs(r);
        p /= s;
        q /= s;
        r /= s;

        if (m == e)
          break;

        a1 = fabs(a.m_matrix[m - 1][m - 2]);
        a2 = fabs(a.m_matrix[m - 2][m - 2]);
        a3 = fabs(a.m_matrix[m][m]);

        if (a1 * (fabs(q) + fabs(r)) <= epsilon * fabs(p) * (a2 + a3))
          break;
      }

      for (i = m + 2; i <= n; i++)
        a.m_matrix[i - 1][i - 3] = 0;

      for (i = m + 3; i <= n; i++)
        a.m_matrix[i - 1][i - 4] = 0;

      /* double QR step */

      for (k = m; k <= n - 1; k++)
      {
        notLast = (k != n - 1);

        if (k != m)
        {
          p = a.m_matrix[k - 1][k - 2];
          q = a.m_matrix[k][k - 2];
          r = notLast ? a.m_matrix[k + 1][k - 2] : 0.0;

          x = fabs(p) + fabs(q) + fabs(r);

          if (x == 0)
            continue;           /* FIXME????? */

          p /= x;
          q /= x;
          r /= x;
        }

        s = sqrt(p * p + q * q + r * r);

        if (p < 0)
          s = -s;

        if (k != m)
          a.m_matrix[k - 1][k - 2] = -s * x;

        else if (e != m)
          a.m_matrix[k - 1][k - 2] *= -1;


        p += s;
        x = p / s;
        y = q / s;
        z = r / s;
        q /= p;
        r /= p;

        /* do row modifications */

        for (j = k; j <= n; j++)
        {
          p = a.m_matrix[k - 1][j - 1] + q * a.m_matrix[k][j - 1];

          if (notLast)
          {
            p += r * a.m_matrix[k + 1][j - 1];
            a.m_matrix[k + 1][j - 1] -= p * z;
          }

          a.m_matrix[k][j - 1] -= p * y;
          a.m_matrix[k - 1][j - 1] -= p * x;
        }

        j = (k + 3 < n) ? (k + 3) : n;

        for (i = e; i <= j; i++)
        {
          p = x * a.m_matrix[i - 1][k - 1] + y * a.m_matrix[i - 1][k];

          if (notLast)
          {
            p += z * a.m_matrix[i - 1][k + 1];
            a.m_matrix[i - 1][k + 1] -= p * r;
          }
          a.m_matrix[i - 1][k] -= p * q;
          a.m_matrix[i - 1][k - 1] -= p;
        }
      }
      goToNextIteration = 1;
    }
    return result;
  };

  /**
   * multiplication operator (Matrix*Matrix)
   * */
  Matrix operator*(const Matrix &mat)
  {
    if (numCols==mat.numRows)
    {
      Matrix res;
      if (hasPolyElements && mat.hasPolyElements)
        res = Matrix(numRows,mat.numCols);
      else
        res = Matrix(numRows,mat.numCols);
      if (!res.hasPolyElements)
      {
        for(int r3=0; r3 < res.numRows; r3++)
          for (int c3=0; c3< res.numCols; c3++)
            for (int c1=0; c1<numCols; c1++)
              res.m_matrix[r3][c3]+=m_matrix[r3][c1]*mat.m_matrix[c1][c3];
      }
      else
      {
        for(int r3=0; r3< res.numRows; r3++)
          for(int c3=0; c3< res.numCols; c3++)
            for(int c1=0; c1<numCols; c1++)
              res.p_matrix[r3][c3]+=(p_matrix[r3][c1]*mat.p_matrix[c1][c3]);
      }
      return res;
    }
    else
    {
      perror("Matrices can not be multiplied. Check number of Rows & Columns");
      exit(EXIT_FAILURE);
    }
  };

  /**
   * minus operator (Matrix-Matrix)
   * */
  Matrix operator-(const Matrix &mat)
  {
    if (numCols==mat.numCols && numRows==mat.numRows)
    {
      Matrix res;
      if (hasPolyElements && mat.hasPolyElements)
        res = Matrix(numRows,numCols);
      else
        res = Matrix(numRows,numCols);

      if (!res.hasPolyElements)
      {
        for (int r=0; r<res.numRows; r++)
          for(int c=0; c<res.numCols; c++)
            res.m_matrix[r][c]=m_matrix[r][c]-mat.m_matrix[r][c];
      }
      else
      {
        for (int r=0; r<res.numRows; r++)
          for (int c=0; c<res.numCols; c++)
            res.p_matrix[r][c]=(p_matrix[r][c]-mat.p_matrix[r][c]);
      }
      return res;
    }
    else
    {
      perror("Matrix Rows & Columns do not match");
      exit(EXIT_FAILURE);
    }
  };

  /**
   * plus operator (Matrix+Matrix)
   * */
  Matrix operator+(const Matrix &mat)
  {
    if (numCols==mat.numCols && numRows==mat.numRows)
    {
      Matrix *res;

      if (hasPolyElements && mat.hasPolyElements)
        res=new Matrix(numRows,numCols);
      else
        res=new Matrix(numRows,numCols);

      if (!res->hasPolyElements)
      {
        for (int r=0; r<res->numRows; r++)
          for(int c=0; c<res->numCols; c++)
            res->m_matrix[r][c]=m_matrix[r][c]+mat.m_matrix[r][c];
      }
      else
      {
        for (int r=0; r<res->numRows; r++)
          for (int c=0; c<res->numCols; c++)
            res->p_matrix[r][c]=(p_matrix[r][c]+mat.p_matrix[r][c]);
      }
      return *res;
    }
    else
    {
      perror("Matrix Rows & Columns do not match");
      exit(EXIT_FAILURE);
    }
  };

  /**
   * Scalar Multiplication operator (Scalar * Matrix)
   * */
  Matrix operator*(double dbl)
  {

    Matrix *res;

    if (hasPolyElements)
      res=new Matrix(numRows,numCols);
    else
      res=new Matrix(numRows,numCols);

    if (!res->hasPolyElements)
    {
      for (int r=0; r<res->numRows; r++)
        for(int c=0; c<res->numCols; c++)
          res->m_matrix[r][c]=m_matrix[r][c]*dbl;
    }
    else
    {
      for (int r=0; r<res->numRows; r++)
        for (int c=0; c<res->numCols; c++)
          res->p_matrix[r][c]=p_matrix[r][c]*dbl;
    }
    return *res;
  };
};

#ifdef JDEBUG
inline std::ostream& operator<<(std::ostream& of, const Matrix &m)
{
  using namespace std;
  for (int i=0;i<m.numRows;++i)
  {
    for (int j=0;j<m.numCols;++j)
      of << m.m_matrix[i][j] << " ";
    of << endl;
  }
  return of;
}
#endif

}

#if defined(_WIN32) || defined(_WIN64)
  #pragma warning(pop)
#endif

