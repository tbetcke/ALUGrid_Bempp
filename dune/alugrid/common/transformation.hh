#ifndef DUNE_ALUGRID_TRANSFORMATION_HH
#define DUNE_ALUGRID_TRANSFORMATION_HH

#include <dune/common/fvector.hh>
#include <dune/common/fmatrix.hh>

namespace Dune
{

  template< class ctype, int dimw >
  struct ALUGridTransformation
  {
    static const int dimension = dimw;

    typedef FieldVector< ctype, dimension > WorldVector;
    typedef FieldMatrix< ctype, dimension, dimension > WorldMatrix;

    ALUGridTransformation ( const WorldMatrix &matrix, const WorldVector &shift )
    : matrix_( matrix ),
      shift_( shift )
    {}

    WorldVector evaluate ( const WorldVector &x ) const
    {
      WorldVector y = shift_;
      matrix_.umv( x, y );
      return y;
    }

    WorldVector evaluateInverse ( const WorldVector &y ) const
    {
      // Note: We assume the matrix to be orthogonal, here
      WorldVector ys = y - shift_;
      WorldVector x;
      matrix_.mtv( ys, x );
      return x;
    }

  private:
    WorldMatrix matrix_;
    WorldVector shift_;
  };

}

#endif // #ifndef DUNE_ALUGRID_TRANSFORMATION_HH
