#ifndef DUNE_ALU2D_BNDPROJECTION_HH
#define DUNE_ALU2D_BNDPROJECTION_HH

#include <dune/alugrid/common/bndprojection.hh>

#include <dune/alugrid/2d/alu2dinclude.hh>

namespace Dune
{

  template< class Grid >
  class ALU2dGridBoundaryProjection
  : public ALU2DSPACE VtxProjection< Grid::dimensionworld,(Grid::elementType == ALU2DSPACE triangle ? 3 : 4) >
  {
    typedef ALU2DSPACE VtxProjection< Grid::dimensionworld,(Grid::elementType == ALU2DSPACE triangle ? 3 : 4) > Base;

  public:
    enum { ncoord = Base::ncoord };

    typedef typename Base::hbndel_t hbndel_t;
    typedef typename Base::helement_t helement_t;

    typedef typename Grid::DuneBoundaryProjectionType DuneBoundaryProjectionType;

    typedef typename DuneBoundaryProjectionType::CoordinateType CoordinateType;

    explicit ALU2dGridBoundaryProjection ( const Grid &grid )
    : grid_( grid )
    {}

    int operator() ( const hbndel_t *hbndel, const double local, double (&global)[ ncoord ] ) const
    {
      return callProjection( grid_.boundaryProjection( hbndel->segmentIndex() ), global );
    }

    int operator() ( const helement_t *helement, const double (&local)[ 2 ], double (&global)[ ncoord ] ) const
    {
      return callProjection( grid_.globalProjection(), global );
    }

  private:
    static int callProjection ( const DuneBoundaryProjectionType *prj, double (&global)[ ncoord ] )
    {
      if( prj )
      {
        CoordinateType x, y;
        for( int i = 0; i < ncoord; ++i )
          x[ i ] = global[ i ];
        y = (*prj)( x );
        for( int i = 0; i < ncoord; ++i )
          global[ i ] = y[ i ];
      }
      return 1;
    }

    const Grid &grid_;
  };

} // end namespace Dune

#endif // #ifndef DUNE_ALU2D_BNDPROJECTION_HH
