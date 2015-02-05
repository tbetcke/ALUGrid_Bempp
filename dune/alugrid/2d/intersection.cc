#ifndef DUNE_ALU2DGRID_INTERSECTION_CC
#define DUNE_ALU2DGRID_INTERSECTION_CC

#if ! COMPILE_ALU2DGRID_INLINE
#include <config.h>
#endif

#include <dune/alugrid/2d/geometry.hh>
#include <dune/alugrid/2d/intersection.hh>

namespace Dune
{

  // Implementation of ALU2DIntersectionGeometryStorage
  // --------------------------------------------------

  template< class GridImp, class LocalGeometryImpl >
  alu2d_inline
  ALU2DIntersectionGeometryStorage< GridImp, LocalGeometryImpl >
    ::ALU2DIntersectionGeometryStorage ()
  {
    for( int i = 0; i < 4; ++i )
    {
      for( int twist = 0; twist < 2; ++twist )
      {
        if ( i < 3 )
          geoms_[ 0 ][ i ][ twist ].buildLocalGeometry( i, twist, 3 );
        geoms_[ 1 ][ i ][ twist ].buildLocalGeometry( i, twist, 4 );
      }
    }
  }



  // Implementation of ALU2dGridLevelIntersectionIterator
  // ----------------------------------------------------

  template< class GridImp >
  alu2d_inline
  int ALU2dGridLevelIntersectionIterator< GridImp >
    ::getOppositeInFather ( const int nrInChild, const int nrOfChild )
  {
    int ret = (nrInChild==0) ? (2-nrOfChild)
              : ((nrInChild-nrOfChild==2 || nrInChild-nrOfChild==0) ? -1 : 0);
    alugrid_assert (ret >= -1 && ret < 3);
    return ret;
  }


  template< class GridImp >
  alu2d_inline
  int ALU2dGridLevelIntersectionIterator< GridImp >
    ::getOppositeInChild ( const int nrInFather, const int nrOfChild )
  {
    int ret = (nrInFather==0) ? (nrOfChild+1) : ((nrInFather-nrOfChild==1) ? -1 : 0);
    alugrid_assert ( ret >= -1 && ret < 3 );
    return ret;
  }


  template< class GridImp >
  alu2d_inline
  void ALU2dGridLevelIntersectionIterator< GridImp >::addNeighboursToStack ()
  {
    alugrid_assert ( current.index_ < current.nFaces() );

    ThinelementType *neighbor = current.inside()->neighbour( current.index_ );
    alugrid_assert ( neighbor );

    IntersectionInfo info;
    if( neighbor->thinis( ThinelementType::bndel_like ) )
    {
      HBndElType *bndel = (HBndElType *)neighbor;
      if( bndel->type() != HBndElType::periodic )
        return;

      PeriodicBndElType *bndnb = ((PeriodicBndElType *)bndel)->periodic_nb;
      alugrid_assert ( bndnb && bndnb->neighbour( 0 ) && bndnb->neighbour( 0 )->thinis( ThinelementType::element_like ) );
      info.first = (HElementType *)bndnb->neighbour( 0 );
      info.second = bndnb->opposite( 0 );
    }
    else
    {
      alugrid_assert ( neighbor->thinis( ThinelementType::element_like ) );
      info.first = (HElementType *)neighbor;
      info.second = current.inside()->opposite( current.index_ );
    }
    alugrid_assert ( info.first );

    while( info.first->level() > walkLevel_ )
    {
      info.second = getOppositeInFather( info.second, info.first->childNr() );
      alugrid_assert ( (info.second >= 0) && (info.second < current.nFaces()) );
      info.first = info.first->father();
    }

    if( info.first->level() >= walkLevel_ )
    {
      nbStack_.push( info );
      return;
    }

    // why should we go up, here?
    while( info.first && (info.first->level() < walkLevel_ - 1) )
    {
      info.second = getOppositeInFather( info.second, info.first->childNr() );
      alugrid_assert ( (info.second >= 0) && (info.second < current.nFaces()) );
      info.first = info.first->father();
    }

    if( info.first )
    {
      alugrid_assert ( info.first->level() == walkLevel_ - 1 );

      const int opposite = info.second;
      for( info.first = info.first->down(); info.first; info.first = info.first->next() )
      {
        info.second = getOppositeInChild( opposite, info.first->childNr() );
        if( info.second != -1 )
          nbStack_.push( info );
      }
    }
  }


  template< class GridImp >
  alu2d_inline
  void ALU2dGridLevelIntersectionIterator< GridImp >::doIncrement ()
  {
    alugrid_assert ( current.index_ < current.nFaces() );

    this->unsetUp2Date();

    if( nbStack_.empty() )
    {
      ++current.index_;
      if( current.index_ >= current.nFaces() )
      {
        alugrid_assert ( current.index_ == current.nFaces() );
        return;
      }

      addNeighboursToStack();
      // if more then one element in stack we have non-conform intersection
      current.useOutside_ = (nbStack_.size() > 1);

      if( nbStack_.empty() )
        return current.setOutside( 0, -1 );
    }

    setupIntersection();

    alugrid_assert ( !current.outside() || (current.outside()->level() == walkLevel_) );
  }


  template< class GridImp >
  alu2d_inline
  void ALU2dGridLevelIntersectionIterator< GridImp >::setFirstItem ( const HElementType &elem, int wLevel )
  {
    // empty stack first
    while( !nbStack_.empty() )
      nbStack_.pop();

    current.setInside( const_cast< HElementType * >( &elem ) );
    current.index_ = -1;
    current.setOutside( 0, -1 );

    walkLevel_ = wLevel;

    alugrid_assert ( current.inside() );

    increment();
  }


  template< class GridImp >
  alu2d_inline
  void ALU2dGridLevelIntersectionIterator< GridImp >::setupIntersection ()
  {
    alugrid_assert ( !nbStack_.empty() );

    IntersectionInfo &info = nbStack_.top();
    current.setOutside( info.first, info.second );
    nbStack_.pop();
  }



  // Implementation of ALU2dGridLeafIntersectionIterator
  // ---------------------------------------------------

  //! increment iterator
  template< class GridImp >
  alu2d_inline
  void ALU2dGridLeafIntersectionIterator< GridImp >::doIncrement ()
  {
    alugrid_assert ( current.index_ < current.nFaces() );

    this->unsetUp2Date();

    // do we still have neighbours?
    if( !nbStack_.empty() )
      return setupIntersection();

    ++current.index_;
    if( current.index_ >= current.nFaces())
    {
      alugrid_assert ( current.index_ == current.nFaces() );
      return;
    }

    ThinelementType *neighbor = current.inside()->neighbour( current.index_ );
    alugrid_assert ( neighbor );

    if( neighbor->thinis( ThinelementType::bndel_like ) )
    {
      HBndElType *bndel = (HBndElType *)neighbor;
      if( bndel->type() != HBndElType::periodic )
      {
        current.useOutside_ = false;
        return current.setOutside( 0, -1 );
      }

      PeriodicBndElType *bndnb = ((PeriodicBndElType *)bndel)->periodic_nb;
      alugrid_assert ( bndnb && bndnb->neighbour( 0 ) && bndnb->neighbour( 0 )->thinis( ThinelementType::element_like ) );
      current.useOutside_ = !bndnb->leaf();
      if( current.useOutside_ )
      {
        IntersectionInfo info;

        // insert left intersection
        HBndElType *left = bndnb->down();
        alugrid_assert ( left && left->leaf() );
        alugrid_assert ( left->neighbour( 0 ) && left->neighbour( 0 )->thinis( ThinelementType::element_like ) );
        info.first = (HElementType *)left->neighbour( 0 );
        info.second = left->opposite( 0 );
        nbStack_.push( info );

        HBndElType *right = left->next();
        alugrid_assert ( right && right->leaf() );
        alugrid_assert ( right->neighbour( 0 ) && right->neighbour( 0 )->thinis( ThinelementType::element_like ) );
        info.first = (HElementType *)right->neighbour( 0 );
        info.second = right->opposite( 0 );
        nbStack_.push( info );

        setupIntersection();
      }
      else
        current.setOutside( (HElementType *)bndnb->neighbour( 0 ), bndnb->opposite( 0 ) );
    }
    else
    {
      current.useOutside_ = current.inside()->hasHangingNode( current.index_ );
      const int opposite = current.inside()->opposite( current.index_ );
      if( current.useOutside_ )
      {
        IntersectionInfo info;

        // insert left intersection
        ThinelementType *left = current.inside()->getLeftIntersection( current.index_ );
        alugrid_assert ( left && left->thinis( ThinelementType::element_like ) );
        info.first = (HElementType *)left;   // neighbor
        info.second = opposite;              // opposite vertex
        alugrid_assert ( info.first->leaf() );
        nbStack_.push( info );

        // insert right intersection
        ThinelementType *right = current.inside()->getRightIntersection( current.index_ );
        alugrid_assert ( right && right->thinis( ThinelementType::element_like ) );
        info.first = (HElementType *)right;  // neighbor
        info.second = opposite;              // opposite vertex
        alugrid_assert ( info.first->leaf() );
        nbStack_.push( info );

        setupIntersection();
      }
      else
        current.setOutside( (HElementType *)current.inside()->neighbour( current.index_ ), opposite );
    }
  }


  template< class GridImp >
  alu2d_inline
  void ALU2dGridLeafIntersectionIterator< GridImp >::setFirstItem ( const HElementType &elem, int wLevel )
  {
    while( !nbStack_.empty() )
      nbStack_.pop();

    current.setInside( const_cast< HElementType * >( &elem ) );
    current.index_ = -1;
    alugrid_assert ( current.inside() );
    walkLevel_ = wLevel;
    increment();
  }


  template< class GridImp >
  alu2d_inline
  void ALU2dGridLeafIntersectionIterator< GridImp >::setupIntersection ()
  {
    alugrid_assert ( !nbStack_.empty() );

    IntersectionInfo &info = nbStack_.top();
    current.setOutside( info.first, info.second );
    nbStack_.pop();
  }

}


#if COMPILE_ALU2DGRID_LIB
// Template Instantiation
// ----------------------

template class Dune::ALU2DIntersectionGeometryStorage
  < Dune::ALU2dGridGeometry< 1, 2, const Dune::ALU2dGrid< 2, 2, ALU2DSPACE triangle > > >;
template class Dune::ALU2DIntersectionGeometryStorage
  < Dune::ALU2dGridGeometry< 1, 2, const Dune::ALU2dGrid< 2, 2, ALU2DSPACE quadrilateral > > >;

template class Dune::ALU2DIntersectionGeometryStorage
  < Dune::ALU2dGridGeometry< 1, 2, const Dune::ALU2dGrid< 2, 3, ALU2DSPACE triangle > > >;
template class Dune::ALU2DIntersectionGeometryStorage
  < Dune::ALU2dGridGeometry< 1, 2, const Dune::ALU2dGrid< 2, 3, ALU2DSPACE quadrilateral > > >;


template class Dune::ALU2dGridLevelIntersectionIterator< const Dune::ALU2dGrid< 2, 2, ALU2DSPACE triangle > >;
template class Dune::ALU2dGridLevelIntersectionIterator< const Dune::ALU2dGrid< 2, 2, ALU2DSPACE quadrilateral > >;

template class Dune::ALU2dGridLevelIntersectionIterator< const Dune::ALU2dGrid< 2, 3, ALU2DSPACE triangle > >;
template class Dune::ALU2dGridLevelIntersectionIterator< const Dune::ALU2dGrid< 2, 3, ALU2DSPACE quadrilateral > >;


template class Dune::ALU2dGridLeafIntersectionIterator< const Dune::ALU2dGrid< 2, 2, ALU2DSPACE triangle > >;
template class Dune::ALU2dGridLeafIntersectionIterator< const Dune::ALU2dGrid< 2, 2, ALU2DSPACE quadrilateral > >;

template class Dune::ALU2dGridLeafIntersectionIterator< const Dune::ALU2dGrid< 2, 3, ALU2DSPACE triangle > >;
template class Dune::ALU2dGridLeafIntersectionIterator< const Dune::ALU2dGrid< 2, 3, ALU2DSPACE quadrilateral > >;

#endif
#endif
