#ifndef DUNE_ALU2DGRID_INCLUDE_HH
#define DUNE_ALU2DGRID_INCLUDE_HH

#include <dune/common/parallel/collectivecommunication.hh>

#include <dune/alugrid/impl/indexstack.h>
#include <dune/alugrid/impl/projectvertex.h>

#include <dune/alugrid/impl/2d/grid.h>
#include <dune/alugrid/impl/2d/triang.h>
#include <dune/alugrid/impl/2d/listwalk.h>
#include <dune/alugrid/impl/2d/handle.h>
#include <dune/alugrid/impl/2d/vtx_btree.h>
#include <dune/alugrid/impl/2d/vmmap.h>

//////////////////////////////////////////////////////////////////////
// compile imp.cc into lib (1 yes, 0 no)
// if you change this, you'll get what you deserve
//////////////////////////////////////////////////////////////////////
#define COMPILE_ALU2DGRID_LIB 0

#if COMPILE_ALU2DGRID_LIB
  #define COMPILE_ALU2DGRID_INLINE 0
#else
  #define COMPILE_ALU2DGRID_INLINE 1
#endif

#if COMPILE_ALU2DGRID_INLINE
#define alu2d_inline inline
#else
#define alu2d_inline
#endif
/////////////////////////////////////////////////////////////////////

#define ALU2DSPACE ::ALU2DGrid::

namespace ALU2DGrid
{

  enum ElementType { triangle, quadrilateral, mixed };

} // namespace ALU2DGrid


namespace Dune
{

  typedef double alu2d_ctype;


#ifndef DOXYGEN
  // ALU2dImplInterface
  // ------------------

  template< int dim, int dimw, ALU2DSPACE ElementType eltype >
  struct ALU2dImplInterface;

  template< int dimw, ALU2DSPACE ElementType eltype >
  struct ALU2dImplInterface< 0, dimw, eltype >
  {
    typedef typename ALU2DSPACE Hmesh_basic< dimw, (eltype == ALU2DSPACE triangle ? 3 : 4) >::vertex_t Type;
  };

  template< int dimw, ALU2DSPACE ElementType eltype >
  struct ALU2dImplInterface< 1, dimw, eltype >
  {
    typedef typename ALU2DSPACE Hmesh_basic< dimw, (eltype == ALU2DSPACE triangle ? 3 : 4) >::helement_t Type;
  };

  template< int dimw, ALU2DSPACE ElementType eltype >
  struct ALU2dImplInterface< 2, dimw, eltype >
  {
    typedef typename ALU2DSPACE Hmesh_basic< dimw, (eltype == ALU2DSPACE triangle ? 3 : 4) >::helement_t Type;
  };
#endif


  // ALU2dImplTraits
  // ---------------

  template< int dimw, ALU2DSPACE ElementType eltype >
  struct ALU2dImplTraits
  {
    template< int cdim >
    struct Codim
    {
      typedef typename ALU2dImplInterface< 2-cdim, dimw, eltype >::Type InterfaceType;
    };

    typedef ALU2DSPACE Hmesh< dimw, (eltype == ALU2DSPACE triangle ? 3 : 4) > HmeshType;
    typedef ALU2DSPACE Thinelement< dimw, (eltype == ALU2DSPACE triangle ? 3 : 4) > ThinelementType;
    typedef ALU2DSPACE Element< dimw, (eltype == ALU2DSPACE triangle ? 3 : 4) > ElementType;
    typedef typename HmeshType::helement_t HElementType;
    typedef typename HmeshType::hbndel_t HBndElType;
    typedef ALU2DSPACE Bndel_periodic< dimw, (eltype == ALU2DSPACE triangle ? 3 : 4) > PeriodicBndElType;
  };



  // ALU2dGridMarkerVector
  // ---------------------

  class ALU2dGridMarkerVector
  {
    typedef std::vector< int > VectorType;
  public:
    ALU2dGridMarkerVector() : valid_(false) {}

    bool valid() const { return valid_; }

    void invalidate() { valid_ = false; }

    bool isOnElement(int elementIndex, int idx, int codim) const
    {
      return marker_[codim-1][idx] == elementIndex;
    }

    template <class GridType>
    void update (const GridType & grid, int level )
    {
      enum { dim = GridType::dimension };
      static const int dimworld = GridType::dimensionworld;
      static const ALU2DSPACE ElementType eltype = GridType::elementType;

      typedef typename ALU2dImplTraits< dimworld, eltype >::template Codim<0>::InterfaceType ElementType;
      typedef ALU2DSPACE Listwalkptr< ElementType > IteratorType;

      // resize
      for(int i=0; i<2; ++i)
      {
        int s = grid.hierSetSize(i+1);
        if((int) marker_[i].size() < s ) marker_[i].resize(s);

        size_t markerSize = marker_[i].size();
        // reset marker vector to default value
        for(size_t k=0; k<markerSize; ++k) marker_[i][k] = -1;
      }

      IteratorType iter(grid.myGrid(), level);

      for(iter->first(); !iter->done(); iter->next())
      {
        ElementType & elem = iter->getitem();
        int elIdx = elem.getIndex();

        for(int i=0; i<elem.numvertices(); ++i)
        {
          enum { vxCodim = 1 };
          int vxIdx = elem.getVertex(i)->getIndex();
          if( marker_[vxCodim][vxIdx] < 0) marker_[vxCodim][vxIdx] = elIdx;

          enum { edgeCodim = 0 };
          int edgeIdx = elem.edge_idx(i);
          if( marker_[edgeCodim][edgeIdx] < 0) marker_[edgeCodim][edgeIdx] = elIdx;
        }
      }
      valid_ = true;
    }

  private:
    VectorType marker_[2];

    bool valid_;
  };

  class ALU2dGridLeafMarkerVector
  {
    typedef std::vector< int > VectorType;
  public:
    ALU2dGridLeafMarkerVector() : valid_(false) {}

    bool valid() const { return valid_; }

    void invalidate() { valid_ = false; }

    // return true, if edge is visited on given element
    bool isOnElement(int elementIndex, int idx, int codim) const
    {
      alugrid_assert ( valid_ );
      // this marker only works for codim 1, i.e. edges
      alugrid_assert ( codim == 1 );
      return marker_[idx] == elementIndex;
    }

    // this is for the LeafIterator
    template <class GridType>
    void update (const GridType & grid)
    {
      static const int dimworld = GridType::dimensionworld;
      static const ALU2DSPACE ElementType eltype = GridType::elementType;

      typedef typename ALU2dImplTraits< dimworld, eltype >::template Codim<0>::InterfaceType ElementType;
      typedef ALU2DSPACE Listwalkptr< ElementType > IteratorType;

      // resize edge marker
      {
        int s = grid.hierSetSize(1);
        if((int) marker_.size() < s ) marker_.resize(s);

        size_t markerSize = marker_.size();
        // reset marker vector to default value
        for(size_t k=0; k<markerSize; ++k) marker_[k] = -1;
      }

      // resize vertex levels
      {
        int s = grid.hierSetSize(2);
        if((int) vertexLevels_.size() < s ) vertexLevels_.resize(s);

        // initialize with -1
        size_t vxSize = vertexLevels_.size();
        for(size_t k=0; k<vxSize; ++k) vertexLevels_[k] = -1;
      }

      enum { dim = GridType::dimension };
      IteratorType iter(grid.myGrid());

      for(iter->first(); !iter->done(); iter->next())
      {
        ElementType & elem = iter->getitem();
        int elIdx = elem.getIndex();

        int level = elem.level();

        for(int i=0; i<elem.numvertices(); ++i)
        {
          int vxIdx = elem.getVertex(i)->getIndex();

          // set max level to vertices, see Grid docu paper
          if(level > vertexLevels_[vxIdx]) vertexLevels_[vxIdx] = level;

          int edgeIdx = elem.edge_idx(i);
          if( marker_[edgeIdx] < 0) marker_[edgeIdx] = elIdx;
        }
      }
      valid_ = true;
    }

    //! return level of vertex
    int levelOfVertex(const int vxIdx) const
    {
      alugrid_assert ( valid_ );
      alugrid_assert ( vxIdx >= 0 && vxIdx < (int) vertexLevels_.size());
      // if this assertion is thrown, the level has not been initialized
      alugrid_assert ( vertexLevels_[vxIdx] >= 0 );
      return vertexLevels_[vxIdx];
    }

    //! return level of vertex
    bool isValidVertex(const int vxIdx) const
    {
      alugrid_assert ( valid_ );
      alugrid_assert ( vxIdx >= 0 && vxIdx < (int) vertexLevels_.size());
      return (vertexLevels_[vxIdx] >= 0);
    }

  private:
    VectorType marker_;
    VectorType vertexLevels_;

    bool valid_;
  };

  // dummy object stream class
  class ALU2dGridObjectStream
  {
    public:
      class EOFException {} ;
      template <class T>
      void readObject (T &) {}
      void readObject (int) {}
      void readObject (double) {}
      template <class T>
      void writeObject (T &) {}
      void writeObject (int) {}
      void writeObject (double) {}

      template <class T>
      void read (T &) const {}
      template <class T>
      void write (const T &) {}
  };

} // namespace Dune

#endif // #ifndef DUNE_ALU2DGRID_INCLUDE_HH
