#ifndef DUNE_ALU2DGRID_IMP_CC
#define DUNE_ALU2DGRID_IMP_CC

#include <fstream>

namespace Dune
{

  template< int dim, int dimworld, ALU2DSPACE ElementType eltype >
  inline typename ALU2dGrid< dim, dimworld, eltype >::HmeshType*
  ALU2dGrid< dim, dimworld, eltype >::
  createGrid(const std::string& macroTriangFilename,
             const int nrOfHangingNodes,
             std::istream* macroFile )
  {
    if( macroFile )
    {
      return new HmeshType(*macroFile,
                           nrOfHangingNodes, (nrOfHangingNodes == 0) ?
                           ALU2DSPACE Refco::ref_1 : ALU2DSPACE Refco::quart);
    }
    else
    {
      return new HmeshType(checkMacroGridFile(macroTriangFilename),
                           nrOfHangingNodes, (nrOfHangingNodes == 0) ?
                           ALU2DSPACE Refco::ref_1 : ALU2DSPACE Refco::quart);
    }
  }

  //--Grid
  template< int dim, int dimworld, ALU2DSPACE ElementType eltype >
  inline ALU2dGrid< dim, dimworld, eltype >::
  ALU2dGrid(const std::string macroTriangFilename,
            const int nrOfHangingNodes,
            const DuneBoundaryProjectionType* bndPrj,
            const DuneBoundaryProjectionVector* bndVec,
            std::istream* macroFile )
  :
    mygrid_ ( createGrid(macroTriangFilename, nrOfHangingNodes, macroFile ) ),
    factory_( *this ),
    hIndexSet_(*this),
    localIdSet_(*this),
    levelIndexVec_( MAXL, (LevelIndexSetImp *) 0 ) ,
    geomTypes_( dim+1 ),
    leafIndexSet_(0),
    maxLevel_(0),
    refineMarked_ (0),
    coarsenMarked_ (0) ,
    nrOfHangingNodes_( nrOfHangingNodes ),
    sizeCache_(0),
    lockPostAdapt_(false),
    bndPrj_ ( bndPrj ),
    bndVec_ ( bndVec ),
    vertexProjection_( (bndPrj || bndVec) ? new ALUGridBoundaryProjectionType( *this ) : 0 )
  {
    alugrid_assert (mygrid_);

    // this feature is available in ALUGrid-1.15
    if( vertexProjection_ )
      myGrid().setVertexProjection( vertexProjection_ );

    makeGeomTypes();
    updateStatus();
  }

  //! Constructor which constructs an empty ALU2dGrid
  template< int dim, int dimworld, ALU2DSPACE ElementType eltype >
  inline ALU2dGrid< dim, dimworld, eltype >::ALU2dGrid ( int nrOfHangingNodes )
  : mygrid_( 0 ),
    factory_( *this ),
    hIndexSet_(*this),
    localIdSet_(*this),
    levelIndexVec_( MAXL, (LevelIndexSetImp *) 0 ),
    geomTypes_( dim+1 ),
    leafIndexSet_(0),
    maxLevel_(0),
    refineMarked_ (0),
    coarsenMarked_ (0),
    nrOfHangingNodes_(nrOfHangingNodes),
    sizeCache_(0),
    lockPostAdapt_(false),
    bndPrj_ ( 0 ),
    bndVec_ ( 0 ),
    vertexProjection_( 0 )
  {
    makeGeomTypes();
  }

  template< int dim, int dimworld, ALU2DSPACE ElementType eltype >
  inline ALU2dGrid< dim, dimworld, eltype >::ALU2dGrid ( const ALU2dGrid< dim, dimworld, eltype > &g )
  : mygrid_ ( 0 ),
    maxLevel_(0),
    coarsenMarked_(0) , refineMarked_(0),
    hIndexSet_(*this),
    localIdSet_ (*this),
    levelIndexVec_( MAXL, (LevelIndexSetImp *) 0 ),
    leafIndexSet_( 0 )
  {
    DUNE_THROW( GridError, "Do not use copy constructor of ALU2dGrid!" );
  }

  template< int dim, int dimworld, ALU2DSPACE ElementType eltype >
  inline ALU2dGrid< dim, dimworld, eltype >::~ALU2dGrid()
  {
    delete vertexProjection_ ;
    if( bndVec_ )
    {
      const size_t bndSize = bndVec_->size();
      for(size_t i=0; i<bndSize; ++i)
      {
        delete (*bndVec_)[i];
      }
      delete bndVec_; bndVec_ = 0;
    }

    for(size_t i=0; i<levelIndexVec_.size(); ++i)
    {
      delete levelIndexVec_[i]; levelIndexVec_[i] = 0;
    }
    delete leafIndexSet_; leafIndexSet_ = 0;
    delete sizeCache_; sizeCache_ = 0;
    delete mygrid_;
  }

  //! Iterator to first entity of given codim on level
  template< int dim, int dimworld, ALU2DSPACE ElementType eltype >
  inline const char *
  ALU2dGrid< dim, dimworld, eltype >::checkMacroGridFile(const std::string & filename)
  {
    std::ifstream file( filename.c_str() );
    if( !file )
    {
      std::cerr << "Couldn't open file '" << filename <<"' !" << std::endl;
      DUNE_THROW(IOError,"Couldn't open file '" << filename <<"' !");
    }

    const std::string aluid("!Triangles");
    std::string idline;
    std::getline(file,idline);
    std::stringstream idstream(idline);
    std::string id;
    idstream >> id;

    if(id != aluid )
    {
      derr << "\nNon or wrong keyword '" << id << "' found! Expected keyword to be '" << aluid << "'! \n";
      DUNE_THROW(IOError,"Wrong file format!");
      return 0;
    }
    return filename.c_str();
  }



  //! Iterator to first entity of given codim on level
  template< int dim, int dimworld, ALU2DSPACE ElementType eltype >
  template<int cd, PartitionIteratorType pitype>
  inline typename ALU2dGrid< dim, dimworld, eltype >::Traits::template Codim<cd>::template Partition<pitype>::LevelIterator
  ALU2dGrid< dim, dimworld, eltype >::lbegin (int level) const {
    return ALU2dGridLevelIterator<cd, pitype, const ThisType>( factory(), level, ( comm().size() <= 1 && pitype == Ghost_Partition ) );
  }

  //! one past the end on this level
  template< int dim, int dimworld, ALU2DSPACE ElementType eltype >
  template<int cd, PartitionIteratorType pitype>
  inline typename ALU2dGrid< dim, dimworld, eltype >::Traits::template Codim<cd>::template Partition<pitype>::LevelIterator
  ALU2dGrid< dim, dimworld, eltype >::lend (int level) const {
    return ALU2dGridLevelIterator<cd, pitype, const ThisType>( factory(), level, true);
  }

  //! Iterator to first entity of given codim on level
  template< int dim, int dimworld, ALU2DSPACE ElementType eltype >
  template<int cd>
  inline typename ALU2dGrid< dim, dimworld, eltype >::Traits::template Codim<cd>::template Partition<All_Partition>::LevelIterator
  ALU2dGrid< dim, dimworld, eltype >::lbegin (int level) const {
    return ALU2dGridLevelIterator<cd, All_Partition, const ThisType>( factory(), level, false);
  }

  //! one past the end on this level
  template< int dim, int dimworld, ALU2DSPACE ElementType eltype >
  template<int cd>
  inline typename ALU2dGrid< dim, dimworld, eltype >::Traits::template Codim<cd>::template Partition<All_Partition>::LevelIterator
  ALU2dGrid< dim, dimworld, eltype >::lend (int level) const {
    return ALU2dGridLevelIterator<cd, All_Partition, const ThisType>( factory(), level, true);
  }

  //! Iterator to first entity of codim 0 on level
  template< int dim, int dimworld, ALU2DSPACE ElementType eltype >
  inline typename ALU2dGrid< dim, dimworld, eltype >::LevelIteratorType
  ALU2dGrid< dim, dimworld, eltype >::lbegin (int level) const {
    return LevelIteratorImp( factory(), level, false);
  }

  //! last entity of codim 0 on level
  template< int dim, int dimworld, ALU2DSPACE ElementType eltype >
  inline typename ALU2dGrid< dim, dimworld, eltype >::LevelIteratorType
  ALU2dGrid< dim, dimworld, eltype >::lend (int level) const {
    return LevelIteratorImp( factory(), level, true);
  }

  //! General definiton for a leaf iterator
  template< int dim, int dimworld, ALU2DSPACE ElementType eltype >
  template <int codim, PartitionIteratorType pitype>
  inline typename ALU2dGrid< dim, dimworld, eltype >::Traits::template Codim<codim>::template Partition<pitype>::LeafIterator
  ALU2dGrid< dim, dimworld, eltype >::leafbegin() const {
    return ALU2dGridLeafIterator<codim, pitype, const ThisType> ( factory(), ( comm().size() <= 1 && pitype == Ghost_Partition ) );
  }

  //! General definition for an end iterator on leaf level
  template< int dim, int dimworld, ALU2DSPACE ElementType eltype >
  template <int codim, PartitionIteratorType pitype>
  inline typename ALU2dGrid< dim, dimworld, eltype >::Traits::template Codim<codim>::template Partition<pitype>::LeafIterator
  ALU2dGrid< dim, dimworld, eltype >::leafend() const {
    return ALU2dGridLeafIterator<codim, pitype, const ThisType> ( factory(), true);
  }

  //! General definiton for a leaf iterator
  template< int dim, int dimworld, ALU2DSPACE ElementType eltype >
  template <int codim>
  inline typename ALU2dGrid< dim, dimworld, eltype >::Traits::template Codim<codim>::LeafIterator
  ALU2dGrid< dim, dimworld, eltype >::leafbegin() const {
    return ALU2dGridLeafIterator<codim, All_Partition, const ThisType> ( factory(), false);
  }

  //! General definition for an end iterator on leaf level
  template< int dim, int dimworld, ALU2DSPACE ElementType eltype >
  template <int codim>
  inline typename ALU2dGrid< dim, dimworld, eltype >::Traits::template Codim<codim>::LeafIterator
  ALU2dGrid< dim, dimworld, eltype >::leafend() const {
    return ALU2dGridLeafIterator<codim, All_Partition, const ThisType> ( factory(), true);
  }

  //! Iterator to first entity of codim 0 on leaf level (All_Partition)
  template< int dim, int dimworld, ALU2DSPACE ElementType eltype >
  inline typename ALU2dGrid< dim, dimworld, eltype >::LeafIteratorType
  ALU2dGrid< dim, dimworld, eltype >::leafbegin () const {
    return LeafIteratorImp( factory(), false);
  }

  //! one past the end on this leaf level (codim 0 and All_Partition)
  template< int dim, int dimworld, ALU2DSPACE ElementType eltype >
  inline typename ALU2dGrid< dim, dimworld, eltype >::LeafIteratorType
  ALU2dGrid< dim, dimworld, eltype >::leafend () const {
    return LeafIteratorImp( factory(), true);
  }

  //! Return maximum level defined in this grid. Levels are numbered
  //! 0 ... maxLevel with 0 the coarsest level.
  template< int dim, int dimworld, ALU2DSPACE ElementType eltype >
  inline int ALU2dGrid< dim, dimworld, eltype >::maxLevel() const {
    return maxLevel_;
  }

  template< int dim, int dimworld, ALU2DSPACE ElementType eltype >
  inline void ALU2dGrid< dim, dimworld, eltype >::calcMaxlevel()
  {
    maxLevel_ = 0;
    // walk the leaf level and take maximum as maxLevel
    ALU2DSPACE Listwalkptr< HElementType > walk( mesh() );
    for( walk->first(); !walk->done(); walk->next() )
      maxLevel_ = std::max( maxLevel_, walk->getitem().level() );
  }

  template< int dim, int dimworld, ALU2DSPACE ElementType eltype >
  inline void ALU2dGrid< dim, dimworld, eltype >::calcExtras()
  {
    // only in truly parallel runs update this
    //if( comm_.size() > 1 ) rankManager_.update();

    ////////////////////////////////////
    // reset size cache
    ////////////////////////////////////
    if(sizeCache_) delete sizeCache_;
    if ( eltype == ALU2DSPACE mixed )
      DUNE_THROW( NotImplemented, "size for mixed grids" );

    sizeCache_ = new SizeCacheType (*this);
    /////////////////////////////////////

    // place this before update of level index, because
    // marker vector is used by level index set
    for(int i=0; i<MAXL; ++i) marker_[i].invalidate();
    leafMarker_.invalidate();

    ///////////////////////////////////////////////
    // update existing index sets
    ///////////////////////////////////////////////
    const int levelSize = levelIndexVec_.size();
    for(int i=0; i<levelSize; ++i)
    {
      if(levelIndexVec_[i])
        (*(levelIndexVec_[i])).calcNewIndex(
                               this->template lbegin<0>( i ),
                               this->template lend<0>( i ) );
    }

    if(leafIndexSet_)
      leafIndexSet_->calcNewIndex( this->template leafbegin<0>(),
                                   this->template leafend<0>() );
    ////////////////////////////////////////////////
  }

  //! Every time the grid is refined, data should be updated
  template< int dim, int dimworld, ALU2DSPACE ElementType eltype >
  inline void ALU2dGrid< dim, dimworld, eltype >::updateStatus() {
    calcMaxlevel();
    calcExtras();
  }

  //! number of grid entities in the entire grid for given codim
  template< int dim, int dimworld, ALU2DSPACE ElementType eltype >
  inline int ALU2dGrid< dim, dimworld, eltype >::hierSetSize (int cd) const {
    return mesh().indexManagerSize(cd);
  }

  //! number of grid entities per level and codim
  template< int dim, int dimworld, ALU2DSPACE ElementType eltype >
  inline int ALU2dGrid< dim, dimworld, eltype >::size (int level, int cd) const
  {
    alugrid_assert ( sizeCache_ );
    return sizeCache_->size(level,cd);
  }

  //! number of entities per level, codim and geometry type in this process
  template< int dim, int dimworld, ALU2DSPACE ElementType eltype >
  inline int ALU2dGrid< dim, dimworld, eltype >::size (int level, GeometryType type) const
  {
    switch (eltype)
    {
      case ALU2DSPACE triangle:
        return type.isSimplex()? size(level, dim-type.dim()):0;
      case ALU2DSPACE quadrilateral:
        return type.isCube()? size(level, dim-type.dim()):0;
      case ALU2DSPACE mixed:
        DUNE_THROW( NotImplemented, "size method for mixed grids" );
    }
  }

  //! number of leaf entities per codim and geometry type in this process
  template< int dim, int dimworld, ALU2DSPACE ElementType eltype >
  inline int ALU2dGrid< dim, dimworld, eltype >::size (GeometryType type) const
  {
    switch (eltype)
    {
      case ALU2DSPACE triangle:
        return type.isSimplex()? size(dim-type.dim()):0;
      case ALU2DSPACE quadrilateral:
        return type.isCube()? size(dim-type.dim()):0;
      case ALU2DSPACE mixed:
        DUNE_THROW( NotImplemented, "size method for mixed grids" );
    }
  }

  template< int dim, int dimworld, ALU2DSPACE ElementType eltype >
  inline int ALU2dGrid< dim, dimworld, eltype >::size (int codim) const
  {
    alugrid_assert ( sizeCache_ );
    return sizeCache_->size(codim);
  }

  //! refine grid refCount times
  template< int dim, int dimworld, ALU2DSPACE ElementType eltype >
  inline void ALU2dGrid< dim, dimworld, eltype >::globalRefine(int refCount)
  {
    if( refCount <= 0 )
      return;

    for (int j = 0; j < refCount; ++j)
    {
      ALU2DSPACE Listwalkptr< HElementType > walk( mesh() );
      for( walk->first(); !walk->done(); walk->next() )
      {
        HElementType &item = walk->getitem();
        item.ALU2DSPACE Refco_el::mark( ALU2DSPACE Refco::ref );
      }

      mesh().refine();
    }

    //update data
    updateStatus();

    // cleanup markers
    postAdapt();
  }


  // global refine
  template< int dim, int dimworld, ALU2DSPACE ElementType eltype >
  template< class GridImp, class DataHandle >
  inline void ALU2dGrid< dim, dimworld, eltype >
    ::globalRefine ( int refCount, AdaptDataHandleInterface< GridImp, DataHandle > &handle )
  {
    alugrid_assert ( (refCount + maxLevel()) < MAXL );

    for( int count = refCount; count > 0; --count )
    {
      typedef typename Traits::template Codim< 0 >::template Partition< Interior_Partition >::LeafIterator LeafIterator;
      const LeafIterator end = leafend< 0, Interior_Partition >();
      for( LeafIterator it = leafbegin< 0, Interior_Partition >(); it != end; ++it )
        mark( 1 , *it );
      adapt( handle );
    }
  }


  //! returns true if a least one entity was marked for coarseing
  template< int dim, int dimworld, ALU2DSPACE ElementType eltype >
  inline bool ALU2dGrid< dim, dimworld, eltype >::preAdapt () {
    return (coarsenMarked_ > 0);
  }


  template< int dim, int dimworld, ALU2DSPACE ElementType eltype >
  inline void ALU2dGrid< dim, dimworld, eltype >::hierarchicClear ( HElementType *el )
  {
    // clear actual tag
    el->ALU2DSPACE Refco_el::clear();
    // clear refined tag
    el->ALU2DSPACE Refco_el::clearWas();
    // go to children
    for( HElementType *child = el->down(); child; child = child->next() )
    {
      // clear marker for child
      hierarchicClear(child);
    }
  }


  //! clear all entity new markers
  template< int dim, int dimworld, ALU2DSPACE ElementType eltype >
  inline void ALU2dGrid< dim, dimworld, eltype >::postAdapt ()
  {
    // clear refinement markers throughout the grid
    typedef ALU2DSPACE Macro < ElementType > macro_t;

    // get macro element iterator
    ALU2DSPACE Listwalkptr <macro_t> walk(mesh());
    for( walk->first() ; ! walk->done() ; walk->next())
    {
      // get element pointer
      HElementType *el = walk->getitem().operator ->();
      // hierarchically clear all markers
      hierarchicClear( el );
    }

    // reset marked element counters
    coarsenMarked_ = 0;
    refineMarked_  = 0;

    // unmark postAdapt flag
    lockPostAdapt_ = false;
  }


  /**! refine all positive marked leaf entities,
    return true if a least one entity was refined
  */
  template< int dim, int dimworld, ALU2DSPACE ElementType eltype >
  inline bool ALU2dGrid< dim, dimworld, eltype >::adapt ()
  {
    if( lockPostAdapt_ == true )
      DUNE_THROW( InvalidStateException, "Make sure that postAdapt is called after adapt was called and returned true!" );

    if( (refineMarked_ > 0) || (coarsenMarked_ > 0) )
    {
      // refine only will be done if
      // at least one element was marked for refinement
      bool adapted = (refineMarked_) ? true : false;
      mesh().refine();
      mesh().coarse();

      updateStatus();

      // notify that postAdapt must be called
      lockPostAdapt_ = true;

      return adapted;
    }
    else
      return false;
  }


  // --adapt
  template< int dim, int dimworld, ALU2DSPACE ElementType eltype >
  template< class GridImp, class DataHandle >
  inline bool ALU2dGrid< dim, dimworld, eltype >
    ::adapt ( AdaptDataHandleInterface< GridImp, DataHandle > &handle )
  {
    typedef AdaptDataHandleInterface< GridImp, DataHandle > AdaptDataHandle;

    typedef typename EntityObject::ImplementationType EntityImp;
    EntityObject father( EntityImp( factory(), maxLevel()) );
    EntityObject son   ( EntityImp( factory(), maxLevel()) );

    // int defaultChunk = newElementsChunk_;
    // int actChunk     = refineEstimate_ * refineMarked_;

    // guess how many new elements we get
    // int newElements = std::max( actChunk , defaultChunk );

    // reserve memory
    // handle.preAdapt( newElements );

    bool ref = false;
    {
      ALU2DSPACE AdaptRestrictProlong2dImpl< ThisType, AdaptDataHandle >
        rp(*this,
           father,this->getRealImplementation(father),
           son,   this->getRealImplementation(son),
           handle );

      ref = myGrid().duneAdapt(rp); // adapt grid
      //if(rp.maxLevel() >= 0) maxlevel_ = rp.maxLevel();
    }

    if(ref) {
      updateStatus();
    }

    // check whether we have balance
    // handle.postAdapt();

    postAdapt();
    return ref;
  }

  template< int dim, int dimworld, ALU2DSPACE ElementType eltype >
  inline bool ALU2dGrid< dim, dimworld, eltype >::
  mark( int refCount , const typename Traits::template Codim<0>::Entity & en )
  {
    bool marked = this->getRealImplementation(en).mark(refCount);
    if(marked)
    {
      if(refCount > 0) ++refineMarked_;
      if(refCount < 0) ++coarsenMarked_;
    }
    return marked;
  }

  template< int dim, int dimworld, ALU2DSPACE ElementType eltype >
  inline int ALU2dGrid< dim, dimworld, eltype >::getMark(const typename Traits::template Codim<0>::Entity & e ) const
  {
    return this->getRealImplementation(e).getMark();
  }

  template< int dim, int dimworld, ALU2DSPACE ElementType eltype >
  inline void ALU2dGrid< dim, dimworld, eltype >::makeGeomTypes()
  {
    alugrid_assert ( eltype != ALU2DSPACE mixed ); // ?????????
    const GeometryType :: BasicType basic = ( eltype == ALU2DSPACE triangle )?
      GeometryType :: simplex : GeometryType::cube;
    // stored is the dim, where is the codim
    for (int i=0;i<3;++i)
      geomTypes_[i].clear();
    geomTypes_[ 2 ].push_back( GeometryType( basic, 0 ) );
    geomTypes_[ 1 ].push_back( GeometryType( basic, 1 ) );
    geomTypes_[ 0 ].push_back( GeometryType( basic, 2 ) );

    // initialize static storage variables
    ALU2dGridGeometry<0,   dimworld, const ThisType> :: geoProvider();
    ALU2dGridGeometry<1,   dimworld, const ThisType> :: geoProvider();
    ALU2dGridGeometry<dim, dimworld, const ThisType> :: geoProvider();

    // local geometry storage in Intersection
    ALU2dGridIntersectionBase < const ThisType > :: LocalGeometryStorageType :: instance();

    static bool firstCall = true ;
    if( firstCall )
    {
      // unset flag first since local geometries are generated by generating a grid,
      // which would result in an infinite loop
      firstCall = false ;
      // geometryInFather storage
      ALU2dGridEntity<0, dim, const ThisType> :: geometryInFatherImpl( nonConform(), geomTypes_[ 0 ][ 0 ], 0 );
    }
  }

   //! get global id set of grid
  template< int dim, int dimworld, ALU2DSPACE ElementType eltype >
  inline const typename ALU2dGrid< dim, dimworld, eltype >::GlobalIdSet &
  ALU2dGrid< dim, dimworld, eltype >::globalIdSet () const
  {
    return localIdSet();
  }

   //! get global id set of grid
  template< int dim, int dimworld, ALU2DSPACE ElementType eltype >
  inline const typename ALU2dGrid< dim, dimworld, eltype >::LocalIdSet &
  ALU2dGrid< dim, dimworld, eltype >::localIdSet () const
  {
    return localIdSet_;
  }

  //! get hierarchic index set of the grid
  template< int dim, int dimworld, ALU2DSPACE ElementType eltype >
  inline const typename ALU2dGrid< dim, dimworld, eltype >::HierarchicIndexSet &
  ALU2dGrid< dim, dimworld, eltype >::hierarchicIndexSet () const
  {
    return hIndexSet_;
  }

  template< int dim, int dimworld, ALU2DSPACE ElementType eltype >
  inline const typename ALU2dGrid< dim, dimworld, eltype >::Traits::LeafIndexSet &
  ALU2dGrid< dim, dimworld, eltype >::leafIndexSet() const
  {
    if(!leafIndexSet_)
      leafIndexSet_ = new LeafIndexSetImp ( *this,
                                            this->template leafbegin<0>(),
                                            this->template leafend<0>() );
    return *leafIndexSet_;
  }

  template< int dim, int dimworld, ALU2DSPACE ElementType eltype >
  inline const typename ALU2dGrid< dim, dimworld, eltype >::Traits::LevelIndexSet &
  ALU2dGrid< dim, dimworld, eltype >::levelIndexSet( int level ) const
  {
    // check if level fits in vector
    alugrid_assert ( level >= 0 );
    alugrid_assert ( levelIndexVec_.size() == MAXL );
    alugrid_assert ( level < (int) levelIndexVec_.size() );

    if( levelIndexVec_[level] == 0 )
      levelIndexVec_[level] =
        new LevelIndexSetImp ( *this,
                               this->template lbegin<0>(level),
                               this->template lend<0>(level),
                               level );
    return *(levelIndexVec_[level]);
  }


  template< int dim, int dimworld, ALU2DSPACE ElementType eltype >
  inline ALU2dGrid< dim, dimworld, eltype > &
  ALU2dGrid< dim, dimworld, eltype >::operator= ( const ALU2dGrid< dim, dimworld, eltype > &g)
  {
    DUNE_THROW(GridError,"Do not use assignment operator of ALU2dGrid! \n");
    return (*this);
  }

  // private methods that return underlying ALU2D Grid
  template< int dim, int dimworld, ALU2DSPACE ElementType eltype >
  inline typename ALU2dGrid< dim, dimworld, eltype >::HmeshType &
  ALU2dGrid< dim, dimworld, eltype >::myGrid()
  {
    return mesh();
  }

  template< int dim, int dimworld, ALU2DSPACE ElementType eltype >
  inline typename ALU2dGrid< dim, dimworld, eltype >::HmeshType &
  ALU2dGrid< dim, dimworld, eltype >::myGrid() const
  {
    return mesh();
  }

  //! return dummy communication
  template< int dim, int dimworld, ALU2DSPACE ElementType eltype >
  inline const typename ALU2dGrid< dim, dimworld, eltype >::CollectiveCommunicationType &
  ALU2dGrid< dim, dimworld, eltype >::comm() const
  {
    return comm_;
  }

  template< int dim, int dimworld, ALU2DSPACE ElementType eltype >
  inline void ALU2dGrid< dim, dimworld, eltype >::
  backup( std::ostream& stream, const ALU3DSPACE MacroFileHeader::Format format  ) const
  {
    // write grid to stream
    myGrid().storeGrid( stream );
  }

  template< int dim, int dimworld, ALU2DSPACE ElementType eltype >
  inline void ALU2dGrid< dim,dimworld, eltype >::
  restore( std::istream& stream )
  {
    // if grid exists delete first
    if( mygrid_ ) delete mygrid_;

    // create new grid with given istream
    mygrid_ = new HmeshType ( stream,
                              nrOfHangingNodes_,
                              ( nrOfHangingNodes_ == 0) ? ALU2DSPACE Refco::ref_1 : ALU2DSPACE Refco::quart );

    if( ! mygrid_ )
    {
      DUNE_THROW(InvalidStateException,"ALUGrid::restore failed");
    }

    // calculate new maxlevel
    // calculate indices
    updateStatus();

    // cleanup markers
    postAdapt();
  }


  // communicate level data
  template< int dim, int dimworld, ALU2DSPACE ElementType eltype >
  template <class DataHandleImp,class DataType>
  inline void ALU2dGrid< dim, dimworld, eltype >::
  communicate (CommDataHandleIF<DataHandleImp,DataType> & data,
               InterfaceType iftype, CommunicationDirection dir, int level ) const
  {}

  // communicate level data
  template< int dim, int dimworld, ALU2DSPACE ElementType eltype >
  template <class DataHandleImp,class DataType>
  inline void ALU2dGrid< dim, dimworld, eltype >::
  communicate (CommDataHandleIF<DataHandleImp,DataType> & data,
               InterfaceType iftype, CommunicationDirection dir) const
  {}

  // re-balance load of grid
  template< int dim, int dimworld, ALU2DSPACE ElementType eltype >
  inline bool ALU2dGrid< dim, dimworld, eltype >::loadBalance ()
  {
    return false;
  }

  // re-balance load of grid
  template< int dim, int dimworld, ALU2DSPACE ElementType eltype >
  template <class DataHandleImp>
  inline bool ALU2dGrid< dim, dimworld, eltype >::loadBalance (DataHandleImp& data)
  {
    return false;
  }

} // namespace Dune

#endif // #ifndef DUNE_ALU2DGRID_IMP_CC
