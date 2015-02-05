#ifndef __HEADER__LISTWALK
#define __HEADER__LISTWALK

namespace ALU2DGrid
{

  template < class A > class Hier;   // hierarchie-template in grid.h

  template < class A > class Listwalkptr;
  template < class A > class Listagent;
  template < class A > class Listwalk_impl;  // echter iterator
  template < class A > class Listwalk;         // interface-klasse
  template < class A > class Macro;


  template < class A > class Listagency {

    IndexProvider *hdl;

    A * first_list;

    A * last_list;

    int number_list;

    //int nlistwalks;

    void listwalkattach() {} //nlistwalks ++; }

    void listwalkdetach() {} // nlistwalks --; }

    int operator == (const Listagency & a) { return (void *)this == (void *) & a; }

    public :

    void detach(A * a) {

      if(a == last_list)  last_list = a->prev_list;

      if(a == first_list) first_list = a->next_list;

      if(a->prev_list) a->prev_list->next_list = a->next_list;

      if(a->next_list) a->next_list->prev_list = a->prev_list;

      a->next_list = 0;

      a->prev_list = 0;

      a->number_list = 0;

      number_list --;


    }

    public :

    Listagency(IndexProvider *phdl) : hdl(phdl),
    first_list(0), last_list(0), number_list(0) //, nlistwalks(0)
    {}

     ~Listagency() {
        alugrid_assert (! busy());

        A * curr = first_list;

        while(curr) {

          A * next_list = curr->next_list;

          delete curr;

          curr = next_list;
        }

      }

      int renumber() {

        A * curr = first_list;

        int count = 0;

        while(curr) {

          curr->number_list = count ++;

          curr = curr->next_list;

        }

        return count;

      }

      void insert(A * a)
      {
        if (last_list)
          last_list->next_list = a;
        else
          first_list = a;

        a->prev_list = last_list;
        last_list = a;
        a->sethdl(hdl);
        number_list ++;
      }

      int size() const { return number_list; }

      int busy() const { /*return nlistwalks;*/ return 0; }

    friend class Listagent < A >;

    friend class Listwalk_impl < A >;

  };


  // +----------------------------------------------------------------------------+
  // | template Listwalk < . > stellt das interface dar.                          |
  // | Alle abgeleiteten Klassen muessen diese Funktionalitaet implementieren.    |
  // +----------------------------------------------------------------------------+

  template < class A > class Listwalk {

    protected:
      // do not creat instances of interface class
      Listwalk() {}
    public :
      // destructor
      virtual ~Listwalk() {}

      // set iterator to first element
      virtual void first() = 0;

      // set iterator to last element
      virtual void last() = 0;

      // set iterator to next element
      virtual void next() = 0;

      // set iterator to previous element
      virtual void prev() = 0;

      // return done
      virtual int done() const = 0;

      // return size of iterated set
      virtual int size() = 0;

      // return current element
      virtual A & getitem() const = 0;

      // create clone of object
      virtual Listwalk < A >  * clone () const = 0;
  };


  // +----------------------------------------------------------------------------+
  // | template Listwalk_empty < . > stellt den null iterator dar.                |
  // +----------------------------------------------------------------------------+
  template < class A > class Listwalk_empty : public Listwalk< A >
  {
    public :
      // createing empty list walk
      Listwalk_empty() {}

      // copy construtor
      Listwalk_empty(const Listwalk_empty &) {}

      // assign list walk
      Listwalk_empty & operator=(const Listwalk_empty &) { return * this; }

      virtual ~Listwalk_empty() { }

      virtual void first() { }

      virtual void last() { }

      virtual void next() { }

      virtual void prev() { }

      virtual int done() const { return 1; }

      virtual int size() { return 0; }

      virtual A & getitem() const {
        alugrid_assert ( ! done () );
        void * p = (void *)(0);  // schlechter Trick ...
        // in der Tat, sehr schlechter Trick, R.K.
        abort();      // .. und Tsch"uss !
        return *(A *)(p);  // never reached .
      }

      // create clone of object
      virtual Listwalk < A >  * clone () const { return new  Listwalk_empty < A > (*this); }
  };


  // +----------------------------------------------------------------------------+
  // |  das template Alignwalk < ., ., . > erm"oglicht das hintereinanderh"angen  |
  // |  von bestehenden iteratoren, falls es eine gemeinsamen basisklasse gibt,   |
  // |  die durch das dritte argument eingef"uhrt wird.                           |
  // +----------------------------------------------------------------------------+

  template < class A, class B, class C > class Alignwalk : public Listwalk < C > {

    Listwalk < A > & walk1;

    Listwalk < B > & walk2;

    int curr;

    public :

      Alignwalk(Listwalk < A > & a, Listwalk < B > & b) : walk1(a), walk2(b), curr(0) { }

      Alignwalk(const Alignwalk < A , B , C > & org )
        : walk1(org.walk1), walk2(org.walk2), curr(org.curr)
      {}

     ~Alignwalk() { }

      void first() {

        curr = 0;

        walk1.first();

        if(walk1.done()) {

          curr = 1;

          walk2.first();

        }

      }

      void last() {

        curr = 1;

        walk2.last();

        if(walk2.done()) {

          curr = 0;

          walk1.last();

        }

      }

      void next() {

        !curr ? (walk1.next(), (walk1.done() ? (walk2.first(), curr = 1) : 0)) : (walk2.next(), 0);

      }

      void prev() {

        curr ? (walk2.prev(), (walk2.done() ? (walk1.last(), curr = 0) : 0)) : (walk1.prev(), 0);

      }

      int size() {

        return walk1.size() + walk2.size();

      }

      int done() const { return curr ? walk2.done() : 0;  }

      C & getitem() const {

        if(curr) return walk2.getitem();

        else return walk1.getitem();

      }

      virtual Listwalk < C > * clone () const
      {
        return new Alignwalk < A , B , C > (*this);
      }
  };


  // +----------------------------------------------------------------------------+
  // |  f"ur intrusive listen die das template Listagent < . > enthalten, gibt    |
  // |  der folgende aufruf eine iterator instanz: Listwalk_impl < . > (agency), wobei |
  // |  der konstruktor mit einer referenz auf eine Listagency < . > instanz ver- |
  // |  sorgt werden muss. keine leeren listen als default.                       |
  // +----------------------------------------------------------------------------+

  template < class A > class Listwalk_impl : public Listwalk < A > {

    Listagency < A > & agency;

    A * curr;

    public :

      // copy constructor
      Listwalk_impl ( const Listwalk_impl & w)
        : agency(w.agency), curr(w.curr)
      {
        agency.listwalkattach();
      }

      Listwalk_impl ( Listagency < A > & a)
        : agency(a), curr(0)
      {
        agency.listwalkattach();
      }

     ~Listwalk_impl() { agency.listwalkdetach(); }

      void first() { curr = agency.first_list; }

      void last()  { curr = agency.last_list; }

      void next()  { if(curr) curr = curr->Listagent < A > :: next(); }

      void prev()  { if(curr) curr = curr->Listagent < A > :: prev(); }

      int done() const { return (curr == 0) ? 1 : 0; }

      int size() { return agency.number_list; }

      A & getitem() const
      {
        alugrid_assert ( curr );
        return * curr;
      }

      virtual Listwalk < A > * clone() const
      {
        return new Listwalk_impl < A > (*this);
      }
  };


  template < class A > class Macro : public Listagent < Macro < A > > {

    public :

      typedef Hier < A > hier_t;

    private :

      hier_t & el;

      Macro(const Macro & e) : el(e.el) { }

      Macro & operator = (const Macro &) { }

    public :

      void sethdl(IndexProvider *phdl) {el.sethdl(phdl);}

      Macro(hier_t & e) : el(e) { }

     ~Macro() { delete & el; }

      hier_t * operator->() const { return & el; }

      hier_t & operator * () const { return el; }

      hier_t & operator [](int ) const { return el; }

  };


  // +----------------------------------------------------------------------------+
  // | Das Leafwalk template versetzt einen in die Lage, f"ur alle Typen, die     |
  // | eine h < typ > basisklasse haben, und f"ur die eine Liste von Macro < typ >|
  // | vorliegt mit einer Listagency < Macro < typ > > im handle, einen Listen -  |
  // | Iterator zu bekommen, der nur "uber die Bl"atter l"auft.                   |
  // +----------------------------------------------------------------------------+

  template < class A > class Leafwalk : public Listwalk < Hier < A > >
  {
    public :

      typedef Macro < A > macro_t;

      typedef Hier < A > hier_t;

    private :

      enum { max = 100 };

      Listagency < macro_t > & agency;

      Listwalk_impl < macro_t > macro;

      hier_t ** stack;

      int pos;

    public :

      Leafwalk(Listagency < macro_t > & a)
        : agency(a), macro(agency) , pos(0)
      {
        stack = new hier_t * [max];
        alugrid_assert (stack);
      }

      Leafwalk(const Leafwalk & w) : agency(w.agency)
                                   , macro(w.macro)
                                   , pos(w.pos)
      {
        stack = new hier_t * [max];
        alugrid_assert (stack);
        // pos is implicitly set by iterating until w.pos
        for(int p = 0; p <= w.pos; ++ p) stack [p] = w.stack [p];
      }

      ~Leafwalk()
      {
        delete[] stack;
      }

      void first()
      {
        pos = 0;

        macro.first();

        * stack = macro.done() ? 0 : & * macro.getitem();

        for( hier_t * e = stack[pos]; e ? ! e->leaf() : 0; alugrid_assert (pos < max) )

          stack[ ++ pos] = (e = e->down());
      }

      void last()
      {
        macro.last();

        pos = 0;

        hier_t * e = ( * stack = macro.done() ? 0 : & * macro.getitem() );

        for(; e ? ! e->leaf() : 0; ) {

          for( e = e->down(); e->next(); e = e->next());

          stack [ ++ pos] = e;

        }

      }

      void next() {

        for(; pos > 0; pos -- )
        {
          stack[pos] = stack[pos]->next();
          if(stack[pos]) break;
        }

        if(pos == 0)
        {
          macro.next();

          * stack = macro.done() ? 0 : & * macro.getitem();
        }

        for( hier_t * e = stack[pos]; e ? ! e->leaf() : 0; alugrid_assert (pos < max))
          stack[ ++ pos] = (e = e->down());
      }

      void prev() {

        for(; pos > 0; pos -- )
          if(stack[pos - 1]->down() != stack[pos]) break;


        if(pos == 0)
        {
          macro.prev();

          * stack = macro.done() ? 0 : & * macro.getitem();
        }

        else pos --;

        hier_t * e = stack [pos];

        hier_t * u = stack [pos + 1];

        for(; e ? ! e->leaf() : 0; ) {

          for( e = e->down(); ! (e->next() == u || e->next() == 0 ); e = e->next());

          stack[ ++ pos] = e;

        }

      }

      int size()
      {

        Listwalk_impl < macro_t > walk (agency);
        int lsize = 0;

        for( walk.first(); ! walk.done(); walk.next())
          lsize += walk.getitem()->count();

        return lsize;
     }

     int done() const { return macro.done(); }

     hier_t & getitem() const
     {
        alugrid_assert ( stack[pos] );
        return * stack[pos];
     }

     virtual Listwalk< Hier< A > > * clone () const {
       return new Leafwalk< A > (*this);
     }
  };


  // +----------------------------------------------------------------------------+
  // | Das Levelwalk template versetzt einen in die Lage, einen Listen - Iterator |
  // | zu bekommen, der nur "uber die Elemente eines Levels l"auft.               |
  // +----------------------------------------------------------------------------+

  template < class A > class Levelwalk : public Listwalk < Hier < A > > {

    public :

      typedef Macro < A > macro_t;

      typedef Hier < A > hier_t;

    private :

      Listagency < macro_t > & agency;

      Listwalk_impl < macro_t > macro;

      hier_t ** stack;

      int pos;

      int depth;

      int pushdown() {

        for( hier_t * e = stack [pos]; pos < depth; stack[ ++ pos] = (e = e->down()))

          if(e->leaf()) break;

        return pos == depth ? 1 : 0;

      }

      int pullup()
      {
        for(; pos > 0; pos -- )
        {
          stack[pos] = stack[pos]->next();
          if(stack[pos]) break;
        }
        return pos == 0 ? 0 : 1;
      }

    public :

      Levelwalk( Listagency < macro_t > & a, int d)
        : agency(a), macro(agency), pos(0) , depth(d)
      {
        stack = new hier_t * [depth + 1];
        alugrid_assert (stack);
      }

      Levelwalk( const Levelwalk & w )
        : agency(w.agency)
        , macro(w.macro)
        , pos(w.pos)
        , depth(w.depth)
      {
        stack = new hier_t * [depth + 1];
        alugrid_assert (stack);
        // pos is implicitly set by iterating until w.pos
        for(int p = 0; p <= w.pos; ++p) stack [p] = w.stack [p];
      }

      ~Levelwalk()
      {
        delete[] stack;
      }

      void first()
      {
        pos = 0;

        for( macro.first(); ! macro.done(); macro.next())
        {
          * stack = & * macro.getitem();

          pos = 0;

          do {

            if(pushdown()) return; // treffer

          } while(pullup());

        }

      }

      void last ()
      {
        std::cerr << "ERROR: Levelwalk< . >.last() not implemented." << std::endl;
        abort();
      }

      void next ()
      {
        do {
          while ( pullup() )
          {
            if( pushdown() )
              return;
          }

          pos = 0;

          macro.next();

          if(macro.done()) { * stack = 0; return; }

          else {

            * stack = & * macro.getitem();

            if(pushdown()) return;

          }
        } while(depth);
      }

      void prev ()
      {
        std::cerr << "Levelwalk< . >.prev() not implemented." << std::endl;
        abort();
      }

      int size()
      {

        Listwalk_impl < macro_t > walk (agency);

        int lsize = 0;

        for( walk.first(); ! walk.done(); walk.next())
          lsize += walk.getitem()->count(depth);

        return lsize;
      }

      int done() const { return macro.done(); }

      hier_t & getitem() const
      {
        alugrid_assert ( stack[pos] );
        alugrid_assert ( stack[pos]->level() == depth);
        return * stack[pos];
      }

      macro_t & getmacro() const
      {
        alugrid_assert ( stack[pos] );
        alugrid_assert (stack[pos]->level() == depth);
        alugrid_assert (depth==0);
        return macro.getitem();
      }

      virtual Listwalk< Hier< A > > * clone () const
      {
        return new Levelwalk< A > (*this);
      }
  };


  /*
   * SubtreeIterator-Instanzen koennen benutz werden um ueber einen Unterbaum
   * einer Hier<t>-Struktur zu iterieren. Dazu benutzt man die stIterator-Methode der
   * Hier<t> Klasse mit dem entsprechenden Knoten als Argument.
   */

  template <class T>
  class SubtreeIterator {
    static const int nbr = 100;
    Hier<T>* stack[nbr];
    Hier<T>* root_;
    bool done_;
    int pos;

  public:
    SubtreeIterator(Hier<T>* root) {
      pos = 0;
      done_ = false;
      stack[pos] = root;
      root_ = root;
    }

    SubtreeIterator(const SubtreeIterator& other) {
      pos = other.pos;
      for( int i=0; i<=pos; i++ )
        stack[pos] = other.stack[pos];
    }

    SubtreeIterator& operator++ () { // prefix
      if( !done_ ) {
        if( stack[pos]->down() ) {
          stack[pos+1] = stack[pos]->down();
          pos++;
        } else if( stack[pos]->next() ) {
          stack[pos+1] = stack[pos]->next();
          pos++;
        } else {
          for(; pos > 0; pos-- ) {
            if( stack[pos-1]->next() == stack[pos] ||
                stack[pos-1]->next() == 0 )
              stack[pos] = 0;
            else {
              stack[pos] = stack[pos-1]->next();
              break;
            }
          }
          if( pos == 0 )
            done_ = true;
        }
      }
      return *this;
    }

    void first() {
      pos = 0;
      done_ = false;
      stack [pos] = root_;
    }
    void next() {
      ++(*this);
    }
    int done() {
      return (!done_ ? 0 : 1);
    }

    void first_leaf() {
      pos = 0;
      done_ = false;
      while( !done_ && !(stack[pos]->leaf()) )
        ++(*this);
    }

    void next_leaf() {
      while( !done_ && !(stack[pos]->leaf()) )
        ++(*this);
    }

    int done_leaf() {
      return (!done_ ? 0 : 1);
    }

    operator int() const { return (done_ == true ? 0 : 1); }

    Hier<T>* operator-> () const {
      return (done_ ? NULL : stack[pos]);
    }
    Hier<T>& getitem() const {
      alugrid_assert (!done_);
      return *(stack[pos]);
    }

    friend class Hier<T>;
  };


  template < class A > class Listwalkptr {

    Listwalk < A > * walk;

    IndexProvider * hdl;

    A * a;

    // this class is a pointer itself, therefore no creation of pointers
    // allowed
    void * operator new (size_t );
    void operator delete (void *);

    public :

      // create empty list walk
      Listwalkptr()
        : walk ( new Listwalk_empty < A > () ) , hdl(0) , a(0)
      {}

      template <int N,int NV>
      Listwalkptr(Hmesh_basic<N,NV> &h) : hdl(&h) , a(0) { walk = h.walk(a); }

      template <int N,int NV>
      Listwalkptr(Hmesh_basic<N,NV> &h, int level) : hdl(&h) , a(0) { walk = h.walk(a, level); }

      Listwalkptr(const Listwalkptr & p) : hdl(p.hdl) , a(0)
      {
        alugrid_assert ( p.walk );
        walk = p.walk->clone();
      }

      ~Listwalkptr() { delete walk; }

      Listwalkptr & operator = (const Listwalkptr & p)
      {
        delete walk;
        hdl = p.hdl;
        alugrid_assert ( p.walk );
        walk = p.walk->clone();
        return *this;
      }

      Listwalk < A > * operator -> () { return walk; }
      const Listwalk < A > * operator -> () const { return walk; }

      Listwalk < A > & operator * ()  { return * walk; }
      const Listwalk < A > & operator * () const { return * walk; }

      Listwalk < A > & operator [] (int ) { return * walk; }
      const Listwalk < A > & operator [] (int ) const { return * walk; }
  };

} // namespace ALU2DGrid

#endif // #ifndef __HEADER__LISTWALK
