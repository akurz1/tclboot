subroutine digraph_adj_components (adj, lda, nnode, ncomp, comp, dad, order,ier)
!*****************************************************************************80
! Subroutine adapted from
! Burkardt, John. (2020, November). GRAFPACK - Graph Computations 
! [FORTRAN90 library]. Retrieved from
! \url{https://people.sc.fsu.edu/~jburkardt/f_src/grafpack/grafpack.html}
!*****************************************************************************80
!
!! DIGRAPH_ADJ_COMPONENTS finds the strongly connected components of a digraph.
!
!  Discussion:
!
!    A digraph is a directed graph.
!
!    A strongly connected component of a directed graph is the largest
!    set of nodes such that there is a directed path from any node to 
!    any other node in the same component.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    15 April 1999
!
!  Reference:
!
!    K Thulasiraman, M Swamy,
!    Graph Theory and Algorithms,
!    John Wiley, New York, 1992.
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) ADJ(LDA,NNODE), the adjacency information.
!    ADJ(I,J) is 1 if there is a direct link from node I to node J.
!
!    Input, integer ( kind = 4 ) LDA, the leading dimension of ADJ.
!
!    Input, integer ( kind = 4 ) NNODE, the number of nodes.
!
!    Output, integer ( kind = 4 ) NCOMP, the number of strongly connected 
!    components.
!
!    Output, integer ( kind = 4 ) COMP(NNODE), lists the connected component 
!    to which each node belongs.
!
!    Output, integer ( kind = 4 ) DAD(NNODE), the father array for the depth 
!    first search trees.  DAD(I) = 0 means that node I is the root of 
!    one of the trees.  DAD(I) = J means that the search descended
!    from node J to node I.
!
!    Output, integer ( kind = 4 ) ORDER(NNODE), the order in which the nodes 
!    were traversed, from 1 to NNODE.
!
!        ier: output error code  added
!             ier = 0: O.K.
!                   1: Stop Illegal stack reference.'

  implicit none
  
  integer(kind=4), intent(out)     :: ier    ! added error code
  integer ( kind = 4 ) lda
  integer ( kind = 4 ) nnode

  integer ( kind = 4 ) adj(lda,nnode)
  integer ( kind = 4 ) comp(nnode)
  integer ( kind = 4 ) dad(nnode)
  integer ( kind = 4 ) iorder
  integer ( kind = 4 ) lowlink(nnode)
  integer ( kind = 4 ) mark(nnode)
  integer ( kind = 4 ) ncomp
  integer ( kind = 4 ) nstack
  integer ( kind = 4 ) order(nnode)
  integer ( kind = 4 ) point(nnode)
  integer ( kind = 4 ) stack(nnode)
  integer ( kind = 4 ) v
  integer ( kind = 4 ) w
  integer ( kind = 4 ) x
!
!  Initialization.
!
  comp(1:nnode) = 0
  dad(1:nnode) = 0
  order(1:nnode) = 0
  lowlink(1:nnode) = 0
  mark(1:nnode) = 0
  point(1:nnode) = 0

  iorder = 0
  nstack = 0
  ncomp = 0
  
  ier=0 ! added error code
!
!  Select any node V not stored in the stack, that is, with MARK(V) = 0.
!
  do

    v = 0

    do

      v = v + 1

      if ( nnode < v ) then
        adj(1:nnode,1:nnode) = abs ( adj(1:nnode,1:nnode) )
        return
      end if

      if ( mark(v) /= 1 ) then
        exit
      end if

    end do

    iorder = iorder + 1

    order(v) = iorder
    lowlink(v) = iorder
    mark(v) = 1
 
    nstack = nstack + 1
    stack(nstack) = v
    point(v) = 1

30  continue
!
!  Consider each node W.
!
    do w = 1, nnode
!
!  Is there an edge (V,W) and has it not been examined yet?
!
      if ( 0 < adj(v,w) ) then

        adj(v,w) = - adj(v,w)
!
!  Is the node on the other end of the edge undiscovered yet?
!
        if ( mark(w) == 0 ) then

          iorder = iorder + 1
          order(w) = iorder
          lowlink(w) = iorder
          dad(w) = v
          mark(w) = 1

          nstack = nstack + 1
          stack(nstack) = w
          point(w) = 1

          v = w

        else if ( mark(w) == 1 ) then

          if ( order(w) < order(v) .and. point(w) == 1 ) then
            lowlink(v) = min ( lowlink(v), order(w) )
          end if

        end if

        go to 30

      end if

    end do

    if ( lowlink(v) == order(v) ) then

      ncomp = ncomp + 1

      do

        if ( nstack <= 0 ) then
 !         write ( *, '(a)' ) ' '
 !         write ( *, '(a)' ) 'DIGRAPH_ADJ_COMPONENTS - Fatal error!'
 !         write ( *, '(a)' ) '  Illegal stack reference.'
 !         stop
         ier = 1 ! added error code, replaces I/O operation 
          return
        end if

        x = stack(nstack)
        nstack = nstack - 1

        point(x) = 0
        comp(x) = ncomp

        if ( x == v ) then
          exit
        end if

      end do

    end if

    if ( dad(v) /= 0 ) then
      lowlink(dad(v)) = min ( lowlink(dad(v)), lowlink(v) )
      v = dad(v)
      go to 30
    end if

  end do

  return
end
