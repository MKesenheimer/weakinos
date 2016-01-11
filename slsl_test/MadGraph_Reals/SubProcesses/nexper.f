      subroutine nexper(n,a,mtc,even)
c*******************************************************************************
c     Gives all permutations for a group                                       
c http://www.cs.sunysb.edu/~algorith/implement/wilf/distrib/processed/nexper_2.f
c next permutation of {1,...,n}. Ref NW p 59.                                  
c*******************************************************************************
      integer a(n),s,d
      logical mtc,even
      if(mtc)goto 10
      nm3=n-3
      do 1 i=1,n
 1        a(i)=i
      mtc=.true.
 5     even=.true.
      if(n.eq.1)goto 8
 6     if(a(n).ne.1.or.a(1).ne.2+mod(n,2))return
      if(n.le.3)goto 8
      do 7 i=1,nm3
      if(a(i+1).ne.a(i)+1)return
 7     continue
 8      mtc=.false.
      return
 10    if(n.eq.1)goto 27
      if(.not.even)goto 20
      ia=a(1)
      a(1)=a(2)
      a(2)=ia
      even=.false.
      goto 6
 20    s=0
      do 26 i1=2,n
 25       ia=a(i1)
      i=i1-1
      d=0
      do 30 j=1,i
 30       if(a(j).gt.ia) d=d+1
      s=d+s
      if(d.ne.i*mod(s,2)) goto 35
 26    continue
 27     a(1)=0
      goto 8
 35    m=mod(s+1,2)*(n+1)
      do 40 j=1,i
      if(isign(1,a(j)-ia).eq.isign(1,a(j)-m))goto 40
      m=a(j)
      l=j
 40    continue
      a(l)=ia
      a(i1)=m
      even=.true.
      return
      end



      subroutine switchmom(p,p1,ic,nexternal)
      implicit none
      integer nexternal
      integer ic(nexternal)
      real*8 p1(0:3,nexternal),p(0:3,nexternal)
      integer i,j
      do i=1,nexternal
         do j=0,3
            p1(j,ic(i))=p(j,i)
         enddo
      enddo
      end

      subroutine switchlegs(legs,legs1,ic,nexternal)
      implicit none
      integer nexternal
      integer ic(nexternal)
      integer legs1(nexternal),legs(nexternal)
      integer i,j
      do i=1,nexternal
         legs1(ic(i))=legs(i)
      enddo
      end

      subroutine switchborns(born,born1,bornjk,bornjk1,bmunu,bmunu1,
     &     ic,npart)
      implicit none
      integer npart
      integer ic(npart)
      double precision born(npart),born1(npart),bmunu(0:3,0:3,npart),
     &    bmunu1(0:3,0:3,npart),bornjk(npart,npart),bornjk1(npart,npart)
      integer i,j,k,mu,nu
      do j=1,npart
         born1(ic(j))=born(j)
         do mu=0,3
            do nu=0,3
               bmunu1(mu,nu,ic(j))=bmunu(mu,nu,j)
            enddo
         enddo
         do k=1,npart
            bornjk1(ic(j),ic(k))=bornjk(j,k)
         enddo
      enddo
      end

      subroutine switchcolor(col,col1,ic,nexternal)
      implicit none
      integer nexternal
      integer ic(nexternal)
      integer col1(2,nexternal),col(2,nexternal)
      integer i,j,k
      do i=1,2
         do j=1,nexternal
            col1(i,ic(j))=col(i,j)
         enddo
      enddo
      end
