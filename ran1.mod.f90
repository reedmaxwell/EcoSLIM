module ran_mod
contains
     function ran1(idum)
        implicit none  !note after use statement
        real*8 ran1
        integer, intent(inout), optional :: idum
        real*8 r(97),rm1,rm2
        integer, parameter :: m1=259200,ia1=7141,ic1=54773
        integer, parameter :: m2=134456,ia2=8121,ic2=28411
        integer, parameter :: m3=243000,ia3=4561,ic3=51349
        integer j
        integer iff,ix1,ix2,ix3
        data iff /0/
        !$OMP THREADPRIVATE(iff,ix1,ix2,ix3,j,r,rm1,rm2)
        save iff,ix1,ix2,ix3,j,r,rm1,rm2
        if(present(idum))then
          if (idum<0.or.iff.eq.0)then
            rm1=1.0/m1
            rm2=1.0/m2
            iff=1
            ix1=mod(ic1-idum,m1)
            ix1=mod(ia1*ix1+ic1,m1)
            ix2=mod(ix1,m2)
            ix1=mod(ia1*ix1+ic1,m1)
            ix3=mod(ix1,m3)
            do j=1,97
                ix1=mod(ia1*ix1+ic1,m1)
                ix2=mod(ia2*ix2+ic2,m2)
                r(j)=(real(ix1)+real(ix2)*rm2)*rm1
            enddo
            idum=1
          endif
        endif
        ix1=mod(ia1*ix1+ic1,m1)
        ix2=mod(ia2*ix2+ic2,m2)
        ix3=mod(ia3*ix3+ic3,m3)
        j=1+(97*ix3)/m3
        if(j>97.or.j<1)then
            write(*,*)' error in ran1 j=',j
            stop
        endif
        ran1=r(j)
        r(j)=(real(ix1)+real(ix2)*rm2)*rm1
        return
     end function ran1
end module ran_mod
