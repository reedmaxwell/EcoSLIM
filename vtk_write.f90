! Copyright 2010 Reed M. Maxwell
!
! This file is part of EcoSLIM; Originally in SLIM-Fast.
!
!    SLIM-Fast is free software: you can redistribute it and/or modify
!    it under the terms of the GNU General Public License as published by
!    the Free Software Foundation, either version 3 of the License.
!
!    SLIM-Fast is distributed in the hope that it will be useful,
!    but WITHOUT ANY WARRANTY; without even the implied warranty of
!    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
!    GNU General Public License for more details.
!
!    You should have received a copy of the GNU General Public License
!    along with SLIM-Fast in /src/gpl.txt.  If not, see <http://www.gnu.org/licenses/>.

SUBROUTINE vtk_write(time,x,conc_header,ixlim,iylim,izlim,icycle,n_constituents,Pnts,vtk_file)
real*8    :: time
REAL*8    :: x(:,:,:,:)
CHARACTER (LEN=20)     :: conc_header(:)
INTEGER*4 :: ixlim
INTEGER*4 :: iylim
INTEGER*4 :: izlim
REAL*8                 :: dx
REAL*8                 :: dy
REAL*8                 :: dz(izlim)
REAL*8                 :: Pnts(:,:)
INTEGER                :: icycle
INTEGER                :: n_constituents
CHARACTER (LEN=200)    :: vtk_file

INTEGER*4 i,j,k, ijk, debug,l, nxyz,nxyzp1
CHARACTER*1 lf
CHARACTER*12 num1, num2, num3
CHARACTER*8 ctime
real*8 number

debug = 0
!
!      Open File
!
Write(ctime,'(i8.8)') icycle
OPEN(15,FILE=trim(vtk_file)//'.'//ctime//'.vtk',FORM='unformatted',  &
    access='stream',convert='BIG_ENDIAN')
!
!      Write header info
!
lf = char(10) ! line feed character
Write(15) "# vtk DataFile Version 3.0"//lf
Write(15) "EcoSLIM Output"//lf
Write(15) "BINARY"//lf
Write(15) "DATASET STRUCTURED_GRID"//lf
!Write(15) "FIELD FieldData 1"//lf
!write(15) "TIME 1 1 double"//lf
!write(ctime, '(f8.3)') time
!write(15) ctime
!write(15) time
!write(15) lf
!Write(15) "CYCLE 1 1 int"//lf
!Write(num1,'(i12)') icycle
!Write(15) num1
!write(15) icycle
!Write(15) lf

Write(num1,'(i12)') ixlim+1
Write(num2,'(i12)') iylim+1
Write(num3,'(i12)') izlim+1

write(15) "DIMENSIONS "//num1//" "//num2//" "//num3//lf
nxyzp1 = (ixlim+1)*(iylim+1)*(izlim+1)
write(num1, '(i12)') nxyzp1
Write(15) "POINTS "//num1//" float"//lf
write(15) ((real(Pnts(i,j),kind=4), j=1,3), i=1,nxyzp1)  ! This forces the expected write order
!write(15) lf
nxyz = (ixlim)*(iylim)*(izlim)
write(num1, '(i12)') nxyz
Write(15) "CELL_DATA "//num1//lf

do l = 1, n_constituents
write(15) "SCALARS "//trim(conc_header(l))//" float"//lf
Write(15) "LOOKUP_TABLE default"//lf

write(15) (((real(X(l,i,j,k),kind=4), i=1,ixlim), j=1,iylim), k=1,izlim)

!DO  k=1,izlim
!  DO  j=1,iylim
!    DO  i=1,ixlim
!        WRITE(15) x(l,i,j,k)
!    END DO
!  END DO
!END DO
!write(15) lf
end do ! l, n_constits
CLOSE(15)


RETURN
END SUBROUTINE vtk_write
