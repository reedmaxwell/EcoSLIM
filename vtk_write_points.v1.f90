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

SUBROUTINE vtk_write_points(P,np_active, np,icycle,vtk_file)
REAL*8    :: P(:,:)
INTEGER                :: icycle
INTEGER*4              :: np_active
INTEGER*4              :: np
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
Write(15) "# vtk DataFile Version 2.0"//lf
Write(15) "EcoSLIM Points Output"//lf
Write(15) "BINARY"//lf
Write(15) "DATASET POLYDATA"//lf

write(num1, '(i12)') np_active
Write(15) "POINTS "//num1//" FLOAT"//lf
write(15) ((real(P(j,i),kind=4), i=1,3), j=1,np_active)  ! This forces the expected write order
!do j =1, np_active
!do i =1 , 3
!  write(15) real(P(j,i),kind=4)
!end do  !i
!end do !!j
          write(15) lf
write(15) "POINT_DATA "//num1//lf
write(15) "SCALARS Time float"//lf
Write(15) "LOOKUP_TABLE default"//lf
write(15) (real(P(j,4),kind=4), j=1,np_active)
write(15) lf
write(15) "SCALARS Mass float"//lf
Write(15) "LOOKUP_TABLE default"//lf
write(15) (real(P(j,6),kind=4), j=1,np_active)
write(15) lf
write(15) "SCALARS Source float"//lf
Write(15) "LOOKUP_TABLE default"//lf
write(15) (real(P(j,7),kind=4), j=1,np_active)
write(15) lf
!write(15) "SCALARS Delta float"//lf
!Write(15) "LOOKUP_TABLE default"//lf
!write(15) (real(P(j,9),kind=4), j=1,np_active)
!write(15) lf
CLOSE(15)


RETURN
END SUBROUTINE vtk_write_points
