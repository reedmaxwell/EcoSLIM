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

SUBROUTINE vtk_write_points(P,np_active, np,icycle,vtk_file, dx, dy,nx,ny,maxZ, dem)
REAL*8                 :: P(:,:)
INTEGER                :: icycle
INTEGER*4              :: np_active
INTEGER*4              :: np
REAL*8                 :: dx
REAL*8                 :: dy
INTEGER*4              :: nx
INTEGER*4              :: ny
REAL*8                 :: maxZ
REAL*8                 :: DEM(:,:)
CHARACTER (LEN=200)    :: vtk_file

INTEGER*4 i,j,k, ijk, debug,l, nxyz,nxyzp1,Px, Py
CHARACTER*1 lf
CHARACTER*12 num1, num2, num3
CHARACTER*8 ctime
real*8 number,X,Clocx, Clocy, Dxl, Dxu, Dyl, Dyu

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
!write(15) ((real(P(j,i),kind=4), i=1,3), j=1,np_active)  ! This forces the expected write order
do j =1, np_active
  write(15) real(P(j,1:2),kind=4)
  ! find integer cell location
  Px = floor(P(j,1) / dx)
  Py = floor(P(j,2) / dy)
  ! Find each particle's factional cell location
  Clocx = (P(j,1) - float(Px)*dx)  / dx
  Clocy = (P(j,2) - float(Py)*dy)  / dy
  ! X is local adjustment for DEM
  ! linearly interpolate particle location in DEM cell
  Dxl = DEM(Px+1,Py+1)
  Dxu = DEM(Px+2,Py+1)
  Dyl = DEM(Px+1,Py+1)
  Dyu = DEM(Px+1,Py+2)
  ! correct for edges
  if (Px == nx -1) Dxu = Dxl
  if (Py == ny -1) Dyu = Dyl
  X = ( ((1.0d0-Clocx)*Dxl &
        + Dxu*Clocx)  +    &
       ((1.0d0-Clocy)*Dyl &
                + Dyu*Clocy) ) / 2.0D0
  ! write new location
  write(15) real(P(j,3)+X -maxZ,kind=4)
end do !!j
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
