      parameter(nx1=600,ny1=300,nx=nx1,ny=ny1,imt=nx1,jmt=ny1,km=40)
c      parameter(nyear1=2013,nyear2=2013,nday=365,nmon=12)
      parameter(nyear1=2009,nyear2=2009,nday=365,nmon=12)

      dimension heff(imt,jmt),aiday(imt,jmt),snow(imt,jmt)
      dimension ossw(imt,jmt),hiday(imt,jmt)
      dimension uice(imt,jmt),vice(imt,jmt),ssh(imt,jmt)
      dimension uo(imt,jmt,km),vo(imt,jmt,km),vdc(imt,jmt,km)
      dimension zoo1(imt,jmt,km),to(imt,jmt,km)
      dimension zoo2(imt,jmt,km),wo(imt,jmt,km)
      dimension zoo3(imt,jmt,km),nitrat(imt,jmt,km)
      dimension flagel(imt,jmt,km),diatom(imt,jmt,km)
      dimension flageli(imt,jmt,km),diatomi(imt,jmt,km)
      dimension clon(imt,jmt),clat(imt,jmt),kmt(imt,jmt)
      dimension ulat(imt,jmt),ulon(imt,jmt),HTN(imt,jmt),HTE(imt,jmt)
     &,HUS(imt,jmt),HUW(imt,jmt),angle(imt,jmt),dxt(imt,jmt)
     &,dyt(imt,jmt)
      dimension dz(km),ZDZ (KM),ZDZZ(KM+1)
      character *80 fopen(5),f1,f2,f3,f4,f5,f6,f7,f8,f9,f10,f11,f12,f13
     &,f14
      character *4 cyear(1900:2100),cyear1(1900:2100)
      character *12 char
      integer SLEN

c      f4='uiday_600_300.H'   ! daily ice velocity (m/s)
c      f5='snowday_600_300.H' ! daily snow depth in water equivalent (m)

      f2='to_600_300.H'        ! daily ocean temperature T (C)
      f3='uo_600_300.H'        ! daily ocean velocity U,V (cm/s)
      f4='wo_600_300.H'        ! daily ocean velocity W (m/s)
      f5='flagel_600_300.H'    ! daily flagellates in ocean (mmol-N/m**3)
      f6='diatom_600_300.H'    ! daily diatom in ocean (mmol-N/m**3)
      f7='zoo1_600_300.H'      ! daily microzooplankton (mmol-N/m**3)
      f8='hiday_600_300.H'     ! daily ice thickness (m)
      f9='aiday_600_300.H'     ! daily ice concentration (fraction)
      f10='osswday_600_300.H'  ! daily shortwave radiation (W/m*m)
      f11='zoo2_600_300.H'     ! daily copepod in ocean (mmol-N/m**3)
      f12='zoo3_600_300.H'     ! daily predatory zooplankton in ocean (mmol-N/m**3)
      f13='nitrat_600_300.H'   ! daily nitrate in ocean (mmol-N/m**3)
      f14='vdcday_600_300.H'   ! daily ocean vertical diffusivity (cm**2/s)

cccc-notes---
c Covert ossw to PAR by multiplying 0.43
cccc

cccc-notes---
c model is on B grid, all data except W are at the center of each ocean level
c W is centered at a scalar grid cell at the base of each ocean level
cccc

c read lon/lat for scalar fields (like ice thickness, ocean T & biology)
c the lon/lat is also for W
      open(20,file='grid.dat.rot')
      read(20,'(10f8.2)') ((clon(i,j),i=1,nx1),j=1,ny1)
      read(20,'(10f8.2)') ((clat(i,j),i=1,nx1),j=1,ny1)
      close(20)

c read lon/lat for vector fields (like ice & ocean veclocities, except W)
      open(24,file='grid.dat.pop')
        read(24,'(10f8.2)') ((ulat(i,j),i=1,nx),j=1,ny)
        read(24,'(10f8.2)') ((ulon(i,j),i=1,nx),j=1,ny)
c HTN, HTE are lengths of the northern and eastern sides of a scaler grid cell in km, HTN*HTE is the area of a scaler grid cell in km**2
c HUS, HUW are lengths of the southern and western sides of a vector grid cell in km
        read(24,'(10f8.2)') ((HTN  (i,j),i=1,nx),j=1,ny) 
        read(24,'(10f8.2)') ((HTE  (i,j),i=1,nx),j=1,ny)
        read(24,'(10f8.2)') ((HUS  (i,j),i=1,nx),j=1,ny)
        read(24,'(10f8.2)') ((HUW  (i,j),i=1,nx),j=1,ny)
c angle is the angle between latitude line and  grid cell x-coordinate line, needed for vector rotation for plotting the vectors in spherical lat-lon coordinate system
c** Do not use the angle variable any more because the rotation has already be made **
        read(24,'(10f8.2)') ((angle(i,j),i=1,nx),j=1,ny)
      close(24)

c read model grid mask with ocean levels; ocean: levels > 0, land: levels = 0
      open(20,file='levels_40_t_aBering1')
      read(20,'(600i2)') kmt
      close(20)

c read ocean level thickness dz(k) in cm
      open(20,file='dz.dta40')
      do k=1,km
        read(20,*)dz(k)
        dz(k)=dz(k)*0.01
      end do
      ZDZ(1)=DZ(1)
      do k=2,km
        ZDZ(K)=ZDZ(K-1)+DZ(K)
      end do
      do k=1,km
        zdzz(k)=zdz(k)-0.5*dz(k)
      end do
      close(20)

      do 999 iyear=nyear1,nyear2

      write(unit=cyear(iyear),fmt='(i4)') iyear
      i=slen(f2)
      open(2,file=f2(1:i)//cyear(iyear)
     &,access='direct',form='unformatted',recl=nx1*ny1*4
     &,status='old')

      i=slen(f3)
      open(3,file=f3(1:i)//cyear(iyear)
     &,access='direct',form='unformatted',recl=nx1*ny1*4
     &,status='old')

      i=slen(f4)
      open(4,file=f4(1:i)//cyear(iyear)
     &,access='direct',form='unformatted',recl=nx1*ny1*4
     &,status='old')
 
      i=slen(f5)
      open(5,file=f5(1:i)//cyear(iyear)
     &,access='direct',form='unformatted',recl=nx1*ny1*4
     &,status='old')

      i=slen(f6)
      open(61,file=f6(1:i)//cyear(iyear)
     &,access='direct',form='unformatted',recl=nx1*ny1*4
     &,status='old')
 
      i=slen(f7)
      open(7,file=f7(1:i)//cyear(iyear)
     &,access='direct',form='unformatted',recl=nx1*ny1*4
     &,status='old')
 
      i=slen(f8)
      open(8,file=f8(1:i)//cyear(iyear)
     &,access='direct',form='unformatted',recl=nx1*ny1*4
     &,status='old')
 
      i=slen(f9)
      open(9,file=f9(1:i)//cyear(iyear)
     &,access='direct',form='unformatted',recl=nx1*ny1*4
     &,status='old')
 
      i=slen(f10)
      open(10,file=f10(1:i)//cyear(iyear)
     &,access='direct',form='unformatted',recl=nx1*ny1*4
     &,status='old')
 
      i=slen(f11)
      open(11,file=f11(1:i)//cyear(iyear)
     &,access='direct',form='unformatted',recl=nx1*ny1*4
     &,status='old')
 
      i=slen(f12)
      open(12,file=f12(1:i)//cyear(iyear)
     &,access='direct',form='unformatted',recl=nx1*ny1*4
     &,status='old')

      i=slen(f13)
      open(13,file=f13(1:i)//cyear(iyear)
     &,access='direct',form='unformatted',recl=nx1*ny1*4
     &,status='old')
 
      i=slen(f14)
      open(14,file=f14(1:i)//cyear(iyear)
     &,access='direct',form='unformatted',recl=nx1*ny1*4
     &,status='old')
 
      do iday=1,nday

        do k=1,km
          read(2)((to(i,j,k),i=1,nx1),j=1,ny1)
        end do

        do k=1,km
          read(3)((uo(i,j,k),i=1,nx1),j=1,ny1)
        end do
        do k=1,km
          read(3)((vo(i,j,k),i=1,nx1),j=1,ny1)
        end do

        do k=1,km
          read(4)((wo(i,j,k),i=1,nx1),j=1,ny1)
        end do

        do k=1,km
          read(5)((flagel(i,j,k),i=1,nx1),j=1,ny1) ! flagellates
        end do
        do k=1,km
          read(61)((diatom(i,j,k),i=1,nx1),j=1,ny1) ! diatom
        end do
        do k=1,km
          read(7)((zoo1(i,j,k),i=1,nx1),j=1,ny1) ! microzooplankton
        end do

        read(8)((hiday(i,j),i=1,nx1),j=1,ny1) ! ice thickness
        read(9)((aiday(i,j),i=1,nx1),j=1,ny1) ! ice concentration
        read(10)((ossw(i,j),i=1,nx1),j=1,ny1) ! shortwave radiation

        do k=1,km
          read(11)((zoo2(i,j,k),i=1,nx1),j=1,ny1) ! copepod
        end do
        do k=1,km
          read(12)((zoo3(i,j,k),i=1,nx1),j=1,ny1) ! predatory zoop
        end do
        do k=1,km
          read(13)((nitrat(i,j,k),i=1,nx1),j=1,ny1) ! nitrate
        end do
        do k=1,km
          read(14)((vdc(i,j,k),i=1,nx1),j=1,ny1) ! vertical diffusivity
        end do

        i=nx/2
        j=ny/2
        write(*,'(2i6, 2f12.3, 2f8.2)') iyear, iday
     &, clon(i,j), clat(i,j), to(i,j,1), flagel(i,j,1)

      end do ! iday

      close(2)
      close(3)
      close(4)
      close(5)
      close(61)
      close(7)
      close(8)
      close(9)
      close(10)
      close(11)
      close(12)
      close(13)
      close(14)
999   continue
      
      stop  
      end

      INTEGER FUNCTION slen (string)
C ---
C --- this function computes the length of a character string less
C --- trailing blanks
C --- slen > 0, length of string less trailing blanks
C ---      = 0, character string is blank
C ---
      CHARACTER*(*) string
      CHARACTER*1 cblank
      INTEGER i
      DATA cblank/' '/
C ---
      DO 50 i = LEN(string), 1, -1
         IF (string(i:i) .NE. ' ')  GO TO 100
50    CONTINUE
      i = 0
100   CONTINUE
      slen = i
      RETURN
      END
