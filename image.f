c
c
c     ###################################################
c     ##  COPYRIGHT (C)  1990  by  Jay William Ponder  ##
c     ##              All Rights Reserved              ##
c     ###################################################
c
c     ################################################################
c     ##                                                            ##
c     ##  subroutine image  --  compute the minimum image distance  ##
c     ##                                                            ##
c     ################################################################
c
c
c     "image" takes the components of pairwise distance between
c     two points in a periodic box and converts to the components
c     of the minimum image distance
c
c
      subroutine image (xr,yr,zr)
      implicit none
      include 'sizes.i'
      include 'boxes.i'
      include 'cell.i'
      real*8 xr,yr,zr
      real*8 xf,yf,zf
c
c
c     for orthogonal lattice, find the desired image directly
c
      if (orthogonal) then
         do while (abs(xr) .gt. xcell2)
            xr = xr - sign(xcell,xr)
         end do
         do while (abs(yr) .gt. ycell2)
            yr = yr - sign(ycell,yr)
         end do
         do while (abs(zr) .gt. zcell2)
            zr = zr - sign(zcell,zr)
         end do
c
c     for monoclinic lattice, convert "xr" and "zr" to
c     fractional coordinates, find desired image and then
c     translate fractional coordinates back to Cartesian
c
      else if (monoclinic) then
         zf = zr / beta_sin
         xf = xr - zf*beta_cos
         do while (abs(xf) .gt. xcell2)
            xf = xf - sign(xcell,xf)
         end do
         do while (abs(yr) .gt. ycell2)
            yr = yr - sign(ycell,yr)
         end do
         do while (abs(zf) .gt. zcell2)
            zf = zf - sign(zcell,zf)
         end do
         xr = xf + zf*beta_cos
         zr = zf * beta_sin
c
c     for triclinic lattice, convert pairwise components to
c     fractional coordinates, find desired image and then
c     translate fractional coordinates back to Cartesian
c
      else if (triclinic) then
         zf = zr / gamma_term
         yf = (yr - zf*beta_term) / gamma_sin
         xf = xr - yf*gamma_cos - zf*beta_cos
         do while (abs(xf) .gt. xcell2)
            xf = xf - sign(xcell,xf)
         end do
         do while (abs(yf) .gt. ycell2)
            yf = yf - sign(ycell,yf)
         end do
         do while (abs(zf) .gt. zcell2)
            zf = zf - sign(zcell,zf)
         end do
         xr = xf + yf*gamma_cos + zf*beta_cos
         yr = yf*gamma_sin + zf*beta_term
         zr = zf * gamma_term
c
c     for truncated octahedron, use orthogonal box equations,
c     then perform extra tests to remove corner pieces
c
      else if (octahedron) then
         do while (abs(xr) .gt. xbox2)
            xr = xr - sign(xbox,xr)
         end do
         do while (abs(yr) .gt. ybox2)
            yr = yr - sign(ybox,yr)
         end do
         do while (abs(zr) .gt. zbox2)
            zr = zr - sign(zbox,zr)
         end do
         if (abs(xr)+abs(yr)+abs(zr) .gt. box34) then
            xr = xr - sign(xbox2,xr)
            yr = yr - sign(ybox2,yr)
            zr = zr - sign(zbox2,zr)
         end if
      end if
      return
      end
c
c
c     ###############################################################
c     ##                                                           ##
c     ##  subroutine imager  --  replicate minimum image distance  ##
c     ##                                                           ##
c     ###############################################################
c
c
c     "imager" takes the components of pairwise distance between
c     two points in the same or neighboring periodic boxes and
c     converts to the components of the minimum image distance
c
c
      subroutine imager (xr,yr,zr,i)
      implicit none
      include 'sizes.i'
      include 'boxes.i'
      include 'cell.i'
      integer i
      real*8 xr,yr,zr
      real*8 xf,yf,zf
      real*8 xsize,ysize,zsize
      real*8 xsize2,ysize2,zsize2
      real*8 xmove,ymove,zmove
c
c
c     set dimensions for either single box or replicated cell
c
      if (i .ge. 0) then
         xsize = xcell
         ysize = ycell
         zsize = zcell
         xsize2 = xcell2
         ysize2 = ycell2
         zsize2 = zcell2
      else
         xsize = xbox
         ysize = ybox
         zsize = zbox
         xsize2 = xbox2
         ysize2 = ybox2
         zsize2 = zbox2
      end if
c
c     compute the distance to translate along each cell axis
c
      if (i .le. 0) then
         xmove = 0.0d0
         ymove = 0.0d0
         zmove = 0.0d0
      else
         xmove = icell(1,i) * xbox
         ymove = icell(2,i) * ybox
         zmove = icell(3,i) * zbox
      end if
c
c     for orthogonal lattice, find the desired image directly
c
      if (orthogonal) then
         xr = xr + xmove
         do while (abs(xr) .gt. xsize2)
            xr = xr - sign(xsize,xr)
         end do
         yr = yr + ymove
         do while (abs(yr) .gt. ysize2)
            yr = yr - sign(ysize,yr)
         end do
         zr = zr + zmove
         do while (abs(zr) .gt. zsize2)
            zr = zr - sign(zsize,zr)
         end do
c
c     for monoclinic lattice, convert "xr" and "zr" to
c     fractional coordinates, find desired image and then
c     translate fractional coordinates back to Cartesian
c
      else if (monoclinic) then
         zf = zr / beta_sin
         xf = xr - zf*beta_cos
         xf = xf + xmove
         do while (abs(xf) .gt. xsize2)
            xf = xf - sign(xsize,xf)
         end do
         yr = yr + ymove
         do while (abs(yr) .gt. ysize2)
            yr = yr - sign(ysize,yr)
         end do
         zf = zf + zmove
         do while (abs(zf) .gt. zsize2)
            zf = zf - sign(zsize,zf)
         end do
         xr = xf + zf*beta_cos
         zr = zf * beta_sin
c
c     for triclinic lattice, convert pairwise components to
c     fractional coordinates, find desired image and then
c     translate fractional coordinates back to Cartesian
c
      else if (triclinic) then
         zf = zr / gamma_term
         yf = (yr - zf*beta_term) / gamma_sin
         xf = xr - yf*gamma_cos - zf*beta_cos
         xf = xf + xmove
         do while (abs(xf) .gt. xsize2)
            xf = xf - sign(xsize,xf)
         end do
         yf = yf + ymove
         do while (abs(yf) .gt. ysize2)
            yf = yf - sign(ysize,yf)
         end do
         zf = zf + zmove
         do while (abs(zf) .gt. zsize2)
            zf = zf - sign(zsize,zf)
         end do
         xr = xf + yf*gamma_cos + zf*beta_cos
         yr = yf*gamma_sin + zf*beta_term
         zr = zf * gamma_term
c
c     for truncated octahedron, use orthogonal box equations,
c     then perform extra tests to remove corner pieces
c
      else if (octahedron) then
         do while (abs(xr) .gt. xbox2)
            xr = xr - sign(xbox,xr)
         end do
         do while (abs(yr) .gt. ybox2)
            yr = yr - sign(ybox,yr)
         end do
         do while (abs(zr) .gt. zbox2)
            zr = zr - sign(zbox,zr)
         end do
         if (abs(xr)+abs(yr)+abs(zr) .gt. box34) then
            xr = xr - sign(xbox2,xr)
            yr = yr - sign(ybox2,yr)
            zr = zr - sign(zbox2,zr)
         end if
      end if
      return
      end
