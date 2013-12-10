c
c
c     ###################################################
c     ##  COPYRIGHT (C)  1992  by  Jay William Ponder  ##
c     ##              All Rights Reserved              ##
c     ###################################################
c
c     ################################################################
c     ##                                                            ##
c     ##  polar.i  --  polarizabilities and induced dipole moments  ##
c     ##                                                            ##
c     ################################################################
c
c
c     polarity  dipole polarizability for each multipole site (Ang**3)
c     thole     Thole polarizability damping value for each site
c     pdamp     value of polarizability scale factor for each site
c     uind      induced dipole components at each multipole site
c     uinp      induced dipoles in field used for energy interactions
c     uinds     GK or PB induced dipoles at each multipole site
c     uinps     induced dipoles in field used for GK or PB energy
c     npolar    total number of polarizable sites in the system
c
c
      integer npolar
      real*8 polarity
      real*8 thole,pdamp
      common /polar/ polarity(maxatm),thole(maxatm),pdamp(maxatm),
     &               npolar
