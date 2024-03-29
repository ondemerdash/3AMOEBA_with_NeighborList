c
c     ###############################################################
c     ##                                                           ##
c     ##               Subroutine empole1c_3b                      ##
c     ##  Liam O-Suilleabhain, Omar Demerdash, Teresa Head-Gordon  ##
c     ##                 Spring 2013                               ##
c     ###############################################################
c
c
c     "empole1c_3b" calculates the atomic multipole and dipole
c     polarizability interaction energy using the 3-body approximation
c
c
      subroutine empole1c_3b
      implicit none
      include 'sizes.i'
      include 'energi.i'
c      include 'atmtyp.i'
c      include 'atoms.i'
      include 'molcul.i'
c      include 'polar.i'
      include 'deriv.i'
c      include 'combo.i'
      include 'mpole.i'
      include 'virial.i'
      include 'openmp.i'
c      real*8, allocatable :: field(:,:)
c      real*8, allocatable :: fieldp(:,:)
c      real*8  ep3b,dep3b(3,npole)
      real*8  ep1moli1,ep1moli2,ep2moli12,ep1moli3
      real*8  dep1moli1(3,30),dep1moli2(3,30),dep1moli3(3,30)
      real*8  dep2moli12(3,30)
      integer i,ii,j,l1,i1,i2,i3
c      real*8 xcm1,ycm1,zcm1,xcm2,ycm2,zcm2
c      real*8 xcm3,ycm3,zcm3,xcm,ycm,zcm
c      real*8 M,M1,M2,M3,r_123
c      real*8 virep3b(3,3),eptemp,deptemp(3,30)

      real*8 eptemp,deptemp(3,npole)
c      real*8 eptemp,deptemp(3,30)
      real*8 vir1moli1(3,3),vir2moli12(3,3),vir1moli3(3,3)
      real*8 vir1moli2(3,3),virtemp(3,3)
      integer pnum(30),npole3b,moli1,moli2,moli3
      integer pnum2(npole),nmol2,np1,np2,np3 
      real*8 ep3bt,dep3bt(3,npole),virep3bt(3,3)
      integer omp_get_num_threads,nthread2,tid
      integer omp_get_thread_num,triplecount
      real*8 findr3_minimage2func,r_123
      real*8 findr2_minimage2func,ewaldcut3b
      external findr3_minimage2func
      external findr2_minimage2func
      
      triplecount=0
c      ep3b = 0.0d0
      ep3bt=0.0d0 
c      allocate (field(3,npole))
c      allocate (fieldp(3,npole))

      do i = 1, npole
        do j = 1, 3
c          field(j,i) = 0.0d0
c          fieldp(j,i) = 0.0d0
c          dep3b(j,i) = 0.0d0
          dep3bt(j,i) = 0.0d0
c          dep1moli1(j,i)=0.0d0
c          dep2moli12(j,i)=0.0d0
c          dep1moli2(j,i)=0.0d0
        end do
      end do

c        do i=1,npole
c          do j=1,3
c            dep1moli1(j,i)=0.0d0
c          end do 
c        end do   

      do i=1,3
         do j=1,3
c           virep3b(i,j)=0.0d0
           virep3bt(i,j)=0.0d0
c           vir1moli1(i,j)=0.0d0
c           vir2moli12(i,j)=0.0d0
c           vir1moli2(i,j)=0.0d0
         end do
      end do 
c
c    Usual 2-body permanent electrostatics but without polarization
c
c      call empole1c_3b_Perm(field,fieldp)
      call empole1c_3b_Perm
c
c    Begin 2-body and 3-body approximation
c
      nmol2=nmol
      moli1=0
      moli2=0
      moli3=0
ccccc !$OMP PARALLEL DO default(private) shared(nmol2,imol)
ccccc !$OMP& reduction(+:ep3bt,virep3bt,dep3bt) schedule(dynamic)

!$OMP PARALLEL default(private) shared(nmol,imol,ep3bt,
!$OMP& virep3bt,dep3bt)
!$OMP DO reduction(+:ep3bt,virep3bt,dep3bt)
!$OMP& schedule(dynamic)
      do moli1 = 1, nmol
c
c    1-body polarization
c
        npole3b=3
        pnum(1)=imol(1,moli1)
        pnum(2)=imol(1,moli1)+1
        pnum(3)=imol(2,moli1)

c        call ewaldcof(aewald3b,30.0d0) 
c        jmax3b=int(30.0d0*aewald3b*xbox*aewald3b/pi)
c        kmax3b=int(30.0d0*aewald3b*ybox*aewald3b/pi)
c        lmax3b=int(30.0d0*aewald3b*zbox*aewald3b/pi)
        ewaldcut3b=25.0d0
        call empole1c_3b_Polar(ewaldcut3b,npole3b,pnum,eptemp,
     &   deptemp,virtemp)
        ep3bt = ep3bt + eptemp
        ep1moli1 = eptemp
        do l1 = 1, npole3b
          i = pnum(l1)
          do j = 1, 3
            dep3bt(j,i) = dep3bt(j,i)+deptemp(j,i)
            dep1moli1(j,l1) = deptemp(j,i)
          end do
        end do
        do i=1,3
           do j=1,3
              virep3bt(i,j) = virep3bt(i,j)+virtemp(i,j)
              vir1moli1(i,j) = virtemp(i,j)
           end do
        end do
 
        do moli2 = moli1+1, nmol

          npole3b=6
          np1=3
          np2=6
          pnum(1)=imol(1,moli1)
          pnum(2)=imol(1,moli1)+1
          pnum(3)=imol(2,moli1)
          pnum(4)=imol(1,moli2)
          pnum(5)=imol(1,moli2)+1
          pnum(6)=imol(2,moli2)

          r_123=findr2_minimage2func(pnum,np1,np2)

         if(r_123 .lt. 25.0d0) then
          npole3b=3
          pnum(1)=imol(1,moli1)
          pnum(2)=imol(1,moli1)+1
          pnum(3)=imol(2,moli1)
          ep3bt = ep3bt - ep1moli1 
          do l1 = 1, npole3b
            i = pnum(l1)
            do j = 1, 3
              dep3bt(j,i) = dep3bt(j,i)-dep1moli1(j,l1) 
            end do
          end do
          do i=1,3
             do j=1,3
               virep3bt(i,j)=virep3bt(i,j)-vir1moli1(i,j)
             end do
          end do

          npole3b=6
          pnum(4)=imol(1,moli2)
          pnum(5)=imol(1,moli2)+1
          pnum(6)=imol(2,moli2)

        ewaldcut3b=r_123+10.0d0
        call empole1c_3b_Polar(ewaldcut3b,npole3b,pnum,eptemp,
     &   deptemp,virtemp)
          ep3bt = ep3bt + eptemp
          ep2moli12=eptemp
          do l1 = 1, npole3b
            i = pnum(l1)
            do j = 1, 3
              dep3bt(j,i) = dep3bt(j,i)+deptemp(j,i)
              dep2moli12(j,l1)=deptemp(j,i)
            end do
          end do
          do i=1,3
             do j=1,3
                virep3bt(i,j)=virep3bt(i,j)+virtemp(i,j)
                vir2moli12(i,j)=virtemp(i,j)
             end do
          end do

          npole3b=3
          pnum(1)=imol(1,moli2)
          pnum(2)=imol(1,moli2)+1
          pnum(3)=imol(2,moli2)

        ewaldcut3b=25.0d0
        call empole1c_3b_Polar(ewaldcut3b,npole3b,pnum,eptemp,
     &      deptemp,virtemp)
          ep3bt = ep3bt - eptemp
          ep1moli2 = eptemp
          do l1 = 1, npole3b
            i = pnum(l1)
            do j = 1, 3
              dep3bt(j,i) = dep3bt(j,i)-deptemp(j,i)
              dep1moli2(j,l1) = deptemp(j,i)
            end do
          end do
          do i=1,3
             do j=1,3
               virep3bt(i,j)=virep3bt(i,j)-virtemp(i,j)
               vir1moli2(i,j)=virtemp(i,j)
             end do
          end do 
         else
          npole3b=6
          pnum(1)=imol(1,moli1)
          pnum(2)=imol(1,moli1)+1
          pnum(3)=imol(2,moli1)
          pnum(4)=imol(1,moli2)
          pnum(5)=imol(1,moli2)+1
          pnum(6)=imol(2,moli2)
          ep2moli12=0.0d0
          do l1 = 1, npole3b
            i = pnum(l1)
            do j = 1, 3
              dep2moli12(j,l1)=0.0d0
            end do
          end do
          do i=1,3
             do j=1,3
                vir2moli12(i,j)=0.0d0
             end do
          end do
          npole3b=3
          pnum(1)=imol(1,moli2)
          pnum(2)=imol(1,moli2)+1
          pnum(3)=imol(2,moli2)
        ewaldcut3b=25.0d0
        call empole1c_3b_Polar(ewaldcut3b,npole3b,pnum,eptemp,
     &   deptemp,virtemp)
          ep1moli2 = eptemp
          do l1 = 1, npole3b
            i = pnum(l1)
            do j = 1, 3
              dep1moli2(j,l1) = deptemp(j,i)
            end do
          end do
          do i=1,3
             do j=1,3
               vir1moli2(i,j)=virtemp(i,j)
             end do
          end do
         end if

          do moli3 = moli2+1, nmol
            npole3b=9
            np1=3
            np2=6
            np3=9
            pnum(1)=imol(1,moli1)
            pnum(2)=imol(1,moli1)+1
            pnum(3)=imol(2,moli1)
            pnum(4)=imol(1,moli2)
            pnum(5)=imol(1,moli2)+1
            pnum(6)=imol(2,moli2)
            pnum(7)=imol(1,moli3)
            pnum(8)=imol(1,moli3)+1
            pnum(9)=imol(2,moli3)
            r_123=findr3_minimage2func(pnum,np1,np2,np3)
           if(r_123 .lt. 25.0d0) then
            triplecount=triplecount+1
            npole3b=9
            np1=3
            np2=6
            np3=9
            pnum(1)=imol(1,moli1)
            pnum(2)=imol(1,moli1)+1
            pnum(3)=imol(2,moli1)
            pnum(4)=imol(1,moli2)
            pnum(5)=imol(1,moli2)+1
            pnum(6)=imol(2,moli2)
            pnum(7)=imol(1,moli3)
            pnum(8)=imol(1,moli3)+1
            pnum(9)=imol(2,moli3)
            ewaldcut3b=r_123+10.0d0
        call empole1c_3b_Polar(ewaldcut3b,npole3b,pnum,eptemp,
     &      deptemp,virtemp)
            ep3bt = ep3bt + eptemp
            do l1 = 1, npole3b
              i = pnum(l1)
              do j = 1, 3
                dep3bt(j,i) = dep3bt(j,i)+deptemp(j,i)
              end do
            end do
            do i=1,3
               do j=1,3
                 virep3bt(i,j)=virep3bt(i,j)+virtemp(i,j) 
               end do
            end do

            npole3b=3
            pnum(1)=imol(1,moli3)
            pnum(2)=imol(1,moli3)+1
            pnum(3)=imol(2,moli3)
            ewaldcut3b=25.0d0
            call empole1c_3b_Polar(ewaldcut3b,npole3b,pnum,eptemp,
     &      deptemp,virtemp)
            ep1moli3=eptemp

            do l1 = 1, npole3b
              i = pnum(l1)
              do j = 1, 3
               dep1moli3(j,l1)=deptemp(j,i)
              end do
            end do

            do i=1,3
              do j=1,3
               vir1moli3(i,j)=virtemp(i,j)
              end do
            end do


            if(ep2moli12.ne.0.0d0) then
              pnum(1)=imol(1,moli1)
              pnum(2)=imol(1,moli1)+1
              pnum(3)=imol(2,moli1)
              pnum(4)=imol(1,moli2)
              pnum(5)=imol(1,moli2)+1
              pnum(6)=imol(2,moli2)

              npole3b=6
              ep3bt = ep3bt - ep2moli12
              do l1 = 1, npole3b
                i = pnum(l1)
                do j = 1, 3
                 dep3bt(j,i) = dep3bt(j,i)-dep2moli12(j,l1)
                end do
              end do
              do i=1,3
                do j=1,3
                 virep3bt(i,j)=virep3bt(i,j)-vir2moli12(i,j)
                end do 
              end do             
              npole3b=3
              ep3bt = ep3bt + ep1moli1
              do l1 = 1, npole3b
                 i = pnum(l1)
                 do j = 1, 3
                  dep3bt(j,i) = dep3bt(j,i)+dep1moli1(j,l1)
                 end do
              end do
              do i=1,3
                 do j=1,3
                 virep3bt(i,j)=virep3bt(i,j)+vir1moli1(i,j)
                 end do
              end do
              pnum(1)=imol(1,moli2)
              pnum(2)=imol(1,moli2)+1
              pnum(3)=imol(2,moli2)
              npole3b=3
              ep3bt = ep3bt + ep1moli2
              do l1 = 1, npole3b
                i = pnum(l1)
                do j = 1, 3
                dep3bt(j,i) = dep3bt(j,i)+dep1moli2(j,l1)
                end do
              end do
              do i=1,3
                do j=1,3
                 virep3bt(i,j)=virep3bt(i,j)+vir1moli2(i,j)
                end do
              end do
            end if
            npole3b=6
            pnum(1)=imol(1,moli1)
            pnum(2)=imol(1,moli1)+1
            pnum(3)=imol(2,moli1)
            pnum(4)=imol(1,moli3)
            pnum(5)=imol(1,moli3)+1
            pnum(6)=imol(2,moli3)
            r_123=findr2_minimage2func(pnum,np1,np2)
            if(r_123 .lt. 25.0d0) then
              npole3b=6
              pnum(1)=imol(1,moli1)
              pnum(2)=imol(1,moli1)+1
              pnum(3)=imol(2,moli1)
              pnum(4)=imol(1,moli3)
              pnum(5)=imol(1,moli3)+1
              pnum(6)=imol(2,moli3)
              ewaldcut3b=r_123+10.0d0
              call empole1c_3b_Polar(ewaldcut3b,npole3b,pnum,eptemp,
     &        deptemp,virtemp)
              ep3bt = ep3bt - eptemp
              do l1 = 1, npole3b
                i = pnum(l1)
                do j = 1, 3
                  dep3bt(j,i) = dep3bt(j,i)-deptemp(j,i)
                end do
              end do
              do i=1,3
                do j=1,3
                 virep3bt(i,j)=virep3bt(i,j)-virtemp(i,j)
                end do
              end do
               npole3b=3
              ep3bt = ep3bt + ep1moli1
              do l1 = 1, npole3b
                i = pnum(l1)
                do j = 1, 3
                 dep3bt(j,i) = dep3bt(j,i)+dep1moli1(j,l1)
                end do
              end do
              do i=1,3
               do j=1,3
                virep3bt(i,j)=virep3bt(i,j)+vir1moli1(i,j)
               end do
              end do
              npole3b=3
              pnum(1)=imol(1,moli3)
              pnum(2)=imol(1,moli3)+1
              pnum(3)=imol(2,moli3)
              ep3bt = ep3bt + ep1moli3
c              ep1moli3=eptemp
              do l1 = 1, npole3b
                i = pnum(l1)
                do j = 1, 3
                  dep3bt(j,i) = dep3bt(j,i)+dep1moli3(j,l1)
c                  dep1moli3(j,l1)=deptemp(j,i)
                end do
              end do
              do i=1,3
               do j=1,3
                 virep3bt(i,j)=virep3bt(i,j)+vir1moli3(i,j)
c                 vir1moli3(i,j)=virtemp(i,j)
               end do
              end do
            end if  

            npole3b=6
            pnum(1)=imol(1,moli2)
            pnum(2)=imol(1,moli2)+1
            pnum(3)=imol(2,moli2)
            pnum(4)=imol(1,moli3)
            pnum(5)=imol(1,moli3)+1
            pnum(6)=imol(2,moli3)
            r_123=findr2_minimage2func(pnum,np1,np2)
            if(r_123 .lt. 25.0d0) then
              npole3b=6
              pnum(1)=imol(1,moli2)
              pnum(2)=imol(1,moli2)+1
              pnum(3)=imol(2,moli2)
              pnum(4)=imol(1,moli3)
              pnum(5)=imol(1,moli3)+1
              pnum(6)=imol(2,moli3)

              ewaldcut3b=r_123+10.0d0
              call empole1c_3b_Polar(ewaldcut3b,npole3b,pnum,eptemp,
     &        deptemp,virtemp)
              ep3bt = ep3bt - eptemp
              do l1 = 1, npole3b
                i = pnum(l1)
                do j = 1, 3
                 dep3bt(j,i) = dep3bt(j,i)-deptemp(j,i)
                end do
              end do
              do i=1,3
                do j=1,3
                 virep3bt(i,j)=virep3bt(i,j)-virtemp(i,j)
                end do
              end do

              npole3b=3
              ep3bt = ep3bt + ep1moli2
              do l1 = 1, npole3b
                i = pnum(l1)
                do j = 1, 3
                 dep3bt(j,i) = dep3bt(j,i)+dep1moli2(j,l1)
                end do
              end do
              do i=1,3
               do j=1,3
                 virep3bt(i,j)=virep3bt(i,j)+vir1moli2(i,j)
               end do
              end do

              npole3b=3
              pnum(1)=imol(1,moli3)
              pnum(2)=imol(1,moli3)+1
              pnum(3)=imol(2,moli3)
              ep3bt = ep3bt + ep1moli3
              do l1 = 1, npole3b
                 i = pnum(l1)
                 do j = 1, 3
                  dep3bt(j,i)=dep3bt(j,i)+dep1moli3(j,l1)
                 end do
              end do
              do i=1,3
                do j=1,3
                 virep3bt(i,j)=virep3bt(i,j)+vir1moli3(i,j)
                end do
              end do
            end if

            npole3b=3
            pnum(1)=imol(1,moli1)
            pnum(2)=imol(1,moli1)+1
            pnum(3)=imol(2,moli1)
            ep3bt=ep3bt-ep1moli1
            do l1 = 1, npole3b
              i = pnum(l1)
              do j = 1, 3
               dep3bt(j,i)=dep3bt(j,i)-dep1moli1(j,l1)
              end do
            end do
            do i=1,3
              do j=1,3
               virep3bt(i,j)=virep3bt(i,j)-vir1moli1(i,j)
              end do
            end do

            npole3b=3
            pnum(1)=imol(1,moli2)
            pnum(2)=imol(1,moli2)+1
            pnum(3)=imol(2,moli2)
            ep3bt=ep3bt-ep1moli2
            do l1 = 1, npole3b
              i = pnum(l1)
              do j = 1, 3
               dep3bt(j,i)=dep3bt(j,i)-dep1moli2(j,l1)
              end do
            end do
            do i=1,3
              do j=1,3
               virep3bt(i,j)=virep3bt(i,j)-vir1moli2(i,j)
              end do
            end do

            npole3b=3
            pnum(1)=imol(1,moli3)
            pnum(2)=imol(1,moli3)+1
            pnum(3)=imol(2,moli3)
            ep3bt=ep3bt-ep1moli3
            do l1 = 1, npole3b
              i = pnum(l1)
              do j = 1, 3
               dep3bt(j,i)=dep3bt(j,i)-dep1moli3(j,l1)
              end do
            end do
            do i=1,3
              do j=1,3
               virep3bt(i,j)=virep3bt(i,j)-vir1moli3(i,j)
              end do
            end do

           end if

          end do
        end do
      end do
!$OMP END DO
!$OMP END PARALLEL


cccccc!$OMP END PARALLEL DO

ccc 100  continue
c      ep3b=ep3b+ep3bt
      ep = ep+ep3bt
      write(6,*)ep
      print*,"Triplecount=",triplecount
      do i = 1, npole
         do j = 1, 3
            dep(j,i) =dep(j,i)+dep3bt(j,i)
         end do
      end do

      do i=1,3
         do j=1,3
           vir(i,j)=vir(i,j)+virep3bt(i,j)
         end do
      end do 
c      deallocate (field)
c      deallocate (fieldp)
      return
      end
c
c     ###############################################################
c     ##                                                           ##
c     ##               Subroutine empole1c_3b                      ##
c     ##  Liam O-Suilleabhain, Omar Demerdash, Teresa Head-Gordon  ##
c     ##                 Spring 2013                               ##
c     ###############################################################
c
c
c     "empole1c_3b" calculates the atomic multipole and dipole
c     polarizability interaction energy using the 3-body approximation
c
c
      subroutine empole1c_3b2
      implicit none
      include 'sizes.i'
      include 'energi.i'
      include 'atoms.i'
      include 'molcul.i'
      include 'polar.i'
      include 'deriv.i'
      include 'combo.i'
      include 'mpole.i'
      real*8, allocatable :: field(:,:)
      real*8, allocatable :: fieldp(:,:)
      real*8  ep3b,dep3b(3,maxatm)
      real*8  ep1(nmol),ep2(nmol,nmol)
      real*8  dep1(3,200,nmol),dep2(3,200,nmol,nmol)
      integer i,ii,j,l1,k,l
      
      ep3b = 0.0d0

      allocate (field(3,npole))
      allocate (fieldp(3,npole))

      do i = 1, npole
        do j = 1, 3
          field(j,i) = 0.0d0
          fieldp(j,i) = 0.0d0
          dep3b(j,i) = 0.0d0
          do k = 1, nmol
           dep1(j,i,k) = 0.0d0
           do l = 1, nmol
             dep2(j,i,k,l) = 0.0d0
           end do
         end do
        end do
      end do
c
c    Usual 2-body permanent electrostatics but without polarization
c
      call empole1c_3b_Perm(field,fieldp)
c
c    1-body polarization
c
      do moli1 = 1, nmol
        npole3b=3
        pnum(1)=imol(1,moli1)
        pnum(2)=imol(1,moli1)+1
        pnum(3)=imol(2,moli1)
        do i = 1, npole
          do j = 1, 3
            field(j,i) = 0.0d0
            fieldp(j,i) = 0.0d0
          end do
        end do
        call field_3b(field,fieldp)
        call empole1c_3b_Polar(field,fieldp)
        ep1(moli1) = ep
        ep3b = ep3b + ep1(moli1)
        do l1 = 1, npole3b
          i = pnum(l1)
          do j = 1, 3
            dep1(j,i,moli1) = dep1(j,i,moli1)+dep(j,i)
            dep3b(j,i) = dep3b(j,i)+dep1(j,i,moli1)
          end do
        end do
      end do
      write(6,*)ep3b
c
c    Subtract of 2-body polarization
c
      do moli1 = 1, nmol-1
        do moli2 = moli1+1, nmol
          npole3b=6
          pnum(1)=imol(1,moli1)
          pnum(2)=imol(1,moli1)+1
          pnum(3)=imol(2,moli1)
          pnum(4)=imol(1,moli2)
          pnum(5)=imol(1,moli2)+1
          pnum(6)=imol(2,moli2)
          do i = 1, npole
            do j = 1, 3
              field(j,i) = 0.0d0
              fieldp(j,i) = 0.0d0
            end do
          end do
          call field_3b(field,fieldp)
          call empole1c_3b_Polar(field,fieldp)
          ep2(moli1,moli2) = ep-ep1(moli1)-ep1(moli2)
          ep3b = ep3b + ep2(moli1,moli2)
          do l1 = 1, npole3b
            i = pnum(l1)
            do j = 1, 3
              dep2(j,i,moli1,moli2) = dep(j,i)-
     &                                dep1(j,i,moli1)-dep1(j,i,moli2)
              dep3b(j,i) = dep3b(j,i)+dep2(j,i,moli1,moli2)
            end do
          end do
        end do
      end do
      write(6,*)ep3b
c
c    Begin 3-body approximation
c
      do moli1 = 1, nmol-2
        do moli2 = moli1+1, nmol-1
          do moli3 = moli2+1, nmol
            npole3b=9
            pnum(1)=imol(1,moli1)
            pnum(2)=imol(1,moli1)+1
            pnum(3)=imol(2,moli1)
            pnum(4)=imol(1,moli2)
            pnum(5)=imol(1,moli2)+1
            pnum(6)=imol(2,moli2)
            pnum(7)=imol(1,moli3)
            pnum(8)=imol(1,moli3)+1
            pnum(9)=imol(2,moli3)
            do i = 1, npole
              do j = 1, 3
                field(j,i) = 0.0d0
                fieldp(j,i) = 0.0d0
              end do
            end do
c
c     compute the fixed electrostatic field for 3-body only
c
            call field_3b(field,fieldp)
c
c     Now do 3-body polarization for this trimer
c
            call empole1c_3b_Polar(field,fieldp)
            ep3b = ep3b + ep - ep2(moli1,moli2)-ep1(moli1)
     &                       - ep2(moli1,moli3)-ep1(moli2)
     &                       - ep2(moli2,moli3)-ep1(moli3)
            do l1 = 1, npole3b
              i = pnum(l1)
              do j = 1, 3
                dep3b(j,i) = dep3b(j,i)+dep(j,i)-dep2(j,i,moli1,moli2)-
     &                      dep2(j,i,moli1,moli3)-dep2(j,i,moli2,moli3)-
     &                   dep1(j,i,moli1)-dep1(j,i,moli2)-dep1(j,i,moli3)
              end do
            end do
          end do
        end do
      end do
      write(6,*)ep3b

      ep = ep3b
      do i = 1, npole
         do j = 1, 3
            dep(j,i) = dep3b(j,i)
         end do
      end do

      deallocate (field)
      deallocate (fieldp)
      return
      end
c
c     ###############################################################
c     ##                                                           ##
c     ##               Subroutine empole1c_3b                      ##
c     ##  Liam O-Suilleabhain, Omar Demerdash, Teresa Head-Gordon  ##
c     ##                 Spring 2013                               ##
c     ###############################################################
c
c
c     "empole1c_3b" calculates the atomic multipole and dipole
c     polarizability interaction energy using the 3-body approximation
c
c
      subroutine empole1c_3b3
      implicit none
      include 'sizes.i'
      include 'energi.i'
      include 'atoms.i'
      include 'molcul.i'
      include 'polar.i'
      include 'deriv.i'
      include 'combo.i'
      include 'mpole.i'
      real*8, allocatable :: field(:,:)
      real*8, allocatable :: fieldp(:,:)
      real*8  ep3b,dep3b(3,maxatm)
      integer i,ii,j,l1
      
      ep3b = 0.0d0

      allocate (field(3,npole))
      allocate (fieldp(3,npole))

      do i = 1, npole
        do j = 1, 3
          field(j,i) = 0.0d0
          fieldp(j,i) = 0.0d0
          dep3b(j,i) = 0.0d0
        end do
      end do
c
c    Usual 2-body permanent electrostatics but without polarization
c
      call empole1c_3b_Perm(field,fieldp)
c
c    1-body polarization
c
      do moli1 = 1, nmol
        npole3b=3
        pnum(1)=imol(1,moli1)
        pnum(2)=imol(1,moli1)+1
        pnum(3)=imol(2,moli1)
        call field_3b(field,fieldp)
        call empole1c_3b_Polar(field,fieldp)
        ep3b = ep3b + ep
        do l1 = 1, npole3b
          i = pnum(l1)
          do j = 1, 3
            dep3b(j,i) = dep3b(j,i)+dep(j,i)
          end do
        end do
      end do
      write(6,*)ep3b
c
c    2-body polarization
c
      do moli1 = 1, nmol-1
        do moli2 = moli1+1, nmol
          npole3b=3
          pnum(1)=imol(1,moli1)
          pnum(2)=imol(1,moli1)+1
          pnum(3)=imol(2,moli1)
          call field_3b(field,fieldp)
          call empole1c_3b_Polar(field,fieldp)
          ep3b = ep3b - ep
          do l1 = 1, npole3b
            i = pnum(l1)
            do j = 1, 3
              dep3b(j,i) = dep3b(j,i)-dep(j,i)
            end do
          end do
          npole3b=6
          pnum(4)=imol(1,moli2)
          pnum(5)=imol(1,moli2)+1
          pnum(6)=imol(2,moli2)
          call field_3b(field,fieldp)
          call empole1c_3b_Polar(field,fieldp)
          ep3b = ep3b + ep
          do l1 = 1, npole3b
            i = pnum(l1)
            do j = 1, 3
              dep3b(j,i) = dep3b(j,i)+dep(j,i)
            end do
          end do
          npole3b=3
          pnum(1)=imol(1,moli2)
          pnum(2)=imol(1,moli2)+1
          pnum(3)=imol(2,moli2)
          call field_3b(field,fieldp)
          call empole1c_3b_Polar(field,fieldp)
          ep3b = ep3b - ep
          do l1 = 1, npole3b
            i = pnum(l1)
            do j = 1, 3
              dep3b(j,i) = dep3b(j,i)-dep(j,i)
            end do
          end do
        end do
      end do
      write(6,*)ep3b
c
c    Begin 3-body approximation
c
      do moli1 = 1, nmol-2
        do moli2 = moli1+1, nmol-1
          do moli3 = moli2+1, nmol
            npole3b=9
            pnum(1)=imol(1,moli1)
            pnum(2)=imol(1,moli1)+1
            pnum(3)=imol(2,moli1)
            pnum(4)=imol(1,moli2)
            pnum(5)=imol(1,moli2)+1
            pnum(6)=imol(2,moli2)
            pnum(7)=imol(1,moli3)
            pnum(8)=imol(1,moli3)+1
            pnum(9)=imol(2,moli3)
            call field_3b(field,fieldp)
            call empole1c_3b_Polar(field,fieldp)
            ep3b = ep3b + ep
            do l1 = 1, npole3b
              i = pnum(l1)
              do j = 1, 3
                dep3b(j,i) = dep3b(j,i)+dep(j,i)
              end do
            end do
            npole3b=6
            pnum(1)=imol(1,moli1)
            pnum(2)=imol(1,moli1)+1
            pnum(3)=imol(2,moli1)
            pnum(4)=imol(1,moli2)
            pnum(5)=imol(1,moli2)+1
            pnum(6)=imol(2,moli2)
            call field_3b(field,fieldp)
            call empole1c_3b_Polar(field,fieldp)
            ep3b = ep3b - ep
            do l1 = 1, npole3b
              i = pnum(l1)
              do j = 1, 3
                dep3b(j,i) = dep3b(j,i)-dep(j,i)
              end do
            end do
            npole3b=6
            pnum(1)=imol(1,moli1)
            pnum(2)=imol(1,moli1)+1
            pnum(3)=imol(2,moli1)
            pnum(4)=imol(1,moli3)
            pnum(5)=imol(1,moli3)+1
            pnum(6)=imol(2,moli3)
            call field_3b(field,fieldp)
            call empole1c_3b_Polar(field,fieldp)
            ep3b = ep3b - ep
            do l1 = 1, npole3b
              i = pnum(l1)
              do j = 1, 3
                dep3b(j,i) = dep3b(j,i)-dep(j,i)
              end do
            end do
            npole3b=6
            pnum(1)=imol(1,moli2)
            pnum(2)=imol(1,moli2)+1
            pnum(3)=imol(2,moli2)
            pnum(4)=imol(1,moli3)
            pnum(5)=imol(1,moli3)+1
            pnum(6)=imol(2,moli3)
            call field_3b(field,fieldp)
            call empole1c_3b_Polar(field,fieldp)
            ep3b = ep3b - ep
            do l1 = 1, npole3b
              i = pnum(l1)
              do j = 1, 3
                dep3b(j,i) = dep3b(j,i)-dep(j,i)
              end do
            end do
            npole3b=3
            pnum(1)=imol(1,moli1)
            pnum(2)=imol(1,moli1)+1
            pnum(3)=imol(2,moli1)
            call field_3b(field,fieldp)
            call empole1c_3b_Polar(field,fieldp)
            ep3b = ep3b + ep
            do l1 = 1, npole3b
              i = pnum(l1)
              do j = 1, 3
                dep3b(j,i) = dep3b(j,i)+dep(j,i)
              end do
            end do
            npole3b=3
            pnum(1)=imol(1,moli2)
            pnum(2)=imol(1,moli2)+1
            pnum(3)=imol(2,moli2)
            call field_3b(field,fieldp)
            call empole1c_3b_Polar(field,fieldp)
            ep3b = ep3b + ep
            do l1 = 1, npole3b
              i = pnum(l1)
              do j = 1, 3
                dep3b(j,i) = dep3b(j,i)+dep(j,i)
              end do
            end do
            npole3b=3
            pnum(1)=imol(1,moli3)
            pnum(2)=imol(1,moli3)+1
            pnum(3)=imol(2,moli3)
            call field_3b(field,fieldp)
            call empole1c_3b_Polar(field,fieldp)
            ep3b = ep3b + ep
            do l1 = 1, npole3b
              i = pnum(l1)
              do j = 1, 3
                dep3b(j,i) = dep3b(j,i)+dep(j,i)
              end do
            end do
          end do
        end do
      end do
      write(6,*)ep3b

      ep = ep3b
      do i = 1, npole
         do j = 1, 3
            dep(j,i) = dep3b(j,i)
         end do
      end do

      deallocate (field)
      deallocate (fieldp)
      return
      end
c
c     ################################################################
c     ##                                                            ##
c     ##              Subroutine empole1c_3b_Perm                   ##
c     ##  Liam O-Suilleabhain, Omar Demerdash, Teresa Head-Gordon   ##
c     ##                 Spring 2013                                ##
c     ##                                                            ##
c     ################################################################
c
c
c     "empole1c_3b_Perm" calculates the multipole energy and derivatives 
c     with respect to Cartesian coordinates using regular Ewald 
c
c
c      subroutine empole1c_3b_Perm(field,fieldp)
      subroutine empole1c_3b_Perm
      implicit none
      include 'sizes.i'
      include 'atoms.i'
      include 'boxes.i'
      include 'combo.i'
      include 'chgpot.i'
      include 'deriv.i'
      include 'energi.i'
      include 'ewald.i'
      include 'inter.i'
      include 'math.i'
      include 'mpole.i'
      include 'polar.i'
      include 'potent.i'
      include 'polpot.i'
      include 'virial.i'
      integer i,j,ii
      real*8 e,ei,eintra
      real*8 f,term,fterm
      real*8 cii,dii,qii,uii
      real*8 xd,yd,zd
      real*8 xu,yu,zu
      real*8 xup,yup,zup
      real*8 xq,yq,zq
      real*8 xv,yv,zv,vterm
      real*8 ci,dix,diy,diz
      real*8 uix,uiy,uiz
      real*8 qixx,qixy,qixz
      real*8 qiyy,qiyz,qizz
      real*8 xdfield,xufield
      real*8 ydfield,yufield
      real*8 zdfield,zufield
      real*8 trq(3),trqi(3)
      real*8 frcx(3),frcy(3),frcz(3)
      real*8 aewald

      aewald=aewaldPerm
c      real*8 field(3,npole)
c      real*8 fieldp(3,npole)
c
c     zero out multipole energy and derivatives
c
      em = 0.0d0
      do i = 1, n
        do j = 1, 3
          dem(j,i) = 0.0d0
        end do
      end do
      if (.not.use_mpole)return
c
c     set the energy unit conversion factor
c
      f = electric / dielec
c
c     check the sign of multipole components at chiral sites
c
      call chkpole
c
c     rotate the multipole components into the global frame
c
      call rotpole
c
c     get the field at each atom
c

c      call field_3b_Perm(field,fieldp)

c
c     compute the reciprocal space part of the Ewald summation
c
      call erecip1_3b_Perm
c
c     compute the real space part of the Ewald summation
c
      call ereal1c_3b_Perm(eintra)
c
c     compute the Ewald self-energy term over all the atoms
c
      term = 2.0d0 * aewald * aewald
      fterm = -f * aewald / sqrtpi
      do i = 1, npole
         ci = rpole(1,i)
         dix = rpole(2,i)
         diy = rpole(3,i)
         diz = rpole(4,i)
         qixx = rpole(5,i)
         qixy = rpole(6,i)
         qixz = rpole(7,i)
         qiyy = rpole(9,i)
         qiyz = rpole(10,i)
         qizz = rpole(13,i)
         cii = ci*ci
         dii = dix*dix + diy*diy + diz*diz
         qii = qixx*qixx + qiyy*qiyy + qizz*qizz
     &            + 2.0d0*(qixy*qixy+qixz*qixz+qiyz*qiyz)
         e = fterm * (cii + term*(dii/3.0d0+2.0d0*term*qii/5.0d0))
         em = em + e
      end do
c
c     compute the cell dipole boundary correction term
c

      if (boundary .eq. 'VACUUM') then
        trqi(1) = 0.0d0
        trqi(2) = 0.0d0
        trqi(3) = 0.0d0
        xd = 0.0d0
        yd = 0.0d0
        zd = 0.0d0
        xu = 0.0d0
        yu = 0.0d0
        zu = 0.0d0
        xup = 0.0d0
        yup = 0.0d0
        zup = 0.0d0
        do i = 1, npole
          ii = ipole(i)
          xd = xd + rpole(2,i) + rpole(1,i)*x(ii)
          yd = yd + rpole(3,i) + rpole(1,i)*y(ii)
          zd = zd + rpole(4,i) + rpole(1,i)*z(ii)
        end do
        term = (2.0d0/3.0d0) * f * (pi/volbox)
        em = em + term*(xd*xd+yd*yd+zd*zd)
        do i = 1, npole
          ii = ipole(i)
          dem(1,ii) = dem(1,ii) + 2.0d0*term*rpole(1,i)*xd
          dem(2,ii) = dem(2,ii) + 2.0d0*term*rpole(1,i)*yd
          dem(3,ii) = dem(3,ii) + 2.0d0*term*rpole(1,i)*zd
        end do
        xdfield = -2.0d0 * term * xd
        ydfield = -2.0d0 * term * yd
        zdfield = -2.0d0 * term * zd
        do i = 1, npole
          trq(1) = rpole(3,i)*zdfield - rpole(4,i)*ydfield
          trq(2) = rpole(4,i)*xdfield - rpole(2,i)*zdfield
          trq(3) = rpole(2,i)*ydfield - rpole(3,i)*xdfield
          call torque_3b_Perm (i,trq,trqi,frcx,frcy,frcz)
        end do
c
c     boundary correction to virial due to overall cell dipole
c
        xd = 0.0d0
        yd = 0.0d0
        zd = 0.0d0
        xq = 0.0d0
        yq = 0.0d0
        zq = 0.0d0
        do i = 1, npole
          ii = ipole(i)
          xd = xd + rpole(2,i)
          yd = yd + rpole(3,i)
          zd = zd + rpole(4,i)
          xq = xq + rpole(1,i)*x(ii)
          yq = yq + rpole(1,i)*y(ii)
          zq = zq + rpole(1,i)*z(ii)
        end do
        xv = xq * (xd+0.5d0*(xu+xup))
        yv = yq * (yd+0.5d0*(yu+yup))
        zv = zq * (zd+0.5d0*(zu+zup))
        vterm = term * (xq*xq + yq*yq + zq*zq + 2.0d0*(xv+yv+zv)
     &                      + xu*xup + yu*yup + zu*zup
     &                      + xd*(xd+xu+xup) + yd*(yd+yu+yup)
     &                      + zd*(zd+zu+zup))
        vir(1,1) = vir(1,1) + 2.0d0*term*(xq*xq+xv) + vterm
        vir(2,1) = vir(2,1) + 2.0d0*term*(xq*yq+xv)
        vir(3,1) = vir(3,1) + 2.0d0*term*(xq*zq+xv)
        vir(1,2) = vir(1,2) + 2.0d0*term*(yq*xq+yv)
        vir(2,2) = vir(2,2) + 2.0d0*term*(yq*yq+yv) + vterm
        vir(3,2) = vir(3,2) + 2.0d0*term*(yq*zq+yv)
        vir(1,3) = vir(1,3) + 2.0d0*term*(zq*xq+zv)
        vir(2,3) = vir(2,3) + 2.0d0*term*(zq*yq+zv)
        vir(3,3) = vir(3,3) + 2.0d0*term*(zq*zq+zv) + vterm
        if (poltyp .eq. 'DIRECT') then
          vterm = term * (xu*xup+yu*yup+zu*zup)
          vir(1,1) = vir(1,1) + vterm
          vir(2,2) = vir(2,2) + vterm
          vir(3,3) = vir(3,3) + vterm
        end if
      end if
c
c     intermolecular energy is total minus intramolecular part
c
      einter = einter + em - eintra
c      do i = 1, n
c          print*,'dem after Perm x', i, dem(1,i)
c          print*,'dem after Perm y', i, dem(2,i)
c          print*,'dem after Perm z', i, dem(3,i) 
c      end do

      return
      end
c
c     ###############################################################
c     ##                                                           ##
c     ##              Subroutine field_3b_Perm                     ##
c     ##  Liam O-Suilleabhain, Omar Demerdash, Teresa Head-Gordon  ##
c     ##                 Spring 2013                               ##
c     ##                                                           ##
c     ###############################################################
c
c
c     "field_3b_Perm" computes the field due to multipoles
c
c
      subroutine field_3b_Perm(field,fieldp)
      implicit none
      include 'sizes.i'
      include 'atoms.i'
      include 'boxes.i'
      include 'combo.i'
      include 'ewald.i'
      include 'inform.i'
      include 'iounit.i'
      include 'math.i'
      include 'mpole.i'
      include 'polar.i'
      include 'polpot.i'
      include 'potent.i'
      include 'units.i'
      include 'uprior.i'
      integer i,j,k,ii
      integer iter,maxiter
      real*8 eps,term
      real*8 epsd,epsp
      real*8 epsold
      real*8 udsum,upsum
      real*8 ucell(3)
      real*8 ucellp(3)
      real*8 field(3,*)
      real*8 fieldp(3,*),aewald
      logical done
c
c     get the reciprical space part of the electrostatic field
c
      call udirect1_Perm(field)
      do i = 1, npole
         do j = 1, 3
            fieldp(j,i) = field(j,i)
         end do
      end do
c
c     get the real space portion of the electrostatic field
c
      call udirect2a(field,fieldp)
c
c     get the self-energy portion of the electrostatic field
c
      aewald=aewaldPerm
      term = (4.0d0/3.0d0) * aewald**3 / sqrtpi
      do i = 1, npole
         do j = 1, 3
            field(j,i) = field(j,i) + term*rpole(j+1,i)
            fieldp(j,i) = fieldp(j,i) + term*rpole(j+1,i)
         end do
      end do
c
c     compute the cell dipole boundary correction to field
c
      if (boundary .eq. 'VACUUM') then
         do i = 1, 3
            ucell(i) = 0.0d0
         end do
         do i = 1, npole
            ii = ipole(i)
            ucell(1) = ucell(1) + rpole(2,i) + rpole(1,i)*x(ii)
            ucell(2) = ucell(2) + rpole(3,i) + rpole(1,i)*y(ii)
            ucell(3) = ucell(3) + rpole(4,i) + rpole(1,i)*z(ii)
         end do
         term = (4.0d0/3.0d0) * pi/volbox
         do i = 1, npole
            do j = 1, 3
               field(j,i) = field(j,i) - term*ucell(j)
               fieldp(j,i) = fieldp(j,i) - term*ucell(j)
            end do
         end do
      end if
      return
      end
c
c     #################################################################
c     ##                                                             ##
c     ##              Subroutine udirect1_Perm                       ##
c     ##  Liam O-Suilleabhain, Omar Demerdash, Teresa Head-Gordon    ##
c     ##                 Spring 2013                                 ##
c     ##                                                             ##
c     #################################################################
c
c
c     "udirect1_Perm" computes the reciprocal space contribution of
c     permanent atomic multipole moments to the electrostatic field
c     via a regular Ewald summation
c
c
      subroutine udirect1_Perm(field)
      implicit none
      include 'sizes.i'
      include 'atoms.i'
      include 'boxes.i'
      include 'ewald.i'
      include 'ewreg.i'
      include 'math.i'
      include 'mpole.i'
      include 'units.i'
      integer i,j,k,l,ii
      integer jmin,jmax
      integer kmin,kmax
      integer lmin,lmax
      real*8 expterm,cut
      real*8 term,fterm
      real*8 xfr,yfr,zfr
      real*8 rj,rk,rl
      real*8 h1,h2,h3,hsq
      real*8 qf,t1,t2
      real*8 ck,dk,qk
      real*8 q1,q2,q3
      real*8 ckr(maxatm)
      real*8 skr(maxatm)
      real*8 cjk(maxatm)
      real*8 sjk(maxatm)
      real*8 cm(maxatm)
      real*8 dm(3,maxatm)
      real*8 qm(9,maxatm)
      real*8 field(3,maxatm)
      real*8 frecip,aewald
      aewald=aewaldPerm
c
c     return if the Ewald coefficient is zero
c
      if (aewald .lt. 1.0d-6)  return
      term = -0.25d0 / aewald**2
      fterm = 8.0d0 * pi / volbox
c
c     set the number of vectors based on box dimensions
c
      frecip = 0.5d0
      cut = 4.0d0 * pi * pi * frecip
      jmin = 0
      kmin = 0
      lmin = 1
      jmax = min(maxvec,int(frecip/recip(1,1)))
      kmax = min(maxvec,int(frecip/recip(2,2)))
      lmax = min(maxvec,int(frecip/recip(3,3)))
c
c     copy the multipole moments into local storage areas
c
      do i = 1, npole
         cm(i) = rpole(1,i)
         dm(1,i) = rpole(2,i)
         dm(2,i) = rpole(3,i)
         dm(3,i) = rpole(4,i)
         qm(1,i) = rpole(5,i)
         qm(2,i) = rpole(6,i)
         qm(3,i) = rpole(7,i)
         qm(4,i) = rpole(8,i)
         qm(5,i) = rpole(9,i)
         qm(6,i) = rpole(10,i)
         qm(7,i) = rpole(11,i)
         qm(8,i) = rpole(12,i)
         qm(9,i) = rpole(13,i)
      end do
c
c     calculate and store the exponential factors
c
      do i = 1, npole
         ii = ipole(i)
         zfr = (z(ii)/gamma_term) / zbox
         yfr = ((y(ii)-zfr*zbox*beta_term)/gamma_sin) / ybox
         xfr = (x(ii)-yfr*ybox*gamma_cos-zfr*zbox*beta_cos) / xbox
         xfr = 2.0d0 * pi * xfr
         yfr = 2.0d0 * pi * yfr
         zfr = 2.0d0 * pi * zfr
         ejc(i,0) = 1.0d0
         ejs(i,0) = 0.0d0
         ekc(i,0) = 1.0d0
         eks(i,0) = 0.0d0
         elc(i,0) = 1.0d0
         els(i,0) = 0.0d0
         ejc(i,1) = cos(xfr)
         ejs(i,1) = sin(xfr)
         ekc(i,1) = cos(yfr)
         eks(i,1) = sin(yfr)
         elc(i,1) = cos(zfr)
         els(i,1) = sin(zfr)
         ekc(i,-1) = ekc(i,1)
         eks(i,-1) = -eks(i,1)
         elc(i,-1) = elc(i,1)
         els(i,-1) = -els(i,1)
         do j = 2, jmax
            ejc(i,j) = ejc(i,j-1)*ejc(i,1) - ejs(i,j-1)*ejs(i,1)
            ejs(i,j) = ejs(i,j-1)*ejc(i,1) + ejc(i,j-1)*ejs(i,1)
         end do
         do j = 2, kmax
            ekc(i,j) = ekc(i,j-1)*ekc(i,1) - eks(i,j-1)*eks(i,1)
            eks(i,j) = eks(i,j-1)*ekc(i,1) + ekc(i,j-1)*eks(i,1)
            ekc(i,-j) = ekc(i,j)
            eks(i,-j) = -eks(i,j)
         end do
         do j = 2, lmax
            elc(i,j) = elc(i,j-1)*elc(i,1) - els(i,j-1)*els(i,1)
            els(i,j) = els(i,j-1)*elc(i,1) + elc(i,j-1)*els(i,1)
            elc(i,-j) = elc(i,j)
            els(i,-j) = -els(i,j)
         end do
      end do
c
c     loop over all k vectors from the reciprocal lattice
c
      do j = jmin, jmax
         rj = 2.0d0 * pi * dble(j)
         do k = kmin, kmax
            rk = 2.0d0 * pi * dble(k)
            do i = 1, npole
               cjk(i) = ejc(i,j)*ekc(i,k) - ejs(i,j)*eks(i,k)
               sjk(i) = ejs(i,j)*ekc(i,k) + ejc(i,j)*eks(i,k)
            end do
            do l = lmin, lmax
               rl = 2.0d0 * pi * dble(l)
               h1 = recip(1,1)*rj
               h2 = recip(2,1)*rj + recip(2,2)*rk
               h3 = recip(3,1)*rj + recip(3,2)*rk + recip(3,3)*rl
               hsq = h1*h1 + h2*h2 + h3*h3
               if (hsq .le. cut) then
                  t1 = 0.0d0
                  t2 = 0.0d0
                  do i = 1, npole
                     ckr(i) = cjk(i)*elc(i,l) - sjk(i)*els(i,l)
                     skr(i) = sjk(i)*elc(i,l) + cjk(i)*els(i,l)
                     ck = cm(i)
                     dk = h1*dm(1,i) + h2*dm(2,i) + h3*dm(3,i)
                     q1 = h1*qm(1,i) + h2*qm(4,i) + h3*qm(7,i)
                     q2 = h1*qm(2,i) + h2*qm(5,i) + h3*qm(8,i)
                     q3 = h1*qm(3,i) + h2*qm(6,i) + h3*qm(9,i)
                     qk = h1*q1 + h2*q2 + h3*q3
                     t1 = t1 + (ck-qk)*skr(i) + dk*ckr(i)
                     t2 = t2 + (ck-qk)*ckr(i) - dk*skr(i)
                  end do
                  expterm = fterm * exp(term*hsq) / hsq
                  if (octahedron) then
                     if (mod(j+k+l,2) .ne. 0)  expterm = 0.0d0
                  end if
                  do i = 1, npole
                     qf = expterm * (skr(i)*t2-ckr(i)*t1)
                     field(1,i) = field(1,i) + h1*qf
                     field(2,i) = field(2,i) + h2*qf
                     field(3,i) = field(3,i) + h3*qf
                  end do
               end if
            end do
            lmin = -lmax
         end do
         kmin = -kmax
      end do
      return
      end
c
c     #################################################################
c     ##                                                             ##
c     ##                 Subroutine erecip1_3b_Perm                  ##
c     ##  Liam O-Suilleabhain, Omar Demerdash, Teresa Head-Gordon    ##
c     ##                 Spring 2013                                 ##
c     ##                                                             ##
c     #################################################################
c
c
c     "erecip1_3b_Perm" evaluates the reciprocal space portion of the regular
c     Ewald summation energy and gradient due to atomic multipole
c     interactions only
c
c     literature reference:
c
c     W. Smith, "Point Multipoles in the Ewald Summation (Revisited)",
c     CCP5 Newsletter, 46, 18-30, 1998  (see http://www.ccp5.org/)
c
c
      subroutine erecip1_3b_Perm
      implicit none
      include 'sizes.i'
      include 'atoms.i'
      include 'boxes.i'
      include 'chgpot.i'
      include 'deriv.i'
      include 'energi.i'
      include 'ewald.i'
      include 'ewreg.i'
      include 'math.i'
      include 'mpole.i'
      include 'polar.i'
      include 'polpot.i'
      include 'units.i'
      include 'virial.i'
      integer i,j,k,l
      integer ii,m1,m2
      integer jmin,jmax
      integer kmin,kmax
      integer lmin,lmax
      real*8 e,ei,etot,f,cut
      real*8 expterm,term,eterm
      real*8 uterm,vterm
      real*8 xfr,yfr,zfr
      real*8 rj,rk,rl
      real*8 h1,h2,h3,hsq
      real*8 ck,dk,qk,uk
      real*8 q1,q2,q3
      real*8 t1,t2,t3,t4
      real*8 ukp,t3p,t4p
      real*8 de,det1,det2
      real*8 dei,det1i,det2i
      real*8 wterm(3,3)
      real*8 t5(3,3),t6(3,3)
      real*8 t5u(3,3),t5p(3,3)
      real*8 t6u(3,3),t6p(3,3)
      real*8 qt(3,3),dt(3,3)
      real*8 dtu(3,3),dtp(3,3)
      real*8 ckr(maxatm),skr(maxatm)
      real*8 cjk(maxatm),sjk(maxatm)
      real*8 s1(maxatm),s2(maxatm)
      real*8 s3(maxatm),s4(maxatm)
      real*8 s3p(maxatm),s4p(maxatm)
      real*8 cm(maxatm),dm(3,maxatm)
      real*8 qm(9,maxatm),um(3,maxatm)
      real*8 trq(3,maxatm),trqi(3,maxatm)
      real*8 dkx(maxatm),qkx(maxatm)
      real*8 dky(maxatm),qky(maxatm)
      real*8 dkz(maxatm),qkz(maxatm)
      real*8 frecip,aewald
      frecip = 0.5
c
c     return if the Ewald coefficient is zero
c
      aewald=aewaldPerm

      if (aewald .lt. 1.0d-6)  return
      f = electric / dielec
      term = -0.25d0 / aewald**2
      eterm = 4.0d0 * pi * f / volbox
c
c     set the number of vectors based on box dimensions
c
      cut = 4.0d0 * pi * pi * frecip
      jmin = 0
      kmin = 0
      lmin = 1
c      jmax = min(maxvec,int(frecip/recip(1,1)))
c      kmax = min(maxvec,int(frecip/recip(2,2)))
c      lmax = min(maxvec,int(frecip/recip(3,3)))
      jmax=min(maxvec,jmaxPerm)
      kmax=min(maxvec,kmaxPerm)
      lmax=min(maxvec,lmaxPerm)
c      print*,"jmax,kmax,lmax in erecip1_3b",jmax3b,kmax3b,lmax3b
      print*,"jmaxPerm,kmaxPerm,lmaxPerm",jmaxPerm,kmaxPerm,lmaxPerm

      print*,"jmax,kmax,lmax in erecip1_3b_Perm",jmax,kmax,lmax
c
c     copy the multipole moments into local storage areas
c
      do i = 1, npole
         cm(i) = rpole(1,i)
         dm(1,i) = rpole(2,i)
         dm(2,i) = rpole(3,i)
         dm(3,i) = rpole(4,i)
         qm(1,i) = rpole(5,i)
         qm(2,i) = rpole(6,i)
         qm(3,i) = rpole(7,i)
         qm(4,i) = rpole(8,i)
         qm(5,i) = rpole(9,i)
         qm(6,i) = rpole(10,i)
         qm(7,i) = rpole(11,i)
         qm(8,i) = rpole(12,i)
         qm(9,i) = rpole(13,i)
      end do
c
c     zero out local accumulation arrays for derivatives
c
      do i = 1, n
         trq(1,i) = 0.0d0
         trq(2,i) = 0.0d0
         trq(3,i) = 0.0d0
         trqi(1,i) = 0.0d0
         trqi(2,i) = 0.0d0
         trqi(3,i) = 0.0d0
      end do
c
c     calculate and store the exponential factors
c
      do i = 1, npole
         ii = ipole(i)
         zfr = (z(ii)/gamma_term) / zbox
         yfr = ((y(ii)-zfr*zbox*beta_term)/gamma_sin) / ybox
         xfr = (x(ii)-yfr*ybox*gamma_cos-zfr*zbox*beta_cos) / xbox
         xfr = 2.0d0 * pi * xfr
         yfr = 2.0d0 * pi * yfr
         zfr = 2.0d0 * pi * zfr
         ejc(i,0) = 1.0d0
         ejs(i,0) = 0.0d0
         ekc(i,0) = 1.0d0
         eks(i,0) = 0.0d0
         elc(i,0) = 1.0d0
         els(i,0) = 0.0d0
         ejc(i,1) = cos(xfr)
         ejs(i,1) = sin(xfr)
         ekc(i,1) = cos(yfr)
         eks(i,1) = sin(yfr)
         elc(i,1) = cos(zfr)
         els(i,1) = sin(zfr)
         ekc(i,-1) = ekc(i,1)
         eks(i,-1) = -eks(i,1)
         elc(i,-1) = elc(i,1)
         els(i,-1) = -els(i,1)
         do j = 2, jmax
            ejc(i,j) = ejc(i,j-1)*ejc(i,1) - ejs(i,j-1)*ejs(i,1)
            ejs(i,j) = ejs(i,j-1)*ejc(i,1) + ejc(i,j-1)*ejs(i,1)
         end do
         do j = 2, kmax
            ekc(i,j) = ekc(i,j-1)*ekc(i,1) - eks(i,j-1)*eks(i,1)
            eks(i,j) = eks(i,j-1)*ekc(i,1) + ekc(i,j-1)*eks(i,1)
            ekc(i,-j) = ekc(i,j)
            eks(i,-j) = -eks(i,j)
         end do
         do j = 2, lmax
            elc(i,j) = elc(i,j-1)*elc(i,1) - els(i,j-1)*els(i,1)
            els(i,j) = els(i,j-1)*elc(i,1) + elc(i,j-1)*els(i,1)
            elc(i,-j) = elc(i,j)
            els(i,-j) = -els(i,j)
         end do
      end do
c
c     loop over all k vectors from the reciprocal lattice
c
      do j = jmin, jmax
         rj = 2.0d0 * pi * dble(j)
         do k = kmin, kmax
            rk = 2.0d0 * pi * dble(k)
            do i = 1, npole
               cjk(i) = ejc(i,j)*ekc(i,k) - ejs(i,j)*eks(i,k)
               sjk(i) = ejs(i,j)*ekc(i,k) + ejc(i,j)*eks(i,k)
            end do
            do l = lmin, lmax
               rl = 2.0d0 * pi * dble(l)
               h1 = recip(1,1)*rj
               h2 = recip(2,1)*rj + recip(2,2)*rk
               h3 = recip(3,1)*rj + recip(3,2)*rk + recip(3,3)*rl
               hsq = h1*h1 + h2*h2 + h3*h3
c               if (hsq .le. cut) then
                  t1 = 0.0d0
                  t2 = 0.0d0
                  t3 = 0.0d0
                  t4 = 0.0d0
                  t3p = 0.0d0
                  t4p = 0.0d0
                  do m2 = 1, 3
                     do m1 = 1, 3
                        t5(m1,m2) = 0.0d0
                        t5u(m1,m2) = 0.0d0
                        t5p(m1,m2) = 0.0d0
                        t6(m1,m2) = 0.0d0
                        t6u(m1,m2) = 0.0d0
                        t6p(m1,m2) = 0.0d0
                     end do
                  end do
                  do i = 1, npole
                     ckr(i) = cjk(i)*elc(i,l) - sjk(i)*els(i,l)
                     skr(i) = sjk(i)*elc(i,l) + cjk(i)*els(i,l)
                     ck = cm(i)
                     dk = h1*dm(1,i) + h2*dm(2,i) + h3*dm(3,i)
                     dkx(i) = h3*dm(2,i) - h2*dm(3,i)
                     dky(i) = h1*dm(3,i) - h3*dm(1,i)
                     dkz(i) = h2*dm(1,i) - h1*dm(2,i)
                     q1 = h1*qm(1,i) + h2*qm(4,i) + h3*qm(7,i)
                     q2 = h1*qm(2,i) + h2*qm(5,i) + h3*qm(8,i)
                     q3 = h1*qm(3,i) + h2*qm(6,i) + h3*qm(9,i)
                     qk = h1*q1 + h2*q2 + h3*q3
                     qkx(i) = h3*q2 - h2*q3
                     qky(i) = h1*q3 - h3*q1
                     qkz(i) = h2*q1 - h1*q2
                     s1(i) = (ck-qk)*skr(i) + dk*ckr(i)
                     s2(i) = (ck-qk)*ckr(i) - dk*skr(i)
                     t1 = t1 + s1(i)
                     t2 = t2 + s2(i)
c                     t3 = t3 + s3(i)
c                     t4 = t4 + s4(i)
c
c     terms needed for subsequent virial tensor calculation
c
                     qt(1,1) = h1*(h1*qm(1,i) + h2*qm(4,i) + h3*qm(7,i))
                     qt(2,1) = h1*(h1*qm(2,i) + h2*qm(5,i) + h3*qm(8,i))
                     qt(3,1) = h1*(h1*qm(3,i) + h2*qm(6,i) + h3*qm(9,i))
                     qt(1,2) = h2*(h1*qm(1,i) + h2*qm(4,i) + h3*qm(7,i))
                     qt(2,2) = h2*(h1*qm(2,i) + h2*qm(5,i) + h3*qm(8,i))
                     qt(3,2) = h2*(h1*qm(3,i) + h2*qm(6,i) + h3*qm(9,i))
                     qt(1,3) = h3*(h1*qm(1,i) + h2*qm(4,i) + h3*qm(7,i))
                     qt(2,3) = h3*(h1*qm(2,i) + h2*qm(5,i) + h3*qm(8,i))
                     qt(3,3) = h3*(h1*qm(3,i) + h2*qm(6,i) + h3*qm(9,i))
                     dt(1,1) = h1 * dm(1,i)
                     dt(2,1) = h1 * dm(2,i)
                     dt(3,1) = h1 * dm(3,i)
                     dt(1,2) = h2 * dm(1,i)
                     dt(2,2) = h2 * dm(2,i)
                     dt(3,2) = h2 * dm(3,i)
                     dt(1,3) = h3 * dm(1,i)
                     dt(2,3) = h3 * dm(2,i)
                     dt(3,3) = h3 * dm(3,i)
                     do m2 = 1, 3
                        do m1 = 1, 3
                           t5(m1,m2) = t5(m1,m2) - dt(m1,m2)*ckr(i)
     &                                    + 2.0d0*qt(m1,m2)*skr(i)
                           t6(m1,m2) = t6(m1,m2) + dt(m1,m2)*skr(i)
     &                                    + 2.0d0*qt(m1,m2)*ckr(i)
                        end do
                     end do
                  end do
c
c     get the energy contributions for current reciprocal vector
c
                  expterm = eterm * exp(term*hsq) / hsq
                  if (octahedron) then
                     if (mod(j+k+l,2) .ne. 0)  expterm = 0.0d0
                  end if
                  e = expterm * (t1*t1+t2*t2)
                  etot = e 
                  em = em + e
c
c     get the virial contributions for current reciprocal vector
c
c                  uterm = expterm * (t1*(t1+t3) + t2*(t2+t4))
                  uterm = expterm * (t1*(t1) + t2*(t2))
                  do m2 = 1, 3
                     do m1 = 1, 3
                        wterm(m1,m2) = 2.0d0 * expterm
     &                     * (t1*t5(m1,m2) + t2*t6(m1,m2))
c     &                        + (t3)*t5(m1,m2)
c     &                        + (t4)*t6(m1,m2))
                     end do
                  end do
                  wterm(2,1) = 0.5d0 * (wterm(2,1)+wterm(1,2))
                  wterm(3,1) = 0.5d0 * (wterm(3,1)+wterm(1,3))
                  wterm(3,2) = 0.5d0 * (wterm(3,2)+wterm(2,3))
                  wterm(1,2) = wterm(2,1)
                  wterm(1,3) = wterm(3,1)
                  wterm(2,3) = wterm(3,2)
                  vterm = 2.0d0 * uterm * (1.0d0-term*hsq) / hsq
                  vir(1,1) = vir(1,1) + h1*h1*vterm + wterm(1,1) - uterm
                  vir(2,1) = vir(2,1) + h2*h1*vterm + wterm(2,1)
                  vir(3,1) = vir(3,1) + h3*h1*vterm + wterm(3,1)
                  vir(1,2) = vir(1,2) + h1*h2*vterm + wterm(1,2)
                  vir(2,2) = vir(2,2) + h2*h2*vterm + wterm(2,2) - uterm
                  vir(3,2) = vir(3,2) + h3*h2*vterm + wterm(3,2)
                  vir(1,3) = vir(1,3) + h1*h3*vterm + wterm(1,3)
                  vir(2,3) = vir(2,3) + h2*h3*vterm + wterm(2,3)
                  vir(3,3) = vir(3,3) + h3*h3*vterm + wterm(3,3) - uterm
c
c     get the force contributions for current reciprocal vector
c
                  expterm = 2.0d0 * expterm
                  do i = 1, npole
                     ii = ipole(i)
                     de = expterm * (s2(i)*t1-s1(i)*t2)
                     det1 = expterm * (skr(i)*t2-ckr(i)*t1)
                     det2 = 2.0d0 * expterm * (ckr(i)*t2+skr(i)*t1)
                     dem(1,ii) = dem(1,ii) + h1*de
                     dem(2,ii) = dem(2,ii) + h2*de
                     dem(3,ii) = dem(3,ii) + h3*de
                     trq(1,ii) = trq(1,ii) + dkx(i)*det1 + qkx(i)*det2
                     trq(2,ii) = trq(2,ii) + dky(i)*det1 + qky(i)*det2
                     trq(3,ii) = trq(3,ii) + dkz(i)*det1 + qkz(i)*det2
                  end do
c               end if
            end do
            lmin = -lmax
         end do
         kmin = -kmax
      end do
c
c     convert the torques to forces and increment the totals
c
      call torque_reg (trq,dem)   

      return
      end
c
c     ################################################################
c     ##                                                            ##
c     ##               Subroutine ereal1c_3b_Perm                   ##
c     ##  Liam O-Suilleabhain, Omar Demerdash, Teresa Head-Gordon   ##
c     ##                 Spring 2013                                ##
c     ##                                                            ##
c     ################################################################
c
c
c     "ereal1c_3b_Perm" evaluates the real space portion of the regular Ewald
c     summation energy and gradient due to atomic multipole interactions
c
c
      subroutine ereal1c_3b_Perm(eintra)
      implicit none
      include 'sizes.i'
      include 'atoms.i'
      include 'bound.i'
      include 'boxes.i'
      include 'cell.i'
      include 'chgpot.i'
      include 'couple.i'
      include 'cutoff.i'
      include 'deriv.i'
      include 'energi.i'
      include 'combo.i'
      include 'ewald.i'
      include 'math.i'
      include 'molcul.i'
      include 'mplpot.i'
      include 'mpole.i'
      include 'polar.i'
      include 'polgrp.i'
      include 'polpot.i'
      include 'shunt.i'
      include 'virial.i'
      integer i,j,k,l1,l3
      integer ii,kk,jcell
      integer iax,iay,iaz
      integer kax,kay,kaz
      real*8 e,ei,f,bfac
      real*8 eintra,erfc
      real*8 damp,expdamp
      real*8 pdi,pti,pgamma
      real*8 scale3,scale5
      real*8 scale7
      real*8 temp3,temp5,temp7
      real*8 dsc3,dsc5,dsc7
      real*8 psc3,psc5,psc7
      real*8 usc3,usc5
      real*8 alsq2,alsq2n
      real*8 exp2a,ralpha
      real*8 gfd,gfdr
      real*8 xr,yr,zr
      real*8 xix,yix,zix
      real*8 xiy,yiy,ziy
      real*8 xiz,yiz,ziz
      real*8 xkx,ykx,zkx
      real*8 xky,yky,zky
      real*8 xkz,ykz,zkz
      real*8 r,r2,rr1,rr3
      real*8 rr5,rr7,rr9,rr11
      real*8 erl,erli
      real*8 vxx,vyy,vzz
      real*8 vyx,vzx,vzy
      real*8 frcxi(3),frcxk(3)
      real*8 frcyi(3),frcyk(3)
      real*8 frczi(3),frczk(3)
      real*8 ci,di(3),qi(9)
      real*8 ck,dk(3),qk(9)
      real*8 fridmp(3),findmp(3)
      real*8 ftm2(3),ftm2i(3)
      real*8 ftm2r(3),ftm2ri(3)
      real*8 ttm2(3),ttm3(3)
      real*8 ttm2i(3),ttm3i(3)
      real*8 ttm2r(3),ttm3r(3)
      real*8 ttm2ri(3),ttm3ri(3)
      real*8 fdir(3),dixdk(3)
      real*8 dkxui(3),dixuk(3)
      real*8 dixukp(3),dkxuip(3)
      real*8 uixqkr(3),ukxqir(3)
      real*8 uixqkrp(3),ukxqirp(3)
      real*8 qiuk(3),qkui(3)
      real*8 qiukp(3),qkuip(3)
      real*8 rxqiuk(3),rxqkui(3)
      real*8 rxqiukp(3),rxqkuip(3)
      real*8 qidk(3),qkdi(3)
      real*8 qir(3),qkr(3)
      real*8 qiqkr(3),qkqir(3)
      real*8 qixqk(3),rxqir(3)
      real*8 dixr(3),dkxr(3)
      real*8 dixqkr(3),dkxqir(3)
      real*8 rxqkr(3),qkrxqir(3)
      real*8 rxqikr(3),rxqkir(3)
      real*8 rxqidk(3),rxqkdi(3)
      real*8 ddsc3(3),ddsc5(3)
      real*8 ddsc7(3)
      real*8 bn(0:5)
      real*8 sc(10),gl(0:8)
      real*8 sci(8),scip(8)
      real*8 gli(7),glip(7)
      real*8 gf(7),gfi(6)
      real*8 gfr(7),gfri(6)
      real*8 gti(6),gtri(6)
      real*8, allocatable :: mscale(:)
      real*8, allocatable :: pscale(:)
      real*8, allocatable :: dscale(:)
      real*8, allocatable :: uscale(:)
      real*8 aewald
      character*6 mode
      external erfc
c
c
c     zero out the intramolecular portion of the Ewald energy
c
      aewald=aewaldPerm
      eintra = 0.0d0
      if (npole .eq. 0)  return
c
c     perform dynamic allocation of some local arrays
c
      allocate (mscale(n))
      allocate (pscale(n))
      allocate (dscale(n))
      allocate (uscale(n))
c
c     set arrays needed to scale connected atom interactions
c
      do i = 1, n
         mscale(i) = 1.0d0
         pscale(i) = 1.0d0
         dscale(i) = 1.0d0
         uscale(i) = 1.0d0
      end do
c
c     set conversion factor, cutoff and switching coefficients
c
      f = electric / dielec
      mode = 'EWALD'
      call switch (mode)
c
c     set the permanent multipole and induced dipole values
c
      do i = 1, npole-1
         ii = ipole(i)
         pdi = pdamp(i)
         pti = thole(i)
         ci = rpole(1,i)
         di(1) = rpole(2,i)
         di(2) = rpole(3,i)
         di(3) = rpole(4,i)
         qi(1) = rpole(5,i)
         qi(2) = rpole(6,i)
         qi(3) = rpole(7,i)
         qi(4) = rpole(8,i)
         qi(5) = rpole(9,i)
         qi(6) = rpole(10,i)
         qi(7) = rpole(11,i)
         qi(8) = rpole(12,i)
         qi(9) = rpole(13,i)
c
c     set interaction scaling coefficients for connected atoms
c
         do j = 1, n12(ii)
            mscale(i12(j,ii)) = m2scale
            pscale(i12(j,ii)) = p2scale
         end do
         do j = 1, n13(ii)
            mscale(i13(j,ii)) = m3scale
            pscale(i13(j,ii)) = p3scale
         end do
         do j = 1, n14(ii)
            mscale(i14(j,ii)) = m4scale
            pscale(i14(j,ii)) = p4scale
            do k = 1, np11(ii)
                if (i14(j,ii) .eq. ip11(k,ii))
     &            pscale(i14(j,ii)) = p4scale * p41scale
            end do
         end do
         do j = 1, n15(ii)
            mscale(i15(j,ii)) = m5scale
            pscale(i15(j,ii)) = p5scale
         end do
         do j = 1, np11(ii)
            dscale(ip11(j,ii)) = d1scale
            uscale(ip11(j,ii)) = u1scale
         end do
         do j = 1, np12(ii)
            dscale(ip12(j,ii)) = d2scale
            uscale(ip12(j,ii)) = u2scale
         end do
         do j = 1, np13(ii)
            dscale(ip13(j,ii)) = d3scale
            uscale(ip13(j,ii)) = u3scale
         end do
         do j = 1, np14(ii)
            dscale(ip14(j,ii)) = d4scale
            uscale(ip14(j,ii)) = u4scale
         end do
         do k = i+1, npole
            kk = ipole(k)
            xr = x(kk) - x(ii)
            yr = y(kk) - y(ii)
            zr = z(kk) - z(ii)
            call image (xr,yr,zr)
            r2 = xr*xr + yr*yr + zr*zr
            if (r2 .le. off2) then
               r = sqrt(r2)
               ck = rpole(1,k)
               dk(1) = rpole(2,k)
               dk(2) = rpole(3,k)
               dk(3) = rpole(4,k)
               qk(1) = rpole(5,k)
               qk(2) = rpole(6,k)
               qk(3) = rpole(7,k)
               qk(4) = rpole(8,k)
               qk(5) = rpole(9,k)
               qk(6) = rpole(10,k)
               qk(7) = rpole(11,k)
               qk(8) = rpole(12,k)
               qk(9) = rpole(13,k)
c
c     calculate the real space error function terms
c
               ralpha = aewald * r
               bn(0) = erfc(ralpha) / r
               alsq2 = 2.0d0 * aewald**2
               alsq2n = 0.0d0
               if (aewald .gt. 0.0d0)  alsq2n = 1.0d0 / (sqrtpi*aewald)
               exp2a = exp(-ralpha**2)
               do j = 1, 5
                  bfac = dble(2*j-1)
                  alsq2n = alsq2 * alsq2n
                  bn(j) = (bfac*bn(j-1)+alsq2n*exp2a) / r2
               end do
c
c     apply Thole polarization damping to scale factors
c
               rr1 = 1.0d0 / r
               rr3 = rr1 / r2
               rr5 = 3.0d0 * rr3 / r2
               rr7 = 5.0d0 * rr5 / r2
               rr9 = 7.0d0 * rr7 / r2
               rr11 = 9.0d0 * rr9 / r2
               scale3 = 1.0d0
               scale5 = 1.0d0
               scale7 = 1.0d0
               do j = 1, 3
                  ddsc3(j) = 0.0d0
                  ddsc5(j) = 0.0d0
                  ddsc7(j) = 0.0d0
               end do
               damp = pdi * pdamp(k)
               if (damp .ne. 0.0d0) then
                  pgamma = min(pti,thole(k))
                  damp = -pgamma * (r/damp)**3
                  if (damp .gt. -50.0d0) then
                     expdamp = exp(damp)
                     scale3 = 1.0d0 - expdamp
                     scale5 = 1.0d0 - (1.0d0-damp)*expdamp
                     scale7 = 1.0d0 - (1.0d0-damp+0.6d0*damp**2)
     &                                       *expdamp
                     temp3 = -3.0d0 * damp * expdamp / r2
                     temp5 = -damp
                     temp7 = -0.2d0 - 0.6d0*damp
                     ddsc3(1) = temp3 * xr
                     ddsc3(2) = temp3 * yr
                     ddsc3(3) = temp3 * zr
                     ddsc5(1) = temp5 * ddsc3(1)
                     ddsc5(2) = temp5 * ddsc3(2)
                     ddsc5(3) = temp5 * ddsc3(3)
                     ddsc7(1) = temp7 * ddsc5(1)
                     ddsc7(2) = temp7 * ddsc5(2)
                     ddsc7(3) = temp7 * ddsc5(3)
                  end if
               end if
               dsc3 = 1.0d0 - scale3*dscale(kk)
               dsc5 = 1.0d0 - scale5*dscale(kk)
               dsc7 = 1.0d0 - scale7*dscale(kk)
               psc3 = 1.0d0 - scale3*pscale(kk)
               psc5 = 1.0d0 - scale5*pscale(kk)
               psc7 = 1.0d0 - scale7*pscale(kk)
               usc3 = 1.0d0 - scale3*uscale(kk)
               usc5 = 1.0d0 - scale5*uscale(kk)
c
c     construct necessary auxiliary vectors
c
               dixdk(1) = di(2)*dk(3) - di(3)*dk(2)
               dixdk(2) = di(3)*dk(1) - di(1)*dk(3)
               dixdk(3) = di(1)*dk(2) - di(2)*dk(1)
               dixr(1) = di(2)*zr - di(3)*yr
               dixr(2) = di(3)*xr - di(1)*zr
               dixr(3) = di(1)*yr - di(2)*xr
               dkxr(1) = dk(2)*zr - dk(3)*yr
               dkxr(2) = dk(3)*xr - dk(1)*zr
               dkxr(3) = dk(1)*yr - dk(2)*xr
               qir(1) = qi(1)*xr + qi(4)*yr + qi(7)*zr
               qir(2) = qi(2)*xr + qi(5)*yr + qi(8)*zr
               qir(3) = qi(3)*xr + qi(6)*yr + qi(9)*zr
               qkr(1) = qk(1)*xr + qk(4)*yr + qk(7)*zr
               qkr(2) = qk(2)*xr + qk(5)*yr + qk(8)*zr
               qkr(3) = qk(3)*xr + qk(6)*yr + qk(9)*zr
               qiqkr(1) = qi(1)*qkr(1) + qi(4)*qkr(2) + qi(7)*qkr(3)
               qiqkr(2) = qi(2)*qkr(1) + qi(5)*qkr(2) + qi(8)*qkr(3)
               qiqkr(3) = qi(3)*qkr(1) + qi(6)*qkr(2) + qi(9)*qkr(3)
               qkqir(1) = qk(1)*qir(1) + qk(4)*qir(2) + qk(7)*qir(3)
               qkqir(2) = qk(2)*qir(1) + qk(5)*qir(2) + qk(8)*qir(3)
               qkqir(3) = qk(3)*qir(1) + qk(6)*qir(2) + qk(9)*qir(3)
               qixqk(1) = qi(2)*qk(3) + qi(5)*qk(6) + qi(8)*qk(9)
     &                       - qi(3)*qk(2) - qi(6)*qk(5) - qi(9)*qk(8)
               qixqk(2) = qi(3)*qk(1) + qi(6)*qk(4) + qi(9)*qk(7)
     &                       - qi(1)*qk(3) - qi(4)*qk(6) - qi(7)*qk(9)
               qixqk(3) = qi(1)*qk(2) + qi(4)*qk(5) + qi(7)*qk(8)
     &                       - qi(2)*qk(1) - qi(5)*qk(4) - qi(8)*qk(7)
               rxqir(1) = yr*qir(3) - zr*qir(2)
               rxqir(2) = zr*qir(1) - xr*qir(3)
               rxqir(3) = xr*qir(2) - yr*qir(1)
               rxqkr(1) = yr*qkr(3) - zr*qkr(2)
               rxqkr(2) = zr*qkr(1) - xr*qkr(3)
               rxqkr(3) = xr*qkr(2) - yr*qkr(1)
               rxqikr(1) = yr*qiqkr(3) - zr*qiqkr(2)
               rxqikr(2) = zr*qiqkr(1) - xr*qiqkr(3)
               rxqikr(3) = xr*qiqkr(2) - yr*qiqkr(1)
               rxqkir(1) = yr*qkqir(3) - zr*qkqir(2)
               rxqkir(2) = zr*qkqir(1) - xr*qkqir(3)
               rxqkir(3) = xr*qkqir(2) - yr*qkqir(1)
               qkrxqir(1) = qkr(2)*qir(3) - qkr(3)*qir(2)
               qkrxqir(2) = qkr(3)*qir(1) - qkr(1)*qir(3)
               qkrxqir(3) = qkr(1)*qir(2) - qkr(2)*qir(1)
               qidk(1) = qi(1)*dk(1) + qi(4)*dk(2) + qi(7)*dk(3)
               qidk(2) = qi(2)*dk(1) + qi(5)*dk(2) + qi(8)*dk(3)
               qidk(3) = qi(3)*dk(1) + qi(6)*dk(2) + qi(9)*dk(3)
               qkdi(1) = qk(1)*di(1) + qk(4)*di(2) + qk(7)*di(3)
               qkdi(2) = qk(2)*di(1) + qk(5)*di(2) + qk(8)*di(3)
               qkdi(3) = qk(3)*di(1) + qk(6)*di(2) + qk(9)*di(3)
               dixqkr(1) = di(2)*qkr(3) - di(3)*qkr(2)
               dixqkr(2) = di(3)*qkr(1) - di(1)*qkr(3)
               dixqkr(3) = di(1)*qkr(2) - di(2)*qkr(1)
               dkxqir(1) = dk(2)*qir(3) - dk(3)*qir(2)
               dkxqir(2) = dk(3)*qir(1) - dk(1)*qir(3)
               dkxqir(3) = dk(1)*qir(2) - dk(2)*qir(1)
               rxqidk(1) = yr*qidk(3) - zr*qidk(2)
               rxqidk(2) = zr*qidk(1) - xr*qidk(3)
               rxqidk(3) = xr*qidk(2) - yr*qidk(1)
               rxqkdi(1) = yr*qkdi(3) - zr*qkdi(2)
               rxqkdi(2) = zr*qkdi(1) - xr*qkdi(3)
               rxqkdi(3) = xr*qkdi(2) - yr*qkdi(1)
c
c     calculate the scalar products for permanent components
c
               sc(2) = di(1)*dk(1) + di(2)*dk(2) + di(3)*dk(3)
               sc(3) = di(1)*xr + di(2)*yr + di(3)*zr
               sc(4) = dk(1)*xr + dk(2)*yr + dk(3)*zr
               sc(5) = qir(1)*xr + qir(2)*yr + qir(3)*zr
               sc(6) = qkr(1)*xr + qkr(2)*yr + qkr(3)*zr
               sc(7) = qir(1)*dk(1) + qir(2)*dk(2) + qir(3)*dk(3)
               sc(8) = qkr(1)*di(1) + qkr(2)*di(2) + qkr(3)*di(3)
               sc(9) = qir(1)*qkr(1) + qir(2)*qkr(2) + qir(3)*qkr(3)
               sc(10) = qi(1)*qk(1) + qi(2)*qk(2) + qi(3)*qk(3)
     &                     + qi(4)*qk(4) + qi(5)*qk(5) + qi(6)*qk(6)
     &                     + qi(7)*qk(7) + qi(8)*qk(8) + qi(9)*qk(9)
c
c     calculate the gl functions for permanent components
c
               gl(0) = ci*ck
               gl(1) = ck*sc(3) - ci*sc(4)
               gl(2) = ci*sc(6) + ck*sc(5) - sc(3)*sc(4)
               gl(3) = sc(3)*sc(6) - sc(4)*sc(5)
               gl(4) = sc(5)*sc(6)
               gl(5) = -4.0d0 * sc(9)
               gl(6) = sc(2)
               gl(7) = 2.0d0 * (sc(7)-sc(8))
               gl(8) = 2.0d0 * sc(10)
c
c     compute the energy contributions for this interaction
c
               e = bn(0)*gl(0) + bn(1)*(gl(1)+gl(6))
     &                + bn(2)*(gl(2)+gl(7)+gl(8))
     &                + bn(3)*(gl(3)+gl(5)) + bn(4)*gl(4)
c
c     get the real energy without any screening function
c
               erl = rr1*gl(0) + rr3*(gl(1)+gl(6))
     &                  + rr5*(gl(2)+gl(7)+gl(8))
     &                  + rr7*(gl(3)+gl(5)) + rr9*gl(4)
               e = e - (1.0d0-mscale(kk))*erl
               e = f * e
               em = em + e
c
c     increment the total intramolecular energy; assumes
c     intramolecular distances are less than half of cell
c     length and less than the ewald cutoff
c
               if (molcule(ii) .eq. molcule(kk)) then
                 eintra = eintra + mscale(kk)*erl*f
               end if
c
c     intermediate variables for permanent force terms
c
               gf(1) = bn(1)*gl(0) + bn(2)*(gl(1)+gl(6))
     &                    + bn(3)*(gl(2)+gl(7)+gl(8))
     &                    + bn(4)*(gl(3)+gl(5)) + bn(5)*gl(4)
               gf(2) = -ck*bn(1) + sc(4)*bn(2) - sc(6)*bn(3)
               gf(3) = ci*bn(1) + sc(3)*bn(2) + sc(5)*bn(3)
               gf(4) = 2.0d0 * bn(2)
               gf(5) = 2.0d0 * (-ck*bn(2)+sc(4)*bn(3)-sc(6)*bn(4))
               gf(6) = 2.0d0 * (-ci*bn(2)-sc(3)*bn(3)-sc(5)*bn(4))
               gf(7) = 4.0d0 * bn(3)
               gfr(1) = rr3*gl(0) + rr5*(gl(1)+gl(6))
     &                     + rr7*(gl(2)+gl(7)+gl(8))
     &                     + rr9*(gl(3)+gl(5)) + rr11*gl(4)
               gfr(2) = -ck*rr3 + sc(4)*rr5 - sc(6)*rr7
               gfr(3) = ci*rr3 + sc(3)*rr5 + sc(5)*rr7
               gfr(4) = 2.0d0 * rr5
               gfr(5) = 2.0d0 * (-ck*rr5+sc(4)*rr7-sc(6)*rr9)
               gfr(6) = 2.0d0 * (-ci*rr5-sc(3)*rr7-sc(5)*rr9)
               gfr(7) = 4.0d0 * rr7
c
c     get the permanent force with screening
c
               ftm2(1) = gf(1)*xr + gf(2)*di(1) + gf(3)*dk(1)
     &                      + gf(4)*(qkdi(1)-qidk(1)) + gf(5)*qir(1)
     &                      + gf(6)*qkr(1) + gf(7)*(qiqkr(1)+qkqir(1))
               ftm2(2) = gf(1)*yr + gf(2)*di(2) + gf(3)*dk(2)
     &                      + gf(4)*(qkdi(2)-qidk(2)) + gf(5)*qir(2)
     &                      + gf(6)*qkr(2) + gf(7)*(qiqkr(2)+qkqir(2))
               ftm2(3) = gf(1)*zr + gf(2)*di(3) + gf(3)*dk(3)
     &                      + gf(4)*(qkdi(3)-qidk(3)) + gf(5)*qir(3)
     &                      + gf(6)*qkr(3) + gf(7)*(qiqkr(3)+qkqir(3))
c
c     get the permanent force without screening
c
               ftm2r(1) = gfr(1)*xr + gfr(2)*di(1) + gfr(3)*dk(1)
     &                      + gfr(4)*(qkdi(1)-qidk(1)) + gfr(5)*qir(1)
     &                      + gfr(6)*qkr(1) + gfr(7)*(qiqkr(1)+qkqir(1))
               ftm2r(2) = gfr(1)*yr + gfr(2)*di(2) + gfr(3)*dk(2)
     &                      + gfr(4)*(qkdi(2)-qidk(2)) + gfr(5)*qir(2)
     &                      + gfr(6)*qkr(2) + gfr(7)*(qiqkr(2)+qkqir(2))
               ftm2r(3) = gfr(1)*zr + gfr(2)*di(3) + gfr(3)*dk(3)
     &                      + gfr(4)*(qkdi(3)-qidk(3)) + gfr(5)*qir(3)
     &                      + gfr(6)*qkr(3) + gfr(7)*(qiqkr(3)+qkqir(3))

c
c     account for partially excluded induced interactions
c
c               temp3 = 0.5d0 * rr3 * ((gli(1)+gli(6))*pscale(kk)
c     &                                  +(glip(1)+glip(6))*dscale(kk))
c               temp5 = 0.5d0 * rr5 * ((gli(2)+gli(7))*pscale(kk)
c     &                                  +(glip(2)+glip(7))*dscale(kk))
c               temp7 = 0.5d0 * rr7 * (gli(3)*pscale(kk)
c     &                                  +glip(3)*dscale(kk))
               fridmp(1) = temp3*ddsc3(1) + temp5*ddsc5(1)
     &                        + temp7*ddsc7(1)
               fridmp(2) = temp3*ddsc3(2) + temp5*ddsc5(2)
     &                        + temp7*ddsc7(2)
               fridmp(3) = temp3*ddsc3(3) + temp5*ddsc5(3)
     &                        + temp7*ddsc7(3)
c
c     get the permanent torque with screening
c
               ttm2(1) = -bn(1)*dixdk(1) + gf(2)*dixr(1)
     &           + gf(4)*(dixqkr(1)+dkxqir(1)+rxqidk(1)-2.0d0*qixqk(1))
     &           - gf(5)*rxqir(1) - gf(7)*(rxqikr(1)+qkrxqir(1))
               ttm2(2) = -bn(1)*dixdk(2) + gf(2)*dixr(2)
     &           + gf(4)*(dixqkr(2)+dkxqir(2)+rxqidk(2)-2.0d0*qixqk(2))
     &           - gf(5)*rxqir(2) - gf(7)*(rxqikr(2)+qkrxqir(2))
               ttm2(3) = -bn(1)*dixdk(3) + gf(2)*dixr(3)
     &           + gf(4)*(dixqkr(3)+dkxqir(3)+rxqidk(3)-2.0d0*qixqk(3))
     &           - gf(5)*rxqir(3) - gf(7)*(rxqikr(3)+qkrxqir(3))
               ttm3(1) = bn(1)*dixdk(1) + gf(3)*dkxr(1)
     &           - gf(4)*(dixqkr(1)+dkxqir(1)+rxqkdi(1)-2.0d0*qixqk(1))
     &           - gf(6)*rxqkr(1) - gf(7)*(rxqkir(1)-qkrxqir(1))
               ttm3(2) = bn(1)*dixdk(2) + gf(3)*dkxr(2)
     &           - gf(4)*(dixqkr(2)+dkxqir(2)+rxqkdi(2)-2.0d0*qixqk(2))
     &           - gf(6)*rxqkr(2) - gf(7)*(rxqkir(2)-qkrxqir(2))
               ttm3(3) = bn(1)*dixdk(3) + gf(3)*dkxr(3)
     &           - gf(4)*(dixqkr(3)+dkxqir(3)+rxqkdi(3)-2.0d0*qixqk(3))
     &           - gf(6)*rxqkr(3) - gf(7)*(rxqkir(3)-qkrxqir(3))
c
c     get the permanent torque without screening
c
               ttm2r(1) = -rr3*dixdk(1) + gfr(2)*dixr(1)-gfr(5)*rxqir(1)
     &           + gfr(4)*(dixqkr(1)+dkxqir(1)+rxqidk(1)-2.0d0*qixqk(1))
     &           - gfr(7)*(rxqikr(1)+qkrxqir(1))
               ttm2r(2) = -rr3*dixdk(2) + gfr(2)*dixr(2)-gfr(5)*rxqir(2)
     &           + gfr(4)*(dixqkr(2)+dkxqir(2)+rxqidk(2)-2.0d0*qixqk(2))
     &           - gfr(7)*(rxqikr(2)+qkrxqir(2))
               ttm2r(3) = -rr3*dixdk(3) + gfr(2)*dixr(3)-gfr(5)*rxqir(3)
     &           + gfr(4)*(dixqkr(3)+dkxqir(3)+rxqidk(3)-2.0d0*qixqk(3))
     &           - gfr(7)*(rxqikr(3)+qkrxqir(3))
               ttm3r(1) = rr3*dixdk(1) + gfr(3)*dkxr(1) -gfr(6)*rxqkr(1)
     &           - gfr(4)*(dixqkr(1)+dkxqir(1)+rxqkdi(1)-2.0d0*qixqk(1))
     &           - gfr(7)*(rxqkir(1)-qkrxqir(1))
               ttm3r(2) = rr3*dixdk(2) + gfr(3)*dkxr(2) -gfr(6)*rxqkr(2)
     &           - gfr(4)*(dixqkr(2)+dkxqir(2)+rxqkdi(2)-2.0d0*qixqk(2))
     &           - gfr(7)*(rxqkir(2)-qkrxqir(2))
               ttm3r(3) = rr3*dixdk(3) + gfr(3)*dkxr(3) -gfr(6)*rxqkr(3)
     &           - gfr(4)*(dixqkr(3)+dkxqir(3)+rxqkdi(3)-2.0d0*qixqk(3))
     &           - gfr(7)*(rxqkir(3)-qkrxqir(3))
c
c     handle the case where scaling is used
c
               do j = 1, 3
                  ftm2(j) = f * (ftm2(j)-(1.0d0-mscale(kk))*ftm2r(j))
                  ttm2(j) = f * (ttm2(j)-(1.0d0-mscale(kk))*ttm2r(j))
                  ttm3(j) = f * (ttm3(j)-(1.0d0-mscale(kk))*ttm3r(j))
               end do
c
c     increment gradient due to force and torque on first site
c
               dem(1,ii) = dem(1,ii) + ftm2(1)
               dem(2,ii) = dem(2,ii) + ftm2(2)
               dem(3,ii) = dem(3,ii) + ftm2(3)
               call torque_3b_Perm (i,ttm2,ttm2i,frcxi,frcyi,frczi)
c
c     increment gradient due to force and torque on second site
c
               dem(1,kk) = dem(1,kk) - ftm2(1)
               dem(2,kk) = dem(2,kk) - ftm2(2)
               dem(3,kk) = dem(3,kk) - ftm2(3)
               call torque_3b_Perm (k,ttm3,ttm3i,frcxk,frcyk,frczk)
c
c     increment the internal virial tensor components
c
               iaz = zaxis(i)
               iax = xaxis(i)
               iay = yaxis(i)
               kaz = zaxis(k)
               kax = xaxis(k)
               kay = yaxis(k)
               if (iaz .eq. 0)  iaz = ii
               if (iax .eq. 0)  iax = ii
               if (iay .eq. 0)  iay = ii
               if (kaz .eq. 0)  kaz = kk
               if (kax .eq. 0)  kax = kk
               if (kay .eq. 0)  kay = kk
               xiz = x(iaz) - x(ii)
               yiz = y(iaz) - y(ii)
               ziz = z(iaz) - z(ii)
               xix = x(iax) - x(ii)
               yix = y(iax) - y(ii)
               zix = z(iax) - z(ii)
               xiy = x(iay) - x(ii)
               yiy = y(iay) - y(ii)
               ziy = z(iay) - z(ii)
               xkz = x(kaz) - x(kk)
               ykz = y(kaz) - y(kk)
               zkz = z(kaz) - z(kk)
               xkx = x(kax) - x(kk)
               ykx = y(kax) - y(kk)
               zkx = z(kax) - z(kk)
               xky = x(kay) - x(kk)
               yky = y(kay) - y(kk)
               zky = z(kay) - z(kk)
c      WRONG! NEED frc<xyz>i AND frc<xyz>k
c               vxx = -xr*ftm2(1)+xkx*frcxk(1)+xky*frcyk(1)+xkz*frczk(1)
c               vyx = -yr*ftm2(1)+ykx*frcxk(1)+yky*frcyk(1)+ykz*frczk(1)
c               vzx = -zr*ftm2(1)+zkx*frcxk(1)+zky*frcyk(1)+zkz*frczk(1)
c               vyy = -yr*ftm2(2)+ykx*frcxk(2)+yky*frcyk(2)+ykz*frczk(2)
c               vzy = -zr*ftm2(2)+zkx*frcxk(2)+zky*frcyk(2)+zkz*frczk(2)
c               vzz = -zr*ftm2(3)+zkx*frcxk(3)+zky*frcyk(3)+zkz*frczk(3)
               vxx = -xr*(ftm2(1)) + xix*frcxi(1)
     &                  + xiy*frcyi(1) + xiz*frczi(1) + xkx*frcxk(1)
     &                  + xky*frcyk(1) + xkz*frczk(1)
               vyx = -yr*(ftm2(1)) + yix*frcxi(1)
     &                  + yiy*frcyi(1) + yiz*frczi(1) + ykx*frcxk(1)
     &                  + yky*frcyk(1) + ykz*frczk(1)
               vzx = -zr*(ftm2(1)) + zix*frcxi(1)
     &                  + ziy*frcyi(1) + ziz*frczi(1) + zkx*frcxk(1)
     &                  + zky*frcyk(1) + zkz*frczk(1)
               vyy = -yr*(ftm2(2)) + yix*frcxi(2)
     &                  + yiy*frcyi(2) + yiz*frczi(2) + ykx*frcxk(2)
     &                  + yky*frcyk(2) + ykz*frczk(2)
               vzy = -zr*(ftm2(2)) + zix*frcxi(2)
     &                  + ziy*frcyi(2) + ziz*frczi(2) + zkx*frcxk(2)
     &                  + zky*frcyk(2) + zkz*frczk(2)
               vzz = -zr*(ftm2(3)) + zix*frcxi(3)
     &                  + ziy*frcyi(3) + ziz*frczi(3) + zkx*frcxk(3)
     &                  + zky*frcyk(3) + zkz*frczk(3)
               vir(1,1) = vir(1,1) + vxx
               vir(2,1) = vir(2,1) + vyx
               vir(3,1) = vir(3,1) + vzx
               vir(1,2) = vir(1,2) + vyx
               vir(2,2) = vir(2,2) + vyy
               vir(3,2) = vir(3,2) + vzy
               vir(1,3) = vir(1,3) + vzx
               vir(2,3) = vir(2,3) + vzy
               vir(3,3) = vir(3,3) + vzz
            end if
         end do
c
c     reset interaction scaling coefficients for connected atoms
c
         do j = 1, n12(ii)
            mscale(i12(j,ii)) = 1.0d0
            pscale(i12(j,ii)) = 1.0d0
         end do
         do j = 1, n13(ii)
            mscale(i13(j,ii)) = 1.0d0
            pscale(i13(j,ii)) = 1.0d0
         end do
         do j = 1, n14(ii)
            mscale(i14(j,ii)) = 1.0d0
            pscale(i14(j,ii)) = 1.0d0
         end do
         do j = 1, n15(ii)
            mscale(i15(j,ii)) = 1.0d0
            pscale(i15(j,ii)) = 1.0d0
         end do
         do j = 1, np11(ii)
            dscale(ip11(j,ii)) = 1.0d0
            uscale(ip11(j,ii)) = 1.0d0
         end do
         do j = 1, np12(ii)
            dscale(ip12(j,ii)) = 1.0d0
            uscale(ip12(j,ii)) = 1.0d0
         end do
         do j = 1, np13(ii)
            dscale(ip13(j,ii)) = 1.0d0
            uscale(ip13(j,ii)) = 1.0d0
         end do
         do j = 1, np14(ii)
            dscale(ip14(j,ii)) = 1.0d0
            uscale(ip14(j,ii)) = 1.0d0
         end do
      end do
c
c     set the permanent multipole and induced dipole values
c
      if (use_replica) then
      do i = 1, npole
         ii = ipole(i)
         pdi = pdamp(i)
         pti = thole(i)
         ci = rpole(1,i)
         di(1) = rpole(2,i)
         di(2) = rpole(3,i)
         di(3) = rpole(4,i)
         qi(1) = rpole(5,i)
         qi(2) = rpole(6,i)
         qi(3) = rpole(7,i)
         qi(4) = rpole(8,i)
         qi(5) = rpole(9,i)
         qi(6) = rpole(10,i)
         qi(7) = rpole(11,i)
         qi(8) = rpole(12,i)
         qi(9) = rpole(13,i)
c
c     set interaction scaling coefficients for connected atoms
c
         do j = 1, n12(ii)
            mscale(i12(j,ii)) = m2scale
            pscale(i12(j,ii)) = p2scale
         end do
         do j = 1, n13(ii)
            mscale(i13(j,ii)) = m3scale
            pscale(i13(j,ii)) = p3scale
         end do
         do j = 1, n14(ii)
            mscale(i14(j,ii)) = m4scale
            pscale(i14(j,ii)) = p4scale
            do k = 1, np11(ii)
                if (i14(j,ii) .eq. ip11(k,ii))
     &            pscale(i14(j,ii)) = p4scale * p41scale
            end do
         end do
         do j = 1, n15(ii)
            mscale(i15(j,ii)) = m5scale
            pscale(i15(j,ii)) = p5scale
         end do
         do j = 1, np11(ii)
            dscale(ip11(j,ii)) = d1scale
            uscale(ip11(j,ii)) = u1scale
         end do
         do j = 1, np12(ii)
            dscale(ip12(j,ii)) = d2scale
            uscale(ip12(j,ii)) = u2scale
         end do
         do j = 1, np13(ii)
            dscale(ip13(j,ii)) = d3scale
            uscale(ip13(j,ii)) = u3scale
         end do
         do j = 1, np14(ii)
            dscale(ip14(j,ii)) = d4scale
            uscale(ip14(j,ii)) = u4scale
         end do
         do k = i, npole
            kk = ipole(k)
            do jcell=1,ncell
            xr = x(kk) - x(ii)
            yr = y(kk) - y(ii)
            zr = z(kk) - z(ii)
            call imager (xr,yr,zr,jcell)
            r2 = xr*xr + yr*yr + zr*zr
            if (.not.(use_polymer .and. r2.le.polycut2)) then
              mscale(kk)=1.0d0
              pscale(kk)=1.0d0
              dscale(kk)=1.0d0
              uscale(kk)=1.0d0
            end if
            if (r2 .le. off2) then
               r = sqrt(r2)
               ck = rpole(1,k)
               dk(1) = rpole(2,k)
               dk(2) = rpole(3,k)
               dk(3) = rpole(4,k)
               qk(1) = rpole(5,k)
               qk(2) = rpole(6,k)
               qk(3) = rpole(7,k)
               qk(4) = rpole(8,k)
               qk(5) = rpole(9,k)
               qk(6) = rpole(10,k)
               qk(7) = rpole(11,k)
               qk(8) = rpole(12,k)
               qk(9) = rpole(13,k)
c
c     calculate the real space error function terms
c
               ralpha = aewald * r
               bn(0) = erfc(ralpha) / r
               alsq2 = 2.0d0 * aewald**2
               alsq2n = 0.0d0
               if (aewald .gt. 0.0d0)  alsq2n = 1.0d0 / (sqrtpi*aewald)
               exp2a = exp(-ralpha**2)
               do j = 1, 5
                  bfac = dble(2*j-1)
                  alsq2n = alsq2 * alsq2n
                  bn(j) = (bfac*bn(j-1)+alsq2n*exp2a) / r2
               end do
c
c     apply Thole polarization damping to scale factors
c
               rr1 = 1.0d0 / r
               rr3 = rr1 / r2
               rr5 = 3.0d0 * rr3 / r2
               rr7 = 5.0d0 * rr5 / r2
               rr9 = 7.0d0 * rr7 / r2
               rr11 = 9.0d0 * rr9 / r2
               scale3 = 1.0d0
               scale5 = 1.0d0
               scale7 = 1.0d0
               do j = 1, 3
                  ddsc3(j) = 0.0d0
                  ddsc5(j) = 0.0d0
                  ddsc7(j) = 0.0d0
               end do
               damp = pdi * pdamp(k)
               if (damp .ne. 0.0d0) then
                  pgamma = min(pti,thole(k))
                  damp = -pgamma * (r/damp)**3
                  if (damp .gt. -50.0d0) then
                     expdamp = exp(damp)
                     scale3 = 1.0d0 - expdamp
                     scale5 = 1.0d0 - (1.0d0-damp)*expdamp
                     scale7 = 1.0d0 - (1.0d0-damp+0.6d0*damp**2)
     &                                       *expdamp
                     temp3 = -3.0d0 * damp * expdamp / r2
                     temp5 = -damp
                     temp7 = -0.2d0 - 0.6d0*damp
                     ddsc3(1) = temp3 * xr
                     ddsc3(2) = temp3 * yr
                     ddsc3(3) = temp3 * zr
                     ddsc5(1) = temp5 * ddsc3(1)
                     ddsc5(2) = temp5 * ddsc3(2)
                     ddsc5(3) = temp5 * ddsc3(3)
                     ddsc7(1) = temp7 * ddsc5(1)
                     ddsc7(2) = temp7 * ddsc5(2)
                     ddsc7(3) = temp7 * ddsc5(3)
                  end if
               end if
               dsc3 = 1.0d0 - scale3*dscale(kk)
               dsc5 = 1.0d0 - scale5*dscale(kk)
               dsc7 = 1.0d0 - scale7*dscale(kk)
               psc3 = 1.0d0 - scale3*pscale(kk)
               psc5 = 1.0d0 - scale5*pscale(kk)
               psc7 = 1.0d0 - scale7*pscale(kk)
               usc3 = 1.0d0 - scale3*uscale(kk)
               usc5 = 1.0d0 - scale5*uscale(kk)
c
c     construct necessary auxiliary vectors
c
               dixdk(1) = di(2)*dk(3) - di(3)*dk(2)
               dixdk(2) = di(3)*dk(1) - di(1)*dk(3)
               dixdk(3) = di(1)*dk(2) - di(2)*dk(1)
               dixr(1) = di(2)*zr - di(3)*yr
               dixr(2) = di(3)*xr - di(1)*zr
               dixr(3) = di(1)*yr - di(2)*xr
               dkxr(1) = dk(2)*zr - dk(3)*yr
               dkxr(2) = dk(3)*xr - dk(1)*zr
               dkxr(3) = dk(1)*yr - dk(2)*xr
               qir(1) = qi(1)*xr + qi(4)*yr + qi(7)*zr
               qir(2) = qi(2)*xr + qi(5)*yr + qi(8)*zr
               qir(3) = qi(3)*xr + qi(6)*yr + qi(9)*zr
               qkr(1) = qk(1)*xr + qk(4)*yr + qk(7)*zr
               qkr(2) = qk(2)*xr + qk(5)*yr + qk(8)*zr
               qkr(3) = qk(3)*xr + qk(6)*yr + qk(9)*zr
               qiqkr(1) = qi(1)*qkr(1) + qi(4)*qkr(2) + qi(7)*qkr(3)
               qiqkr(2) = qi(2)*qkr(1) + qi(5)*qkr(2) + qi(8)*qkr(3)
               qiqkr(3) = qi(3)*qkr(1) + qi(6)*qkr(2) + qi(9)*qkr(3)
               qkqir(1) = qk(1)*qir(1) + qk(4)*qir(2) + qk(7)*qir(3)
               qkqir(2) = qk(2)*qir(1) + qk(5)*qir(2) + qk(8)*qir(3)
               qkqir(3) = qk(3)*qir(1) + qk(6)*qir(2) + qk(9)*qir(3)
               qixqk(1) = qi(2)*qk(3) + qi(5)*qk(6) + qi(8)*qk(9)
     &                       - qi(3)*qk(2) - qi(6)*qk(5) - qi(9)*qk(8)
               qixqk(2) = qi(3)*qk(1) + qi(6)*qk(4) + qi(9)*qk(7)
     &                       - qi(1)*qk(3) - qi(4)*qk(6) - qi(7)*qk(9)
               qixqk(3) = qi(1)*qk(2) + qi(4)*qk(5) + qi(7)*qk(8)
     &                       - qi(2)*qk(1) - qi(5)*qk(4) - qi(8)*qk(7)
               rxqir(1) = yr*qir(3) - zr*qir(2)
               rxqir(2) = zr*qir(1) - xr*qir(3)
               rxqir(3) = xr*qir(2) - yr*qir(1)
               rxqkr(1) = yr*qkr(3) - zr*qkr(2)
               rxqkr(2) = zr*qkr(1) - xr*qkr(3)
               rxqkr(3) = xr*qkr(2) - yr*qkr(1)
               rxqikr(1) = yr*qiqkr(3) - zr*qiqkr(2)
               rxqikr(2) = zr*qiqkr(1) - xr*qiqkr(3)
               rxqikr(3) = xr*qiqkr(2) - yr*qiqkr(1)
               rxqkir(1) = yr*qkqir(3) - zr*qkqir(2)
               rxqkir(2) = zr*qkqir(1) - xr*qkqir(3)
               rxqkir(3) = xr*qkqir(2) - yr*qkqir(1)
               qkrxqir(1) = qkr(2)*qir(3) - qkr(3)*qir(2)
               qkrxqir(2) = qkr(3)*qir(1) - qkr(1)*qir(3)
               qkrxqir(3) = qkr(1)*qir(2) - qkr(2)*qir(1)
               qidk(1) = qi(1)*dk(1) + qi(4)*dk(2) + qi(7)*dk(3)
               qidk(2) = qi(2)*dk(1) + qi(5)*dk(2) + qi(8)*dk(3)
               qidk(3) = qi(3)*dk(1) + qi(6)*dk(2) + qi(9)*dk(3)
               qkdi(1) = qk(1)*di(1) + qk(4)*di(2) + qk(7)*di(3)
               qkdi(2) = qk(2)*di(1) + qk(5)*di(2) + qk(8)*di(3)
               qkdi(3) = qk(3)*di(1) + qk(6)*di(2) + qk(9)*di(3)
               dixqkr(1) = di(2)*qkr(3) - di(3)*qkr(2)
               dixqkr(2) = di(3)*qkr(1) - di(1)*qkr(3)
               dixqkr(3) = di(1)*qkr(2) - di(2)*qkr(1)
               dkxqir(1) = dk(2)*qir(3) - dk(3)*qir(2)
               dkxqir(2) = dk(3)*qir(1) - dk(1)*qir(3)
               dkxqir(3) = dk(1)*qir(2) - dk(2)*qir(1)
               rxqidk(1) = yr*qidk(3) - zr*qidk(2)
               rxqidk(2) = zr*qidk(1) - xr*qidk(3)
               rxqidk(3) = xr*qidk(2) - yr*qidk(1)
               rxqkdi(1) = yr*qkdi(3) - zr*qkdi(2)
               rxqkdi(2) = zr*qkdi(1) - xr*qkdi(3)
               rxqkdi(3) = xr*qkdi(2) - yr*qkdi(1)
c
c     calculate the scalar products for permanent components
c
               sc(2) = di(1)*dk(1) + di(2)*dk(2) + di(3)*dk(3)
               sc(3) = di(1)*xr + di(2)*yr + di(3)*zr
               sc(4) = dk(1)*xr + dk(2)*yr + dk(3)*zr
               sc(5) = qir(1)*xr + qir(2)*yr + qir(3)*zr
               sc(6) = qkr(1)*xr + qkr(2)*yr + qkr(3)*zr
               sc(7) = qir(1)*dk(1) + qir(2)*dk(2) + qir(3)*dk(3)
               sc(8) = qkr(1)*di(1) + qkr(2)*di(2) + qkr(3)*di(3)
               sc(9) = qir(1)*qkr(1) + qir(2)*qkr(2) + qir(3)*qkr(3)
               sc(10) = qi(1)*qk(1) + qi(2)*qk(2) + qi(3)*qk(3)
     &                     + qi(4)*qk(4) + qi(5)*qk(5) + qi(6)*qk(6)
     &                     + qi(7)*qk(7) + qi(8)*qk(8) + qi(9)*qk(9)
c
c     calculate the gl functions for permanent components
c
               gl(0) = ci*ck
               gl(1) = ck*sc(3) - ci*sc(4)
               gl(2) = ci*sc(6) + ck*sc(5) - sc(3)*sc(4)
               gl(3) = sc(3)*sc(6) - sc(4)*sc(5)
               gl(4) = sc(5)*sc(6)
               gl(5) = -4.0d0 * sc(9)
               gl(6) = sc(2)
               gl(7) = 2.0d0 * (sc(7)-sc(8))
               gl(8) = 2.0d0 * sc(10)
c
c     compute the energy contributions for this interaction
c
               e = bn(0)*gl(0) + bn(1)*(gl(1)+gl(6))
     &                + bn(2)*(gl(2)+gl(7)+gl(8))
     &                + bn(3)*(gl(3)+gl(5)) + bn(4)*gl(4)
c
c     get the real energy without any screening function
c
               erl = rr1*gl(0) + rr3*(gl(1)+gl(6))
     &                  + rr5*(gl(2)+gl(7)+gl(8))
     &                  + rr7*(gl(3)+gl(5)) + rr9*gl(4)
               if (use_polymer .and. r2.le.polycut2)
     &          e = e - (1.0d0-mscale(kk))*erl
               e = f * e
               if (ii .eq. kk) then
                 e = 0.5d0 * e
               end if
               em = em + e
c
c     increment the total intramolecular energy; assumes
c     intramolecular distances are less than half of cell
c     length and less than the ewald cutoff
c
               if (molcule(ii) .eq. molcule(kk)) then
                  eintra = eintra + mscale(kk)*erl*f
               end if
c
c     intermediate variables for permanent force terms
c
               gf(1) = bn(1)*gl(0) + bn(2)*(gl(1)+gl(6))
     &                    + bn(3)*(gl(2)+gl(7)+gl(8))
     &                    + bn(4)*(gl(3)+gl(5)) + bn(5)*gl(4)
               gf(2) = -ck*bn(1) + sc(4)*bn(2) - sc(6)*bn(3)
               gf(3) = ci*bn(1) + sc(3)*bn(2) + sc(5)*bn(3)
               gf(4) = 2.0d0 * bn(2)
               gf(5) = 2.0d0 * (-ck*bn(2)+sc(4)*bn(3)-sc(6)*bn(4))
               gf(6) = 2.0d0 * (-ci*bn(2)-sc(3)*bn(3)-sc(5)*bn(4))
               gf(7) = 4.0d0 * bn(3)
               gfr(1) = rr3*gl(0) + rr5*(gl(1)+gl(6))
     &                     + rr7*(gl(2)+gl(7)+gl(8))
     &                     + rr9*(gl(3)+gl(5)) + rr11*gl(4)
               gfr(2) = -ck*rr3 + sc(4)*rr5 - sc(6)*rr7
               gfr(3) = ci*rr3 + sc(3)*rr5 + sc(5)*rr7
               gfr(4) = 2.0d0 * rr5
               gfr(5) = 2.0d0 * (-ck*rr5+sc(4)*rr7-sc(6)*rr9)
               gfr(6) = 2.0d0 * (-ci*rr5-sc(3)*rr7-sc(5)*rr9)
               gfr(7) = 4.0d0 * rr7
c
c     get the permanent force with screening
c
               ftm2(1) = gf(1)*xr + gf(2)*di(1) + gf(3)*dk(1)
     &                      + gf(4)*(qkdi(1)-qidk(1)) + gf(5)*qir(1)
     &                      + gf(6)*qkr(1) + gf(7)*(qiqkr(1)+qkqir(1))
               ftm2(2) = gf(1)*yr + gf(2)*di(2) + gf(3)*dk(2)
     &                      + gf(4)*(qkdi(2)-qidk(2)) + gf(5)*qir(2)
     &                      + gf(6)*qkr(2) + gf(7)*(qiqkr(2)+qkqir(2))
               ftm2(3) = gf(1)*zr + gf(2)*di(3) + gf(3)*dk(3)
     &                      + gf(4)*(qkdi(3)-qidk(3)) + gf(5)*qir(3)
     &                      + gf(6)*qkr(3) + gf(7)*(qiqkr(3)+qkqir(3))
c
c     get the permanent force without screening
c
               ftm2r(1) = gfr(1)*xr + gfr(2)*di(1) + gfr(3)*dk(1)
     &                      + gfr(4)*(qkdi(1)-qidk(1)) + gfr(5)*qir(1)
     &                      + gfr(6)*qkr(1) + gfr(7)*(qiqkr(1)+qkqir(1))
               ftm2r(2) = gfr(1)*yr + gfr(2)*di(2) + gfr(3)*dk(2)
     &                      + gfr(4)*(qkdi(2)-qidk(2)) + gfr(5)*qir(2)
     &                      + gfr(6)*qkr(2) + gfr(7)*(qiqkr(2)+qkqir(2))
               ftm2r(3) = gfr(1)*zr + gfr(2)*di(3) + gfr(3)*dk(3)
     &                      + gfr(4)*(qkdi(3)-qidk(3)) + gfr(5)*qir(3)
     &                      + gfr(6)*qkr(3) + gfr(7)*(qiqkr(3)+qkqir(3))
c
c     get the permanent torque with screening
c
               ttm2(1) = -bn(1)*dixdk(1) + gf(2)*dixr(1)
     &           + gf(4)*(dixqkr(1)+dkxqir(1)+rxqidk(1)-2.0d0*qixqk(1))
     &           - gf(5)*rxqir(1) - gf(7)*(rxqikr(1)+qkrxqir(1))
               ttm2(2) = -bn(1)*dixdk(2) + gf(2)*dixr(2)
     &           + gf(4)*(dixqkr(2)+dkxqir(2)+rxqidk(2)-2.0d0*qixqk(2))
     &           - gf(5)*rxqir(2) - gf(7)*(rxqikr(2)+qkrxqir(2))
               ttm2(3) = -bn(1)*dixdk(3) + gf(2)*dixr(3)
     &           + gf(4)*(dixqkr(3)+dkxqir(3)+rxqidk(3)-2.0d0*qixqk(3))
     &           - gf(5)*rxqir(3) - gf(7)*(rxqikr(3)+qkrxqir(3))
               ttm3(1) = bn(1)*dixdk(1) + gf(3)*dkxr(1)
     &           - gf(4)*(dixqkr(1)+dkxqir(1)+rxqkdi(1)-2.0d0*qixqk(1))
     &           - gf(6)*rxqkr(1) - gf(7)*(rxqkir(1)-qkrxqir(1))
               ttm3(2) = bn(1)*dixdk(2) + gf(3)*dkxr(2)
     &           - gf(4)*(dixqkr(2)+dkxqir(2)+rxqkdi(2)-2.0d0*qixqk(2))
     &           - gf(6)*rxqkr(2) - gf(7)*(rxqkir(2)-qkrxqir(2))
               ttm3(3) = bn(1)*dixdk(3) + gf(3)*dkxr(3)
     &           - gf(4)*(dixqkr(3)+dkxqir(3)+rxqkdi(3)-2.0d0*qixqk(3))
     &           - gf(6)*rxqkr(3) - gf(7)*(rxqkir(3)-qkrxqir(3))
c
c     get the permanent torque without screening
c
               ttm2r(1) = -rr3*dixdk(1) + gfr(2)*dixr(1)-gfr(5)*rxqir(1)
     &           + gfr(4)*(dixqkr(1)+dkxqir(1)+rxqidk(1)-2.0d0*qixqk(1))
     &           - gfr(7)*(rxqikr(1)+qkrxqir(1))
               ttm2r(2) = -rr3*dixdk(2) + gfr(2)*dixr(2)-gfr(5)*rxqir(2)
     &           + gfr(4)*(dixqkr(2)+dkxqir(2)+rxqidk(2)-2.0d0*qixqk(2))
     &           - gfr(7)*(rxqikr(2)+qkrxqir(2))
               ttm2r(3) = -rr3*dixdk(3) + gfr(2)*dixr(3)-gfr(5)*rxqir(3)
     &           + gfr(4)*(dixqkr(3)+dkxqir(3)+rxqidk(3)-2.0d0*qixqk(3))
     &           - gfr(7)*(rxqikr(3)+qkrxqir(3))
               ttm3r(1) = rr3*dixdk(1) + gfr(3)*dkxr(1) -gfr(6)*rxqkr(1)
     &           - gfr(4)*(dixqkr(1)+dkxqir(1)+rxqkdi(1)-2.0d0*qixqk(1))
     &           - gfr(7)*(rxqkir(1)-qkrxqir(1))
               ttm3r(2) = rr3*dixdk(2) + gfr(3)*dkxr(2) -gfr(6)*rxqkr(2)
     &           - gfr(4)*(dixqkr(2)+dkxqir(2)+rxqkdi(2)-2.0d0*qixqk(2))
     &           - gfr(7)*(rxqkir(2)-qkrxqir(2))
               ttm3r(3) = rr3*dixdk(3) + gfr(3)*dkxr(3) -gfr(6)*rxqkr(3)
     &           - gfr(4)*(dixqkr(3)+dkxqir(3)+rxqkdi(3)-2.0d0*qixqk(3))
     &           - gfr(7)*(rxqkir(3)-qkrxqir(3))
c
c     handle the case where scaling is used
c
               if (use_polymer .and. r2.le.polycut2) then
               do j = 1, 3
                  ftm2(j) = f * (ftm2(j)-(1.0d0-mscale(kk))*ftm2r(j))
                  ttm2(j) = f * (ttm2(j)-(1.0d0-mscale(kk))*ttm2r(j))
                  ttm3(j) = f * (ttm3(j)-(1.0d0-mscale(kk))*ttm3r(j))
               end do
               else
               do j = 1, 3
                  ftm2(j) = f * ftm2(j)
                  ttm2(j) = f * ttm2(j)
                  ttm3(j) = f * ttm3(j)
               end do
               end if
               if (ii .eq. kk) then
               do j = 1, 3
                  ftm2(j) = 0.5d0 * ftm2(j)
                  ttm2(j) = 0.5d0 * ttm2(j)
                  ttm3(j) = 0.5d0 * ttm3(j)
               end do
               end if    
c
c     increment gradient due to force and torque on first site
c
               dem(1,ii) = dem(1,ii) + ftm2(1)
               dem(2,ii) = dem(2,ii) + ftm2(2)
               dem(3,ii) = dem(3,ii) + ftm2(3)
               call torque_3b_Perm(i,ttm2,ttm2i,frcxi,frcyi,frczi)
c
c     increment gradient due to force and torque on second site
c
               dem(1,kk) = dem(1,kk) - ftm2(1)
               dem(2,kk) = dem(2,kk) - ftm2(2)
               dem(3,kk) = dem(3,kk) - ftm2(3)
               call torque_3b_Perm(k,ttm3,ttm3i,frcxk,frcyk,frczk)
c
c     increment the internal virial tensor components
c
               iaz = zaxis(i)
               iax = xaxis(i)
               iay = yaxis(i)
               kaz = zaxis(k)
               kax = xaxis(k)
               kay = yaxis(k)
               if (iaz .eq. 0)  iaz = ii
               if (iax .eq. 0)  iax = ii
               if (iay .eq. 0)  iay = ii
               if (kaz .eq. 0)  kaz = kk
               if (kax .eq. 0)  kax = kk
               if (kay .eq. 0)  kay = kk
               xiz = x(iaz) - x(ii)
               yiz = y(iaz) - y(ii)
               ziz = z(iaz) - z(ii)
               xix = x(iax) - x(ii)
               yix = y(iax) - y(ii)
               zix = z(iax) - z(ii)
               xiy = x(iay) - x(ii)
               yiy = y(iay) - y(ii)
               ziy = z(iay) - z(ii)
               xkz = x(kaz) - x(kk)
               ykz = y(kaz) - y(kk)
               zkz = z(kaz) - z(kk)
               xkx = x(kax) - x(kk)
               ykx = y(kax) - y(kk)
               zkx = z(kax) - z(kk)
               xky = x(kay) - x(kk)
               yky = y(kay) - y(kk)
               zky = z(kay) - z(kk)
c               vxx = -xr*ftm2(1)+xkx*frcxk(1)+xky*frcyk(1)+xkz*frczk(1)
c               vyx = -yr*ftm2(1)+ykx*frcxk(1)+yky*frcyk(1)+ykz*frczk(1)
c               vzx = -zr*ftm2(1)+zkx*frcxk(1)+zky*frcyk(1)+zkz*frczk(1)
c               vyy = -yr*ftm2(2)+ykx*frcxk(2)+yky*frcyk(2)+ykz*frczk(2)
c               vzy = -zr*ftm2(2)+zkx*frcxk(2)+zky*frcyk(2)+zkz*frczk(2)
c               vzz = -zr*ftm2(3)+zkx*frcxk(3)+zky*frcyk(3)+zkz*frczk(3)
               vxx = -xr*(ftm2(1)) + xix*frcxi(1)
     &                  + xiy*frcyi(1) + xiz*frczi(1) + xkx*frcxk(1)
     &                  + xky*frcyk(1) + xkz*frczk(1)
               vyx = -yr*(ftm2(1)) + yix*frcxi(1)
     &                  + yiy*frcyi(1) + yiz*frczi(1) + ykx*frcxk(1)
     &                  + yky*frcyk(1) + ykz*frczk(1)
               vzx = -zr*(ftm2(1)) + zix*frcxi(1)
     &                  + ziy*frcyi(1) + ziz*frczi(1) + zkx*frcxk(1)
     &                  + zky*frcyk(1) + zkz*frczk(1)
               vyy = -yr*(ftm2(2)) + yix*frcxi(2)
     &                  + yiy*frcyi(2) + yiz*frczi(2) + ykx*frcxk(2)
     &                  + yky*frcyk(2) + ykz*frczk(2)
               vzy = -zr*(ftm2(2)) + zix*frcxi(2)
     &                  + ziy*frcyi(2) + ziz*frczi(2) + zkx*frcxk(2)
     &                  + zky*frcyk(2) + zkz*frczk(2)
               vzz = -zr*(ftm2(3)) + zix*frcxi(3)
     &                  + ziy*frcyi(3) + ziz*frczi(3) + zkx*frcxk(3)
     &                  + zky*frcyk(3) + zkz*frczk(3)
               vir(1,1) = vir(1,1) + vxx
               vir(2,1) = vir(2,1) + vyx
               vir(3,1) = vir(3,1) + vzx
               vir(1,2) = vir(1,2) + vyx
               vir(2,2) = vir(2,2) + vyy
               vir(3,2) = vir(3,2) + vzy
               vir(1,3) = vir(1,3) + vzx
               vir(2,3) = vir(2,3) + vzy
               vir(3,3) = vir(3,3) + vzz
            end if
         end do
         end do
c
c     reset interaction scaling coefficients for connected atoms
c
         do j = 1, n12(ii)
            mscale(i12(j,ii)) = 1.0d0
            pscale(i12(j,ii)) = 1.0d0
         end do
         do j = 1, n13(ii)
            mscale(i13(j,ii)) = 1.0d0
            pscale(i13(j,ii)) = 1.0d0
         end do
         do j = 1, n14(ii)
            mscale(i14(j,ii)) = 1.0d0
            pscale(i14(j,ii)) = 1.0d0
         end do
         do j = 1, n15(ii)
            mscale(i15(j,ii)) = 1.0d0
            pscale(i15(j,ii)) = 1.0d0
         end do
         do j = 1, np11(ii)
            dscale(ip11(j,ii)) = 1.0d0
            uscale(ip11(j,ii)) = 1.0d0
         end do
         do j = 1, np12(ii)
            dscale(ip12(j,ii)) = 1.0d0
            uscale(ip12(j,ii)) = 1.0d0
         end do
         do j = 1, np13(ii)
            dscale(ip13(j,ii)) = 1.0d0
            uscale(ip13(j,ii)) = 1.0d0
         end do
         do j = 1, np14(ii)
            dscale(ip14(j,ii)) = 1.0d0
            uscale(ip14(j,ii)) = 1.0d0
         end do
      end do
      end if
c
c     perform deallocation of some local arrays
c
      deallocate (mscale)
      deallocate (pscale)
      deallocate (dscale)
      deallocate (uscale)
      return
      end
c
c     ###############################################################
c     ##                                                           ##
c     ##                 Subroutine field3b                        ##
c     ##  Liam O-Suilleabhain, Omar Demerdash, Teresa Head-Gordon  ##
c     ##                 Spring 2013                               ##
c     ##                                                           ##
c     ###############################################################
c
c
c     "field_3b" computes the field due to 1-, 2-, or 3-body
c
c
      subroutine field_3b(field,fieldp,npole3b,pnum)
      implicit none
      include 'sizes.i'
      include 'atoms.i'
      include 'boxes.i'
c      include 'combo.i'
      include 'ewald.i'
      include 'inform.i'
      include 'iounit.i'
      include 'math.i'
      include 'mpole.i'
c      include 'polar.i'
      include 'polar2.i'
      include 'polpot.i'
      include 'potent.i'
      include 'units.i'
      include 'uprior.i'
      integer i,j,k,ii,l1
      integer iter,maxiter
      real*8 eps,term
      real*8 epsd,epsp
      real*8 epsold
      real*8 udsum,upsum
      real*8 ucell(3)
      real*8 ucellp(3)
      real*8 field(3,*)
      real*8 fieldp(3,*)
      logical done
      integer npole3b,pnum(*)
      real*8 aewald,aewald3b
      aewald=aewald3b
c
c     get the reciprical space part of the electrostatic field
c
        do i = 1, npole3b
          do j = 1, 3
            field(j,i) = 0.0d0
            fieldp(j,i) = 0.0d0
          end do
        end do
      call udirect1_3b (field,npole3b,pnum)
      do l1 = 1, npole3b
        i=pnum(l1)
        do j = 1, 3
c          fieldp(j,i) = field(j,i)
          fieldp(j,l1)=field(j,l1)
        end do
      end do
c
c     get the real space portion of the electrostatic field
c
      call udirect2a_3b (field,fieldp,npole3b,pnum)
c
c     get the self-energy portion of the electrostatic field
c
      term = (4.0d0/3.0d0) * aewald**3 / sqrtpi
      do l1 = 1, npole3b
        i=pnum(l1)
         do j = 1, 3
c            field(j,i) = field(j,i) + term*rpole(j+1,i)
c            fieldp(j,i) = fieldp(j,i) + term*rpole(j+1,i)
            field(j,l1)=field(j,l1) + term*rpole(j+1,i)
            fieldp(j,l1)=fieldp(j,l1) + term*rpole(j+1,i)
         end do
      end do
      return
      end
c
c
c     #################################################################
c     ##                                                             ##
c     ##                 Subroutine udirect1_3b                      ##
c     ##  Liam O-Suilleabhain, Omar Demerdash, Teresa Head-Gordon    ##
c     ##                 Spring 2013                                 ##
c     ##                                                             ##
c     #################################################################
c
c
c     "udirect1_3b" computes the reciprocal space contribution of
c     the permanent atomic multipole moments to the electrostatic field
c     for use in finding the direct induced dipole moments for 1-, 2-,
c     3-body via a regular Ewald summation
c
c
      subroutine udirect1_3b (field,npole3b,pnum)
      implicit none
      include 'sizes.i'
      include 'atoms.i'
      include 'boxes.i'
      include 'ewald.i'
c      include 'ewreg.i'
      include 'ewreg2.i'       
      include 'math.i'
      include 'mpole.i'
      include 'units.i'
c      include 'molcul.i'
c      include 'combo.i'
      integer l1
      integer i,j,k,l,ii
      integer jmin,jmax
      integer kmin,kmax
      integer lmin,lmax
      real*8 expterm,cut
      real*8 term,fterm
      real*8 xfr,yfr,zfr
      real*8 rj,rk,rl
      real*8 h1,h2,h3,hsq
      real*8 qf,t1,t2
      real*8 ck,dk,qk
      real*8 q1,q2,q3
c      real*8 ckr(maxatm)
c      real*8 skr(maxatm)
c      real*8 cjk(maxatm)
c      real*8 sjk(maxatm)
c      real*8 cm(maxatm)
c      real*8 dm(3,maxatm)
c      real*8 qm(9,maxatm)
      integer npole3b,pnum(*)
      real*8 ckr(npole3b)
      real*8 skr(npole3b)
      real*8 cjk(npole3b)
      real*8 sjk(npole3b)
      real*8 cm(npole3b)
      real*8 dm(3,npole3b)
      real*8 qm(9,npole3b)
      real*8 ejc(npole3b,0:maxvec),ejs(npole3b,0:maxvec)
      real*8 ekc(npole3b,-maxvec:maxvec),eks(npole3b,-maxvec:maxvec)
      real*8 elc(npole3b,-maxvec:maxvec),els(npole3b,-maxvec:maxvec)
      real*8 field(3,*)
      real*8 frecip,aewald,aewald3b

      aewald=aewald3b
c
c     return if the Ewald coefficient is zero
c
      if (aewald .lt. 1.0d-6)  return
      term = -0.25d0 / aewald**2
      fterm = 8.0d0 * pi / volbox
c
c     set the number of vectors based on box dimensions
c
      frecip = 0.5d0
      cut = 4.0d0 * pi * pi * frecip
      jmin = 0
      kmin = 0
      lmin = 1
      jmax = min(maxvec,int(frecip/recip(1,1)))
      kmax = min(maxvec,int(frecip/recip(2,2)))
      lmax = min(maxvec,int(frecip/recip(3,3)))
c
c     copy the multipole moments into local storage areas
c
      do l1 = 1, npole3b
         i = pnum(l1)
         cm(l1) = rpole(1,i)
         dm(1,l1) = rpole(2,i)
         dm(2,l1) = rpole(3,i)
         dm(3,l1) = rpole(4,i)
         qm(1,l1) = rpole(5,i)
         qm(2,l1) = rpole(6,i)
         qm(3,l1) = rpole(7,i)
         qm(4,l1) = rpole(8,i)
         qm(5,l1) = rpole(9,i)
         qm(6,l1) = rpole(10,i)
         qm(7,l1) = rpole(11,i)
         qm(8,l1) = rpole(12,i)
         qm(9,l1) = rpole(13,i)
      end do
c
c     calculate and store the exponential factors
c
      do l1 = 1, npole3b
         i = pnum(l1)
         ii = ipole(i)
         zfr = (z(ii)/gamma_term) / zbox
         yfr = ((y(ii)-zfr*zbox*beta_term)/gamma_sin) / ybox
         xfr = (x(ii)-yfr*ybox*gamma_cos-zfr*zbox*beta_cos) / xbox
         xfr = 2.0d0 * pi * xfr
         yfr = 2.0d0 * pi * yfr
         zfr = 2.0d0 * pi * zfr
c         ejc(i,0) = 1.0d0
c         ejs(i,0) = 0.0d0
c         ekc(i,0) = 1.0d0
c         eks(i,0) = 0.0d0
c         elc(i,0) = 1.0d0
c         els(i,0) = 0.0d0
c         ejc(i,1) = cos(xfr)
c         ejs(i,1) = sin(xfr)
c         ekc(i,1) = cos(yfr)
c         eks(i,1) = sin(yfr)
c         elc(i,1) = cos(zfr)
c         els(i,1) = sin(zfr)
c         ekc(i,-1) = ekc(i,1)
c         eks(i,-1) = -eks(i,1)
c         elc(i,-1) = elc(i,1)
c         els(i,-1) = -els(i,1)
         ejc(l1,0) = 1.0d0
         ejs(l1,0) = 0.0d0
         ekc(l1,0) = 1.0d0
         eks(l1,0) = 0.0d0
         elc(l1,0) = 1.0d0
         els(l1,0) = 0.0d0
         ejc(l1,1) = cos(xfr)
         ejs(l1,1) = sin(xfr)
         ekc(l1,1) = cos(yfr)
         eks(l1,1) = sin(yfr)
         elc(l1,1) = cos(zfr)
         els(l1,1) = sin(zfr)
         ekc(l1,-1) = ekc(l1,1)
         eks(l1,-1) = -eks(l1,1)
         elc(l1,-1) = elc(l1,1)
         els(l1,-1) = -els(l1,1)

c         do j = 2, jmax
c            ejc(i,j) = ejc(i,j-1)*ejc(i,1) - ejs(i,j-1)*ejs(i,1)
c            ejs(i,j) = ejs(i,j-1)*ejc(i,1) + ejc(i,j-1)*ejs(i,1)
c         end do
c         do j = 2, kmax
c            ekc(i,j) = ekc(i,j-1)*ekc(i,1) - eks(i,j-1)*eks(i,1)
c            eks(i,j) = eks(i,j-1)*ekc(i,1) + ekc(i,j-1)*eks(i,1)
c            ekc(i,-j) = ekc(i,j)
c            eks(i,-j) = -eks(i,j)
c         end do
c         do j = 2, lmax
c            elc(i,j) = elc(i,j-1)*elc(i,1) - els(i,j-1)*els(i,1)
c            els(i,j) = els(i,j-1)*elc(i,1) + elc(i,j-1)*els(i,1)
c            elc(i,-j) = elc(i,j)
c            els(i,-j) = -els(i,j)
c         end do
         do j = 2, jmax
            ejc(l1,j) = ejc(l1,j-1)*ejc(l1,1) - ejs(l1,j-1)*ejs(l1,1)
            ejs(l1,j) = ejs(l1,j-1)*ejc(l1,1) + ejc(l1,j-1)*ejs(l1,1)
         end do
         do j = 2, kmax
            ekc(l1,j) = ekc(l1,j-1)*ekc(l1,1) - eks(l1,j-1)*eks(l1,1)
            eks(l1,j) = eks(l1,j-1)*ekc(l1,1) + ekc(l1,j-1)*eks(l1,1)
            ekc(l1,-j) = ekc(l1,j)
            eks(l1,-j) = -eks(l1,j)
         end do
         do j = 2, lmax
            elc(l1,j) = elc(l1,j-1)*elc(l1,1) - els(l1,j-1)*els(l1,1)
            els(l1,j) = els(l1,j-1)*elc(l1,1) + elc(l1,j-1)*els(l1,1)
            elc(l1,-j) = elc(l1,j)
            els(l1,-j) = -els(l1,j)
         end do
      end do
c
c     loop over all k vectors from the reciprocal lattice
c
c
      do j = jmin, jmax
         rj = 2.0d0 * pi * dble(j)
         do k = kmin, kmax
            rk = 2.0d0 * pi * dble(k)
            do l1 = 1, npole3b
               i = pnum(l1)
c               cjk(i) = ejc(i,j)*ekc(i,k) - ejs(i,j)*eks(i,k)
c               sjk(i) = ejs(i,j)*ekc(i,k) + ejc(i,j)*eks(i,k)
               cjk(l1) = ejc(l1,j)*ekc(l1,k) - ejs(l1,j)*eks(l1,k)
               sjk(l1) = ejs(l1,j)*ekc(l1,k) + ejc(l1,j)*eks(l1,k)
            end do
            do l = lmin, lmax
               rl = 2.0d0 * pi * dble(l)
               h1 = recip(1,1)*rj
               h2 = recip(2,1)*rj + recip(2,2)*rk
               h3 = recip(3,1)*rj + recip(3,2)*rk + recip(3,3)*rl
               hsq = h1*h1 + h2*h2 + h3*h3
               if (hsq .le. cut) then
                  t1 = 0.0d0
                  t2 = 0.0d0
                  do l1 = 1, npole3b
                     i = pnum(l1)
c                     ckr(i) = cjk(i)*elc(i,l) - sjk(i)*els(i,l)
c                     skr(i) = sjk(i)*elc(i,l) + cjk(i)*els(i,l)
                     ckr(l1) = cjk(l1)*elc(l1,l) - sjk(l1)*els(l1,l)
                     skr(l1) = sjk(l1)*elc(l1,l) + cjk(l1)*els(l1,l)

c                     ck = cm(i)
c                     dk = h1*dm(1,i) + h2*dm(2,i) + h3*dm(3,i)
c                     q1 = h1*qm(1,i) + h2*qm(4,i) + h3*qm(7,i)
c                     q2 = h1*qm(2,i) + h2*qm(5,i) + h3*qm(8,i)
c                     q3 = h1*qm(3,i) + h2*qm(6,i) + h3*qm(9,i)
                     ck = cm(l1)
                     dk = h1*dm(1,l1) + h2*dm(2,l1) + h3*dm(3,l1)
                     q1 = h1*qm(1,l1) + h2*qm(4,l1) + h3*qm(7,l1)
                     q2 = h1*qm(2,l1) + h2*qm(5,l1) + h3*qm(8,l1)
                     q3 = h1*qm(3,l1) + h2*qm(6,l1) + h3*qm(9,l1)

                     qk = h1*q1 + h2*q2 + h3*q3
c                     t1 = t1 + (ck-qk)*skr(i) + dk*ckr(i)
c                     t2 = t2 + (ck-qk)*ckr(i) - dk*skr(i)
                     t1 = t1 + (ck-qk)*skr(l1) + dk*ckr(l1)
                     t2 = t2 + (ck-qk)*ckr(l1) - dk*skr(l1)

                  end do
                  expterm = fterm * exp(term*hsq) / hsq
                  if (octahedron) then
                     if (mod(j+k+l,2) .ne. 0)  expterm = 0.0d0
                  end if
                  do l1 = 1, npole3b
                     i= pnum(l1)
c                     qf = expterm * (skr(i)*t2-ckr(i)*t1)
c                     field(1,i) = field(1,i) + h1*qf
c                     field(2,i) = field(2,i) + h2*qf
c                     field(3,i) = field(3,i) + h3*qf
                     qf = expterm * (skr(l1)*t2-ckr(l1)*t1)
                     field(1,l1) = field(1,l1) + h1*qf
                     field(2,l1) = field(2,l1) + h2*qf
                     field(3,l1) = field(3,l1) + h3*qf
                  end do
               end if
            end do
            lmin = -lmax
         end do
         kmin = -kmax
      end do
      return
      end
c
c     ################################################################
c     ##                                                            ##
c     ##          Subroutine empole1c_3b_Polar                      ##
c     ##  Liam O-Suilleabhain, Omar Demerdash, Teresa Head-Gordon   ##
c     ##                 Spring 2013                                ##
c     ##                                                            ##
c     ################################################################
c
c
c     "empole1c_3b_Polar" calculates the dipole polarization
c     energy and derivatives with respect to Cartesian coordinates
c     using regular Ewald summation under 3-body approximation
c
c  call empole1c_3b_Polar(npole3b,pnum,eptemp,deptemp,virtemp)
c      subroutine empole1c_3b_Polar(field,fieldp)
c      subroutine empole1c_3b_Polar(npole3b,pnum,eptemp,deptemp,
c     &   virtemp,moli1,moli2,moli3)
      subroutine empole1c_3b_Polar(ewaldcut3b,npole3b,pnum,eptemp,
     &   deptemp,virtemp)
      implicit none
      include 'sizes.i'
      include 'atoms.i'
      include 'boxes.i'
      include 'chgpot.i'
c      include 'deriv.i'
c      include 'energi.i'
      include 'ewald.i'
      include 'inter.i'
      include 'math.i'
      include 'mpole.i'
      include 'polar2.i'
      include 'polpot.i'
c      include 'virial.i'
c      include 'combo.i'
      integer i,j,ii,l1,l3
      integer moli1,moli2,moli3
      real*8 e,ei,eintra
      real*8 f,term,fterm
      real*8 cii,dii,qii,uii
      real*8 xd,yd,zd
      real*8 xu,yu,zu
      real*8 xup,yup,zup
      real*8 xq,yq,zq
      real*8 xv,yv,zv,vterm
      real*8 ci,dix,diy,diz
      real*8 uix,uiy,uiz
      real*8 qixx,qixy,qixz
      real*8 qiyy,qiyz,qizz
      real*8 xdfield,xufield
      real*8 ydfield,yufield
      real*8 zdfield,zufield
      real*8 trq(3),trqi(3)
      real*8 frcx(3),frcy(3),frcz(3)
c      real*8 field(3,*),fieldp(3,*)
      real*8 field(3,npole3b),fieldp(3,npole3b)
      real*8 eptemp,deptemp(3,npole),virtemp(3,3) 
      integer npole3b,pnum(*),pnum2(npole)
      real*8 uind(3,npole3b)
      real*8 uinp(3,npole3b)
      real*8 aewald,ewaldcut3b,aewald3b
      integer jmax3b,kmax3b,lmax3b

        call ewaldcof(aewald3b,ewaldcut3b) 
        jmax3b=int(ewaldcut3b*aewald3b*xbox*aewald3b/pi)
        kmax3b=int(ewaldcut3b*aewald3b*ybox*aewald3b/pi)
        lmax3b=int(ewaldcut3b*aewald3b*zbox*aewald3b/pi)

      aewald=aewald3b
      eptemp = 0.0d0

c   Zero out temporary gradient of polarization energy
      do l1 = 1, npole3b
        i = pnum(l1)
        do j = 1, 3
          deptemp(j,i) = 0.0d0
        end do
      end do

c   Zero out temporary virial
      do i=1,3
         do j=1,3
           virtemp(i,j)=0.0d0
         end do
      end do

c      do l1 =1, npole3b
c        pnum2(pnum(l1))=l1
c      end do

c
c     set the energy unit conversion factor
c
      f = electric / dielec
c
c     compute the induced dipole moment at each atom
c

c      call field_3b(field,fieldp,npole3b,pnum)

c
c     compute the induced dipole moment at each atom
c
      call induce0c_3b_new(aewald3b,ewaldcut3b,jmax3b,kmax3b,lmax3b,
     &     npole3b,pnum,uind,uinp)

c
c     compute the reciprocal space part of the Ewald summation
c

c      call erecip1_3b(npole3b,pnum,pnum2,uind,uinp,
c     &  eptemp,deptemp,virtemp)

      call erecip1_3b(aewald3b,ewaldcut3b,jmax3b,kmax3b,lmax3b,npole3b,
     &  pnum,uind,uinp,eptemp,deptemp,virtemp)

c      subroutine erecip1_3b(npole3b,pnum,pnum2,uind,uinp,
c     &   eptemp,deptemp,virtemp)

c
c     compute the real space part of the Ewald summation
c
c      call ereal1c_3b (eintra)

c      call ereal1c_3b(npole3b,pnum,pnum2,uind,uinp,
c     &   eptemp,deptemp,virtemp,moli1,moli2,moli3)

      call ereal1c_3b(aewald3b,ewaldcut3b,npole3b,pnum,uind,uinp,
     &   eptemp,deptemp,virtemp)

c       subroutine ereal1c_3b(npole3b,pnum,pnum2,uind,uinp,
c     &   eptemp,deptemp,virtemp)


c
c     compute the Ewald self-energy term over all the atoms
c
      term = 2.0d0 * aewald * aewald
      fterm = -f * aewald / sqrtpi
      do l1 = 1, npole3b
         i=pnum(l1)
         dix = rpole(2,i)
         diy = rpole(3,i)
         diz = rpole(4,i)
         uix = uind(1,l1)
         uiy = uind(2,l1)
         uiz = uind(3,l1)
         uii = dix*uix + diy*uiy + diz*uiz
         ei = fterm * term * uii / 3.0d0
         eptemp = eptemp + ei
      end do
c
c     compute the self-energy torque term due to induced dipole
c
      trq(1) = 0.0d0
      trq(2) = 0.0d0
      trq(3) = 0.0d0
      term = (4.0d0/3.0d0) * f * aewald**3 / sqrtpi
      do l1 = 1, npole3b
         i=pnum(l1)
         dix = rpole(2,i)
         diy = rpole(3,i)
         diz = rpole(4,i)
         uix = 0.5d0 * (uind(1,l1)+uinp(1,l1))
         uiy = 0.5d0 * (uind(2,l1)+uinp(2,l1))
         uiz = 0.5d0 * (uind(3,l1)+uinp(3,l1))
         trqi(1) = term * (diy*uiz-diz*uiy)
         trqi(2) = term * (diz*uix-dix*uiz)
         trqi(3) = term * (dix*uiy-diy*uix)
c         call torque_3b (pnum2,deptemp,i,trq,trqi,frcx,frcy,frcz)
         call torque_3b (deptemp,i,trq,trqi,frcx,frcy,frcz)
      end do

c      subroutine torque_3b(pnum2,deptemp,
c     &    i,trq1,trq2,frcx,frcy,frcz)

c
c     compute the cell dipole boundary correction term
c
      if (boundary .eq. 'VACUUM') then
         xd = 0.0d0
         yd = 0.0d0
         zd = 0.0d0
         xu = 0.0d0
         yu = 0.0d0
         zu = 0.0d0
         xup = 0.0d0
         yup = 0.0d0
         zup = 0.0d0
         do l1 = 1, npole3b
            i=pnum(l1)
            ii = ipole(i)
            xd = xd + rpole(2,i) + rpole(1,i)*x(ii)
            yd = yd + rpole(3,i) + rpole(1,i)*y(ii)
            zd = zd + rpole(4,i) + rpole(1,i)*z(ii)
c            xu = xu + uind(1,i)
c            yu = yu + uind(2,i)
c            zu = zu + uind(3,i)
c            xup = xup + uinp(1,i)
c            yup = yup + uinp(2,i)
c            zup = zup + uinp(3,i)
            xu = xu + uind(1,l1)
            yu = yu + uind(2,l1)
            zu = zu + uind(3,l1)
            xup = xup + uinp(1,l1)
            yup = yup + uinp(2,l1)
            zup = zup + uinp(3,l1)

         end do
         term = (2.0d0/3.0d0) * f * (pi/volbox)
         eptemp=eptemp + term*(xd*xu+yd*yu+zd*zu)
         do l1 = 1, npole3b
            i=pnum(l1)
            ii = ipole(i)
c            dep(1,ii) = dep(1,ii) + term*rpole(1,i)*(xu+xup)
c            dep(2,ii) = dep(2,ii) + term*rpole(1,i)*(yu+yup)
c            dep(3,ii) = dep(3,ii) + term*rpole(1,i)*(zu+zup)
            deptemp(1,l1) = deptemp(1,l1) + term*rpole(1,i)*(xu+xup)
            deptemp(2,l1) = deptemp(2,l1) + term*rpole(1,i)*(yu+yup)
            deptemp(3,l1) = deptemp(3,l1) + term*rpole(1,i)*(zu+zup)

         end do
         xdfield = -2.0d0 * term * xd
         ydfield = -2.0d0 * term * yd
         zdfield = -2.0d0 * term * zd
         xufield = -term * (xu+xup)
         yufield = -term * (yu+yup)
         zufield = -term * (zu+zup)
         do l1 = 1, npole3b
            i=pnum(l1)
            trq(1) = rpole(3,i)*zdfield - rpole(4,i)*ydfield
            trq(2) = rpole(4,i)*xdfield - rpole(2,i)*zdfield
            trq(3) = rpole(2,i)*ydfield - rpole(3,i)*xdfield
            trqi(1) = rpole(3,i)*zufield - rpole(4,i)*yufield
            trqi(2) = rpole(4,i)*xufield - rpole(2,i)*zufield
            trqi(3) = rpole(2,i)*yufield - rpole(3,i)*xufield
            call torque_3b (pnum2,deptemp,i,trq,trqi,frcx,frcy,frcz)
         end do
c
c     boundary correction to virial due to overall cell dipole
c
         xd = 0.0d0
         yd = 0.0d0
         zd = 0.0d0
         xq = 0.0d0
         yq = 0.0d0
         zq = 0.0d0
         do l1 = 1, npole3b
            i=pnum(l1)
            ii = ipole(i)
            xd = xd + rpole(2,i)
            yd = yd + rpole(3,i)
            zd = zd + rpole(4,i)
            xq = xq + rpole(1,i)*x(ii)
            yq = yq + rpole(1,i)*y(ii)
            zq = zq + rpole(1,i)*z(ii)
         end do
         xv = xq * (xd+0.5d0*(xu+xup))
         yv = yq * (yd+0.5d0*(yu+yup))
         zv = zq * (zd+0.5d0*(zu+zup))
         vterm = term * (xq*xq + yq*yq + zq*zq + 2.0d0*(xv+yv+zv)
     &                      + xu*xup + yu*yup + zu*zup
     &                      + xd*(xd+xu+xup) + yd*(yd+yu+yup)
     &                      + zd*(zd+zu+zup))
         virtemp(1,1)=virtemp(1,1)+2.0d0*term*(xq*xq+xv) + vterm
         virtemp(2,1)=virtemp(2,1)+2.0d0*term*(xq*yq+xv)
         virtemp(3,1)=virtemp(3,1)+2.0d0*term*(xq*zq+xv)
         virtemp(1,2)=virtemp(1,2) + 2.0d0*term*(yq*xq+yv)
         virtemp(2,2)=virtemp(2,2) + 2.0d0*term*(yq*yq+yv) + vterm
         virtemp(3,2)=virtemp(3,2) + 2.0d0*term*(yq*zq+yv)
         virtemp(1,3)=virtemp(1,3) + 2.0d0*term*(zq*xq+zv)
         virtemp(2,3)=virtemp(2,3) + 2.0d0*term*(zq*yq+zv)
         virtemp(3,3)=virtemp(3,3) + 2.0d0*term*(zq*zq+zv) + vterm
         if (poltyp .eq. 'DIRECT') then
            vterm = term * (xu*xup+yu*yup+zu*zup)
            virtemp(1,1) = virtemp(1,1) + vterm
            virtemp(2,2) = virtemp(2,2) + vterm
            virtemp(3,3) = virtemp(3,3) + vterm
         end if
      end if
c
c     intermolecular energy is total minus intramolecular part
c
c      einter = einter + em + eptemp - eintra
  
c      print*,"End of empole1c_3b_Polar ep=",eptemp,moli1,moli2,moli3
c      do l1=1,npole3b
c         print*,"End of empole1c_3b_Polar l1 Pnum(l1)",l1,pnum(l1)
c      end do
      return
      end
c
c     ###############################################################
c     ##                                                           ##
c     ##              Subroutine induce0c_3b                       ##
c     ##  Liam O-Suilleabhain, Omar Demerdash, Teresa Head-Gordon  ##
c     ##                 Spring 2013                               ##
c     ##                                                           ##
c     ###############################################################
c
c
c     "induce0c_3b" computes the induced dipole moments at polarizable
c     sites using regular ewald and matrix inversion for the 3-body
c
c
      subroutine induce0c_3b(field,fieldp,npole3b,pnum,uind,uinp)
      implicit none
      include 'sizes.i'
c      include 'atoms.i'
      include 'boxes.i'
c      include 'combo.i'
      include 'ewald.i'
      include 'inform.i'
      include 'iounit.i'
      include 'math.i'
c      include 'mpole.i'
      include 'polar2.i'
      include 'polpot.i'
      include 'potent.i'
      include 'units.i'
c      include 'uprior.i'
c      include 'molcul.i'
c      include 'pme.i'
      integer i,j,k,ii,i1,k1
      integer l1,l3,i2,j2,k2
c      integer iter,maxiter
c      real*8 eps,term
c      real*8 epsd,epsp
c      real*8 epsold
c      real*8 udsum,upsum
      real*8 ucell(3)
      real*8 ucellp(3)
      real*8 field(3,*)
      real*8 fieldp(3,*)      
c      real*8 liam(3,maxatm)
      integer npole3b,pnum(*)
      real*8 M_recip(3*npole3b,3*npole3b)
      real*8 M_real(3*npole3b,3*npole3b)
      real*8 M_tot(3*npole3b,3*npole3b)
      real timediff,timearray(2)
      logical done, use_liam
      real*8 uind(3,*)
      real*8 uinp(3,*)
      real*8 term,aewald,aewald3b
      aewald=aewald3b
     
c
c     zero out the induced dipole and the field at each site
c
      do i = 1, npole3b
         do j = 1, 3
            uind(j,i) = 0.0d0
            uinp(j,i) = 0.0d0
         end do
      end do
      do l1 = 1, 3*npole3b
         do l3 = 1, 3*npole3b
            M_recip(l1,l3) = 0.0d0
            M_real(l1,l3) = 0.0d0
            M_tot(l1,l3) = 0.0d0
         end do
      end do
      if (.not. use_polar)  return

      call umutual_ewrp_3b (M_recip,npole3b,pnum) 
      call umutual_ewrl_3b (M_real,npole3b,pnum) 

      do l1 = 1, 3*npole3b-1
        do l3 = l1+1, 3*npole3b
           M_tot(l1,l3) = - (M_real(l1,l3) + M_recip(l1,l3))
           M_tot(l3,l1) = - (M_real(l3,l1) + M_recip(l3,l1))
        end do
      end do

      term = (4.0d0/3.0d0) * aewald**3 / sqrtpi
      do l1 = 1, npole3b
        i = pnum(l1)
        do j = 1, 3
          i1 = 3*(l1-1)+j
          M_tot(i1,i1) = - (M_real(i1,i1) + M_recip(i1,i1)) - 
     &                      term + 1.0d0/polarity(i)
        end do
      end do

      call invert(3*npole3b,M_tot)

      do l1 = 1, npole3b
        i = pnum(l1)
        do i1 = 1, 3
          i2 = 3*(l1-1) + i1
          do l3 = 1, npole3b
             k = pnum(l3)
             k2 = 3*(l3-1)
c             uind(i1,i) = uind(i1,i)+M_tot(i2,k2+1)*field(1,k)+
c     &                               M_tot(i2,k2+2)*field(2,k)+
c     &                               M_tot(i2,k2+3)*field(3,k)
c             uind(i1,i) = uind(i1,i)+M_tot(i2,k2+1)*field(1,l3)+
c     &                               M_tot(i2,k2+2)*field(2,l3)+
c     &                               M_tot(i2,k2+3)*field(3,l3)

             uind(i1,l1) = uind(i1,l1)+M_tot(i2,k2+1)*field(1,l3)+
     &                               M_tot(i2,k2+2)*field(2,l3)+
     &                               M_tot(i2,k2+3)*field(3,l3)

          end do
        end do
      end do
      do l1 = 1, npole3b
c         i = pnum(l1)
         do j = 1, 3
            uinp(j,l1)=uind(j,l1)
         end do
      end do
      if (boundary .eq. 'VACUUM') then
         do i = 1, 3
            ucell(i) = 0.0d0
            ucellp(i) = 0.0d0
         end do
         do l1 = 1, npole3b
           i=pnum(l1)
           do j=1,3
c            ucell(j) = ucell(j) + uind(j,i)
c            ucellp(j) = ucellp(j) + uinp(j,i)
            ucell(j) = ucell(j) + uind(j,l1)
            ucellp(j) = ucellp(j) + uinp(j,l1)
           end do
         end do
         term = (4.0d0/3.0d0) * pi/volbox
         do l1 = 1, npole3b
           i=pnum(l1)
            do j = 1, 3
               field(j,i) = field(j,i) - term*ucell(j)
               fieldp(j,i) = fieldp(j,i) - term*ucellp(j)
            end do
         end do
      end if
 
      return
      end
c
c     #################################################################
c     ##                                                             ##
c     ##              Subroutine umutual_ewrp3b                      ##
c     ##  Liam O-Suilleabhain, Omar Demerdash, Teresa Head-Gordon    ##
c     ##                 Spring 2013                                 ##
c     ##                                                             ##
c     #################################################################
c
c
c     "umutual_ewrp3b" computes the reciprocal space contribution of the
c     induced atomic dipole moments to the electrostatic field for
c     use in matrix inversion of induced dipole moments via a
c     regular Ewald summation
c
c
      subroutine umutual_ewrp_3b (M_recip,npole3b,pnum)
      implicit none
      include 'sizes.i'
      include 'atoms.i'
      include 'boxes.i'
      include 'ewald.i'
c      include 'ewreg.i'
      include 'ewreg2.i'
      include 'math.i'
      include 'mpole.i'
c      include 'polar.i'
      include 'polar2.i'
      include 'units.i'
c      include 'combo.i'
      integer i,j,k,l,ii
      integer jmin,jmax
      integer kmin,kmax
      integer lmin,lmax
      integer l1,l3,i1,k1,i2,j2
      real*8 expterm,cut
      real*8 term,fterm
      real*8 xfr,yfr,zfr
      real*8 rj,rk,rl
      real*8 hsq
      real*8 duk,puk
      real*8 dqf,pqf
c      real*8 ckr(maxatm)
c      real*8 skr(maxatm)
c      real*8 cjk(maxatm)
c      real*8 sjk(maxatm)
c      real*8 field(3,maxatm)
c      real*8 fieldp(3,maxatm)
      integer npole3b,pnum(*)
      real*8 frecip
      real*8 h(3)
      integer a,a1,a2,a3
      integer b,b1,b2,b3
      real*8 T1xx,T1xy,T1xz,T1yx,T1yy,T1yz,T1zx,T1zy,T1zz
      real*8 T2xx,T2xy,T2xz,T2yx,T2yy,T2yz,T2zx,T2zy,T2zz
      real*8 T1_rec(3*npole3b,3*npole3b)
      real*8 T2_rec(3*npole3b,3*npole3b)
      real*8 M_recip(3*npole3b,3*npole3b)
      real timediff,timearray(2)
      real*8 ckr(npole3b)
      real*8 skr(npole3b)
      real*8 cjk(npole3b)
      real*8 sjk(npole3b)
      real*8 ejc(npole3b,0:maxvec)
      real*8 ejs(npole3b,0:maxvec)
      real*8 ekc(npole3b,-maxvec:maxvec)
      real*8 eks(npole3b,-maxvec:maxvec)
      real*8 elc(npole3b,-maxvec:maxvec)
      real*8 els(npole3b,-maxvec:maxvec)
      real*8 aewald,aewald3b

      aewald=aewald3b
      frecip = 0.5d0
      do l1 = 1, 3*npole3b
         do l3 = 1, 3*npole3b
            M_recip(l1,l3) = 0.0d0
            T1_rec(l1,l3) = 0.0d0
            T2_rec(l1,l3) = 0.0d0
         end do
      end do
c
c     return if the Ewald coefficient is zero
c
      if (aewald .lt. 1.0d-6)  return
      term = -0.25d0 / aewald**2
      fterm = 8.0d0 * pi / volbox
c
c     set the number of vectors based on box dimensions
c
      cut = 4.0d0 * pi * pi * frecip
      jmin = 0
      kmin = 0
      lmin = 1
      jmax = min(maxvec,int(frecip/recip(1,1)))
      kmax = min(maxvec,int(frecip/recip(2,2)))
      lmax = min(maxvec,int(frecip/recip(3,3)))
c
c     calculate and store the exponential factors
c
      do l1 = 1, npole3b
         i = pnum(l1)
         ii = ipole(i)
         zfr = (z(ii)/gamma_term) / zbox
         yfr = ((y(ii)-zfr*zbox*beta_term)/gamma_sin) / ybox
         xfr = (x(ii)-yfr*ybox*gamma_cos-zfr*zbox*beta_cos) / xbox
         xfr = 2.0d0 * pi * xfr
         yfr = 2.0d0 * pi * yfr
         zfr = 2.0d0 * pi * zfr
c         ejc(i,0) = 1.0d0
c         ejs(i,0) = 0.0d0
c         ekc(i,0) = 1.0d0
c         eks(i,0) = 0.0d0
c         elc(i,0) = 1.0d0
c         els(i,0) = 0.0d0
c         ejc(i,1) = cos(xfr)
c         ejs(i,1) = sin(xfr)
c         ekc(i,1) = cos(yfr)
c         eks(i,1) = sin(yfr)
c         elc(i,1) = cos(zfr)
c         els(i,1) = sin(zfr)
c         ekc(i,-1) = ekc(i,1)
c         eks(i,-1) = -eks(i,1)
c         elc(i,-1) = elc(i,1)
c         els(i,-1) = -els(i,1)
         ejc(l1,0) = 1.0d0
         ejs(l1,0) = 0.0d0
         ekc(l1,0) = 1.0d0
         eks(l1,0) = 0.0d0
         elc(l1,0) = 1.0d0
         els(l1,0) = 0.0d0
         ejc(l1,1) = cos(xfr)
         ejs(l1,1) = sin(xfr)
         ekc(l1,1) = cos(yfr)
         eks(l1,1) = sin(yfr)
         elc(l1,1) = cos(zfr)
         els(l1,1) = sin(zfr)
         ekc(l1,-1) = ekc(l1,1)
         eks(l1,-1) = -eks(l1,1)
         elc(l1,-1) = elc(l1,1)
         els(l1,-1) = -els(l1,1)

c         do j = 2, jmax
c            ejc(i,j) = ejc(i,j-1)*ejc(i,1) - ejs(i,j-1)*ejs(i,1)
c            ejs(i,j) = ejs(i,j-1)*ejc(i,1) + ejc(i,j-1)*ejs(i,1)
c         end do
c         do j = 2, kmax
c            ekc(i,j) = ekc(i,j-1)*ekc(i,1) - eks(i,j-1)*eks(i,1)
c            eks(i,j) = eks(i,j-1)*ekc(i,1) + ekc(i,j-1)*eks(i,1)
c            ekc(i,-j) = ekc(i,j)
c            eks(i,-j) = -eks(i,j)
c         end do
c         do j = 2, lmax
c            elc(i,j) = elc(i,j-1)*elc(i,1) - els(i,j-1)*els(i,1)
c            els(i,j) = els(i,j-1)*elc(i,1) + elc(i,j-1)*els(i,1)
c            elc(i,-j) = elc(i,j)
c            els(i,-j) = -els(i,j)
c         end do
         do j = 2, jmax
            ejc(l1,j) = ejc(l1,j-1)*ejc(l1,1) - ejs(l1,j-1)*ejs(l1,1)
            ejs(l1,j) = ejs(l1,j-1)*ejc(l1,1) + ejc(l1,j-1)*ejs(l1,1)
         end do
         do j = 2, kmax
            ekc(l1,j) = ekc(l1,j-1)*ekc(l1,1) - eks(l1,j-1)*eks(l1,1)
            eks(l1,j) = eks(l1,j-1)*ekc(l1,1) + ekc(l1,j-1)*eks(l1,1)
            ekc(l1,-j) = ekc(l1,j)
            eks(l1,-j) = -eks(l1,j)
         end do
         do j = 2, lmax
            elc(l1,j) = elc(l1,j-1)*elc(l1,1) - els(l1,j-1)*els(l1,1)
            els(l1,j) = els(l1,j-1)*elc(l1,1) + elc(l1,j-1)*els(l1,1)
            elc(l1,-j) = elc(l1,j)
            els(l1,-j) = -els(l1,j)
         end do
      end do
c
c     loop over all k vectors from the reciprocal lattice
c
       do j = jmin, jmax
         rj = 2.0d0 * pi * dble(j)
         do k = kmin, kmax
            rk = 2.0d0 * pi * dble(k)
            do l1 = 1, npole3b
               i = pnum(l1)
c               cjk(i) = ejc(i,j)*ekc(i,k) - ejs(i,j)*eks(i,k)
c               sjk(i) = ejs(i,j)*ekc(i,k) + ejc(i,j)*eks(i,k)
               cjk(l1) = ejc(l1,j)*ekc(l1,k) - ejs(l1,j)*eks(l1,k)
               sjk(l1) = ejs(l1,j)*ekc(l1,k) + ejc(l1,j)*eks(l1,k)
            end do
            do l = lmin, lmax
               rl = 2.0d0 * pi * dble(l)
               h(1) = recip(1,1)*rj
               h(2) = recip(2,1)*rj + recip(2,2)*rk
               h(3) = recip(3,1)*rj + recip(3,2)*rk + recip(3,3)*rl
               hsq = h(1)*h(1) + h(2)*h(2) + h(3)*h(3)
               if (hsq .le. cut) then
                  do l1 = 1, npole3b
                     i = pnum(l1)
                     ckr(l1) = cjk(l1)*elc(l1,l) - sjk(l1)*els(l1,l)
                     skr(l1) = sjk(l1)*elc(l1,l) + cjk(l1)*els(l1,l)

                  end do
                  expterm = fterm * exp(term*hsq) / hsq
                  if (octahedron) then
                     if (mod(j+k+l,2) .ne. 0)  expterm = 0.0d0
                  end if
c
c Matrix of reciprocal stuff
c
                  do l1 = 1, npole3b 
                     a = pnum(l1)
                     a1 = 3*(l1-1)
                     do l3 = l1, npole3b
                        b = pnum(l3)
                        b1 = 3*(l3-1)
c                        T1zz = ckr(a)*ckr(b)*h(3)*h(3)
                        T1xx = ckr(l1)*ckr(l3)*h(1)*h(1)
                        T1xy = ckr(l1)*ckr(l3)*h(1)*h(2)
                        T1xz = ckr(l1)*ckr(l3)*h(1)*h(3)
                        T1yy = ckr(l1)*ckr(l3)*h(2)*h(2)
                        T1yz = ckr(l1)*ckr(l3)*h(2)*h(3)
                        T1zz = ckr(l1)*ckr(l3)*h(3)*h(3)

                        T1yx = T1xy
                        T1zx = T1xz
                        T1zy = T1yz

c                        T2xx = -skr(a)*skr(b)*h(1)*h(1)
c                        T2xy = -skr(a)*skr(b)*h(1)*h(2)
c                        T2xz = -skr(a)*skr(b)*h(1)*h(3)
c                        T2yy = -skr(a)*skr(b)*h(2)*h(2)
c                        T2yz = -skr(a)*skr(b)*h(2)*h(3)
c                        T2zz = -skr(a)*skr(b)*h(3)*h(3)

                        T2xx = -skr(l1)*skr(l3)*h(1)*h(1)
                        T2xy = -skr(l1)*skr(l3)*h(1)*h(2)
                        T2xz = -skr(l1)*skr(l3)*h(1)*h(3)
                        T2yy = -skr(l1)*skr(l3)*h(2)*h(2)
                        T2yz = -skr(l1)*skr(l3)*h(2)*h(3)
                        T2zz = -skr(l1)*skr(l3)*h(3)*h(3)

                        T2yx = T2xy
                        T2zx = T2xz
                        T2zy = T2yz

                        M_recip(a1+1,b1+1) =  M_recip(a1+1,b1+1) +
     &                   expterm*(T2xx - T1xx)
                        M_recip(a1+1,b1+2) =  M_recip(a1+1,b1+2) +
     &                   expterm*(T2xy - T1xy)
                        M_recip(a1+1,b1+3) =  M_recip(a1+1,b1+3) +
     &                   expterm*(T2xz - T1xz)
                        M_recip(a1+2,b1+1) =  M_recip(a1+2,b1+1) +
     &                   expterm*(T2yx - T1yx)
                        M_recip(a1+2,b1+2) =  M_recip(a1+2,b1+2) +
     &                   expterm*(T2yy - T1yy)
                        M_recip(a1+2,b1+3) =  M_recip(a1+2,b1+3) +
     &                   expterm*(T2yz - T1yz)
                        M_recip(a1+3,b1+1) =  M_recip(a1+3,b1+1) +
     &                   expterm*(T2zx - T1zx)
                        M_recip(a1+3,b1+2) =  M_recip(a1+3,b1+2) +
     &                   expterm*(T2zy - T1zy)
                        M_recip(a1+3,b1+3) =  M_recip(a1+3,b1+3) +
     &                   expterm*(T2zz - T1zz)
 
                     end do
                  end do
               end if
            end do
            lmin = -lmax
         end do
         kmin = -kmax
      end do
      do l1 = 1 , 3*npole3b-1
         do l3 = l1+1, 3*npole3b
            M_recip(l3,l1) = M_recip(l1,l3)
         end do
      end do
      return
      end
c
c     ##################################################################
c     ##                                                              ##
c     ##              Subroutine umutual_ewrl_3b                      ##
c     ##  Liam O-Suilleabhain, Omar Demerdash, Teresa Head-Gordon     ##
c     ##                 Spring 2013                                  ##
c     ##                                                              ##
c     ##################################################################
c
c
c     "umutual_ewrl_3b" computes the real space contribution of the induced
c     atomic dipole moments to the field for the 3-body apprrox
c
c
      subroutine umutual_ewrl_3b (M_real,npole3b,pnum)
      implicit none
      include 'sizes.i'
      include 'atoms.i'
      include 'boxes.i'
      include 'bound.i'
      include 'cell.i'
      include 'couple.i'
      include 'ewald.i'
      include 'math.i'
      include 'mpole.i'
c      include 'polar.i'
      include 'polar2.i'
      include 'polgrp.i'
      include 'polpot.i'
      include 'shunt.i'
      include 'units.i'
c      include 'molcul.i'
c      include 'combo.i'
      integer i,j,k,m
      integer ii,kk
      real*8 xr,yr,zr
      real*8 r,r2,rr3,rr5
      real*8 erfc,bfac,exp2a
      real*8 duir,dukr
      real*8 puir,pukr
      real*8 duix,duiy,duiz
      real*8 puix,puiy,puiz
      real*8 dukx,duky,dukz
      real*8 pukx,puky,pukz
      real*8 ralpha
      real*8 alsq2,alsq2n
      real*8 pdi,pti,pgamma
      real*8 damp,expdamp
      real*8 scale3,scale5
      real*8 bn(0:2)
      real*8 fimd(3),fkmd(3)
      real*8 fimp(3),fkmp(3)
      real*8 fid(3),fkd(3)
      real*8 fip(3),fkp(3)
      real*8 field(3,npole3b)
      real*8, allocatable :: dscale(:)
      real*8 term
      character*6 mode
      external erfc
      integer l1,l3,i2,k2
      integer npole3b,pnum(*)
      real*8 M_real(3*npole3b,3*npole3b)
      real*8 Txx,Txy,Txz,Tyx
      real*8 Tyy,Tyz,Tzx,Tzy,Tzz
      real*8 fid_xx,fid_xy,fid_xz
      real*8 fid_yx,fid_yy,fid_yz
      real*8 fid_zx,fid_zy,fid_zz
      real*8 fimd_xx,fimd_xy,fimd_xz
      real*8 fimd_yx,fimd_yy,fimd_yz
      real*8 fimd_zx,fimd_zy,fimd_zz
      real*8 Mk_x,Mk_y,Mk_z,fkx,fky,fkz
      real*8 Mk0_x,Mk0_y,Mk0_z
      real*8 aewald,aewald3b
  
      aewald=aewald3b
c
c     check for multipoles and set cutoff coefficients
c
      if (npole3b .eq. 0)  return
      mode = 'EWALD'
      call switch (mode)
c
c     perform dynamic allocation of some local arrays
c
      allocate (dscale(n))
c
c     set array needed to scale connected atom interactions
c
      do i = 1, n
         dscale(i) = 1.0d0
      end do
c
c     compute the real space portion of the Ewald summation
c
      do l1 = 1, npole3b-1
         i = pnum(l1)
         i2 = 3*(l1-1)
         ii = ipole(i)
         pdi = pdamp(i)
         pti = thole(i)
c         duix = uind(1,i)
c         duiy = uind(2,i)
c         duiz = uind(3,i)
         do j = 1, np11(ii)
            dscale(ip11(j,ii)) = u1scale
         end do
         do j = 1, np12(ii)
            dscale(ip12(j,ii)) = u2scale
         end do
         do j = 1, np13(ii)
            dscale(ip13(j,ii)) = u3scale
         end do
         do j = 1, np14(ii)
            dscale(ip14(j,ii)) = u4scale
         end do
         do l3 = l1+1, npole3b
            k = pnum(l3)
            k2 = 3*(l3-1)
            kk = ipole(k)
            xr = x(kk) - x(ii)
            yr = y(kk) - y(ii)
            zr = z(kk) - z(ii)
            call image (xr,yr,zr)
            r2 = xr*xr + yr* yr + zr*zr
            if (r2 .le. cut2) then
               r = sqrt(r2)
c
c     calculate the error function damping terms
c
               ralpha = aewald * r
               bn(0) = erfc(ralpha) / r
               alsq2 = 2.0d0 * aewald**2
               alsq2n = 0.0d0
               if (aewald .gt. 0.0d0)
     &            alsq2n = 1.0d0 / (sqrtpi*aewald)
                  exp2a = exp(-ralpha**2)
               do j = 1, 2
                  bfac = dble(j+j-1)
                  alsq2n = alsq2 * alsq2n
                  bn(j) = (bfac*bn(j-1)+alsq2n*exp2a) / r2
               end do
c
c     compute the error function scaled and unscaled terms
c
               scale3 = dscale(kk)
               scale5 = dscale(kk)
               damp = pdi * pdamp(k)
               if (damp .ne. 0.0d0) then
                  pgamma = min(pti,thole(k))
                  damp = -pgamma * (r/damp)**3
                  if (damp .gt. -50.0d0) then
                     expdamp = exp(damp)
                     scale3 = scale3 * (1.0d0-expdamp)
                     scale5 = scale5 * (1.0d0-(1.0d0-damp)*expdamp)
                  end if
               end if
               rr3 = (1.0d0-scale3) / (r*r2)
               rr5 = 3.0d0 * (1.0d0-scale5) / (r*r2*r2)
               
               fid_xx = -rr3 + xr*xr*rr5
               fid_xy =  xr*yr*rr5
               fid_xz =  xr*zr*rr5
               fid_yx =  yr*xr*rr5
               fid_yy = -rr3 + yr*yr*rr5
               fid_yz =  yr*zr*rr5
               fid_zx =  zr*xr*rr5
               fid_zy =  zr*yr*rr5
               fid_zz = -rr3 + zr*zr*rr5

               fimd_xx = -bn(1) + xr*xr*bn(2)
               fimd_xy =  xr*yr*bn(2)
               fimd_xz =  xr*zr*bn(2)
               fimd_yx =  yr*xr*bn(2)
               fimd_yy = -bn(1) + yr*yr*bn(2)
               fimd_yz =  yr*zr*bn(2)
               fimd_zx =  zr*xr*bn(2)
               fimd_zy =  zr*yr*bn(2)
               fimd_zz = -bn(1) + zr*zr*bn(2)

               Txx = (fimd_xx - fid_xx)
               Txy = (fimd_xy - fid_xy)
               Txz = (fimd_xz - fid_xz)
               Tyx = (fimd_yx - fid_yx)
               Tyy = (fimd_yy - fid_yy)
               Tyz = (fimd_yz - fid_yz)
               Tzx = (fimd_zx - fid_zx)
               Tzy = (fimd_zy - fid_zy)
               Tzz = (fimd_zz - fid_zz)
               
               M_real(i2+1,k2+1) = Txx
               M_real(i2+1,k2+2) = Txy
               M_real(i2+1,k2+3) = Txz
               M_real(i2+2,k2+1) = Tyx
               M_real(i2+2,k2+2) = Tyy
               M_real(i2+2,k2+3) = Tyz
               M_real(i2+3,k2+1) = Tzx
               M_real(i2+3,k2+2) = Tzy
               M_real(i2+3,k2+3) = Tzz

               M_real(k2+1,i2+1) = Txx
               M_real(k2+1,i2+2) = Txy
               M_real(k2+1,i2+3) = Txz
               M_real(k2+2,i2+1) = Tyx
               M_real(k2+2,i2+2) = Tyy
               M_real(k2+2,i2+3) = Tyz
               M_real(k2+3,i2+1) = Tzx
               M_real(k2+3,i2+2) = Tzy
               M_real(k2+3,i2+3) = Tzz

            end if
         end do
c
c     reset interaction scaling coefficients for connected atoms
c
         do j = 1, np11(ii)
            dscale(ip11(j,ii)) = 1.0d0
         end do
         do j = 1, np12(ii)
            dscale(ip12(j,ii)) = 1.0d0
         end do
         do j = 1, np13(ii)
            dscale(ip13(j,ii)) = 1.0d0
         end do
         do j = 1, np14(ii)

            dscale(ip14(j,ii)) = 1.0d0
         end do
      end do
c
c     periodic boundary for lareg cutoffs via replicates method
c
      if (use_replica) then
        do l1 = 1, npole3b
          i = pnum(l1)
          i2 = 3*(l1-1)
          ii = ipole(i)
          pdi = pdamp(i)
          pti = thole(i)
c          duix = uind(1,i)
c          duiy = uind(2,i)
c          duiz = uind(3,i)
          do j = 1, np11(ii)
             dscale(ip11(j,ii)) = u1scale
          end do
          do j = 1, np12(ii)
             dscale(ip12(j,ii)) = u2scale
          end do
          do j = 1, np13(ii)
             dscale(ip13(j,ii)) = u3scale
          end do
          do j = 1, np14(ii)
             dscale(ip14(j,ii)) = u4scale
          end do
          do l3 = l1, npole3b
            k = pnum(l3)
            k2 = 3*(l3-1)
            kk = ipole(k)
            do m = 1,ncell
              xr = x(kk) - x(ii)
              yr = y(kk) - y(ii)
              zr = z(kk) - z(ii)
              call imager (xr,yr,zr,m)
              r2 = xr*xr + yr* yr + zr*zr
              if (r2 .le. cut2) then
                 r = sqrt(r2)
c
c     calculate the error function damping terms
c
                 ralpha = aewald * r
                 bn(0) = erfc(ralpha) / r
                 alsq2 = 2.0d0 * aewald**2
                 alsq2n = 0.0d0
                 if (aewald .gt. 0.0d0)
     &              alsq2n = 1.0d0 / (sqrtpi*aewald)
                    exp2a = exp(-ralpha**2)
                 do j = 1, 2
                    bfac = dble(j+j-1)
                    alsq2n = alsq2 * alsq2n
                    bn(j) = (bfac*bn(j-1)+alsq2n*exp2a) / r2
                 end do
c
c     compute the error function scaled and unscaled terms
c
               scale3 = 1.0d0
               scale5 = 1.0d0
               damp = pdi * pdamp(k)
               if (damp .ne. 0.0d0) then
                  pgamma = min(pti,thole(k))
                  damp = -pgamma * (r/damp)**3
                  if (damp .gt. -50.0d0) then
                     expdamp = exp(damp)
                     scale3 = scale3 * (1.0d0-expdamp)
                     scale5 = scale5 * (1.0d0-(1.0d0-damp)*expdamp)
                  end if
               end if
               if (use_polymer) then
                 if (r2 .le. polycut2) then
                   scale3 = scale3 * dscale(kk)
                   scale5 = scale5 *dscale(kk)
                 end if
               end if
               rr3 = (1.0d0-scale3) / (r*r2)
               rr5 = 3.0d0 * (1.0d0-scale5) / (r*r2*r2)
               
               fid_xx = -rr3 + xr*xr*rr5
               fid_xy =  xr*yr*rr5
               fid_xz =  xr*zr*rr5
               fid_yx =  yr*xr*rr5
               fid_yy = -rr3 + yr*yr*rr5
               fid_yz =  yr*zr*rr5
               fid_zx =  zr*xr*rr5
               fid_zy =  zr*yr*rr5
               fid_zz = -rr3 + zr*zr*rr5

               fimd_xx = -bn(1) + xr*xr*bn(2)
               fimd_xy =  xr*yr*bn(2)
               fimd_xz =  xr*zr*bn(2)
               fimd_yx =  yr*xr*bn(2)
               fimd_yy = -bn(1) + yr*yr*bn(2)
               fimd_yz =  yr*zr*bn(2)
               fimd_zx =  zr*xr*bn(2)
               fimd_zy =  zr*yr*bn(2)
               fimd_zz = -bn(1) + zr*zr*bn(2)

               Txx = (fimd_xx - fid_xx)
               Txy = (fimd_xy - fid_xy)
               Txz = (fimd_xz - fid_xz)
               Tyx = (fimd_yx - fid_yx)
               Tyy = (fimd_yy - fid_yy)
               Tyz = (fimd_yz - fid_yz)
               Tzx = (fimd_zx - fid_zx)
               Tzy = (fimd_zy - fid_zy)
               Tzz = (fimd_zz - fid_zz)
               
               M_real(i2+1,k2+1) = M_real(i2+1,k2+1)+Txx
               M_real(i2+1,k2+2) = M_real(i2+1,k2+2)+Txy
               M_real(i2+1,k2+3) = M_real(i2+1,k2+3)+Txz
               M_real(i2+2,k2+1) = M_real(i2+2,k2+1)+Tyx
               M_real(i2+2,k2+2) = M_real(i2+2,k2+2)+Tyy
               M_real(i2+2,k2+3) = M_real(i2+2,k2+3)+Tyz
               M_real(i2+3,k2+1) = M_real(i2+3,k2+1)+Tzx
               M_real(i2+3,k2+2) = M_real(i2+3,k2+2)+Tzy
               M_real(i2+3,k2+3) = M_real(i2+3,k2+3)+Tzz

               if (ii.ne.kk) then
                 M_real(k2+1,i2+1) = M_real(k2+1,i2+1)+Txx
                 M_real(k2+1,i2+2) = M_real(k2+1,i2+2)+Txy
                 M_real(k2+1,i2+3) = M_real(k2+1,i2+3)+Txz
                 M_real(k2+2,i2+1) = M_real(k2+2,i2+1)+Tyx
                 M_real(k2+2,i2+2) = M_real(k2+2,i2+2)+Tyy
                 M_real(k2+2,i2+3) = M_real(k2+2,i2+3)+Tyz
                 M_real(k2+3,i2+1) = M_real(k2+3,i2+1)+Tzx
                 M_real(k2+3,i2+2) = M_real(k2+3,i2+2)+Tzy
                 M_real(k2+3,i2+3) = M_real(k2+3,i2+3)+Tzz
               end if

            end if
           end do
         end do
c
c     reset interaction scaling coefficients for connected atoms
c
         do j = 1, np11(ii)
            dscale(ip11(j,ii)) = 1.0d0
         end do
         do j = 1, np12(ii)
            dscale(ip12(j,ii)) = 1.0d0
         end do
         do j = 1, np13(ii)
            dscale(ip13(j,ii)) = 1.0d0
         end do
         do j = 1, np14(ii)
            dscale(ip14(j,ii)) = 1.0d0
         end do
      end do
      end if
c
c     perform deallocation of some local arrays
c
      deallocate (dscale)
      return
      end
c
c     #################################################################
c     ##                                                             ##
c     ##  subroutine erecip1  --  mpole Ewald recip energy & derivs  ##
c     ##  Liam O-Suilleabhain, Omar Demerdash, Teresa Head-Gordon    ##
c     ##                 Spring 2013                                 ##
c     ##                                                             ##
c     #################################################################
c
c
c     "erecip1_3b" evaluates the reciprocal space portion of the regular
c     Ewald summation energy and gradient due to atomic multipole
c     interactions under 3-body
c
c     literature reference:
c
c     W. Smith, "Point Multipoles in the Ewald Summation (Revisited)",
c     CCP5 Newsletter, 46, 18-30, 1998  (see http://www.ccp5.org/)
c
c
c      subroutine erecip1_3b(npole3b,pnum,pnum2,uind,uinp,
c     &   eptemp,deptemp,virtemp)
      subroutine erecip1_3b(aewald3b,ewaldcut3b,jmax3b,kmax3b,lmax3b,
     &   npole3b,pnum,uind,uinp,eptemp,deptemp,virtemp)
      implicit none
      include 'sizes.i'
      include 'atoms.i'
      include 'boxes.i'
      include 'chgpot.i'
c      include 'deriv.i'
c      include 'energi.i'
      include 'ewald.i'
c      include 'ewreg.i'
      include 'ewreg2.i'
      include 'math.i'
      include 'mpole.i'
c      include 'polar.i'
      include 'polar2.i'
      include 'polpot.i'
      include 'units.i'
      include 'virial.i'
c      include 'combo.i'
      integer i,j,k,l,l1
      integer ii,m1,m2
      integer jmin,jmax
      integer kmin,kmax
      integer lmin,lmax
      real*8 e,ei,etot,f,cut
      real*8 expterm,term,eterm
      real*8 uterm,vterm
      real*8 xfr,yfr,zfr
      real*8 rj,rk,rl
      real*8 h1,h2,h3,hsq
      real*8 ck,dk,qk,uk
      real*8 q1,q2,q3
      real*8 t1,t2,t3,t4
      real*8 ukp,t3p,t4p
      real*8 de,det1,det2
      real*8 dei,det1i,det2i
      real*8 wterm(3,3)
      real*8 t5(3,3),t6(3,3)
      real*8 t5u(3,3),t5p(3,3)
      real*8 t6u(3,3),t6p(3,3)
      real*8 qt(3,3),dt(3,3)
      real*8 dtu(3,3),dtp(3,3)
c      real*8 ckr(maxatm),skr(maxatm)
c      real*8 cjk(maxatm),sjk(maxatm)
c      real*8 s1(maxatm),s2(maxatm)
c      real*8 s3(maxatm),s4(maxatm)
c      real*8 s3p(maxatm),s4p(maxatm)
c      real*8 cm(maxatm),dm(3,maxatm)
c      real*8 qm(9,maxatm),um(3,maxatm)
c      real*8 trq(3,maxatm),trqi(3,maxatm)
c      real*8 dkx(maxatm),qkx(maxatm)
c      real*8 dky(maxatm),qky(maxatm)
c      real*8 dkz(maxatm),qkz(maxatm)
c      integer npole3b,pnum(*),pnum2(*)
      integer npole3b,pnum(*)
      real*8 ckr(npole3b),skr(npole3b)
      real*8 cjk(npole3b),sjk(npole3b)
      real*8 s1(npole3b),s2(npole3b)
      real*8 s3(npole3b),s4(npole3b)
      real*8 s3p(npole3b),s4p(npole3b)
      real*8 cm(npole3b),dm(3,npole3b)
      real*8 qm(9,npole3b),um(3,npole3b)
c      real*8 trq(3,npole3b),trqi(3,npole3b)
      real*8 trq(3,npole),trqi(3,npole)
      real*8 dkx(npole3b),qkx(npole3b)
      real*8 dky(npole3b),qky(npole3b)
      real*8 dkz(npole3b),qkz(npole3b)
      real*8 ejc(npole3b,0:maxvec)
      real*8 ejs(npole3b,0:maxvec)
      real*8 ekc(npole3b,-maxvec:maxvec)
      real*8 eks(npole3b,-maxvec:maxvec)
      real*8 elc(npole3b,-maxvec:maxvec)
      real*8 els(npole3b,-maxvec:maxvec)
      real*8 uind(3,*),eptemp,deptemp(3,*),virtemp(3,3)
      real*8 frecip,uinp(3,*)
      real*8 aewald,aewald3b,ewaldcut3b
      integer jmax3b,kmax3b,lmax3b
c
c     return if the Ewald coefficient is zero
c
      aewald=aewald3b

      if (aewald .lt. 1.0d-6)  return
      do l1 = 1, npole3b
         i = pnum(l1)
         do j = 1, 3
c            dep(j,i) = 0.0d0
            uinp(j,l1) = uind(j,l1)
         end do
      end do

      frecip = 0.5d0
      f = electric / dielec
      term = -0.25d0 / aewald**2
      eterm = 4.0d0 * pi * f / volbox
c
c     set the number of vectors based on box dimensions
c
      cut = 4.0d0 * pi * pi * frecip
      jmin = 0
      kmin = 0
      lmin = 1
c      jmax = min(maxvec,int(frecip/recip(1,1)))
c      kmax = min(maxvec,int(frecip/recip(2,2)))
c      lmax = min(maxvec,int(frecip/recip(3,3)))
      jmax=min(maxvec,jmax3b)
      kmax=min(maxvec,kmax3b)
      lmax=min(maxvec,lmax3b)
      
      print*,"jmax,kmax,lmax in erecip1_3b",jmax3b,kmax3b,lmax3b
c
c     copy the multipole moments into local storage areas
c
      do l1 = 1, npole3b
         i = pnum(l1)
c         cm(i) = rpole(1,i)
c         dm(1,i) = rpole(2,i)
c         dm(2,i) = rpole(3,i)
c         dm(3,i) = rpole(4,i)
c         qm(1,i) = rpole(5,i)
c         qm(2,i) = rpole(6,i)
c         qm(3,i) = rpole(7,i)
c         qm(4,i) = rpole(8,i)
c         qm(5,i) = rpole(9,i)
c         qm(6,i) = rpole(10,i)
c         qm(7,i) = rpole(11,i)
c         qm(8,i) = rpole(12,i)
c         qm(9,i) = rpole(13,i)
c         um(1,i) = uind(1,i)
c         um(2,i) = uind(2,i)
c         um(3,i) = uind(3,i)
         cm(l1) = rpole(1,i)
         dm(1,l1) = rpole(2,i)
         dm(2,l1) = rpole(3,i)
         dm(3,l1) = rpole(4,i)
         qm(1,l1) = rpole(5,i)
         qm(2,l1) = rpole(6,i)
         qm(3,l1) = rpole(7,i)
         qm(4,l1) = rpole(8,i)
         qm(5,l1) = rpole(9,i)
         qm(6,l1) = rpole(10,i)
         qm(7,l1) = rpole(11,i)
         qm(8,l1) = rpole(12,i)
         qm(9,l1) = rpole(13,i)
         um(1,l1) = uind(1,l1)
         um(2,l1) = uind(2,l1)
         um(3,l1) = uind(3,l1)
      end do
c
c     zero out local accumulation arrays for derivatives
c
      do l1 = 1, npole3b
         i = pnum(l1)
c         trq(1,i) = 0.0d0
c         trq(2,i) = 0.0d0
c         trq(3,i) = 0.0d0
         trqi(1,i) = 0.0d0
         trqi(2,i) = 0.0d0
         trqi(3,i) = 0.0d0

c         trq(1,l1) = 0.0d0
c         trq(2,l1) = 0.0d0
c         trq(3,l1) = 0.0d0
c         trqi(1,l1) = 0.0d0
c         trqi(2,l1) = 0.0d0
c         trqi(3,l1) = 0.0d0
      end do
c
c     calculate and store the exponential factors
c
      do l1 = 1, npole3b
         i = pnum(l1)
         ii = ipole(i)
         zfr = (z(ii)/gamma_term) / zbox
         yfr = ((y(ii)-zfr*zbox*beta_term)/gamma_sin) / ybox
         xfr = (x(ii)-yfr*ybox*gamma_cos-zfr*zbox*beta_cos) / xbox
         xfr = 2.0d0 * pi * xfr
         yfr = 2.0d0 * pi * yfr
         zfr = 2.0d0 * pi * zfr
c         ejc(i,0) = 1.0d0
c         ejs(i,0) = 0.0d0
c         ekc(i,0) = 1.0d0
c         eks(i,0) = 0.0d0
c         elc(i,0) = 1.0d0
c         els(i,0) = 0.0d0
c         ejc(i,1) = cos(xfr)
c         ejs(i,1) = sin(xfr)
c         ekc(i,1) = cos(yfr)
c         eks(i,1) = sin(yfr)
c         elc(i,1) = cos(zfr)
c         els(i,1) = sin(zfr)
c         ekc(i,-1) = ekc(i,1)
c         eks(i,-1) = -eks(i,1)
c         elc(i,-1) = elc(i,1)
c         els(i,-1) = -els(i,1)
c         do j = 2, jmax
c            ejc(i,j) = ejc(i,j-1)*ejc(i,1) - ejs(i,j-1)*ejs(i,1)
c            ejs(i,j) = ejs(i,j-1)*ejc(i,1) + ejc(i,j-1)*ejs(i,1)
c         end do
c         do j = 2, kmax
c            ekc(i,j) = ekc(i,j-1)*ekc(i,1) - eks(i,j-1)*eks(i,1)
c            eks(i,j) = eks(i,j-1)*ekc(i,1) + ekc(i,j-1)*eks(i,1)
c            ekc(i,-j) = ekc(i,j)
c            eks(i,-j) = -eks(i,j)
c         end do
c         do j = 2, lmax
c            elc(i,j) = elc(i,j-1)*elc(i,1) - els(i,j-1)*els(i,1)
c            els(i,j) = els(i,j-1)*elc(i,1) + elc(i,j-1)*els(i,1)
c            elc(i,-j) = elc(i,j)
c            els(i,-j) = -els(i,j)
c         end do
         ejc(l1,0) = 1.0d0
         ejs(l1,0) = 0.0d0
         ekc(l1,0) = 1.0d0
         eks(l1,0) = 0.0d0
         elc(l1,0) = 1.0d0
         els(l1,0) = 0.0d0
         ejc(l1,1) = cos(xfr)
         ejs(l1,1) = sin(xfr)
         ekc(l1,1) = cos(yfr)
         eks(l1,1) = sin(yfr)
         elc(l1,1) = cos(zfr)
         els(l1,1) = sin(zfr)
         ekc(l1,-1) = ekc(l1,1)
         eks(l1,-1) = -eks(l1,1)
         elc(l1,-1) = elc(l1,1)
         els(l1,-1) = -els(l1,1)
         do j = 2, jmax
            ejc(l1,j) = ejc(l1,j-1)*ejc(l1,1) - ejs(l1,j-1)*ejs(l1,1)
            ejs(l1,j) = ejs(l1,j-1)*ejc(l1,1) + ejc(l1,j-1)*ejs(l1,1)
         end do
         do j = 2, kmax
            ekc(l1,j) = ekc(l1,j-1)*ekc(l1,1) - eks(l1,j-1)*eks(l1,1)
            eks(l1,j) = eks(l1,j-1)*ekc(l1,1) + ekc(l1,j-1)*eks(l1,1)
            ekc(l1,-j) = ekc(l1,j)
            eks(l1,-j) = -eks(l1,j)
         end do
         do j = 2, lmax
            elc(l1,j) = elc(l1,j-1)*elc(l1,1) - els(l1,j-1)*els(l1,1)
            els(l1,j) = els(l1,j-1)*elc(l1,1) + elc(l1,j-1)*els(l1,1)
            elc(l1,-j) = elc(l1,j)
            els(l1,-j) = -els(l1,j)
         end do
      end do
c
c     loop over all k vectors from the reciprocal lattice
c
      do j = jmin, jmax
         rj = 2.0d0 * pi * dble(j)
         do k = kmin, kmax
            rk = 2.0d0 * pi * dble(k)
            do l1 = 1, npole3b
               i = pnum(l1)
c               cjk(i) = ejc(i,j)*ekc(i,k) - ejs(i,j)*eks(i,k)
c               sjk(i) = ejs(i,j)*ekc(i,k) + ejc(i,j)*eks(i,k)
               cjk(l1) = ejc(l1,j)*ekc(l1,k) - ejs(l1,j)*eks(l1,k)
               sjk(l1) = ejs(l1,j)*ekc(l1,k) + ejc(l1,j)*eks(l1,k)
            end do
            do l = lmin, lmax
               rl = 2.0d0 * pi * dble(l)
               h1 = recip(1,1)*rj
               h2 = recip(2,1)*rj + recip(2,2)*rk
               h3 = recip(3,1)*rj + recip(3,2)*rk + recip(3,3)*rl
               hsq = h1*h1 + h2*h2 + h3*h3
c               if (hsq .le. cut) then
                  t1 = 0.0d0
                  t2 = 0.0d0
                  t3 = 0.0d0
                  t4 = 0.0d0
                  t3p = 0.0d0
                  t4p = 0.0d0
                  do m2 = 1, 3
                     do m1 = 1, 3
                        t5(m1,m2) = 0.0d0
                        t5u(m1,m2) = 0.0d0
                        t5p(m1,m2) = 0.0d0
                        t6(m1,m2) = 0.0d0
                        t6u(m1,m2) = 0.0d0
                        t6p(m1,m2) = 0.0d0
                     end do
                  end do
                  do l1 = 1, npole3b
                     i = pnum(l1)
                     ckr(l1) = cjk(l1)*elc(l1,l) - sjk(l1)*els(l1,l)
                     skr(l1) = sjk(l1)*elc(l1,l) + cjk(l1)*els(l1,l)
                     ck = cm(l1)
                     dk = h1*dm(1,l1) + h2*dm(2,l1) + h3*dm(3,l1)
                     dkx(l1) = h3*dm(2,l1) - h2*dm(3,l1)
                     dky(l1) = h1*dm(3,l1) - h3*dm(1,l1)
                     dkz(l1) = h2*dm(1,l1) - h1*dm(2,l1)
                     q1 = h1*qm(1,l1) + h2*qm(4,l1) + h3*qm(7,l1)
                     q2 = h1*qm(2,l1) + h2*qm(5,l1) + h3*qm(8,l1)
                     q3 = h1*qm(3,l1) + h2*qm(6,l1) + h3*qm(9,l1)
                     qk = h1*q1 + h2*q2 + h3*q3
                     qkx(l1) = h3*q2 - h2*q3
                     qky(l1) = h1*q3 - h3*q1
                     qkz(l1) = h2*q1 - h1*q2
                     uk = h1*uind(1,l1) + h2*uind(2,l1) + h3*uind(3,l1)
                     ukp = h1*uinp(1,l1) + h2*uinp(2,l1) + h3*uinp(3,l1)
                     s1(l1) = (ck-qk)*skr(l1) + dk*ckr(l1)
                     s2(l1) = (ck-qk)*ckr(l1) - dk*skr(l1)
                     s3(l1) = uk * ckr(l1)
                     s4(l1) = -uk * skr(l1)
                     s3p(l1) = ukp * ckr(l1)
                     s4p(l1) = -ukp * skr(l1)
                     t1 = t1 + s1(l1)
                     t2 = t2 + s2(l1)
                     t3 = t3 + s3(l1)
                     t4 = t4 + s4(l1)
                     t3p = t3p + s3p(l1)
                     t4p = t4p + s4p(l1)
c
c     terms needed for subsequent virial tensor calculation
c
                  qt(1,1)=h1*(h1*qm(1,l1)+h2*qm(4,l1)+h3*qm(7,l1))
                  qt(2,1)=h1*(h1*qm(2,l1)+h2*qm(5,l1)+h3*qm(8,l1))
                  qt(3,1) = h1*(h1*qm(3,l1)+h2*qm(6,l1)+h3*qm(9,l1))
                  qt(1,2) = h2*(h1*qm(1,l1)+h2*qm(4,l1)+h3*qm(7,l1))
                  qt(2,2) = h2*(h1*qm(2,l1)+h2*qm(5,l1)+h3*qm(8,l1))
                  qt(3,2) = h2*(h1*qm(3,l1)+h2*qm(6,l1)+h3*qm(9,l1))
                  qt(1,3) = h3*(h1*qm(1,l1)+h2*qm(4,l1)+h3*qm(7,l1))
                  qt(2,3) = h3*(h1*qm(2,l1)+h2*qm(5,l1)+h3*qm(8,l1))
                  qt(3,3) = h3*(h1*qm(3,l1)+h2*qm(6,l1)+h3*qm(9,l1))
                     dt(1,1) = h1 * dm(1,l1)
                     dt(2,1) = h1 * dm(2,l1)
                     dt(3,1) = h1 * dm(3,l1)
                     dt(1,2) = h2 * dm(1,l1)
                     dt(2,2) = h2 * dm(2,l1)
                     dt(3,2) = h2 * dm(3,l1)
                     dt(1,3) = h3 * dm(1,l1)
                     dt(2,3) = h3 * dm(2,l1)
                     dt(3,3) = h3 * dm(3,l1)
                     dtu(1,1) = h1 * uind(1,l1)
                     dtu(2,1) = h1 * uind(2,l1)
                     dtu(3,1) = h1 * uind(3,l1)
                     dtu(1,2) = h2 * uind(1,l1)
                     dtu(2,2) = h2 * uind(2,l1)
                     dtu(3,2) = h2 * uind(3,l1)
                     dtu(1,3) = h3 * uind(1,l1)
                     dtu(2,3) = h3 * uind(2,l1)
                     dtu(3,3) = h3 * uind(3,l1)
                     dtp(1,1) = h1 * uinp(1,l1)
                     dtp(2,1) = h1 * uinp(2,l1)
                     dtp(3,1) = h1 * uinp(3,l1)
                     dtp(1,2) = h2 * uinp(1,l1)
                     dtp(2,2) = h2 * uinp(2,l1)
                     dtp(3,2) = h2 * uinp(3,l1)
                     dtp(1,3) = h3 * uinp(1,l1)
                     dtp(2,3) = h3 * uinp(2,l1)
                     dtp(3,3) = h3 * uinp(3,l1)
                     do m2 = 1, 3
                        do m1 = 1, 3
                           t5(m1,m2) = t5(m1,m2) - dt(m1,m2)*ckr(l1)
     &                                    + 2.0d0*qt(m1,m2)*skr(l1)
                           t5u(m1,m2) = t5u(m1,m2) - dtu(m1,m2)*ckr(l1)
                           t5p(m1,m2) = t5p(m1,m2) - dtp(m1,m2)*ckr(l1)
                           t6(m1,m2) = t6(m1,m2) + dt(m1,m2)*skr(l1)
     &                                    + 2.0d0*qt(m1,m2)*ckr(l1)
                           t6u(m1,m2) = t6u(m1,m2) + dtu(m1,m2)*skr(l1)
                           t6p(m1,m2) = t6p(m1,m2) + dtp(m1,m2)*skr(l1)
                        end do
                     end do
                  end do
c
c     get the energy contributions for current reciprocal vector
c
                  expterm = eterm * exp(term*hsq) / hsq
                  if (octahedron) then
                     if (mod(j+k+l,2) .ne. 0)  expterm = 0.0d0
                  end if
                  ei = expterm * (t1*t3+t2*t4)
                  etot = e + ei
                  eptemp = eptemp + ei
c
c     get the virial contributions for current reciprocal vector
c 
c     OLD uterm
c                  uterm = expterm * (t1*(t1+t3+t3p) + t3*t3p
c     &                                 + t2*(t2+t4+t4p) + t4*t4p)
c     NEW uterm w/ just the permanent elec part removed (already done)
                  uterm = expterm * (t1*(t3+t3p) + t3*t3p
     &                                 + t2*(t4+t4p) + t4*t4p)
                  do m2 = 1, 3
                     do m1 = 1, 3
c                        wterm(m1,m2) = 2.0d0 * expterm
c     &                     * (t1*t5(m1,m2) + t2*t6(m1,m2)
c     &                        + 0.5d0*(t1*(t5u(m1,m2)+t5p(m1,m2))
c     &                        + t2*(t6u(m1,m2)+t6p(m1,m2))
c     &                        + (t3+t3p)*t5(m1,m2)
c     &                        + t3*t5p(m1,m2) + t3p*t5u(m1,m2)
c     &                        + (t4+t4p)*t6(m1,m2)
c     &                        + t4*t6p(m1,m2) + t4p*t6u(m1,m2)))

c     NEW wterm w/ just the permanent elec part removed (already done)
                        wterm(m1,m2) = 2.0d0 * expterm
     &                     * (
     &                        0.5d0*(t1*(t5u(m1,m2)+t5p(m1,m2))
     &                        + t2*(t6u(m1,m2)+t6p(m1,m2))
     &                        + (t3+t3p)*t5(m1,m2)
     &                        + t3*t5p(m1,m2) + t3p*t5u(m1,m2)
     &                        + (t4+t4p)*t6(m1,m2)
     &                        + t4*t6p(m1,m2) + t4p*t6u(m1,m2)))

                     end do
                  end do
                  if (poltyp .eq. 'DIRECT') then
                     uterm = uterm - expterm*(t3*t3p+t4*t4p)
                     do m2 = 1, 3
                        do m1 = 1, 3
                           wterm(m1,m2) = wterm(m1,m2)
     &                        - expterm*(t3*t5p(m1,m2)+t3p*t5u(m1,m2)
     &                                  +t4*t6p(m1,m2)+t4p*t6u(m1,m2))
                        end do
                     end do
                  end if
                  wterm(2,1) = 0.5d0 * (wterm(2,1)+wterm(1,2))
                  wterm(3,1) = 0.5d0 * (wterm(3,1)+wterm(1,3))
                  wterm(3,2) = 0.5d0 * (wterm(3,2)+wterm(2,3))
                  wterm(1,2) = wterm(2,1)
                  wterm(1,3) = wterm(3,1)
                  wterm(2,3) = wterm(3,2)
                  vterm = 2.0d0 * uterm * (1.0d0-term*hsq) / hsq
c                  vir(1,1) = vir(1,1) + h1*h1*vterm + wterm(1,1) - uterm
c                  vir(2,1) = vir(2,1) + h2*h1*vterm + wterm(2,1)
c                  vir(3,1) = vir(3,1) + h3*h1*vterm + wterm(3,1)
c                  vir(1,2) = vir(1,2) + h1*h2*vterm + wterm(1,2)
c                  vir(2,2) = vir(2,2) + h2*h2*vterm + wterm(2,2) - uterm
c                  vir(3,2) = vir(3,2) + h3*h2*vterm + wterm(3,2)
c                  vir(1,3) = vir(1,3) + h1*h3*vterm + wterm(1,3)
c                  vir(2,3) = vir(2,3) + h2*h3*vterm + wterm(2,3)
c                  vir(3,3) = vir(3,3) + h3*h3*vterm + wterm(3,3) - uterm

                  virtemp(1,1) = virtemp(1,1) + h1*h1*vterm + wterm(1,1)
     &                            - uterm
                  virtemp(2,1) = virtemp(2,1) + h2*h1*vterm + wterm(2,1)
                  virtemp(3,1) = virtemp(3,1) + h3*h1*vterm + wterm(3,1)
                  virtemp(1,2) = virtemp(1,2) + h1*h2*vterm + wterm(1,2)
                  virtemp(2,2) = virtemp(2,2) + h2*h2*vterm + wterm(2,2)
     &                            - uterm
                  virtemp(3,2) = virtemp(3,2) + h3*h2*vterm + wterm(3,2)
                  virtemp(1,3) = virtemp(1,3) + h1*h3*vterm + wterm(1,3)
                  virtemp(2,3) = virtemp(2,3) + h2*h3*vterm + wterm(2,3)
                  virtemp(3,3) = virtemp(3,3) + h3*h3*vterm + wterm(3,3)
     &                             - uterm

c
c     get the force contributions for current reciprocal vector
c
                  expterm = 2.0d0 * expterm
                  do l1 = 1, npole3b
                     i = pnum(l1)
                     ii = ipole(i)
                     dei = 0.5d0 * expterm * ((s4(l1)+s4p(l1))*t1
     &                                       -(s3(l1)+s3p(l1))*t2
     &                                +s2(l1)*(t3+t3p)-s1(l1)*(t4+t4p))
                     if (poltyp .eq. 'MUTUAL') then
                         dei = dei + 0.5d0 * expterm
     &                            * (s4p(l1)*t3+s4(l1)*t3p
     &                              -s3p(l1)*t4-s3(l1)*t4p)
                     end if
c                     det1 = expterm * (skr(l1)*t2-ckr(l1)*t1)
c                     det2 = 2.0d0 * expterm * (ckr(l1)*t2+skr(l1)*t1)
                     det1i = 0.5d0 * expterm * (skr(l1)*(t4+t4p)
     &                                         -ckr(l1)*(t3+t3p))
                     det2i =expterm*(ckr(l1)*(t4+t4p)+skr(l1)*(t3+t3p))
                     deptemp(1,i) = deptemp(1,i) + h1*dei
                     deptemp(2,i) = deptemp(2,i) + h2*dei
                     deptemp(3,i) = deptemp(3,i) + h3*dei
c                     trq(1,l1) = trq(1,l1) + dkx(l1)*det1 + qkx(l1)*det2
c                     trq(2,l1) = trq(2,l1) + dky(l1)*det1 + qky(l1)*det2
c                     trq(3,l1) = trq(3,l1) + dkz(l1)*det1 + qkz(l1)*det2
                     trqi(1,i) = trqi(1,i) + dkx(l1)*det1i
     &                               + qkx(l1)*det2i
                     trqi(2,i) = trqi(2,i) + dky(l1)*det1i
     &                               + qky(l1)*det2i
                     trqi(3,i) = trqi(3,i) + dkz(l1)*det1i
     &                               + qkz(l1)*det2i
                  end do
c               end if
            end do
            lmin = -lmax
         end do
         kmin = -kmax
      end do
c
c     convert the torques to forces and increment the totals
c
c      call torque_reg3b (trqi,deptemp,pnum,pnum2,npole3b)
      call torque_reg3b (trqi,deptemp,pnum,npole3b)
ccc subroutine torque_reg3b (trq,derivs,pnum,pnum2,npole3b)

c      print*,'ep after erecip1_3b=',ep,moli1,moli2,moli3      

      return
      end
c
c     ################################################################
c     ##                                                            ##
c     ##             Subroutine ereal1c_3b                          ##
c     ##  Liam O-Suilleabhain, Omar Demerdash, Teresa Head-Gordon   ##
c     ##                 Spring 2013                                ##
c     ##                                                            ##
c     ################################################################
c
c
c     "ereal1c_3b" evaluates the real space portion of the regular Ewald
c     summation energy and gradient due to dipole polarizability
c
c      subroutine ereal1c_3b(npole3b,pnum,pnum2,uind,uinp,
c     &   eptemp,deptemp,virtemp,moli1,moli2,moli3)
c      subroutine ereal1c_3b(npole3b,pnum,uind,uinp,
c     &   eptemp,deptemp,virtemp,moli1,moli2,moli3)
      subroutine ereal1c_3b(aewald3b,ewaldcut3b,npole3b,pnum,uind,uinp,
     &   eptemp,deptemp,virtemp)
      implicit none
      include 'sizes.i'
      include 'atoms.i'
      include 'bound.i'
      include 'boxes.i'
      include 'cell.i'
      include 'chgpot.i'
      include 'couple.i'
      include 'cutoff.i'
c      include 'deriv.i'
c      include 'energi.i'
      include 'ewald.i'
      include 'math.i'
c      include 'molcul.i'
      include 'mplpot.i'
      include 'mpole.i'
c      include 'polar.i'
      include 'polar2.i'
      include 'polgrp.i'
      include 'polpot.i'
      include 'shunt.i'
c      include 'virial.i'
c      include 'combo.i'
      integer moli1,moli2,moli3
      integer i,j,k,l1,l3,l2,l,ll
      integer ii,kk,jcell,jj
      integer iax,iay,iaz
      integer kax,kay,kaz
      real*8 e,ei,f,bfac
      real*8 eintra,erfc
      real*8 damp,expdamp
      real*8 pdi,pti,pgamma
      real*8 scale3,scale5
      real*8 scale7
      real*8 temp3,temp5,temp7
      real*8 dsc3,dsc5,dsc7
      real*8 psc3,psc5,psc7
      real*8 usc3,usc5
      real*8 alsq2,alsq2n
      real*8 exp2a,ralpha
      real*8 gfd,gfdr
      real*8 xr,yr,zr
      real*8 xix,yix,zix
      real*8 xiy,yiy,ziy
      real*8 xiz,yiz,ziz
      real*8 xkx,ykx,zkx
      real*8 xky,yky,zky
      real*8 xkz,ykz,zkz
      real*8 r,r2,rr1,rr3
      real*8 rr5,rr7,rr9,rr11
      real*8 erl,erli
      real*8 vxx,vyy,vzz
      real*8 vyx,vzx,vzy
      real*8 frcxi(3),frcxk(3)
      real*8 frcyi(3),frcyk(3)
      real*8 frczi(3),frczk(3)
      real*8 ci,di(3),qi(9)
      real*8 ck,dk(3),qk(9)
      real*8 fridmp(3),findmp(3)
      real*8 ftm2(3),ftm2i(3)
      real*8 ftm2r(3),ftm2ri(3)
      real*8 ttm2(3),ttm3(3)
      real*8 ttm2i(3),ttm3i(3)
      real*8 ttm2r(3),ttm3r(3)
      real*8 ttm2ri(3),ttm3ri(3)
      real*8 fdir(3),dixdk(3)
      real*8 dkxui(3),dixuk(3)
      real*8 dixukp(3),dkxuip(3)
      real*8 uixqkr(3),ukxqir(3)
      real*8 uixqkrp(3),ukxqirp(3)
      real*8 qiuk(3),qkui(3)
      real*8 qiukp(3),qkuip(3)
      real*8 rxqiuk(3),rxqkui(3)
      real*8 rxqiukp(3),rxqkuip(3)
      real*8 qidk(3),qkdi(3)
      real*8 qir(3),qkr(3)
      real*8 qiqkr(3),qkqir(3)
      real*8 qixqk(3),rxqir(3)
      real*8 dixr(3),dkxr(3)
      real*8 dixqkr(3),dkxqir(3)
      real*8 rxqkr(3),qkrxqir(3)
      real*8 rxqikr(3),rxqkir(3)
      real*8 rxqidk(3),rxqkdi(3)
      real*8 ddsc3(3),ddsc5(3)
      real*8 ddsc7(3)
      real*8 bn(0:5)
      real*8 sc(10),gl(0:8)
      real*8 sci(8),scip(8)
      real*8 gli(7),glip(7)
      real*8 gf(7),gfi(6)
      real*8 gfr(7),gfri(6)
      real*8 gti(6),gtri(6)
      real*8, allocatable :: mscale(:)
      real*8, allocatable :: pscale(:)
      real*8, allocatable :: dscale(:)
      real*8, allocatable :: uscale(:)
      real*8 uind(3,*),eptemp,deptemp(3,*)
      real*8 virtemp(3,3),uinp(3,*),aewald
      real*8 aewald3b,ewaldcut3b,ewaldcut3b_2
      integer npole3b,pnum(*)
c      integer pnum2(*)
      character*6 mode
      logical flag
      external erfc
c
c
c     zero out the intramolecular portion of the Ewald energy
c
      aewald=aewald3b
      eintra = 0.0d0
      if (npole3b .eq. 0)  return
c
c     perform dynamic allocation of some local arrays
c
      allocate (mscale(n))
      allocate (pscale(n))
      allocate (dscale(n))
      allocate (uscale(n))
c
c     set arrays needed to scale connected atom interactions
c
      do i = 1, n
         mscale(i) = 1.0d0
         pscale(i) = 1.0d0
         dscale(i) = 1.0d0
         uscale(i) = 1.0d0
      end do
c
c     set conversion factor, cutoff and switching coefficients
c
      f = electric / dielec
c      mode = 'EWALD3B'
c      call switch (mode)
      ewaldcut3b_2=ewaldcut3b*ewaldcut3b
c
c     set the permanent multipole and induced dipole values
c

c      do l1 =1, npole3b 
c        pnum2(pnum(l1))=l1
c      end do     

      do l1 = 1, npole3b-1
         i = pnum(l1)
         ii = ipole(i)
         pdi = pdamp(i)
         pti = thole(i)
         ci = rpole(1,i)
         di(1) = rpole(2,i)
         di(2) = rpole(3,i)
         di(3) = rpole(4,i)
         qi(1) = rpole(5,i)
         qi(2) = rpole(6,i)
         qi(3) = rpole(7,i)
         qi(4) = rpole(8,i)
         qi(5) = rpole(9,i)
         qi(6) = rpole(10,i)
         qi(7) = rpole(11,i)
         qi(8) = rpole(12,i)
         qi(9) = rpole(13,i)
c
c     set interaction scaling coefficients for connected atoms
c
         do j = 1, n12(ii)
            mscale(i12(j,ii)) = m2scale
            pscale(i12(j,ii)) = p2scale
         end do
         do j = 1, n13(ii)
            mscale(i13(j,ii)) = m3scale
            pscale(i13(j,ii)) = p3scale
         end do
         do j = 1, n14(ii)
            mscale(i14(j,ii)) = m4scale
            pscale(i14(j,ii)) = p4scale
            do k = 1, np11(ii)
                if (i14(j,ii) .eq. ip11(k,ii))
     &            pscale(i14(j,ii)) = p4scale * p41scale
            end do
         end do
         do j = 1, n15(ii)
            mscale(i15(j,ii)) = m5scale
            pscale(i15(j,ii)) = p5scale
         end do
         do j = 1, np11(ii)
            dscale(ip11(j,ii)) = d1scale
            uscale(ip11(j,ii)) = u1scale
         end do
         do j = 1, np12(ii)
            dscale(ip12(j,ii)) = d2scale
            uscale(ip12(j,ii)) = u2scale
         end do
         do j = 1, np13(ii)
            dscale(ip13(j,ii)) = d3scale
            uscale(ip13(j,ii)) = u3scale
         end do
         do j = 1, np14(ii)
            dscale(ip14(j,ii)) = d4scale
            uscale(ip14(j,ii)) = u4scale
         end do
         do l2 = l1+1, npole3b
            k=pnum(l2)
            kk = ipole(k)
c            flag=.true.
c            do l2=1,npole3b
c              l=pnum(l2)
c              if ((k.eq.l).and.(l1.gt.l2)) then
c                flag=.false.
c                go to 100
c              end if
c            end do
c 100        continue
c            if (flag) then
            xr = x(kk) - x(ii)
            yr = y(kk) - y(ii)
            zr = z(kk) - z(ii)
            call image (xr,yr,zr)
            r2 = xr*xr + yr*yr + zr*zr
c            if (r2 .le. off2) then
            if (r2 .le. ewaldcut3b_2) then
               r = sqrt(r2)
               ck = rpole(1,k)
               dk(1) = rpole(2,k)
               dk(2) = rpole(3,k)
               dk(3) = rpole(4,k)
               qk(1) = rpole(5,k)
               qk(2) = rpole(6,k)
               qk(3) = rpole(7,k)
               qk(4) = rpole(8,k)
               qk(5) = rpole(9,k)
               qk(6) = rpole(10,k)
               qk(7) = rpole(11,k)
               qk(8) = rpole(12,k)
               qk(9) = rpole(13,k)
c
c     calculate the real space error function terms
c
               ralpha = aewald * r
               bn(0) = erfc(ralpha) / r
               alsq2 = 2.0d0 * aewald**2
               alsq2n = 0.0d0
               if (aewald .gt. 0.0d0)  alsq2n = 1.0d0 / (sqrtpi*aewald)
               exp2a = exp(-ralpha**2)
               do j = 1, 5
                  bfac = dble(2*j-1)
                  alsq2n = alsq2 * alsq2n
                  bn(j) = (bfac*bn(j-1)+alsq2n*exp2a) / r2
               end do
c
c     apply Thole polarization damping to scale factors
c
               rr1 = 1.0d0 / r
               rr3 = rr1 / r2
               rr5 = 3.0d0 * rr3 / r2
               rr7 = 5.0d0 * rr5 / r2
               rr9 = 7.0d0 * rr7 / r2
               rr11 = 9.0d0 * rr9 / r2
               scale3 = 1.0d0
               scale5 = 1.0d0
               scale7 = 1.0d0
               do j = 1, 3
                  ddsc3(j) = 0.0d0
                  ddsc5(j) = 0.0d0
                  ddsc7(j) = 0.0d0
               end do
               damp = pdi * pdamp(k)
               if (damp .ne. 0.0d0) then
                  pgamma = min(pti,thole(k))
                  damp = -pgamma * (r/damp)**3
                  if (damp .gt. -50.0d0) then
                     expdamp = exp(damp)
                     scale3 = 1.0d0 - expdamp
                     scale5 = 1.0d0 - (1.0d0-damp)*expdamp
                     scale7 = 1.0d0 - (1.0d0-damp+0.6d0*damp**2)
     &                                       *expdamp
                     temp3 = -3.0d0 * damp * expdamp / r2
                     temp5 = -damp
                     temp7 = -0.2d0 - 0.6d0*damp
                     ddsc3(1) = temp3 * xr
                     ddsc3(2) = temp3 * yr
                     ddsc3(3) = temp3 * zr
                     ddsc5(1) = temp5 * ddsc3(1)
                     ddsc5(2) = temp5 * ddsc3(2)
                     ddsc5(3) = temp5 * ddsc3(3)
                     ddsc7(1) = temp7 * ddsc5(1)
                     ddsc7(2) = temp7 * ddsc5(2)
                     ddsc7(3) = temp7 * ddsc5(3)
                  end if
               end if
               dsc3 = 1.0d0 - scale3*dscale(kk)
               dsc5 = 1.0d0 - scale5*dscale(kk)
               dsc7 = 1.0d0 - scale7*dscale(kk)
               psc3 = 1.0d0 - scale3*pscale(kk)
               psc5 = 1.0d0 - scale5*pscale(kk)
               psc7 = 1.0d0 - scale7*pscale(kk)
               usc3 = 1.0d0 - scale3*uscale(kk)
               usc5 = 1.0d0 - scale5*uscale(kk)
c
c     construct necessary auxiliary vectors
c

c               dixdk(1) = di(2)*dk(3) - di(3)*dk(2)
c               dixdk(2) = di(3)*dk(1) - di(1)*dk(3)
c               dixdk(3) = di(1)*dk(2) - di(2)*dk(1)
               dixuk(1) = di(2)*uind(3,l2) - di(3)*uind(2,l2)
               dixuk(2) = di(3)*uind(1,l2) - di(1)*uind(3,l2)
               dixuk(3) = di(1)*uind(2,l2) - di(2)*uind(1,l2)
               dkxui(1) = dk(2)*uind(3,l1) - dk(3)*uind(2,l1)
               dkxui(2) = dk(3)*uind(1,l1) - dk(1)*uind(3,l1)
               dkxui(3) = dk(1)*uind(2,l1) - dk(2)*uind(1,l1)
               dixukp(1) = di(2)*uinp(3,l2) - di(3)*uinp(2,l2)
               dixukp(2) = di(3)*uinp(1,l2) - di(1)*uinp(3,l2)
               dixukp(3) = di(1)*uinp(2,l2) - di(2)*uinp(1,l2)
               dkxuip(1) = dk(2)*uinp(3,l1) - dk(3)*uinp(2,l1)
               dkxuip(2) = dk(3)*uinp(1,l1) - dk(1)*uinp(3,l1)
               dkxuip(3) = dk(1)*uinp(2,l1) - dk(2)*uinp(1,l1)
               dixr(1) = di(2)*zr - di(3)*yr
               dixr(2) = di(3)*xr - di(1)*zr
               dixr(3) = di(1)*yr - di(2)*xr
               dkxr(1) = dk(2)*zr - dk(3)*yr
               dkxr(2) = dk(3)*xr - dk(1)*zr
               dkxr(3) = dk(1)*yr - dk(2)*xr
               qir(1) = qi(1)*xr + qi(4)*yr + qi(7)*zr
               qir(2) = qi(2)*xr + qi(5)*yr + qi(8)*zr
               qir(3) = qi(3)*xr + qi(6)*yr + qi(9)*zr
               qkr(1) = qk(1)*xr + qk(4)*yr + qk(7)*zr
               qkr(2) = qk(2)*xr + qk(5)*yr + qk(8)*zr
               qkr(3) = qk(3)*xr + qk(6)*yr + qk(9)*zr
c               qiqkr(1) = qi(1)*qkr(1) + qi(4)*qkr(2) + qi(7)*qkr(3)
c               qiqkr(2) = qi(2)*qkr(1) + qi(5)*qkr(2) + qi(8)*qkr(3)
c               qiqkr(3) = qi(3)*qkr(1) + qi(6)*qkr(2) + qi(9)*qkr(3)
c               qkqir(1) = qk(1)*qir(1) + qk(4)*qir(2) + qk(7)*qir(3)
c               qkqir(2) = qk(2)*qir(1) + qk(5)*qir(2) + qk(8)*qir(3)
c               qkqir(3) = qk(3)*qir(1) + qk(6)*qir(2) + qk(9)*qir(3)
c               qixqk(1) = qi(2)*qk(3) + qi(5)*qk(6) + qi(8)*qk(9)
c     &                       - qi(3)*qk(2) - qi(6)*qk(5) - qi(9)*qk(8)
c               qixqk(2) = qi(3)*qk(1) + qi(6)*qk(4) + qi(9)*qk(7)
c     &                       - qi(1)*qk(3) - qi(4)*qk(6) - qi(7)*qk(9)
c               qixqk(3) = qi(1)*qk(2) + qi(4)*qk(5) + qi(7)*qk(8)
c     &                       - qi(2)*qk(1) - qi(5)*qk(4) - qi(8)*qk(7)
               rxqir(1) = yr*qir(3) - zr*qir(2)
               rxqir(2) = zr*qir(1) - xr*qir(3)
               rxqir(3) = xr*qir(2) - yr*qir(1)
               rxqkr(1) = yr*qkr(3) - zr*qkr(2)
               rxqkr(2) = zr*qkr(1) - xr*qkr(3)
               rxqkr(3) = xr*qkr(2) - yr*qkr(1)
c               rxqikr(1) = yr*qiqkr(3) - zr*qiqkr(2)
c               rxqikr(2) = zr*qiqkr(1) - xr*qiqkr(3)
c               rxqikr(3) = xr*qiqkr(2) - yr*qiqkr(1)
c               rxqkir(1) = yr*qkqir(3) - zr*qkqir(2)
c               rxqkir(2) = zr*qkqir(1) - xr*qkqir(3)
c               rxqkir(3) = xr*qkqir(2) - yr*qkqir(1)
c               qkrxqir(1) = qkr(2)*qir(3) - qkr(3)*qir(2)
c               qkrxqir(2) = qkr(3)*qir(1) - qkr(1)*qir(3)
c               qkrxqir(3) = qkr(1)*qir(2) - qkr(2)*qir(1)
c               qidk(1) = qi(1)*dk(1) + qi(4)*dk(2) + qi(7)*dk(3)
c               qidk(2) = qi(2)*dk(1) + qi(5)*dk(2) + qi(8)*dk(3)
c               qidk(3) = qi(3)*dk(1) + qi(6)*dk(2) + qi(9)*dk(3)
c               qkdi(1) = qk(1)*di(1) + qk(4)*di(2) + qk(7)*di(3)
c               qkdi(2) = qk(2)*di(1) + qk(5)*di(2) + qk(8)*di(3)
c               qkdi(3) = qk(3)*di(1) + qk(6)*di(2) + qk(9)*di(3)
               qiuk(1) = qi(1)*uind(1,l2) + qi(4)*uind(2,l2)
     &                      + qi(7)*uind(3,l2)
               qiuk(2) = qi(2)*uind(1,l2) + qi(5)*uind(2,l2)
     &                      + qi(8)*uind(3,l2)
               qiuk(3) = qi(3)*uind(1,l2) + qi(6)*uind(2,l2)
     &                      + qi(9)*uind(3,l2)
               qkui(1) = qk(1)*uind(1,l1) + qk(4)*uind(2,l1)
     &                      + qk(7)*uind(3,l1)
               qkui(2) = qk(2)*uind(1,l1) + qk(5)*uind(2,l1)
     &                      + qk(8)*uind(3,l1)
               qkui(3) = qk(3)*uind(1,l1) + qk(6)*uind(2,l1)
     &                      + qk(9)*uind(3,l1)
               qiukp(1) = qi(1)*uinp(1,l2) + qi(4)*uinp(2,l2)
     &                       + qi(7)*uinp(3,l2)
               qiukp(2) = qi(2)*uinp(1,l2) + qi(5)*uinp(2,l2)
     &                       + qi(8)*uinp(3,l2)
               qiukp(3) = qi(3)*uinp(1,l2) + qi(6)*uinp(2,l2)
     &                       + qi(9)*uinp(3,l2)
               qkuip(1) = qk(1)*uinp(1,l1) + qk(4)*uinp(2,l1)
     &                       + qk(7)*uinp(3,l1)
               qkuip(2) = qk(2)*uinp(1,l1) + qk(5)*uinp(2,l1)
     &                       + qk(8)*uinp(3,l1)
               qkuip(3) = qk(3)*uinp(1,l1) + qk(6)*uinp(2,l1)
     &                       + qk(9)*uinp(3,l1)
c               dixqkr(1) = di(2)*qkr(3) - di(3)*qkr(2)
c               dixqkr(2) = di(3)*qkr(1) - di(1)*qkr(3)
c               dixqkr(3) = di(1)*qkr(2) - di(2)*qkr(1)
c               dkxqir(1) = dk(2)*qir(3) - dk(3)*qir(2)
c               dkxqir(2) = dk(3)*qir(1) - dk(1)*qir(3)
c               dkxqir(3) = dk(1)*qir(2) - dk(2)*qir(1)
               uixqkr(1) = uind(2,l1)*qkr(3) - uind(3,l1)*qkr(2)
               uixqkr(2) = uind(3,l1)*qkr(1) - uind(1,l1)*qkr(3)
               uixqkr(3) = uind(1,l1)*qkr(2) - uind(2,l1)*qkr(1)
               ukxqir(1) = uind(2,l2)*qir(3) - uind(3,l2)*qir(2)
               ukxqir(2) = uind(3,l2)*qir(1) - uind(1,l2)*qir(3)
               ukxqir(3) = uind(1,l2)*qir(2) - uind(2,l2)*qir(1)
               uixqkrp(1) = uinp(2,l1)*qkr(3) - uinp(3,l1)*qkr(2)
               uixqkrp(2) = uinp(3,l1)*qkr(1) - uinp(1,l1)*qkr(3)
               uixqkrp(3) = uinp(1,l1)*qkr(2) - uinp(2,l1)*qkr(1)
               ukxqirp(1) = uinp(2,l2)*qir(3) - uinp(3,l2)*qir(2)
               ukxqirp(2) = uinp(3,l2)*qir(1) - uinp(1,l2)*qir(3)
               ukxqirp(3) = uinp(1,l2)*qir(2) - uinp(2,l2)*qir(1)
c               rxqidk(1) = yr*qidk(3) - zr*qidk(2)
c               rxqidk(2) = zr*qidk(1) - xr*qidk(3)
c               rxqidk(3) = xr*qidk(2) - yr*qidk(1)
c               rxqkdi(1) = yr*qkdi(3) - zr*qkdi(2)
c               rxqkdi(2) = zr*qkdi(1) - xr*qkdi(3)
c               rxqkdi(3) = xr*qkdi(2) - yr*qkdi(1)
               rxqiuk(1) = yr*qiuk(3) - zr*qiuk(2)
               rxqiuk(2) = zr*qiuk(1) - xr*qiuk(3)
               rxqiuk(3) = xr*qiuk(2) - yr*qiuk(1)
               rxqkui(1) = yr*qkui(3) - zr*qkui(2)
               rxqkui(2) = zr*qkui(1) - xr*qkui(3)
               rxqkui(3) = xr*qkui(2) - yr*qkui(1)
               rxqiukp(1) = yr*qiukp(3) - zr*qiukp(2)
               rxqiukp(2) = zr*qiukp(1) - xr*qiukp(3)
               rxqiukp(3) = xr*qiukp(2) - yr*qiukp(1)
               rxqkuip(1) = yr*qkuip(3) - zr*qkuip(2)
               rxqkuip(2) = zr*qkuip(1) - xr*qkuip(3)
               rxqkuip(3) = xr*qkuip(2) - yr*qkuip(1)
c
c     calculate the scalar products for permanent components
c
               sc(2) = di(1)*dk(1) + di(2)*dk(2) + di(3)*dk(3)
               sc(3) = di(1)*xr + di(2)*yr + di(3)*zr
               sc(4) = dk(1)*xr + dk(2)*yr + dk(3)*zr
               sc(5) = qir(1)*xr + qir(2)*yr + qir(3)*zr
               sc(6) = qkr(1)*xr + qkr(2)*yr + qkr(3)*zr
               sc(7) = qir(1)*dk(1) + qir(2)*dk(2) + qir(3)*dk(3)
               sc(8) = qkr(1)*di(1) + qkr(2)*di(2) + qkr(3)*di(3)
               sc(9) = qir(1)*qkr(1) + qir(2)*qkr(2) + qir(3)*qkr(3)
               sc(10) = qi(1)*qk(1) + qi(2)*qk(2) + qi(3)*qk(3)
     &                     + qi(4)*qk(4) + qi(5)*qk(5) + qi(6)*qk(6)
     &                     + qi(7)*qk(7) + qi(8)*qk(8) + qi(9)*qk(9)
c
c     calculate the scalar products for induced components
c
               sci(1) = uind(1,l1)*dk(1) + uind(2,l1)*dk(2)
     &                     + uind(3,l1)*dk(3) + di(1)*uind(1,l2)
     &                     + di(2)*uind(2,l2) + di(3)*uind(3,l2)
               sci(2) = uind(1,l1)*uind(1,l2) + uind(2,l1)*uind(2,l2)
     &                     + uind(3,l1)*uind(3,l2)
               sci(3) = uind(1,l1)*xr + uind(2,l1)*yr + uind(3,l1)*zr
               sci(4) = uind(1,l2)*xr + uind(2,l2)*yr + uind(3,l2)*zr
               sci(7) = qir(1)*uind(1,l2) + qir(2)*uind(2,l2)
     &                     + qir(3)*uind(3,l2)
               sci(8) = qkr(1)*uind(1,l1) + qkr(2)*uind(2,l1)
     &                     + qkr(3)*uind(3,l1)
               scip(1) = uinp(1,l1)*dk(1) + uinp(2,l1)*dk(2)
     &                      + uinp(3,l1)*dk(3) + di(1)*uinp(1,l2)
     &                      + di(2)*uinp(2,l2) + di(3)*uinp(3,l2)
               scip(2) = uind(1,l1)*uinp(1,l2)+uind(2,l1)*uinp(2,l2)
     &                   + uind(3,l1)*uinp(3,l2)+uinp(1,l1)*uind(1,l2)
     &                   + uinp(2,l1)*uind(2,l2)+uinp(3,l1)*uind(3,l2)
               scip(3) = uinp(1,l1)*xr + uinp(2,l1)*yr + uinp(3,l1)*zr
               scip(4) = uinp(1,l2)*xr + uinp(2,l2)*yr + uinp(3,l2)*zr
               scip(7) = qir(1)*uinp(1,l2) + qir(2)*uinp(2,l2)
     &                      + qir(3)*uinp(3,l2)
               scip(8) = qkr(1)*uinp(1,l1) + qkr(2)*uinp(2,l1)
     &                      + qkr(3)*uinp(3,l1)
c
c     calculate the gl functions for permanent components
c

c               gl(0) = ci*ck
c               gl(1) = ck*sc(3) - ci*sc(4)
c               gl(2) = ci*sc(6) + ck*sc(5) - sc(3)*sc(4)
c               gl(3) = sc(3)*sc(6) - sc(4)*sc(5)
c               gl(4) = sc(5)*sc(6)
c               gl(5) = -4.0d0 * sc(9)
c               gl(6) = sc(2)
c               gl(7) = 2.0d0 * (sc(7)-sc(8))
c               gl(8) = 2.0d0 * sc(10)

c
c     calculate the gl functions for induced components
c
               gli(1) = ck*sci(3) - ci*sci(4)
               gli(2) = -sc(3)*sci(4) - sci(3)*sc(4)
               gli(3) = sci(3)*sc(6) - sci(4)*sc(5)
               gli(6) = sci(1)
               gli(7) = 2.0d0 * (sci(7)-sci(8))
               glip(1) = ck*scip(3) - ci*scip(4)
               glip(2) = -sc(3)*scip(4) - scip(3)*sc(4)
               glip(3) = scip(3)*sc(6) - scip(4)*sc(5)
               glip(6) = scip(1)
               glip(7) = 2.0d0 * (scip(7)-scip(8))
c
c     compute the energy contributions for this interaction
c
               ei = 0.5d0 * (bn(1)*(gli(1)+gli(6))
     &                      + bn(2)*(gli(2)+gli(7)) + bn(3)*gli(3))
c
c     get the real energy without any screening function
c
               erli = 0.5d0*(rr3*(gli(1)+gli(6))*psc3
     &                   + rr5*(gli(2)+gli(7))*psc5
     &                   + rr7*gli(3)*psc7)
               ei = ei - erli
               ei = f * ei
               eptemp = eptemp + ei
c
c     increment the total intramolecular energy; assumes
c     intramolecular distances are less than half of cell
c     length and less than the ewald cutoff
c
c               if (molcule(ii) .eq. molcule(kk)) then
c                  eintra = eintra + 0.5d0*pscale(kk)
c     &                        * (rr3*(gli(1)+gli(6))*scale3
c     &                              + rr5*(gli(2)+gli(7))*scale5
c     &                              + rr7*gli(3)*scale7)
c               end if
c
c     intermediate variables for permanent force terms
c

c               gf(1) = bn(1)*gl(0) + bn(2)*(gl(1)+gl(6))
c     &                    + bn(3)*(gl(2)+gl(7)+gl(8))
c     &                    + bn(4)*(gl(3)+gl(5)) + bn(5)*gl(4)
c               gf(2) = -ck*bn(1) + sc(4)*bn(2) - sc(6)*bn(3)
c               gf(3) = ci*bn(1) + sc(3)*bn(2) + sc(5)*bn(3)
c               gf(4) = 2.0d0 * bn(2)
c               gf(5) = 2.0d0 * (-ck*bn(2)+sc(4)*bn(3)-sc(6)*bn(4))
c               gf(6) = 2.0d0 * (-ci*bn(2)-sc(3)*bn(3)-sc(5)*bn(4))
c               gf(7) = 4.0d0 * bn(3)
c               gfr(1) = rr3*gl(0) + rr5*(gl(1)+gl(6))
c     &                     + rr7*(gl(2)+gl(7)+gl(8))
c     &                     + rr9*(gl(3)+gl(5)) + rr11*gl(4)
c               gfr(2) = -ck*rr3 + sc(4)*rr5 - sc(6)*rr7
c               gfr(3) = ci*rr3 + sc(3)*rr5 + sc(5)*rr7
c               gfr(4) = 2.0d0 * rr5
c               gfr(5) = 2.0d0 * (-ck*rr5+sc(4)*rr7-sc(6)*rr9)
c               gfr(6) = 2.0d0 * (-ci*rr5-sc(3)*rr7-sc(5)*rr9)
c               gfr(7) = 4.0d0 * rr7

c
c     intermediate variables for induced force terms
c
               gfi(1) = 0.5d0*bn(2)*(gli(1)+glip(1)+gli(6)+glip(6))
     &                     + 0.5d0*bn(2)*scip(2)
     &                     + 0.5d0*bn(3)*(gli(2)+glip(2)+gli(7)+glip(7))
     &                     - 0.5d0*bn(3)*(sci(3)*scip(4)+scip(3)*sci(4))
     &                     + 0.5d0*bn(4)*(gli(3)+glip(3))
               gfi(2) = -ck*bn(1) + sc(4)*bn(2) - sc(6)*bn(3)
               gfi(3) = ci*bn(1) + sc(3)*bn(2) + sc(5)*bn(3)
               gfi(4) = 2.0d0 * bn(2)
               gfi(5) = bn(3) * (sci(4)+scip(4))
               gfi(6) = -bn(3) * (sci(3)+scip(3))
               gfri(1) = 0.5d0*rr5*((gli(1)+gli(6))*psc3
     &                            + (glip(1)+glip(6))*dsc3
     &                            + scip(2)*usc3)
     &                 + 0.5d0*rr7*((gli(7)+gli(2))*psc5
     &                            + (glip(7)+glip(2))*dsc5
     &                     - (sci(3)*scip(4)+scip(3)*sci(4))*usc5)
     &                 + 0.5d0*rr9*(gli(3)*psc7+glip(3)*dsc7)
               gfri(2) = -rr3*ck + rr5*sc(4) - rr7*sc(6)
               gfri(3) = rr3*ci + rr5*sc(3) + rr7*sc(5)
               gfri(4) = 2.0d0 * rr5
               gfri(5) = rr7 * (sci(4)*psc7+scip(4)*dsc7)
               gfri(6) = -rr7 * (sci(3)*psc7+scip(3)*dsc7)
c
c     get the induced force with screening
c
               ftm2i(1) = gfi(1)*xr + 0.5d0*
     &             (gfi(2)*(uind(1,l1)+uinp(1,l1))
     &            + bn(2)*(sci(4)*uinp(1,l1)+scip(4)*uind(1,l1))
     &            + gfi(3)*(uind(1,l2)+uinp(1,l2))
     &            + bn(2)*(sci(3)*uinp(1,l2)+scip(3)*uind(1,l2))
     &            + (sci(4)+scip(4))*bn(2)*di(1)
     &            + (sci(3)+scip(3))*bn(2)*dk(1)
     &            + gfi(4)*(qkui(1)+qkuip(1)-qiuk(1)-qiukp(1)))
     &            + gfi(5)*qir(1) + gfi(6)*qkr(1)
               ftm2i(2) = gfi(1)*yr + 0.5d0*
     &             (gfi(2)*(uind(2,l1)+uinp(2,l1))
     &            + bn(2)*(sci(4)*uinp(2,l1)+scip(4)*uind(2,l1))
     &            + gfi(3)*(uind(2,l2)+uinp(2,l2))
     &            + bn(2)*(sci(3)*uinp(2,l2)+scip(3)*uind(2,l2))
     &            + (sci(4)+scip(4))*bn(2)*di(2)
     &            + (sci(3)+scip(3))*bn(2)*dk(2)
     &            + gfi(4)*(qkui(2)+qkuip(2)-qiuk(2)-qiukp(2)))
     &            + gfi(5)*qir(2) + gfi(6)*qkr(2)
               ftm2i(3) = gfi(1)*zr + 0.5d0*
     &             (gfi(2)*(uind(3,l1)+uinp(3,l1))
     &            + bn(2)*(sci(4)*uinp(3,l1)+scip(4)*uind(3,l1))
     &            + gfi(3)*(uind(3,l2)+uinp(3,l2))
     &            + bn(2)*(sci(3)*uinp(3,l2)+scip(3)*uind(3,l2))
     &            + (sci(4)+scip(4))*bn(2)*di(3)
     &            + (sci(3)+scip(3))*bn(2)*dk(3)
     &            + gfi(4)*(qkui(3)+qkuip(3)-qiuk(3)-qiukp(3)))
     &            + gfi(5)*qir(3) + gfi(6)*qkr(3)
c
c     get the induced force without screening
c
               ftm2ri(1) = gfri(1)*xr + 0.5d0*
     &           (- rr3*ck*(uind(1,l1)*psc3+uinp(1,l1)*dsc3)
     &            + rr5*sc(4)*(uind(1,l1)*psc5+uinp(1,l1)*dsc5)
     &            - rr7*sc(6)*(uind(1,l1)*psc7+uinp(1,l1)*dsc7))
     &            + (rr3*ci*(uind(1,l2)*psc3+uinp(1,l2)*dsc3)
     &            + rr5*sc(3)*(uind(1,l2)*psc5+uinp(1,l2)*dsc5)
     &            + rr7*sc(5)*(uind(1,l2)*psc7+uinp(1,l2)*dsc7))*0.5d0
     &            + rr5*usc5*(sci(4)*uinp(1,l1)+scip(4)*uind(1,l1)
     &            + sci(3)*uinp(1,l2)+scip(3)*uind(1,l2))*0.5d0
     &            + 0.5d0*(sci(4)*psc5+scip(4)*dsc5)*rr5*di(1)
     &            + 0.5d0*(sci(3)*psc5+scip(3)*dsc5)*rr5*dk(1)
     &            + 0.5d0*gfri(4)*((qkui(1)-qiuk(1))*psc5
     &            + (qkuip(1)-qiukp(1))*dsc5)
     &            + gfri(5)*qir(1) + gfri(6)*qkr(1)
               ftm2ri(2) = gfri(1)*yr + 0.5d0*
     &           (- rr3*ck*(uind(2,l1)*psc3+uinp(2,l1)*dsc3)
     &            + rr5*sc(4)*(uind(2,l1)*psc5+uinp(2,l1)*dsc5)
     &            - rr7*sc(6)*(uind(2,l1)*psc7+uinp(2,l1)*dsc7))
     &            + (rr3*ci*(uind(2,l2)*psc3+uinp(2,l2)*dsc3)
     &            + rr5*sc(3)*(uind(2,l2)*psc5+uinp(2,l2)*dsc5)
     &            + rr7*sc(5)*(uind(2,l2)*psc7+uinp(2,l2)*dsc7))*0.5d0
     &            + rr5*usc5*(sci(4)*uinp(2,l1)+scip(4)*uind(2,l1)
     &            + sci(3)*uinp(2,l2)+scip(3)*uind(2,l2))*0.5d0
     &            + 0.5d0*(sci(4)*psc5+scip(4)*dsc5)*rr5*di(2)
     &            + 0.5d0*(sci(3)*psc5+scip(3)*dsc5)*rr5*dk(2)
     &            + 0.5d0*gfri(4)*((qkui(2)-qiuk(2))*psc5
     &            + (qkuip(2)-qiukp(2))*dsc5)
     &            + gfri(5)*qir(2) + gfri(6)*qkr(2)
               ftm2ri(3) = gfri(1)*zr + 0.5d0*
     &           (- rr3*ck*(uind(3,l1)*psc3+uinp(3,l1)*dsc3)
     &            + rr5*sc(4)*(uind(3,l1)*psc5+uinp(3,l1)*dsc5)
     &            - rr7*sc(6)*(uind(3,l1)*psc7+uinp(3,l1)*dsc7))
     &            + (rr3*ci*(uind(3,l2)*psc3+uinp(3,l2)*dsc3)
     &            + rr5*sc(3)*(uind(3,l2)*psc5+uinp(3,l2)*dsc5)
     &            + rr7*sc(5)*(uind(3,l2)*psc7+uinp(3,l2)*dsc7))*0.5d0
     &            + rr5*usc5*(sci(4)*uinp(3,l1)+scip(4)*uind(3,l1)
     &            + sci(3)*uinp(3,l2)+scip(3)*uind(3,l2))*0.5d0
     &            + 0.5d0*(sci(4)*psc5+scip(4)*dsc5)*rr5*di(3)
     &            + 0.5d0*(sci(3)*psc5+scip(3)*dsc5)*rr5*dk(3)
     &            + 0.5d0*gfri(4)*((qkui(3)-qiuk(3))*psc5
     &            + (qkuip(3)-qiukp(3))*dsc5)
     &            + gfri(5)*qir(3) + gfri(6)*qkr(3)
c
c     account for partially excluded induced interactions
c
               temp3 = 0.5d0 * rr3 * ((gli(1)+gli(6))*pscale(kk)
     &                                  +(glip(1)+glip(6))*dscale(kk))
               temp5 = 0.5d0 * rr5 * ((gli(2)+gli(7))*pscale(kk)
     &                                  +(glip(2)+glip(7))*dscale(kk))
               temp7 = 0.5d0 * rr7 * (gli(3)*pscale(kk)
     &                                  +glip(3)*dscale(kk))
               fridmp(1) = temp3*ddsc3(1) + temp5*ddsc5(1)
     &                        + temp7*ddsc7(1)
               fridmp(2) = temp3*ddsc3(2) + temp5*ddsc5(2)
     &                        + temp7*ddsc7(2)
               fridmp(3) = temp3*ddsc3(3) + temp5*ddsc5(3)
     &                        + temp7*ddsc7(3)
c
c     find some scaling terms for induced-induced force
c
               temp3 = 0.5d0 * rr3 * uscale(kk) * scip(2)
               temp5 = -0.5d0 * rr5 * uscale(kk)
     &                    * (sci(3)*scip(4)+scip(3)*sci(4))
               findmp(1) = temp3*ddsc3(1) + temp5*ddsc5(1)
               findmp(2) = temp3*ddsc3(2) + temp5*ddsc5(2)
               findmp(3) = temp3*ddsc3(3) + temp5*ddsc5(3)
c
c     modify the forces for partially excluded interactions
c
               ftm2i(1) = ftm2i(1) - fridmp(1) - findmp(1)
               ftm2i(2) = ftm2i(2) - fridmp(2) - findmp(2)
               ftm2i(3) = ftm2i(3) - fridmp(3) - findmp(3)
c
c     correction to convert mutual to direct polarization force
c

c               if (poltyp .eq. 'DIRECT') then
c                  gfd = 0.5d0 * (bn(2)*scip(2)
c     &                     - bn(3)*(scip(3)*sci(4)+sci(3)*scip(4)))
c                  gfdr = 0.5d0 * (rr5*scip(2)*usc3
c     &                     - rr7*(scip(3)*sci(4)
c     &                           +sci(3)*scip(4))*usc5)
c
c                  ftm2i(1) = ftm2i(1) - gfd*xr - 0.5d0*bn(2)*
c     &                          (sci(4)*uinp(1,l1)+scip(4)*uind(1,l1)
c     &                          +sci(3)*uinp(1,l2)+scip(3)*uind(1,l2))
c                  ftm2i(2) = ftm2i(2) - gfd*yr - 0.5d0*bn(2)*
c     &                          (sci(4)*uinp(2,l1)+scip(4)*uind(2,l1)
c     &                          +sci(3)*uinp(2,l2)+scip(3)*uind(2,l2))
c                  ftm2i(3) = ftm2i(3) - gfd*zr - 0.5d0*bn(2)*
c     &                          (sci(4)*uinp(3,l1)+scip(4)*uind(3,l1)
c     &                          +sci(3)*uinp(3,l2)+scip(3)*uind(3,l2))
c                  fdir(1) = gfdr*xr + 0.5d0*usc5*rr5*
c     &                         (sci(4)*uinp(1,l1)+scip(4)*uind(1,l1)
c     &                        + sci(3)*uinp(1,l2)+scip(3)*uind(1,l2))
c                  fdir(2) = gfdr*yr + 0.5d0*usc5*rr5*
c     &                         (sci(4)*uinp(2,l1)+scip(4)*uind(2,l1)
c     &                        + sci(3)*uinp(2,l2)+scip(3)*uind(2,l2))
c                  fdir(3) = gfdr*zr + 0.5d0*usc5*rr5*
c     &                         (sci(4)*uinp(3,l1)+scip(4)*uind(3,l1)
c     &                        + sci(3)*uinp(3,l2)+scip(3)*uind(3,l2))
c                  ftm2i(1) = ftm2i(1) + fdir(1) + findmp(1)
c                  ftm2i(2) = ftm2i(2) + fdir(2) + findmp(2)
c                  ftm2i(3) = ftm2i(3) + fdir(3) + findmp(3)
c               end if

c
c     intermediate variables for induced torque terms
c
               gti(2) = 0.5d0 * bn(2) * (sci(4)+scip(4))
               gti(3) = 0.5d0 * bn(2) * (sci(3)+scip(3))
               gti(4) = gfi(4)
               gti(5) = gfi(5)
               gti(6) = gfi(6)
               gtri(2) = 0.5d0 * rr5 * (sci(4)*psc5+scip(4)*dsc5)
               gtri(3) = 0.5d0 * rr5 * (sci(3)*psc5+scip(3)*dsc5)
               gtri(4) = gfri(4)
               gtri(5) = gfri(5)
               gtri(6) = gfri(6)
c
c     get the induced torque with screening
c
               ttm2i(1) = -bn(1)*(dixuk(1)+dixukp(1))*0.5d0
     &           + gti(2)*dixr(1) + gti(4)*(ukxqir(1)+rxqiuk(1)
     &           + ukxqirp(1)+rxqiukp(1))*0.5d0 - gti(5)*rxqir(1)
               ttm2i(2) = -bn(1)*(dixuk(2)+dixukp(2))*0.5d0
     &           + gti(2)*dixr(2) + gti(4)*(ukxqir(2)+rxqiuk(2)
     &           + ukxqirp(2)+rxqiukp(2))*0.5d0 - gti(5)*rxqir(2)
               ttm2i(3) = -bn(1)*(dixuk(3)+dixukp(3))*0.5d0
     &           + gti(2)*dixr(3) + gti(4)*(ukxqir(3)+rxqiuk(3)
     &           + ukxqirp(3)+rxqiukp(3))*0.5d0 - gti(5)*rxqir(3)
               ttm3i(1) = -bn(1)*(dkxui(1)+dkxuip(1))*0.5d0
     &           + gti(3)*dkxr(1) - gti(4)*(uixqkr(1)+rxqkui(1)
     &           + uixqkrp(1)+rxqkuip(1))*0.5d0 - gti(6)*rxqkr(1)
               ttm3i(2) = -bn(1)*(dkxui(2)+dkxuip(2))*0.5d0
     &           + gti(3)*dkxr(2) - gti(4)*(uixqkr(2)+rxqkui(2)
     &           + uixqkrp(2)+rxqkuip(2))*0.5d0 - gti(6)*rxqkr(2)
               ttm3i(3) = -bn(1)*(dkxui(3)+dkxuip(3))*0.5d0
     &           + gti(3)*dkxr(3) - gti(4)*(uixqkr(3)+rxqkui(3)
     &           + uixqkrp(3)+rxqkuip(3))*0.5d0 - gti(6)*rxqkr(3)
c
c     get the induced torque without screening
c
               ttm2ri(1) = -rr3*(dixuk(1)*psc3+dixukp(1)*dsc3)*0.5d0
     &           + gtri(2)*dixr(1) + gtri(4)*((ukxqir(1)+rxqiuk(1))*psc5
     &           +(ukxqirp(1)+rxqiukp(1))*dsc5)*0.5d0 - gtri(5)*rxqir(1)
               ttm2ri(2) = -rr3*(dixuk(2)*psc3+dixukp(2)*dsc3)*0.5d0
     &           + gtri(2)*dixr(2) + gtri(4)*((ukxqir(2)+rxqiuk(2))*psc5
     &           +(ukxqirp(2)+rxqiukp(2))*dsc5)*0.5d0 - gtri(5)*rxqir(2)
               ttm2ri(3) = -rr3*(dixuk(3)*psc3+dixukp(3)*dsc3)*0.5d0
     &           + gtri(2)*dixr(3) + gtri(4)*((ukxqir(3)+rxqiuk(3))*psc5
     &           +(ukxqirp(3)+rxqiukp(3))*dsc5)*0.5d0 - gtri(5)*rxqir(3)
               ttm3ri(1) = -rr3*(dkxui(1)*psc3+dkxuip(1)*dsc3)*0.5d0
     &           + gtri(3)*dkxr(1) - gtri(4)*((uixqkr(1)+rxqkui(1))*psc5
     &           +(uixqkrp(1)+rxqkuip(1))*dsc5)*0.5d0 - gtri(6)*rxqkr(1)
               ttm3ri(2) = -rr3*(dkxui(2)*psc3+dkxuip(2)*dsc3)*0.5d0
     &           + gtri(3)*dkxr(2) - gtri(4)*((uixqkr(2)+rxqkui(2))*psc5
     &           +(uixqkrp(2)+rxqkuip(2))*dsc5)*0.5d0 - gtri(6)*rxqkr(2)
               ttm3ri(3) = -rr3*(dkxui(3)*psc3+dkxuip(3)*dsc3)*0.5d0
     &           + gtri(3)*dkxr(3) - gtri(4)*((uixqkr(3)+rxqkui(3))*psc5
     &           +(uixqkrp(3)+rxqkuip(3))*dsc5)*0.5d0 - gtri(6)*rxqkr(3)
c
c     handle the case where scaling is used
c
               do j = 1, 3
                  ftm2i(j) = f * (ftm2i(j)-ftm2ri(j))
                  ttm2i(j) = f * (ttm2i(j)-ttm2ri(j))
                  ttm3i(j) = f * (ttm3i(j)-ttm3ri(j))
               end do
c
c     increment gradient due to force and torque on first site
c
c               deptemp(1,l1) = deptemp(1,l1) + ftm2i(1)
c               deptemp(2,l1) = deptemp(2,l1) + ftm2i(2)
c               deptemp(3,l1) = deptemp(3,l1) + ftm2i(3)
c               call torque_3b (pnum2,deptemp,i,
c     &          ttm2,ttm2i,frcxi,frcyi,frczi)

               deptemp(1,i) = deptemp(1,i) + ftm2i(1)
               deptemp(2,i) = deptemp(2,i) + ftm2i(2)
               deptemp(3,i) = deptemp(3,i) + ftm2i(3)
               call torque_3b (deptemp,i,
     &          ttm2,ttm2i,frcxi,frcyi,frczi)

c
c     increment gradient due to force and torque on second site
c
c               deptemp(1,l2) = deptemp(1,l2) - ftm2i(1)
c               deptemp(2,l2) = deptemp(2,l2) - ftm2i(2)
c               deptemp(3,l2) = deptemp(3,l2) - ftm2i(3)
c               call torque_3b (pnum2,deptemp,k,
c     &            ttm3,ttm3i,frcxk,frcyk,frczk)

               deptemp(1,k) = deptemp(1,k) - ftm2i(1)
               deptemp(2,k) = deptemp(2,k) - ftm2i(2)
               deptemp(3,k) = deptemp(3,k) - ftm2i(3)
               call torque_3b (deptemp,k,
     &            ttm3,ttm3i,frcxk,frcyk,frczk)

c      subroutine torque_3b(pnum2,deptemp,
c     &    i,trq1,trq2,frcx,frcy,frcz)

c
c     increment the internal virial tensor components
c
               iaz = zaxis(i)
               iax = xaxis(i)
               iay = yaxis(i)
               kaz = zaxis(k)
               kax = xaxis(k)
               kay = yaxis(k)
               if (iaz .eq. 0)  iaz = ii
               if (iax .eq. 0)  iax = ii
               if (iay .eq. 0)  iay = ii
               if (kaz .eq. 0)  kaz = kk
               if (kax .eq. 0)  kax = kk
               if (kay .eq. 0)  kay = kk
               xiz = x(iaz) - x(ii)
               yiz = y(iaz) - y(ii)
               ziz = z(iaz) - z(ii)
               xix = x(iax) - x(ii)
               yix = y(iax) - y(ii)
               zix = z(iax) - z(ii)
               xiy = x(iay) - x(ii)
               yiy = y(iay) - y(ii)
               ziy = z(iay) - z(ii)
               xkz = x(kaz) - x(kk)
               ykz = y(kaz) - y(kk)
               zkz = z(kaz) - z(kk)
               xkx = x(kax) - x(kk)
               ykx = y(kax) - y(kk)
               zkx = z(kax) - z(kk)
               xky = x(kay) - x(kk)
               yky = y(kay) - y(kk)
               zky = z(kay) - z(kk)
c               vxx = -xr*ftm2i(1)+xix*frcxi(1)+xiy*frcyi(1)+xiz*frczi(1)
c               vyx = -yr*ftm2i(1)+yix*frcxi(1)+yiy*frcyi(1)+yiz*frczi(1)
c               vzx = -zr*ftm2i(1)+zix*frcxi(1)+ziy*frcyi(1)+ziz*frczi(1)
c               vyy = -yr*ftm2i(2)+yix*frcxi(2)+yiy*frcyi(2)+yiz*frczi(2)
c               vzy = -zr*ftm2i(2)+zix*frcxi(2)+ziy*frcyi(2)+ziz*frczi(2)
c               vzz = -zr*ftm2i(3)+zix*frcxi(3)+ziy*frcyi(3)+ziz*frczi(3)
c               vir(1,1) = vir(1,1) + vxx
c               vir(2,1) = vir(2,1) + vyx
c               vir(3,1) = vir(3,1) + vzx
c               vir(1,2) = vir(1,2) + vyx
c               vir(2,2) = vir(2,2) + vyy
c               vir(3,2) = vir(3,2) + vzy
c               vir(1,3) = vir(1,3) + vzx
c               vir(2,3) = vir(2,3) + vzy
c               vir(3,3) = vir(3,3) + vzz
               vxx = -xr*(ftm2i(1)) + xix*frcxi(1)
     &                  + xiy*frcyi(1) + xiz*frczi(1) + xkx*frcxk(1)
     &                  + xky*frcyk(1) + xkz*frczk(1)
               vyx = -yr*(ftm2i(1)) + yix*frcxi(1)
     &                  + yiy*frcyi(1) + yiz*frczi(1) + ykx*frcxk(1)
     &                  + yky*frcyk(1) + ykz*frczk(1)
               vzx = -zr*(ftm2i(1)) + zix*frcxi(1)
     &                  + ziy*frcyi(1) + ziz*frczi(1) + zkx*frcxk(1)
     &                  + zky*frcyk(1) + zkz*frczk(1)
               vyy = -yr*(ftm2i(2)) + yix*frcxi(2)
     &                  + yiy*frcyi(2) + yiz*frczi(2) + ykx*frcxk(2)
     &                  + yky*frcyk(2) + ykz*frczk(2)
               vzy = -zr*(ftm2i(2)) + zix*frcxi(2)
     &                  + ziy*frcyi(2) + ziz*frczi(2) + zkx*frcxk(2)
     &                  + zky*frcyk(2) + zkz*frczk(2)
               vzz = -zr*(ftm2i(3)) + zix*frcxi(3)
     &                  + ziy*frcyi(3) + ziz*frczi(3) + zkx*frcxk(3)
     &                  + zky*frcyk(3) + zkz*frczk(3)

               virtemp(1,1) = virtemp(1,1) + vxx
               virtemp(2,1) = virtemp(2,1) + vyx
               virtemp(3,1) = virtemp(3,1) + vzx
               virtemp(1,2) = virtemp(1,2) + vyx
               virtemp(2,2) = virtemp(2,2) + vyy
               virtemp(3,2) = virtemp(3,2) + vzy
               virtemp(1,3) = virtemp(1,3) + vzx
               virtemp(2,3) = virtemp(2,3) + vzy
               virtemp(3,3) = virtemp(3,3) + vzz

            end if
c          end if
         end do
c
c     reset interaction scaling coefficients for connected atoms
c
         do j = 1, n12(ii)
            mscale(i12(j,ii)) = 1.0d0
            pscale(i12(j,ii)) = 1.0d0
         end do
         do j = 1, n13(ii)
            mscale(i13(j,ii)) = 1.0d0
            pscale(i13(j,ii)) = 1.0d0
         end do
         do j = 1, n14(ii)
            mscale(i14(j,ii)) = 1.0d0
            pscale(i14(j,ii)) = 1.0d0
         end do
         do j = 1, n15(ii)
            mscale(i15(j,ii)) = 1.0d0
            pscale(i15(j,ii)) = 1.0d0
         end do
         do j = 1, np11(ii)
            dscale(ip11(j,ii)) = 1.0d0
            uscale(ip11(j,ii)) = 1.0d0
         end do
         do j = 1, np12(ii)
            dscale(ip12(j,ii)) = 1.0d0
            uscale(ip12(j,ii)) = 1.0d0
         end do
         do j = 1, np13(ii)
            dscale(ip13(j,ii)) = 1.0d0
            uscale(ip13(j,ii)) = 1.0d0
         end do
         do j = 1, np14(ii)
            dscale(ip14(j,ii)) = 1.0d0
            uscale(ip14(j,ii)) = 1.0d0
         end do
      end do
c
c     dor periodic boundary conditions with large cutoffs
c
c      print*,'ep b4 use_replica ereal1c_3b=',ep,moli1,moli2,moli3

      if (use_replica) then
c
c     calculate interactions with other unit cells
c
c      print*,'After use_replica!',moli1,moli2,moli3
c      do l1 = 1, npole3b
c         print*,"After use_replica l1 Pnum(l1)",l1,pnum(l1)
c      end do 
      do l1 = 1, npole3b
         i = pnum(l1)
         ii = ipole(i)
         pdi = pdamp(i)
         pti = thole(i)
         ci = rpole(1,i)
         di(1) = rpole(2,i)
         di(2) = rpole(3,i)
         di(3) = rpole(4,i)
         qi(1) = rpole(5,i)
         qi(2) = rpole(6,i)
         qi(3) = rpole(7,i)
         qi(4) = rpole(8,i)
         qi(5) = rpole(9,i)
         qi(6) = rpole(10,i)
         qi(7) = rpole(11,i)
         qi(8) = rpole(12,i)
         qi(9) = rpole(13,i)
c
c     set interaction scaling coefficients for connected atoms
c
         do j = 1, n12(ii)
            mscale(i12(j,ii)) = m2scale
            pscale(i12(j,ii)) = p2scale
         end do
         do j = 1, n13(ii)
            mscale(i13(j,ii)) = m3scale
            pscale(i13(j,ii)) = p3scale
         end do
         do j = 1, n14(ii)
            mscale(i14(j,ii)) = m4scale
            pscale(i14(j,ii)) = p4scale
            do k = 1, np11(ii)
                if (i14(j,ii) .eq. ip11(k,ii))
     &            pscale(i14(j,ii)) = p4scale * p41scale
            end do
         end do
         do j = 1, n15(ii)
            mscale(i15(j,ii)) = m5scale
            pscale(i15(j,ii)) = p5scale
         end do
         do j = 1, np11(ii)
            dscale(ip11(j,ii)) = d1scale
            uscale(ip11(j,ii)) = u1scale
         end do
         do j = 1, np12(ii)
            dscale(ip12(j,ii)) = d2scale
            uscale(ip12(j,ii)) = u2scale
         end do
         do j = 1, np13(ii)
            dscale(ip13(j,ii)) = d3scale
            uscale(ip13(j,ii)) = u3scale
         end do
         do j = 1, np14(ii)
            dscale(ip14(j,ii)) = d4scale
            uscale(ip14(j,ii)) = u4scale
         end do
         do l2 = l1, npole3b
            k=pnum(l2)
            kk = ipole(k)
c            flag=.true.
c            do l2=1,npole3b
c              l=pnum(l2)
c              if (k.eq.l.and.i.lt.l) then
c                flag=.false.
c                go to 200
c              end if
c            end do
c 200        continue
c            if (flag) then
            do jcell=1,ncell
            xr = x(kk) - x(ii)
            yr = y(kk) - y(ii)
            zr = z(kk) - z(ii)
            call imager (xr,yr,zr,jcell)
            r2 = xr*xr + yr*yr + zr*zr
            if (.not.(use_polymer .and. r2.le.polycut2)) then
              mscale(kk)=1.0d0
              pscale(kk)=1.0d0
              dscale(kk)=1.0d0
              uscale(kk)=1.0d0
            end if
c            if (r2 .le. off2) then
            if (r2 .le. ewaldcut3b_2) then
               r = sqrt(r2)
               ck = rpole(1,k)
               dk(1) = rpole(2,k)
               dk(2) = rpole(3,k)
               dk(3) = rpole(4,k)
               qk(1) = rpole(5,k)
               qk(2) = rpole(6,k)
               qk(3) = rpole(7,k)
               qk(4) = rpole(8,k)
               qk(5) = rpole(9,k)
               qk(6) = rpole(10,k)
               qk(7) = rpole(11,k)
               qk(8) = rpole(12,k)
               qk(9) = rpole(13,k)
c
c     calculate the real space error function terms
c
               ralpha = aewald * r
               bn(0) = erfc(ralpha) / r
               alsq2 = 2.0d0 * aewald**2
               alsq2n = 0.0d0
               if (aewald .gt. 0.0d0)  alsq2n = 1.0d0 / (sqrtpi*aewald)
               exp2a = exp(-ralpha**2)
               do j = 1, 5
                  bfac = dble(2*j-1)
                  alsq2n = alsq2 * alsq2n
                  bn(j) = (bfac*bn(j-1)+alsq2n*exp2a) / r2
               end do
c
c     apply Thole polarization damping to scale factors
c
               rr1 = 1.0d0 / r
               rr3 = rr1 / r2
               rr5 = 3.0d0 * rr3 / r2
               rr7 = 5.0d0 * rr5 / r2
               rr9 = 7.0d0 * rr7 / r2
               rr11 = 9.0d0 * rr9 / r2
               scale3 = 1.0d0
               scale5 = 1.0d0
               scale7 = 1.0d0
               do j = 1, 3
                  ddsc3(j) = 0.0d0
                  ddsc5(j) = 0.0d0
                  ddsc7(j) = 0.0d0
               end do
               damp = pdi * pdamp(k)
               if (damp .ne. 0.0d0) then
                  pgamma = min(pti,thole(k))
                  damp = -pgamma * (r/damp)**3
                  if (damp .gt. -50.0d0) then
                     expdamp = exp(damp)
                     scale3 = 1.0d0 - expdamp
                     scale5 = 1.0d0 - (1.0d0-damp)*expdamp
                     scale7 = 1.0d0 - (1.0d0-damp+0.6d0*damp**2)
     &                                       *expdamp
                     temp3 = -3.0d0 * damp * expdamp / r2
                     temp5 = -damp
                     temp7 = -0.2d0 - 0.6d0*damp
                     ddsc3(1) = temp3 * xr
                     ddsc3(2) = temp3 * yr
                     ddsc3(3) = temp3 * zr
                     ddsc5(1) = temp5 * ddsc3(1)
                     ddsc5(2) = temp5 * ddsc3(2)
                     ddsc5(3) = temp5 * ddsc3(3)
                     ddsc7(1) = temp7 * ddsc5(1)
                     ddsc7(2) = temp7 * ddsc5(2)
                     ddsc7(3) = temp7 * ddsc5(3)
                  end if
               end if
               dsc3 = 1.0d0 - scale3*dscale(kk)
               dsc5 = 1.0d0 - scale5*dscale(kk)
               dsc7 = 1.0d0 - scale7*dscale(kk)
               psc3 = 1.0d0 - scale3*pscale(kk)
               psc5 = 1.0d0 - scale5*pscale(kk)
               psc7 = 1.0d0 - scale7*pscale(kk)
               usc3 = 1.0d0 - scale3*uscale(kk)
               usc5 = 1.0d0 - scale5*uscale(kk)
c
c     construct necessary auxiliary vectors
c

c               dixdk(1) = di(2)*dk(3) - di(3)*dk(2)
c               dixdk(2) = di(3)*dk(1) - di(1)*dk(3)
c               dixdk(3) = di(1)*dk(2) - di(2)*dk(1)
               dixuk(1) = di(2)*uind(3,l2) - di(3)*uind(2,l2)
               dixuk(2) = di(3)*uind(1,l2) - di(1)*uind(3,l2)
               dixuk(3) = di(1)*uind(2,l2) - di(2)*uind(1,l2)
               dkxui(1) = dk(2)*uind(3,l1) - dk(3)*uind(2,l1)
               dkxui(2) = dk(3)*uind(1,l1) - dk(1)*uind(3,l1)
               dkxui(3) = dk(1)*uind(2,l1) - dk(2)*uind(1,l1)
               dixukp(1) = di(2)*uinp(3,l2) - di(3)*uinp(2,l2)
               dixukp(2) = di(3)*uinp(1,l2) - di(1)*uinp(3,l2)
               dixukp(3) = di(1)*uinp(2,l2) - di(2)*uinp(1,l2)
               dkxuip(1) = dk(2)*uinp(3,l1) - dk(3)*uinp(2,l1)
               dkxuip(2) = dk(3)*uinp(1,l1) - dk(1)*uinp(3,l1)
               dkxuip(3) = dk(1)*uinp(2,l1) - dk(2)*uinp(1,l1)
               dixr(1) = di(2)*zr - di(3)*yr
               dixr(2) = di(3)*xr - di(1)*zr
               dixr(3) = di(1)*yr - di(2)*xr
               dkxr(1) = dk(2)*zr - dk(3)*yr
               dkxr(2) = dk(3)*xr - dk(1)*zr
               dkxr(3) = dk(1)*yr - dk(2)*xr
               qir(1) = qi(1)*xr + qi(4)*yr + qi(7)*zr
               qir(2) = qi(2)*xr + qi(5)*yr + qi(8)*zr
               qir(3) = qi(3)*xr + qi(6)*yr + qi(9)*zr
               qkr(1) = qk(1)*xr + qk(4)*yr + qk(7)*zr
               qkr(2) = qk(2)*xr + qk(5)*yr + qk(8)*zr
               qkr(3) = qk(3)*xr + qk(6)*yr + qk(9)*zr
c               qiqkr(1) = qi(1)*qkr(1) + qi(4)*qkr(2) + qi(7)*qkr(3)
c               qiqkr(2) = qi(2)*qkr(1) + qi(5)*qkr(2) + qi(8)*qkr(3)
c               qiqkr(3) = qi(3)*qkr(1) + qi(6)*qkr(2) + qi(9)*qkr(3)
c               qkqir(1) = qk(1)*qir(1) + qk(4)*qir(2) + qk(7)*qir(3)
c               qkqir(2) = qk(2)*qir(1) + qk(5)*qir(2) + qk(8)*qir(3)
c               qkqir(3) = qk(3)*qir(1) + qk(6)*qir(2) + qk(9)*qir(3)
c               qixqk(1) = qi(2)*qk(3) + qi(5)*qk(6) + qi(8)*qk(9)
c     &                       - qi(3)*qk(2) - qi(6)*qk(5) - qi(9)*qk(8)
c               qixqk(2) = qi(3)*qk(1) + qi(6)*qk(4) + qi(9)*qk(7)
c     &                       - qi(1)*qk(3) - qi(4)*qk(6) - qi(7)*qk(9)
c               qixqk(3) = qi(1)*qk(2) + qi(4)*qk(5) + qi(7)*qk(8)
c     &                       - qi(2)*qk(1) - qi(5)*qk(4) - qi(8)*qk(7)
               rxqir(1) = yr*qir(3) - zr*qir(2)
               rxqir(2) = zr*qir(1) - xr*qir(3)
               rxqir(3) = xr*qir(2) - yr*qir(1)
               rxqkr(1) = yr*qkr(3) - zr*qkr(2)
               rxqkr(2) = zr*qkr(1) - xr*qkr(3)
               rxqkr(3) = xr*qkr(2) - yr*qkr(1)
c               rxqikr(1) = yr*qiqkr(3) - zr*qiqkr(2)
c               rxqikr(2) = zr*qiqkr(1) - xr*qiqkr(3)
c               rxqikr(3) = xr*qiqkr(2) - yr*qiqkr(1)
c               rxqkir(1) = yr*qkqir(3) - zr*qkqir(2)
c               rxqkir(2) = zr*qkqir(1) - xr*qkqir(3)
c               rxqkir(3) = xr*qkqir(2) - yr*qkqir(1)
c               qkrxqir(1) = qkr(2)*qir(3) - qkr(3)*qir(2)
c               qkrxqir(2) = qkr(3)*qir(1) - qkr(1)*qir(3)
c               qkrxqir(3) = qkr(1)*qir(2) - qkr(2)*qir(1)
c               qidk(1) = qi(1)*dk(1) + qi(4)*dk(2) + qi(7)*dk(3)
c               qidk(2) = qi(2)*dk(1) + qi(5)*dk(2) + qi(8)*dk(3)
c               qidk(3) = qi(3)*dk(1) + qi(6)*dk(2) + qi(9)*dk(3)
c               qkdi(1) = qk(1)*di(1) + qk(4)*di(2) + qk(7)*di(3)
c               qkdi(2) = qk(2)*di(1) + qk(5)*di(2) + qk(8)*di(3)
c               qkdi(3) = qk(3)*di(1) + qk(6)*di(2) + qk(9)*di(3)
               qiuk(1) = qi(1)*uind(1,l2) + qi(4)*uind(2,l2)
     &                      + qi(7)*uind(3,l2)
               qiuk(2) = qi(2)*uind(1,l2) + qi(5)*uind(2,l2)
     &                      + qi(8)*uind(3,l2)
               qiuk(3) = qi(3)*uind(1,l2) + qi(6)*uind(2,l2)
     &                      + qi(9)*uind(3,l2)
               qkui(1) = qk(1)*uind(1,l1) + qk(4)*uind(2,l1)
     &                      + qk(7)*uind(3,l1)
               qkui(2) = qk(2)*uind(1,l1) + qk(5)*uind(2,l1)
     &                      + qk(8)*uind(3,l1)
               qkui(3) = qk(3)*uind(1,l1) + qk(6)*uind(2,l1)
     &                      + qk(9)*uind(3,l1)
               qiukp(1) = qi(1)*uinp(1,l2) + qi(4)*uinp(2,l2)
     &                       + qi(7)*uinp(3,l2)
               qiukp(2) = qi(2)*uinp(1,l2) + qi(5)*uinp(2,l2)
     &                       + qi(8)*uinp(3,l2)
               qiukp(3) = qi(3)*uinp(1,l2) + qi(6)*uinp(2,l2)
     &                       + qi(9)*uinp(3,l2)
               qkuip(1) = qk(1)*uinp(1,l1) + qk(4)*uinp(2,l1)
     &                       + qk(7)*uinp(3,l1)
               qkuip(2) = qk(2)*uinp(1,l1) + qk(5)*uinp(2,l1)
     &                       + qk(8)*uinp(3,l1)
               qkuip(3) = qk(3)*uinp(1,l1) + qk(6)*uinp(2,l1)
     &                       + qk(9)*uinp(3,l1)
c               dixqkr(1) = di(2)*qkr(3) - di(3)*qkr(2)
c               dixqkr(2) = di(3)*qkr(1) - di(1)*qkr(3)
c               dixqkr(3) = di(1)*qkr(2) - di(2)*qkr(1)
c               dkxqir(1) = dk(2)*qir(3) - dk(3)*qir(2)
c               dkxqir(2) = dk(3)*qir(1) - dk(1)*qir(3)
c               dkxqir(3) = dk(1)*qir(2) - dk(2)*qir(1)
               uixqkr(1) = uind(2,l1)*qkr(3) - uind(3,l1)*qkr(2)
               uixqkr(2) = uind(3,l1)*qkr(1) - uind(1,l1)*qkr(3)
               uixqkr(3) = uind(1,l1)*qkr(2) - uind(2,l1)*qkr(1)
               ukxqir(1) = uind(2,l2)*qir(3) - uind(3,l2)*qir(2)
               ukxqir(2) = uind(3,l2)*qir(1) - uind(1,l2)*qir(3)
               ukxqir(3) = uind(1,l2)*qir(2) - uind(2,l2)*qir(1)
               uixqkrp(1) = uinp(2,l1)*qkr(3) - uinp(3,l1)*qkr(2)
               uixqkrp(2) = uinp(3,l1)*qkr(1) - uinp(1,l1)*qkr(3)
               uixqkrp(3) = uinp(1,l1)*qkr(2) - uinp(2,l1)*qkr(1)
               ukxqirp(1) = uinp(2,l2)*qir(3) - uinp(3,l2)*qir(2)
               ukxqirp(2) = uinp(3,l2)*qir(1) - uinp(1,l2)*qir(3)
               ukxqirp(3) = uinp(1,l2)*qir(2) - uinp(2,l2)*qir(1)
c               rxqidk(1) = yr*qidk(3) - zr*qidk(2)
c               rxqidk(2) = zr*qidk(1) - xr*qidk(3)
c               rxqidk(3) = xr*qidk(2) - yr*qidk(1)
c               rxqkdi(1) = yr*qkdi(3) - zr*qkdi(2)
c               rxqkdi(2) = zr*qkdi(1) - xr*qkdi(3)
c               rxqkdi(3) = xr*qkdi(2) - yr*qkdi(1)
               rxqiuk(1) = yr*qiuk(3) - zr*qiuk(2)
               rxqiuk(2) = zr*qiuk(1) - xr*qiuk(3)
               rxqiuk(3) = xr*qiuk(2) - yr*qiuk(1)
               rxqkui(1) = yr*qkui(3) - zr*qkui(2)
               rxqkui(2) = zr*qkui(1) - xr*qkui(3)
               rxqkui(3) = xr*qkui(2) - yr*qkui(1)
               rxqiukp(1) = yr*qiukp(3) - zr*qiukp(2)
               rxqiukp(2) = zr*qiukp(1) - xr*qiukp(3)
               rxqiukp(3) = xr*qiukp(2) - yr*qiukp(1)
               rxqkuip(1) = yr*qkuip(3) - zr*qkuip(2)
               rxqkuip(2) = zr*qkuip(1) - xr*qkuip(3)
               rxqkuip(3) = xr*qkuip(2) - yr*qkuip(1)
c
c     calculate the scalar products for permanent components
c
               sc(2) = di(1)*dk(1) + di(2)*dk(2) + di(3)*dk(3)
               sc(3) = di(1)*xr + di(2)*yr + di(3)*zr
               sc(4) = dk(1)*xr + dk(2)*yr + dk(3)*zr
               sc(5) = qir(1)*xr + qir(2)*yr + qir(3)*zr
               sc(6) = qkr(1)*xr + qkr(2)*yr + qkr(3)*zr
               sc(7) = qir(1)*dk(1) + qir(2)*dk(2) + qir(3)*dk(3)
               sc(8) = qkr(1)*di(1) + qkr(2)*di(2) + qkr(3)*di(3)
               sc(9) = qir(1)*qkr(1) + qir(2)*qkr(2) + qir(3)*qkr(3)
               sc(10) = qi(1)*qk(1) + qi(2)*qk(2) + qi(3)*qk(3)
     &                     + qi(4)*qk(4) + qi(5)*qk(5) + qi(6)*qk(6)
     &                     + qi(7)*qk(7) + qi(8)*qk(8) + qi(9)*qk(9)
c
c     calculate the scalar products for induced components
c
               sci(1) = uind(1,l1)*dk(1) + uind(2,l1)*dk(2)
     &                     + uind(3,l1)*dk(3) + di(1)*uind(1,l2)
     &                     + di(2)*uind(2,l2) + di(3)*uind(3,l2)
               sci(2) = uind(1,l1)*uind(1,l2) + uind(2,l1)*uind(2,l2)
     &                     + uind(3,l1)*uind(3,l2)
               sci(3) = uind(1,l1)*xr + uind(2,l1)*yr + uind(3,l1)*zr
               sci(4) = uind(1,l2)*xr + uind(2,l2)*yr + uind(3,l2)*zr
               sci(7) = qir(1)*uind(1,l2) + qir(2)*uind(2,l2)
     &                     + qir(3)*uind(3,l2)
               sci(8) = qkr(1)*uind(1,l1) + qkr(2)*uind(2,l1)
     &                     + qkr(3)*uind(3,l1)
               scip(1) = uinp(1,l1)*dk(1) + uinp(2,l1)*dk(2)
     &                      + uinp(3,l1)*dk(3) + di(1)*uinp(1,l2)
     &                      + di(2)*uinp(2,l2) + di(3)*uinp(3,l2)
               scip(2) = uind(1,l1)*uinp(1,l2)+uind(2,l1)*uinp(2,l2)
     &                   + uind(3,l1)*uinp(3,l2)+uinp(1,l1)*uind(1,l2)
     &                   + uinp(2,l1)*uind(2,l2)+uinp(3,l1)*uind(3,l2)
               scip(3) = uinp(1,l1)*xr + uinp(2,l1)*yr + uinp(3,l1)*zr
               scip(4) = uinp(1,l2)*xr + uinp(2,l2)*yr + uinp(3,l2)*zr
               scip(7) = qir(1)*uinp(1,l2) + qir(2)*uinp(2,l2)
     &                      + qir(3)*uinp(3,l2)
               scip(8) = qkr(1)*uinp(1,l1) + qkr(2)*uinp(2,l1)
     &                      + qkr(3)*uinp(3,l1)
c
c     calculate the gl functions for induced components
c
               gli(1) = ck*sci(3) - ci*sci(4)
               gli(2) = -sc(3)*sci(4) - sci(3)*sc(4)
               gli(3) = sci(3)*sc(6) - sci(4)*sc(5)
               gli(6) = sci(1)
               gli(7) = 2.0d0 * (sci(7)-sci(8))
               glip(1) = ck*scip(3) - ci*scip(4)
               glip(2) = -sc(3)*scip(4) - scip(3)*sc(4)
               glip(3) = scip(3)*sc(6) - scip(4)*sc(5)
               glip(6) = scip(1)
               glip(7) = 2.0d0 * (scip(7)-scip(8))
c
c     compute the energy contributions for this interaction
c
               ei = 0.5d0 * (bn(1)*(gli(1)+gli(6))
     &                      + bn(2)*(gli(2)+gli(7)) + bn(3)*gli(3))
c
c     get the real energy without any screening function
c
               erli = 0.5d0*(rr3*(gli(1)+gli(6))*psc3
     &                   + rr5*(gli(2)+gli(7))*psc5
     &                   + rr7*gli(3)*psc7)
               ei = ei - erli
               ei = f * ei
               if (ii .eq. kk) then
                 ei = 0.5d0 * ei
               end if
               eptemp = eptemp + ei
c
c     increment the total intramolecular energy; assumes
c     intramolecular distances are less than half of cell
c     length and less than the ewald cutoff
c
c               if (molcule(ii) .eq. molcule(kk)) then
c                  eintra = eintra + 0.5d0*pscale(kk)
c     &                        * (rr3*(gli(1)+gli(6))*scale3
c     &                              + rr5*(gli(2)+gli(7))*scale5
c     &                              + rr7*gli(3)*scale7)
c               end if
c
c     intermediate variables for induced force terms
c
               gfi(1) = 0.5d0*bn(2)*(gli(1)+glip(1)+gli(6)+glip(6))
     &                     + 0.5d0*bn(2)*scip(2)
     &                     + 0.5d0*bn(3)*(gli(2)+glip(2)+gli(7)+glip(7))
     &                     - 0.5d0*bn(3)*(sci(3)*scip(4)+scip(3)*sci(4))
     &                     + 0.5d0*bn(4)*(gli(3)+glip(3))
               gfi(2) = -ck*bn(1) + sc(4)*bn(2) - sc(6)*bn(3)
               gfi(3) = ci*bn(1) + sc(3)*bn(2) + sc(5)*bn(3)
               gfi(4) = 2.0d0 * bn(2)
               gfi(5) = bn(3) * (sci(4)+scip(4))
               gfi(6) = -bn(3) * (sci(3)+scip(3))
               gfri(1) = 0.5d0*rr5*((gli(1)+gli(6))*psc3
     &                            + (glip(1)+glip(6))*dsc3
     &                            + scip(2)*usc3)
     &                 + 0.5d0*rr7*((gli(7)+gli(2))*psc5
     &                            + (glip(7)+glip(2))*dsc5
     &                     - (sci(3)*scip(4)+scip(3)*sci(4))*usc5)
     &                 + 0.5d0*rr9*(gli(3)*psc7+glip(3)*dsc7)
               gfri(2) = -rr3*ck + rr5*sc(4) - rr7*sc(6)
               gfri(3) = rr3*ci + rr5*sc(3) + rr7*sc(5)
               gfri(4) = 2.0d0 * rr5
               gfri(5) = rr7 * (sci(4)*psc7+scip(4)*dsc7)
               gfri(6) = -rr7 * (sci(3)*psc7+scip(3)*dsc7)
c
c     get the induced force with screening
c
               ftm2i(1) = gfi(1)*xr + 0.5d0*
     &             (gfi(2)*(uind(1,l1)+uinp(1,l1))
     &            + bn(2)*(sci(4)*uinp(1,l1)+scip(4)*uind(1,l1))
     &            + gfi(3)*(uind(1,l2)+uinp(1,l2))
     &            + bn(2)*(sci(3)*uinp(1,l2)+scip(3)*uind(1,l2))
     &            + (sci(4)+scip(4))*bn(2)*di(1)
     &            + (sci(3)+scip(3))*bn(2)*dk(1)
     &            + gfi(4)*(qkui(1)+qkuip(1)-qiuk(1)-qiukp(1)))
     &            + gfi(5)*qir(1) + gfi(6)*qkr(1)
               ftm2i(2) = gfi(1)*yr + 0.5d0*
     &             (gfi(2)*(uind(2,l1)+uinp(2,l1))
     &            + bn(2)*(sci(4)*uinp(2,l1)+scip(4)*uind(2,l1))
     &            + gfi(3)*(uind(2,l2)+uinp(2,l2))
     &            + bn(2)*(sci(3)*uinp(2,l2)+scip(3)*uind(2,l2))
     &            + (sci(4)+scip(4))*bn(2)*di(2)
     &            + (sci(3)+scip(3))*bn(2)*dk(2)
     &            + gfi(4)*(qkui(2)+qkuip(2)-qiuk(2)-qiukp(2)))
     &            + gfi(5)*qir(2) + gfi(6)*qkr(2)
               ftm2i(3) = gfi(1)*zr + 0.5d0*
     &             (gfi(2)*(uind(3,l1)+uinp(3,l1))
     &            + bn(2)*(sci(4)*uinp(3,l1)+scip(4)*uind(3,l1))
     &            + gfi(3)*(uind(3,l2)+uinp(3,l2))
     &            + bn(2)*(sci(3)*uinp(3,l2)+scip(3)*uind(3,l2))
     &            + (sci(4)+scip(4))*bn(2)*di(3)
     &            + (sci(3)+scip(3))*bn(2)*dk(3)
     &            + gfi(4)*(qkui(3)+qkuip(3)-qiuk(3)-qiukp(3)))
     &            + gfi(5)*qir(3) + gfi(6)*qkr(3)
c
c     get the induced force without screening
c
               ftm2ri(1) = gfri(1)*xr + 0.5d0*
     &           (- rr3*ck*(uind(1,l1)*psc3+uinp(1,l1)*dsc3)
     &            + rr5*sc(4)*(uind(1,l1)*psc5+uinp(1,l1)*dsc5)
     &            - rr7*sc(6)*(uind(1,l1)*psc7+uinp(1,l1)*dsc7))
     &            + (rr3*ci*(uind(1,l2)*psc3+uinp(1,l2)*dsc3)
     &            + rr5*sc(3)*(uind(1,l2)*psc5+uinp(1,l2)*dsc5)
     &            + rr7*sc(5)*(uind(1,l2)*psc7+uinp(1,l2)*dsc7))*0.5d0
     &            + rr5*usc5*(sci(4)*uinp(1,l1)+scip(4)*uind(1,l1)
     &            + sci(3)*uinp(1,l2)+scip(3)*uind(1,l2))*0.5d0
     &            + 0.5d0*(sci(4)*psc5+scip(4)*dsc5)*rr5*di(1)
     &            + 0.5d0*(sci(3)*psc5+scip(3)*dsc5)*rr5*dk(1)
     &            + 0.5d0*gfri(4)*((qkui(1)-qiuk(1))*psc5
     &            + (qkuip(1)-qiukp(1))*dsc5)
     &            + gfri(5)*qir(1) + gfri(6)*qkr(1)
               ftm2ri(2) = gfri(1)*yr + 0.5d0*
     &           (- rr3*ck*(uind(2,l1)*psc3+uinp(2,l1)*dsc3)
     &            + rr5*sc(4)*(uind(2,l1)*psc5+uinp(2,l1)*dsc5)
     &            - rr7*sc(6)*(uind(2,l1)*psc7+uinp(2,l1)*dsc7))
     &            + (rr3*ci*(uind(2,l2)*psc3+uinp(2,l2)*dsc3)
     &            + rr5*sc(3)*(uind(2,l2)*psc5+uinp(2,l2)*dsc5)
     &            + rr7*sc(5)*(uind(2,l2)*psc7+uinp(2,l2)*dsc7))*0.5d0
     &            + rr5*usc5*(sci(4)*uinp(2,l1)+scip(4)*uind(2,l1)
     &            + sci(3)*uinp(2,l2)+scip(3)*uind(2,l2))*0.5d0
     &            + 0.5d0*(sci(4)*psc5+scip(4)*dsc5)*rr5*di(2)
     &            + 0.5d0*(sci(3)*psc5+scip(3)*dsc5)*rr5*dk(2)
     &            + 0.5d0*gfri(4)*((qkui(2)-qiuk(2))*psc5
     &            + (qkuip(2)-qiukp(2))*dsc5)
     &            + gfri(5)*qir(2) + gfri(6)*qkr(2)
               ftm2ri(3) = gfri(1)*zr + 0.5d0*
     &           (- rr3*ck*(uind(3,l1)*psc3+uinp(3,l1)*dsc3)
     &            + rr5*sc(4)*(uind(3,l1)*psc5+uinp(3,l1)*dsc5)
     &            - rr7*sc(6)*(uind(3,l1)*psc7+uinp(3,l1)*dsc7))
     &            + (rr3*ci*(uind(3,l2)*psc3+uinp(3,l2)*dsc3)
     &            + rr5*sc(3)*(uind(3,l2)*psc5+uinp(3,l2)*dsc5)
     &            + rr7*sc(5)*(uind(3,l2)*psc7+uinp(3,l2)*dsc7))*0.5d0
     &            + rr5*usc5*(sci(4)*uinp(3,l1)+scip(4)*uind(3,l1)
     &            + sci(3)*uinp(3,l2)+scip(3)*uind(3,l2))*0.5d0
     &            + 0.5d0*(sci(4)*psc5+scip(4)*dsc5)*rr5*di(3)
     &            + 0.5d0*(sci(3)*psc5+scip(3)*dsc5)*rr5*dk(3)
     &            + 0.5d0*gfri(4)*((qkui(3)-qiuk(3))*psc5
     &            + (qkuip(3)-qiukp(3))*dsc5)
     &            + gfri(5)*qir(3) + gfri(6)*qkr(3)
c
c     account for partially excluded induced interactions
c
               temp3 = 0.5d0 * rr3 * ((gli(1)+gli(6))*pscale(kk)
     &                                  +(glip(1)+glip(6))*dscale(kk))
               temp5 = 0.5d0 * rr5 * ((gli(2)+gli(7))*pscale(kk)
     &                                  +(glip(2)+glip(7))*dscale(kk))
               temp7 = 0.5d0 * rr7 * (gli(3)*pscale(kk)
     &                                  +glip(3)*dscale(kk))
               fridmp(1) = temp3*ddsc3(1) + temp5*ddsc5(1)
     &                        + temp7*ddsc7(1)
               fridmp(2) = temp3*ddsc3(2) + temp5*ddsc5(2)
     &                        + temp7*ddsc7(2)
               fridmp(3) = temp3*ddsc3(3) + temp5*ddsc5(3)
     &                        + temp7*ddsc7(3)
c
c     find some scaling terms for induced-induced force
c
               temp3 = 0.5d0 * rr3 * uscale(kk) * scip(2)
               temp5 = -0.5d0 * rr5 * uscale(kk)
     &                    * (sci(3)*scip(4)+scip(3)*sci(4))
               findmp(1) = temp3*ddsc3(1) + temp5*ddsc5(1)
               findmp(2) = temp3*ddsc3(2) + temp5*ddsc5(2)
               findmp(3) = temp3*ddsc3(3) + temp5*ddsc5(3)
c
c     modify the forces for partially excluded interactions
c
               ftm2i(1) = ftm2i(1) - fridmp(1) - findmp(1)
               ftm2i(2) = ftm2i(2) - fridmp(2) - findmp(2)
               ftm2i(3) = ftm2i(3) - fridmp(3) - findmp(3)
c
c     correction to convert mutual to direct polarization force
c

c               if (poltyp .eq. 'DIRECT') then
c                  gfd = 0.5d0 * (bn(2)*scip(2)
c     &                     - bn(3)*(scip(3)*sci(4)+sci(3)*scip(4)))
c                  gfdr = 0.5d0 * (rr5*scip(2)*usc3
c     &                     - rr7*(scip(3)*sci(4)
c     &                           +sci(3)*scip(4))*usc5)
c
c                  ftm2i(1) = ftm2i(1) - gfd*xr - 0.5d0*bn(2)*
c     &                          (sci(4)*uinp(1,l1)+scip(4)*uind(1,l1)
c     &                          +sci(3)*uinp(1,l2)+scip(3)*uind(1,l2))
c                  ftm2i(2) = ftm2i(2) - gfd*yr - 0.5d0*bn(2)*
c     &                          (sci(4)*uinp(2,l1)+scip(4)*uind(2,l1)
c     &                          +sci(3)*uinp(2,l2)+scip(3)*uind(2,l2))
c                  ftm2i(3) = ftm2i(3) - gfd*zr - 0.5d0*bn(2)*
c     &                          (sci(4)*uinp(3,l1)+scip(4)*uind(3,l1)
c     &                          +sci(3)*uinp(3,l2)+scip(3)*uind(3,l2))
c                  fdir(1) = gfdr*xr + 0.5d0*usc5*rr5*
c     &                         (sci(4)*uinp(1,l1)+scip(4)*uind(1,l1)
c     &                        + sci(3)*uinp(1,l2)+scip(3)*uind(1,l2))
c                  fdir(2) = gfdr*yr + 0.5d0*usc5*rr5*
c     &                         (sci(4)*uinp(2,l1)+scip(4)*uind(2,l1)
c     &                        + sci(3)*uinp(2,l2)+scip(3)*uind(2,l2))
c                  fdir(3) = gfdr*zr + 0.5d0*usc5*rr5*
c     &                         (sci(4)*uinp(3,l1)+scip(4)*uind(3,l1)
c     &                        + sci(3)*uinp(3,l2)+scip(3)*uind(3,l2))
c                  ftm2i(1) = ftm2i(1) + fdir(1) + findmp(1)
c                  ftm2i(2) = ftm2i(2) + fdir(2) + findmp(2)
c                  ftm2i(3) = ftm2i(3) + fdir(3) + findmp(3)
c               end if

c
c     intermediate variables for induced torque terms
c
               gti(2) = 0.5d0 * bn(2) * (sci(4)+scip(4))
               gti(3) = 0.5d0 * bn(2) * (sci(3)+scip(3))
               gti(4) = gfi(4)
               gti(5) = gfi(5)
               gti(6) = gfi(6)
               gtri(2) = 0.5d0 * rr5 * (sci(4)*psc5+scip(4)*dsc5)
               gtri(3) = 0.5d0 * rr5 * (sci(3)*psc5+scip(3)*dsc5)
               gtri(4) = gfri(4)
               gtri(5) = gfri(5)
               gtri(6) = gfri(6)
c
c     get the induced torque with screening
c
               ttm2i(1) = -bn(1)*(dixuk(1)+dixukp(1))*0.5d0
     &           + gti(2)*dixr(1) + gti(4)*(ukxqir(1)+rxqiuk(1)
     &           + ukxqirp(1)+rxqiukp(1))*0.5d0 - gti(5)*rxqir(1)
               ttm2i(2) = -bn(1)*(dixuk(2)+dixukp(2))*0.5d0
     &           + gti(2)*dixr(2) + gti(4)*(ukxqir(2)+rxqiuk(2)
     &           + ukxqirp(2)+rxqiukp(2))*0.5d0 - gti(5)*rxqir(2)
               ttm2i(3) = -bn(1)*(dixuk(3)+dixukp(3))*0.5d0
     &           + gti(2)*dixr(3) + gti(4)*(ukxqir(3)+rxqiuk(3)
     &           + ukxqirp(3)+rxqiukp(3))*0.5d0 - gti(5)*rxqir(3)
               ttm3i(1) = -bn(1)*(dkxui(1)+dkxuip(1))*0.5d0
     &           + gti(3)*dkxr(1) - gti(4)*(uixqkr(1)+rxqkui(1)
     &           + uixqkrp(1)+rxqkuip(1))*0.5d0 - gti(6)*rxqkr(1)
               ttm3i(2) = -bn(1)*(dkxui(2)+dkxuip(2))*0.5d0
     &           + gti(3)*dkxr(2) - gti(4)*(uixqkr(2)+rxqkui(2)
     &           + uixqkrp(2)+rxqkuip(2))*0.5d0 - gti(6)*rxqkr(2)
               ttm3i(3) = -bn(1)*(dkxui(3)+dkxuip(3))*0.5d0
     &           + gti(3)*dkxr(3) - gti(4)*(uixqkr(3)+rxqkui(3)
     &           + uixqkrp(3)+rxqkuip(3))*0.5d0 - gti(6)*rxqkr(3)
c
c     get the induced torque without screening
c
               ttm2ri(1) = -rr3*(dixuk(1)*psc3+dixukp(1)*dsc3)*0.5d0
     &           + gtri(2)*dixr(1) + gtri(4)*((ukxqir(1)+rxqiuk(1))*psc5
     &           +(ukxqirp(1)+rxqiukp(1))*dsc5)*0.5d0 - gtri(5)*rxqir(1)
               ttm2ri(2) = -rr3*(dixuk(2)*psc3+dixukp(2)*dsc3)*0.5d0
     &           + gtri(2)*dixr(2) + gtri(4)*((ukxqir(2)+rxqiuk(2))*psc5
     &           +(ukxqirp(2)+rxqiukp(2))*dsc5)*0.5d0 - gtri(5)*rxqir(2)
               ttm2ri(3) = -rr3*(dixuk(3)*psc3+dixukp(3)*dsc3)*0.5d0
     &           + gtri(2)*dixr(3) + gtri(4)*((ukxqir(3)+rxqiuk(3))*psc5
     &           +(ukxqirp(3)+rxqiukp(3))*dsc5)*0.5d0 - gtri(5)*rxqir(3)
               ttm3ri(1) = -rr3*(dkxui(1)*psc3+dkxuip(1)*dsc3)*0.5d0
     &           + gtri(3)*dkxr(1) - gtri(4)*((uixqkr(1)+rxqkui(1))*psc5
     &           +(uixqkrp(1)+rxqkuip(1))*dsc5)*0.5d0 - gtri(6)*rxqkr(1)
               ttm3ri(2) = -rr3*(dkxui(2)*psc3+dkxuip(2)*dsc3)*0.5d0
     &           + gtri(3)*dkxr(2) - gtri(4)*((uixqkr(2)+rxqkui(2))*psc5
     &           +(uixqkrp(2)+rxqkuip(2))*dsc5)*0.5d0 - gtri(6)*rxqkr(2)
               ttm3ri(3) = -rr3*(dkxui(3)*psc3+dkxuip(3)*dsc3)*0.5d0
     &           + gtri(3)*dkxr(3) - gtri(4)*((uixqkr(3)+rxqkui(3))*psc5
     &           +(uixqkrp(3)+rxqkuip(3))*dsc5)*0.5d0 - gtri(6)*rxqkr(3)
c
c     handle the case where scaling is used
c
              if (use_polymer .and. r2.le.polycut2) then
               do j = 1, 3
                  ftm2i(j) = f * (ftm2i(j)-ftm2ri(j))
                  ttm2i(j) = f * (ttm2i(j)-ttm2ri(j))
                  ttm3i(j) = f * (ttm3i(j)-ttm3ri(j))
               end do
              else
               do j = 1, 3
                  ftm2i(j) = f * (ftm2i(j)-ftm2ri(j))
                  ttm2i(j) = f * (ttm2i(j)-ttm2ri(j))
                  ttm3i(j) = f * (ttm3i(j)-ttm3ri(j))
               end do
              end if
               if (ii .eq. kk) then
                  do j = 1, 3
                     ftm2i(j) = 0.5d0 * ftm2i(j)
                     ttm2i(j) = 0.5d0 * ttm2i(j)
                     ttm3i(j) = 0.5d0 * ttm3i(j)
                  end do
               end if

c
c     increment gradient due to force and torque on first site
c
c               deptemp(1,l1) = deptemp(1,l1) + ftm2i(1)
c               deptemp(2,l1) = deptemp(2,l1) + ftm2i(2)
c               deptemp(3,l1) = deptemp(3,l1) + ftm2i(3)
c               call torque_3b (pnum2,deptemp,i,
c     &           ttm2,ttm2i,frcxi,frcyi,frczi)

               deptemp(1,i) = deptemp(1,i) + ftm2i(1)
               deptemp(2,i) = deptemp(2,i) + ftm2i(2)
               deptemp(3,i) = deptemp(3,i) + ftm2i(3)
               call torque_3b (deptemp,i,
     &           ttm2,ttm2i,frcxi,frcyi,frczi)

c
c     increment gradient due to force and torque on second site
c
c               deptemp(1,l2) = deptemp(1,l2) - ftm2i(1)
c               deptemp(2,l2) = deptemp(2,l2) - ftm2i(2)
c               deptemp(3,l2) = deptemp(3,l2) - ftm2i(3)
c               call torque_3b (pnum2,deptemp,k,
c     &            ttm3,ttm3i,frcxk,frcyk,frczk)

               deptemp(1,k) = deptemp(1,k) - ftm2i(1)
               deptemp(2,k) = deptemp(2,k) - ftm2i(2)
               deptemp(3,k) = deptemp(3,k) - ftm2i(3)
               call torque_3b (deptemp,k,
     &            ttm3,ttm3i,frcxk,frcyk,frczk)

c      subroutine torque_3b(pnum2,deptemp,
c     &    i,trq1,trq2,frcx,frcy,frcz)
c
c     increment the internal virial tensor components
c
               iaz = zaxis(i)
               iax = xaxis(i)
               iay = yaxis(i)
               kaz = zaxis(k)
               kax = xaxis(k)
               kay = yaxis(k)
               if (iaz .eq. 0)  iaz = ii
               if (iax .eq. 0)  iax = ii
               if (iay .eq. 0)  iay = ii
               if (kaz .eq. 0)  kaz = kk
               if (kax .eq. 0)  kax = kk
               if (kay .eq. 0)  kay = kk
               xiz = x(iaz) - x(ii)
               yiz = y(iaz) - y(ii)
               ziz = z(iaz) - z(ii)
               xix = x(iax) - x(ii)
               yix = y(iax) - y(ii)
               zix = z(iax) - z(ii)
               xiy = x(iay) - x(ii)
               yiy = y(iay) - y(ii)
               ziy = z(iay) - z(ii)
               xkz = x(kaz) - x(kk)
               ykz = y(kaz) - y(kk)
               zkz = z(kaz) - z(kk)
               xkx = x(kax) - x(kk)
               ykx = y(kax) - y(kk)
               zkx = z(kax) - z(kk)
               xky = x(kay) - x(kk)
               yky = y(kay) - y(kk)
               zky = z(kay) - z(kk)
c               vxx = -xr*ftm2i(1)+xix*frcxi(1)+xiy*frcyi(1)+xiz*frczi(1)
c               vyx = -yr*ftm2i(1)+yix*frcxi(1)+yiy*frcyi(1)+yiz*frczi(1)
c               vzx = -zr*ftm2i(1)+zix*frcxi(1)+ziy*frcyi(1)+ziz*frczi(1)
c               vyy = -yr*ftm2i(2)+yix*frcxi(2)+yiy*frcyi(2)+yiz*frczi(2)
c               vzy = -zr*ftm2i(2)+zix*frcxi(2)+ziy*frcyi(2)+ziz*frczi(2)
c               vzz = -zr*ftm2i(3)+zix*frcxi(3)+ziy*frcyi(3)+ziz*frczi(3)
c               vir(1,1) = vir(1,1) + vxx
c               vir(2,1) = vir(2,1) + vyx
c               vir(3,1) = vir(3,1) + vzx
c               vir(1,2) = vir(1,2) + vyx
c               vir(2,2) = vir(2,2) + vyy
c               vir(3,2) = vir(3,2) + vzy
c               vir(1,3) = vir(1,3) + vzx
c               vir(2,3) = vir(2,3) + vzy
c               vir(3,3) = vir(3,3) + vzz
               vxx = -xr*(ftm2i(1)) + xix*frcxi(1)
     &                  + xiy*frcyi(1) + xiz*frczi(1) + xkx*frcxk(1)
     &                  + xky*frcyk(1) + xkz*frczk(1)
               vyx = -yr*(ftm2i(1)) + yix*frcxi(1)
     &                  + yiy*frcyi(1) + yiz*frczi(1) + ykx*frcxk(1)
     &                  + yky*frcyk(1) + ykz*frczk(1)
               vzx = -zr*(ftm2i(1)) + zix*frcxi(1)
     &                  + ziy*frcyi(1) + ziz*frczi(1) + zkx*frcxk(1)
     &                  + zky*frcyk(1) + zkz*frczk(1)
               vyy = -yr*(ftm2i(2)) + yix*frcxi(2)
     &                  + yiy*frcyi(2) + yiz*frczi(2) + ykx*frcxk(2)
     &                  + yky*frcyk(2) + ykz*frczk(2)
               vzy = -zr*(ftm2i(2)) + zix*frcxi(2)
     &                  + ziy*frcyi(2) + ziz*frczi(2) + zkx*frcxk(2)
     &                  + zky*frcyk(2) + zkz*frczk(2)
               vzz = -zr*(ftm2i(3)) + zix*frcxi(3)
     &                  + ziy*frcyi(3) + ziz*frczi(3) + zkx*frcxk(3)
     &                  + zky*frcyk(3) + zkz*frczk(3)

               virtemp(1,1) = virtemp(1,1) + vxx
               virtemp(2,1) = virtemp(2,1) + vyx
               virtemp(3,1) = virtemp(3,1) + vzx
               virtemp(1,2) = virtemp(1,2) + vyx
               virtemp(2,2) = virtemp(2,2) + vyy
               virtemp(3,2) = virtemp(3,2) + vzy
               virtemp(1,3) = virtemp(1,3) + vzx
               virtemp(2,3) = virtemp(2,3) + vzy
               virtemp(3,3) = virtemp(3,3) + vzz

            end if
         end do
c         end if
         end do
c
c     reset interaction scaling coefficients for connected atoms
c
         do j = 1, n12(ii)
            mscale(i12(j,ii)) = 1.0d0
            pscale(i12(j,ii)) = 1.0d0
         end do
         do j = 1, n13(ii)
            mscale(i13(j,ii)) = 1.0d0
            pscale(i13(j,ii)) = 1.0d0
         end do
         do j = 1, n14(ii)
            mscale(i14(j,ii)) = 1.0d0
            pscale(i14(j,ii)) = 1.0d0
         end do
         do j = 1, n15(ii)
            mscale(i15(j,ii)) = 1.0d0
            pscale(i15(j,ii)) = 1.0d0
         end do
         do j = 1, np11(ii)
            dscale(ip11(j,ii)) = 1.0d0
            uscale(ip11(j,ii)) = 1.0d0
         end do
         do j = 1, np12(ii)
            dscale(ip12(j,ii)) = 1.0d0
            uscale(ip12(j,ii)) = 1.0d0
         end do
         do j = 1, np13(ii)
            dscale(ip13(j,ii)) = 1.0d0
            uscale(ip13(j,ii)) = 1.0d0
         end do
         do j = 1, np14(ii)
            dscale(ip14(j,ii)) = 1.0d0
            uscale(ip14(j,ii)) = 1.0d0
         end do
      end do
      end if
c
c     for periodic boundary conditions with large cutoffs
c     neighbors must be found by the replicates method
c
c
c     perform deallocation of some local arrays
c
      deallocate (mscale)
      deallocate (pscale)
      deallocate (dscale)
      deallocate (uscale)
c      print*,'ep after ereal1c_3b=',eptemp,moli1,moli2,moli3
c      do l1 = 1, npole3b
c         print*,"After ereal1c_3b l1 Pnum(l1)",l1,pnum(l1)
c      end do

      return
      end

c
c     ################################################################
c     ##                                                            ##
c     ##             Subroutine ereal1c_3b                          ##
c     ##  Liam O-Suilleabhain, Omar Demerdash, Teresa Head-Gordon   ##
c     ##                 Spring 2013                                ##
c     ##                                                            ##
c     ################################################################
c
c
c     "ereal1c_3b" evaluates the real space portion of the regular Ewald
c     summation energy and gradient due to dipole polarizability
c
      subroutine ereal1c_3b2 (eintra)
      implicit none
      include 'sizes.i'
      include 'action.i'
      include 'analyz.i'
      include 'atmtyp.i'
      include 'atoms.i'
      include 'bound.i'
      include 'boxes.i'
      include 'cell.i'
      include 'chgpot.i'
      include 'couple.i'
      include 'cutoff.i'
      include 'deriv.i'
      include 'energi.i'
      include 'ewald.i'
      include 'math.i'
      include 'molcul.i'
      include 'mplpot.i'
      include 'mpole.i'
      include 'polar.i'
      include 'polgrp.i'
      include 'polpot.i'
      include 'shunt.i'
      include 'potent.i'
      include 'virial.i'
      include 'combo.i'
      integer i,j,k,l1,l3
      integer ii,kk,jcell
      integer iax,iay,iaz
      integer kax,kay,kaz
      real*8 e,ei,f,bfac
      real*8 eintra,erfc
      real*8 damp,expdamp
      real*8 pdi,pti,pgamma
      real*8 scale3,scale5
      real*8 scale7
      real*8 temp3,temp5,temp7
      real*8 dsc3,dsc5,dsc7
      real*8 psc3,psc5,psc7
      real*8 usc3,usc5
      real*8 alsq2,alsq2n
      real*8 exp2a,ralpha
      real*8 gfd,gfdr
      real*8 xr,yr,zr
      real*8 xix,yix,zix
      real*8 xiy,yiy,ziy
      real*8 xiz,yiz,ziz
      real*8 xkx,ykx,zkx
      real*8 xky,yky,zky
      real*8 xkz,ykz,zkz
      real*8 r,r2,rr1,rr3
      real*8 rr5,rr7,rr9,rr11
      real*8 erl,erli
      real*8 vxx,vyy,vzz
      real*8 vyx,vzx,vzy
      real*8 frcxi(3),frcxk(3)
      real*8 frcyi(3),frcyk(3)
      real*8 frczi(3),frczk(3)
      real*8 ci,di(3),qi(9)
      real*8 ck,dk(3),qk(9)
      real*8 fridmp(3),findmp(3)
      real*8 ftm2(3),ftm2i(3)
      real*8 ftm2r(3),ftm2ri(3)
      real*8 ttm2(3),ttm3(3)
      real*8 ttm2i(3),ttm3i(3)
      real*8 ttm2r(3),ttm3r(3)
      real*8 ttm2ri(3),ttm3ri(3)
      real*8 fdir(3),dixdk(3)
      real*8 dkxui(3),dixuk(3)
      real*8 dixukp(3),dkxuip(3)
      real*8 uixqkr(3),ukxqir(3)
      real*8 uixqkrp(3),ukxqirp(3)
      real*8 qiuk(3),qkui(3)
      real*8 qiukp(3),qkuip(3)
      real*8 rxqiuk(3),rxqkui(3)
      real*8 rxqiukp(3),rxqkuip(3)
      real*8 qidk(3),qkdi(3)
      real*8 qir(3),qkr(3)
      real*8 qiqkr(3),qkqir(3)
      real*8 qixqk(3),rxqir(3)
      real*8 dixr(3),dkxr(3)
      real*8 dixqkr(3),dkxqir(3)
      real*8 rxqkr(3),qkrxqir(3)
      real*8 rxqikr(3),rxqkir(3)
      real*8 rxqidk(3),rxqkdi(3)
      real*8 ddsc3(3),ddsc5(3)
      real*8 ddsc7(3)
      real*8 bn(0:5)
      real*8 sc(10),gl(0:8)
      real*8 sci(8),scip(8)
      real*8 gli(7),glip(7)
      real*8 gf(7),gfi(6)
      real*8 gfr(7),gfri(6)
      real*8 gti(6),gtri(6)
      real*8, allocatable :: mscale(:)
      real*8, allocatable :: pscale(:)
      real*8, allocatable :: dscale(:)
      real*8, allocatable :: uscale(:)
      real*8 aewald,aewald3b

      character*6 mode
      external erfc
    
      aewald=aewald3b
c
c
c     zero out the intramolecular portion of the Ewald energy
c
      eintra = 0.0d0
      if (npole .eq. 0)  return
c
c     perform dynamic allocation of some local arrays
c
      allocate (mscale(n))
      allocate (pscale(n))
      allocate (dscale(n))
      allocate (uscale(n))
c
c     set arrays needed to scale connected atom interactions
c
      do i = 1, n
         mscale(i) = 1.0d0
         pscale(i) = 1.0d0
         dscale(i) = 1.0d0
         uscale(i) = 1.0d0
      end do
c
c     set conversion factor, cutoff and switching coefficients
c
      f = electric / dielec
      mode = 'EWALD'
      call switch (mode)
c
c     set the permanent multipole and induced dipole values
c
      do l1 = 1, npole3b-1
         i = pnum(l1)
         ii = ipole(i)
         pdi = pdamp(i)
         pti = thole(i)
         ci = rpole(1,i)
         di(1) = rpole(2,i)
         di(2) = rpole(3,i)
         di(3) = rpole(4,i)
         qi(1) = rpole(5,i)
         qi(2) = rpole(6,i)
         qi(3) = rpole(7,i)
         qi(4) = rpole(8,i)
         qi(5) = rpole(9,i)
         qi(6) = rpole(10,i)
         qi(7) = rpole(11,i)
         qi(8) = rpole(12,i)
         qi(9) = rpole(13,i)
c
c     set interaction scaling coefficients for connected atoms
c
         do j = 1, n12(ii)
            mscale(i12(j,ii)) = m2scale
            pscale(i12(j,ii)) = p2scale
         end do
         do j = 1, n13(ii)
            mscale(i13(j,ii)) = m3scale
            pscale(i13(j,ii)) = p3scale
         end do
         do j = 1, n14(ii)
            mscale(i14(j,ii)) = m4scale
            pscale(i14(j,ii)) = p4scale
            do k = 1, np11(ii)
                if (i14(j,ii) .eq. ip11(k,ii))
     &            pscale(i14(j,ii)) = p4scale * p41scale
            end do
         end do
         do j = 1, n15(ii)
            mscale(i15(j,ii)) = m5scale
            pscale(i15(j,ii)) = p5scale
         end do
         do j = 1, np11(ii)
            dscale(ip11(j,ii)) = d1scale
            uscale(ip11(j,ii)) = u1scale
         end do
         do j = 1, np12(ii)
            dscale(ip12(j,ii)) = d2scale
            uscale(ip12(j,ii)) = u2scale
         end do
         do j = 1, np13(ii)
            dscale(ip13(j,ii)) = d3scale
            uscale(ip13(j,ii)) = u3scale
         end do
         do j = 1, np14(ii)
            dscale(ip14(j,ii)) = d4scale
            uscale(ip14(j,ii)) = u4scale
         end do
         do l3 = l1+1, npole3b
            k = pnum(l3)
            kk = ipole(k)
            xr = x(kk) - x(ii)
            yr = y(kk) - y(ii)
            zr = z(kk) - z(ii)
            call image (xr,yr,zr)
            r2 = xr*xr + yr*yr + zr*zr
            if (r2 .le. off2) then
               r = sqrt(r2)
               ck = rpole(1,k)
               dk(1) = rpole(2,k)
               dk(2) = rpole(3,k)
               dk(3) = rpole(4,k)
               qk(1) = rpole(5,k)
               qk(2) = rpole(6,k)
               qk(3) = rpole(7,k)
               qk(4) = rpole(8,k)
               qk(5) = rpole(9,k)
               qk(6) = rpole(10,k)
               qk(7) = rpole(11,k)
               qk(8) = rpole(12,k)
               qk(9) = rpole(13,k)
c
c     calculate the real space error function terms
c
               ralpha = aewald * r
               bn(0) = erfc(ralpha) / r
               alsq2 = 2.0d0 * aewald**2
               alsq2n = 0.0d0
               if (aewald .gt. 0.0d0)  alsq2n = 1.0d0 / (sqrtpi*aewald)
               exp2a = exp(-ralpha**2)
               do j = 1, 5
                  bfac = dble(2*j-1)
                  alsq2n = alsq2 * alsq2n
                  bn(j) = (bfac*bn(j-1)+alsq2n*exp2a) / r2
               end do
c
c     apply Thole polarization damping to scale factors
c
               rr1 = 1.0d0 / r
               rr3 = rr1 / r2
               rr5 = 3.0d0 * rr3 / r2
               rr7 = 5.0d0 * rr5 / r2
               rr9 = 7.0d0 * rr7 / r2
               rr11 = 9.0d0 * rr9 / r2
               scale3 = 1.0d0
               scale5 = 1.0d0
               scale7 = 1.0d0
               do j = 1, 3
                  ddsc3(j) = 0.0d0
                  ddsc5(j) = 0.0d0
                  ddsc7(j) = 0.0d0
               end do
               damp = pdi * pdamp(k)
               if (damp .ne. 0.0d0) then
                  pgamma = min(pti,thole(k))
                  damp = -pgamma * (r/damp)**3
                  if (damp .gt. -50.0d0) then
                     expdamp = exp(damp)
                     scale3 = 1.0d0 - expdamp
                     scale5 = 1.0d0 - (1.0d0-damp)*expdamp
                     scale7 = 1.0d0 - (1.0d0-damp+0.6d0*damp**2)
     &                                       *expdamp
                     temp3 = -3.0d0 * damp * expdamp / r2
                     temp5 = -damp
                     temp7 = -0.2d0 - 0.6d0*damp
                     ddsc3(1) = temp3 * xr
                     ddsc3(2) = temp3 * yr
                     ddsc3(3) = temp3 * zr
                     ddsc5(1) = temp5 * ddsc3(1)
                     ddsc5(2) = temp5 * ddsc3(2)
                     ddsc5(3) = temp5 * ddsc3(3)
                     ddsc7(1) = temp7 * ddsc5(1)
                     ddsc7(2) = temp7 * ddsc5(2)
                     ddsc7(3) = temp7 * ddsc5(3)
                  end if
               end if
               dsc3 = 1.0d0 - scale3*dscale(kk)
               dsc5 = 1.0d0 - scale5*dscale(kk)
               dsc7 = 1.0d0 - scale7*dscale(kk)
               psc3 = 1.0d0 - scale3*pscale(kk)
               psc5 = 1.0d0 - scale5*pscale(kk)
               psc7 = 1.0d0 - scale7*pscale(kk)
               usc3 = 1.0d0 - scale3*uscale(kk)
               usc5 = 1.0d0 - scale5*uscale(kk)
c
c     construct necessary auxiliary vectors
c
               dixdk(1) = di(2)*dk(3) - di(3)*dk(2)
               dixdk(2) = di(3)*dk(1) - di(1)*dk(3)
               dixdk(3) = di(1)*dk(2) - di(2)*dk(1)
               dixuk(1) = di(2)*uind(3,k) - di(3)*uind(2,k)
               dixuk(2) = di(3)*uind(1,k) - di(1)*uind(3,k)
               dixuk(3) = di(1)*uind(2,k) - di(2)*uind(1,k)
               dkxui(1) = dk(2)*uind(3,i) - dk(3)*uind(2,i)
               dkxui(2) = dk(3)*uind(1,i) - dk(1)*uind(3,i)
               dkxui(3) = dk(1)*uind(2,i) - dk(2)*uind(1,i)
               dixukp(1) = di(2)*uinp(3,k) - di(3)*uinp(2,k)
               dixukp(2) = di(3)*uinp(1,k) - di(1)*uinp(3,k)
               dixukp(3) = di(1)*uinp(2,k) - di(2)*uinp(1,k)
               dkxuip(1) = dk(2)*uinp(3,i) - dk(3)*uinp(2,i)
               dkxuip(2) = dk(3)*uinp(1,i) - dk(1)*uinp(3,i)
               dkxuip(3) = dk(1)*uinp(2,i) - dk(2)*uinp(1,i)
               dixr(1) = di(2)*zr - di(3)*yr
               dixr(2) = di(3)*xr - di(1)*zr
               dixr(3) = di(1)*yr - di(2)*xr
               dkxr(1) = dk(2)*zr - dk(3)*yr
               dkxr(2) = dk(3)*xr - dk(1)*zr
               dkxr(3) = dk(1)*yr - dk(2)*xr
               qir(1) = qi(1)*xr + qi(4)*yr + qi(7)*zr
               qir(2) = qi(2)*xr + qi(5)*yr + qi(8)*zr
               qir(3) = qi(3)*xr + qi(6)*yr + qi(9)*zr
               qkr(1) = qk(1)*xr + qk(4)*yr + qk(7)*zr
               qkr(2) = qk(2)*xr + qk(5)*yr + qk(8)*zr
               qkr(3) = qk(3)*xr + qk(6)*yr + qk(9)*zr
               qiqkr(1) = qi(1)*qkr(1) + qi(4)*qkr(2) + qi(7)*qkr(3)
               qiqkr(2) = qi(2)*qkr(1) + qi(5)*qkr(2) + qi(8)*qkr(3)
               qiqkr(3) = qi(3)*qkr(1) + qi(6)*qkr(2) + qi(9)*qkr(3)
               qkqir(1) = qk(1)*qir(1) + qk(4)*qir(2) + qk(7)*qir(3)
               qkqir(2) = qk(2)*qir(1) + qk(5)*qir(2) + qk(8)*qir(3)
               qkqir(3) = qk(3)*qir(1) + qk(6)*qir(2) + qk(9)*qir(3)
               qixqk(1) = qi(2)*qk(3) + qi(5)*qk(6) + qi(8)*qk(9)
     &                       - qi(3)*qk(2) - qi(6)*qk(5) - qi(9)*qk(8)
               qixqk(2) = qi(3)*qk(1) + qi(6)*qk(4) + qi(9)*qk(7)
     &                       - qi(1)*qk(3) - qi(4)*qk(6) - qi(7)*qk(9)
               qixqk(3) = qi(1)*qk(2) + qi(4)*qk(5) + qi(7)*qk(8)
     &                       - qi(2)*qk(1) - qi(5)*qk(4) - qi(8)*qk(7)
               rxqir(1) = yr*qir(3) - zr*qir(2)
               rxqir(2) = zr*qir(1) - xr*qir(3)
               rxqir(3) = xr*qir(2) - yr*qir(1)
               rxqkr(1) = yr*qkr(3) - zr*qkr(2)
               rxqkr(2) = zr*qkr(1) - xr*qkr(3)
               rxqkr(3) = xr*qkr(2) - yr*qkr(1)
               rxqikr(1) = yr*qiqkr(3) - zr*qiqkr(2)
               rxqikr(2) = zr*qiqkr(1) - xr*qiqkr(3)
               rxqikr(3) = xr*qiqkr(2) - yr*qiqkr(1)
               rxqkir(1) = yr*qkqir(3) - zr*qkqir(2)
               rxqkir(2) = zr*qkqir(1) - xr*qkqir(3)
               rxqkir(3) = xr*qkqir(2) - yr*qkqir(1)
               qkrxqir(1) = qkr(2)*qir(3) - qkr(3)*qir(2)
               qkrxqir(2) = qkr(3)*qir(1) - qkr(1)*qir(3)
               qkrxqir(3) = qkr(1)*qir(2) - qkr(2)*qir(1)
               qidk(1) = qi(1)*dk(1) + qi(4)*dk(2) + qi(7)*dk(3)
               qidk(2) = qi(2)*dk(1) + qi(5)*dk(2) + qi(8)*dk(3)
               qidk(3) = qi(3)*dk(1) + qi(6)*dk(2) + qi(9)*dk(3)
               qkdi(1) = qk(1)*di(1) + qk(4)*di(2) + qk(7)*di(3)
               qkdi(2) = qk(2)*di(1) + qk(5)*di(2) + qk(8)*di(3)
               qkdi(3) = qk(3)*di(1) + qk(6)*di(2) + qk(9)*di(3)
               qiuk(1) = qi(1)*uind(1,k) + qi(4)*uind(2,k)
     &                      + qi(7)*uind(3,k)
               qiuk(2) = qi(2)*uind(1,k) + qi(5)*uind(2,k)
     &                      + qi(8)*uind(3,k)
               qiuk(3) = qi(3)*uind(1,k) + qi(6)*uind(2,k)
     &                      + qi(9)*uind(3,k)
               qkui(1) = qk(1)*uind(1,i) + qk(4)*uind(2,i)
     &                      + qk(7)*uind(3,i)
               qkui(2) = qk(2)*uind(1,i) + qk(5)*uind(2,i)
     &                      + qk(8)*uind(3,i)
               qkui(3) = qk(3)*uind(1,i) + qk(6)*uind(2,i)
     &                      + qk(9)*uind(3,i)
               qiukp(1) = qi(1)*uinp(1,k) + qi(4)*uinp(2,k)
     &                       + qi(7)*uinp(3,k)
               qiukp(2) = qi(2)*uinp(1,k) + qi(5)*uinp(2,k)
     &                       + qi(8)*uinp(3,k)
               qiukp(3) = qi(3)*uinp(1,k) + qi(6)*uinp(2,k)
     &                       + qi(9)*uinp(3,k)
               qkuip(1) = qk(1)*uinp(1,i) + qk(4)*uinp(2,i)
     &                       + qk(7)*uinp(3,i)
               qkuip(2) = qk(2)*uinp(1,i) + qk(5)*uinp(2,i)
     &                       + qk(8)*uinp(3,i)
               qkuip(3) = qk(3)*uinp(1,i) + qk(6)*uinp(2,i)
     &                       + qk(9)*uinp(3,i)
               dixqkr(1) = di(2)*qkr(3) - di(3)*qkr(2)
               dixqkr(2) = di(3)*qkr(1) - di(1)*qkr(3)
               dixqkr(3) = di(1)*qkr(2) - di(2)*qkr(1)
               dkxqir(1) = dk(2)*qir(3) - dk(3)*qir(2)
               dkxqir(2) = dk(3)*qir(1) - dk(1)*qir(3)
               dkxqir(3) = dk(1)*qir(2) - dk(2)*qir(1)
               uixqkr(1) = uind(2,i)*qkr(3) - uind(3,i)*qkr(2)
               uixqkr(2) = uind(3,i)*qkr(1) - uind(1,i)*qkr(3)
               uixqkr(3) = uind(1,i)*qkr(2) - uind(2,i)*qkr(1)
               ukxqir(1) = uind(2,k)*qir(3) - uind(3,k)*qir(2)
               ukxqir(2) = uind(3,k)*qir(1) - uind(1,k)*qir(3)
               ukxqir(3) = uind(1,k)*qir(2) - uind(2,k)*qir(1)
               uixqkrp(1) = uinp(2,i)*qkr(3) - uinp(3,i)*qkr(2)
               uixqkrp(2) = uinp(3,i)*qkr(1) - uinp(1,i)*qkr(3)
               uixqkrp(3) = uinp(1,i)*qkr(2) - uinp(2,i)*qkr(1)
               ukxqirp(1) = uinp(2,k)*qir(3) - uinp(3,k)*qir(2)
               ukxqirp(2) = uinp(3,k)*qir(1) - uinp(1,k)*qir(3)
               ukxqirp(3) = uinp(1,k)*qir(2) - uinp(2,k)*qir(1)
               rxqidk(1) = yr*qidk(3) - zr*qidk(2)
               rxqidk(2) = zr*qidk(1) - xr*qidk(3)
               rxqidk(3) = xr*qidk(2) - yr*qidk(1)
               rxqkdi(1) = yr*qkdi(3) - zr*qkdi(2)
               rxqkdi(2) = zr*qkdi(1) - xr*qkdi(3)
               rxqkdi(3) = xr*qkdi(2) - yr*qkdi(1)
               rxqiuk(1) = yr*qiuk(3) - zr*qiuk(2)
               rxqiuk(2) = zr*qiuk(1) - xr*qiuk(3)
               rxqiuk(3) = xr*qiuk(2) - yr*qiuk(1)
               rxqkui(1) = yr*qkui(3) - zr*qkui(2)
               rxqkui(2) = zr*qkui(1) - xr*qkui(3)
               rxqkui(3) = xr*qkui(2) - yr*qkui(1)
               rxqiukp(1) = yr*qiukp(3) - zr*qiukp(2)
               rxqiukp(2) = zr*qiukp(1) - xr*qiukp(3)
               rxqiukp(3) = xr*qiukp(2) - yr*qiukp(1)
               rxqkuip(1) = yr*qkuip(3) - zr*qkuip(2)
               rxqkuip(2) = zr*qkuip(1) - xr*qkuip(3)
               rxqkuip(3) = xr*qkuip(2) - yr*qkuip(1)
c
c     calculate the scalar products for permanent components
c
               sc(2) = di(1)*dk(1) + di(2)*dk(2) + di(3)*dk(3)
               sc(3) = di(1)*xr + di(2)*yr + di(3)*zr
               sc(4) = dk(1)*xr + dk(2)*yr + dk(3)*zr
               sc(5) = qir(1)*xr + qir(2)*yr + qir(3)*zr
               sc(6) = qkr(1)*xr + qkr(2)*yr + qkr(3)*zr
               sc(7) = qir(1)*dk(1) + qir(2)*dk(2) + qir(3)*dk(3)
               sc(8) = qkr(1)*di(1) + qkr(2)*di(2) + qkr(3)*di(3)
               sc(9) = qir(1)*qkr(1) + qir(2)*qkr(2) + qir(3)*qkr(3)
               sc(10) = qi(1)*qk(1) + qi(2)*qk(2) + qi(3)*qk(3)
     &                     + qi(4)*qk(4) + qi(5)*qk(5) + qi(6)*qk(6)
     &                     + qi(7)*qk(7) + qi(8)*qk(8) + qi(9)*qk(9)
c
c     calculate the scalar products for induced components
c
               sci(1) = uind(1,i)*dk(1) + uind(2,i)*dk(2)
     &                     + uind(3,i)*dk(3) + di(1)*uind(1,k)
     &                     + di(2)*uind(2,k) + di(3)*uind(3,k)
               sci(2) = uind(1,i)*uind(1,k) + uind(2,i)*uind(2,k)
     &                     + uind(3,i)*uind(3,k)
               sci(3) = uind(1,i)*xr + uind(2,i)*yr + uind(3,i)*zr
               sci(4) = uind(1,k)*xr + uind(2,k)*yr + uind(3,k)*zr
               sci(7) = qir(1)*uind(1,k) + qir(2)*uind(2,k)
     &                     + qir(3)*uind(3,k)
               sci(8) = qkr(1)*uind(1,i) + qkr(2)*uind(2,i)
     &                     + qkr(3)*uind(3,i)
               scip(1) = uinp(1,i)*dk(1) + uinp(2,i)*dk(2)
     &                      + uinp(3,i)*dk(3) + di(1)*uinp(1,k)
     &                      + di(2)*uinp(2,k) + di(3)*uinp(3,k)
               scip(2) = uind(1,i)*uinp(1,k)+uind(2,i)*uinp(2,k)
     &                      + uind(3,i)*uinp(3,k)+uinp(1,i)*uind(1,k)
     &                      + uinp(2,i)*uind(2,k)+uinp(3,i)*uind(3,k)
               scip(3) = uinp(1,i)*xr + uinp(2,i)*yr + uinp(3,i)*zr
               scip(4) = uinp(1,k)*xr + uinp(2,k)*yr + uinp(3,k)*zr
               scip(7) = qir(1)*uinp(1,k) + qir(2)*uinp(2,k)
     &                      + qir(3)*uinp(3,k)
               scip(8) = qkr(1)*uinp(1,i) + qkr(2)*uinp(2,i)
     &                      + qkr(3)*uinp(3,i)
c
c     calculate the gl functions for permanent components
c
               gl(0) = ci*ck
               gl(1) = ck*sc(3) - ci*sc(4)
               gl(2) = ci*sc(6) + ck*sc(5) - sc(3)*sc(4)
               gl(3) = sc(3)*sc(6) - sc(4)*sc(5)
               gl(4) = sc(5)*sc(6)
               gl(5) = -4.0d0 * sc(9)
               gl(6) = sc(2)
               gl(7) = 2.0d0 * (sc(7)-sc(8))
               gl(8) = 2.0d0 * sc(10)
c
c     calculate the gl functions for induced components
c
               gli(1) = ck*sci(3) - ci*sci(4)
               gli(2) = -sc(3)*sci(4) - sci(3)*sc(4)
               gli(3) = sci(3)*sc(6) - sci(4)*sc(5)
               gli(6) = sci(1)
               gli(7) = 2.0d0 * (sci(7)-sci(8))
               glip(1) = ck*scip(3) - ci*scip(4)
               glip(2) = -sc(3)*scip(4) - scip(3)*sc(4)
               glip(3) = scip(3)*sc(6) - scip(4)*sc(5)
               glip(6) = scip(1)
               glip(7) = 2.0d0 * (scip(7)-scip(8))
c
c     compute the energy contributions for this interaction
c
               ei = 0.5d0 * (bn(1)*(gli(1)+gli(6))
     &                      + bn(2)*(gli(2)+gli(7)) + bn(3)*gli(3))
c
c     get the real energy without any screening function
c
               erli = 0.5d0*(rr3*(gli(1)+gli(6))*psc3
     &                   + rr5*(gli(2)+gli(7))*psc5
     &                   + rr7*gli(3)*psc7)
               ei = ei - erli
               ei = f * ei
               ep = ep + ei
c
c     increment the total intramolecular energy; assumes
c     intramolecular distances are less than half of cell
c     length and less than the ewald cutoff
c
               if (molcule(ii) .eq. molcule(kk)) then
c                  eintra = eintra + 0.5d0*pscale(kk)
c     &                        * (rr3*(gli(1)+gli(6))*scale3
c     &                              + rr5*(gli(2)+gli(7))*scale5
c     &                              + rr7*gli(3)*scale7)
               end if
c
c     intermediate variables for induced force terms
c
               gfi(1) = 0.5d0*bn(2)*(gli(1)+glip(1)+gli(6)+glip(6))
     &                     + 0.5d0*bn(2)*scip(2)
     &                     + 0.5d0*bn(3)*(gli(2)+glip(2)+gli(7)+glip(7))
     &                     - 0.5d0*bn(3)*(sci(3)*scip(4)+scip(3)*sci(4))
     &                     + 0.5d0*bn(4)*(gli(3)+glip(3))
               gfi(2) = -ck*bn(1) + sc(4)*bn(2) - sc(6)*bn(3)
               gfi(3) = ci*bn(1) + sc(3)*bn(2) + sc(5)*bn(3)
               gfi(4) = 2.0d0 * bn(2)
               gfi(5) = bn(3) * (sci(4)+scip(4))
               gfi(6) = -bn(3) * (sci(3)+scip(3))
               gfri(1) = 0.5d0*rr5*((gli(1)+gli(6))*psc3
     &                            + (glip(1)+glip(6))*dsc3
     &                            + scip(2)*usc3)
     &                 + 0.5d0*rr7*((gli(7)+gli(2))*psc5
     &                            + (glip(7)+glip(2))*dsc5
     &                     - (sci(3)*scip(4)+scip(3)*sci(4))*usc5)
     &                 + 0.5d0*rr9*(gli(3)*psc7+glip(3)*dsc7)
               gfri(2) = -rr3*ck + rr5*sc(4) - rr7*sc(6)
               gfri(3) = rr3*ci + rr5*sc(3) + rr7*sc(5)
               gfri(4) = 2.0d0 * rr5
               gfri(5) = rr7 * (sci(4)*psc7+scip(4)*dsc7)
               gfri(6) = -rr7 * (sci(3)*psc7+scip(3)*dsc7)
c
c     get the induced force with screening
c
               ftm2i(1) = gfi(1)*xr + 0.5d0*
     &             (gfi(2)*(uind(1,i)+uinp(1,i))
     &            + bn(2)*(sci(4)*uinp(1,i)+scip(4)*uind(1,i))
     &            + gfi(3)*(uind(1,k)+uinp(1,k))
     &            + bn(2)*(sci(3)*uinp(1,k)+scip(3)*uind(1,k))
     &            + (sci(4)+scip(4))*bn(2)*di(1)
     &            + (sci(3)+scip(3))*bn(2)*dk(1)
     &            + gfi(4)*(qkui(1)+qkuip(1)-qiuk(1)-qiukp(1)))
     &            + gfi(5)*qir(1) + gfi(6)*qkr(1)
               ftm2i(2) = gfi(1)*yr + 0.5d0*
     &             (gfi(2)*(uind(2,i)+uinp(2,i))
     &            + bn(2)*(sci(4)*uinp(2,i)+scip(4)*uind(2,i))
     &            + gfi(3)*(uind(2,k)+uinp(2,k))
     &            + bn(2)*(sci(3)*uinp(2,k)+scip(3)*uind(2,k))
     &            + (sci(4)+scip(4))*bn(2)*di(2)
     &            + (sci(3)+scip(3))*bn(2)*dk(2)
     &            + gfi(4)*(qkui(2)+qkuip(2)-qiuk(2)-qiukp(2)))
     &            + gfi(5)*qir(2) + gfi(6)*qkr(2)
               ftm2i(3) = gfi(1)*zr + 0.5d0*
     &             (gfi(2)*(uind(3,i)+uinp(3,i))
     &            + bn(2)*(sci(4)*uinp(3,i)+scip(4)*uind(3,i))
     &            + gfi(3)*(uind(3,k)+uinp(3,k))
     &            + bn(2)*(sci(3)*uinp(3,k)+scip(3)*uind(3,k))
     &            + (sci(4)+scip(4))*bn(2)*di(3)
     &            + (sci(3)+scip(3))*bn(2)*dk(3)
     &            + gfi(4)*(qkui(3)+qkuip(3)-qiuk(3)-qiukp(3)))
     &            + gfi(5)*qir(3) + gfi(6)*qkr(3)
c
c     get the induced force without screening
c
               ftm2ri(1) = gfri(1)*xr + 0.5d0*
     &           (- rr3*ck*(uind(1,i)*psc3+uinp(1,i)*dsc3)
     &            + rr5*sc(4)*(uind(1,i)*psc5+uinp(1,i)*dsc5)
     &            - rr7*sc(6)*(uind(1,i)*psc7+uinp(1,i)*dsc7))
     &            + (rr3*ci*(uind(1,k)*psc3+uinp(1,k)*dsc3)
     &            + rr5*sc(3)*(uind(1,k)*psc5+uinp(1,k)*dsc5)
     &            + rr7*sc(5)*(uind(1,k)*psc7+uinp(1,k)*dsc7))*0.5d0
     &            + rr5*usc5*(sci(4)*uinp(1,i)+scip(4)*uind(1,i)
     &            + sci(3)*uinp(1,k)+scip(3)*uind(1,k))*0.5d0
     &            + 0.5d0*(sci(4)*psc5+scip(4)*dsc5)*rr5*di(1)
     &            + 0.5d0*(sci(3)*psc5+scip(3)*dsc5)*rr5*dk(1)
     &            + 0.5d0*gfri(4)*((qkui(1)-qiuk(1))*psc5
     &            + (qkuip(1)-qiukp(1))*dsc5)
     &            + gfri(5)*qir(1) + gfri(6)*qkr(1)
               ftm2ri(2) = gfri(1)*yr + 0.5d0*
     &           (- rr3*ck*(uind(2,i)*psc3+uinp(2,i)*dsc3)
     &            + rr5*sc(4)*(uind(2,i)*psc5+uinp(2,i)*dsc5)
     &            - rr7*sc(6)*(uind(2,i)*psc7+uinp(2,i)*dsc7))
     &            + (rr3*ci*(uind(2,k)*psc3+uinp(2,k)*dsc3)
     &            + rr5*sc(3)*(uind(2,k)*psc5+uinp(2,k)*dsc5)
     &            + rr7*sc(5)*(uind(2,k)*psc7+uinp(2,k)*dsc7))*0.5d0
     &            + rr5*usc5*(sci(4)*uinp(2,i)+scip(4)*uind(2,i)
     &            + sci(3)*uinp(2,k)+scip(3)*uind(2,k))*0.5d0
     &            + 0.5d0*(sci(4)*psc5+scip(4)*dsc5)*rr5*di(2)
     &            + 0.5d0*(sci(3)*psc5+scip(3)*dsc5)*rr5*dk(2)
     &            + 0.5d0*gfri(4)*((qkui(2)-qiuk(2))*psc5
     &            + (qkuip(2)-qiukp(2))*dsc5)
     &            + gfri(5)*qir(2) + gfri(6)*qkr(2)
               ftm2ri(3) = gfri(1)*zr + 0.5d0*
     &           (- rr3*ck*(uind(3,i)*psc3+uinp(3,i)*dsc3)
     &            + rr5*sc(4)*(uind(3,i)*psc5+uinp(3,i)*dsc5)
     &            - rr7*sc(6)*(uind(3,i)*psc7+uinp(3,i)*dsc7))
     &            + (rr3*ci*(uind(3,k)*psc3+uinp(3,k)*dsc3)
     &            + rr5*sc(3)*(uind(3,k)*psc5+uinp(3,k)*dsc5)
     &            + rr7*sc(5)*(uind(3,k)*psc7+uinp(3,k)*dsc7))*0.5d0
     &            + rr5*usc5*(sci(4)*uinp(3,i)+scip(4)*uind(3,i)
     &            + sci(3)*uinp(3,k)+scip(3)*uind(3,k))*0.5d0
     &            + 0.5d0*(sci(4)*psc5+scip(4)*dsc5)*rr5*di(3)
     &            + 0.5d0*(sci(3)*psc5+scip(3)*dsc5)*rr5*dk(3)
     &            + 0.5d0*gfri(4)*((qkui(3)-qiuk(3))*psc5
     &            + (qkuip(3)-qiukp(3))*dsc5)
     &            + gfri(5)*qir(3) + gfri(6)*qkr(3)
c
c     account for partially excluded induced interactions
c
               temp3 = 0.5d0 * rr3 * ((gli(1)+gli(6))*pscale(kk)
     &                                  +(glip(1)+glip(6))*dscale(kk))
               temp5 = 0.5d0 * rr5 * ((gli(2)+gli(7))*pscale(kk)
     &                                  +(glip(2)+glip(7))*dscale(kk))
               temp7 = 0.5d0 * rr7 * (gli(3)*pscale(kk)
     &                                  +glip(3)*dscale(kk))
               fridmp(1) = temp3*ddsc3(1) + temp5*ddsc5(1)
     &                        + temp7*ddsc7(1)
               fridmp(2) = temp3*ddsc3(2) + temp5*ddsc5(2)
     &                        + temp7*ddsc7(2)
               fridmp(3) = temp3*ddsc3(3) + temp5*ddsc5(3)
     &                        + temp7*ddsc7(3)
c
c     find some scaling terms for induced-induced force
c
               temp3 = 0.5d0 * rr3 * uscale(kk) * scip(2)
               temp5 = -0.5d0 * rr5 * uscale(kk)
     &                    * (sci(3)*scip(4)+scip(3)*sci(4))
               findmp(1) = temp3*ddsc3(1) + temp5*ddsc5(1)
               findmp(2) = temp3*ddsc3(2) + temp5*ddsc5(2)
               findmp(3) = temp3*ddsc3(3) + temp5*ddsc5(3)
c
c     modify the forces for partially excluded interactions
c
               ftm2i(1) = ftm2i(1) - fridmp(1) - findmp(1)
               ftm2i(2) = ftm2i(2) - fridmp(2) - findmp(2)
               ftm2i(3) = ftm2i(3) - fridmp(3) - findmp(3)
c
c     correction to convert mutual to direct polarization force
c
               if (poltyp .eq. 'DIRECT') then
                  gfd = 0.5d0 * (bn(2)*scip(2)
     &                     - bn(3)*(scip(3)*sci(4)+sci(3)*scip(4)))
                  gfdr = 0.5d0 * (rr5*scip(2)*usc3
     &                     - rr7*(scip(3)*sci(4)
     &                           +sci(3)*scip(4))*usc5)

                  ftm2i(1) = ftm2i(1) - gfd*xr - 0.5d0*bn(2)*
     &                          (sci(4)*uinp(1,i)+scip(4)*uind(1,i)
     &                          +sci(3)*uinp(1,k)+scip(3)*uind(1,k))
                  ftm2i(2) = ftm2i(2) - gfd*yr - 0.5d0*bn(2)*
     &                          (sci(4)*uinp(2,i)+scip(4)*uind(2,i)
     &                          +sci(3)*uinp(2,k)+scip(3)*uind(2,k))
                  ftm2i(3) = ftm2i(3) - gfd*zr - 0.5d0*bn(2)*
     &                          (sci(4)*uinp(3,i)+scip(4)*uind(3,i)
     &                          +sci(3)*uinp(3,k)+scip(3)*uind(3,k))
                  fdir(1) = gfdr*xr + 0.5d0*usc5*rr5*
     &                         (sci(4)*uinp(1,i)+scip(4)*uind(1,i)
     &                        + sci(3)*uinp(1,k)+scip(3)*uind(1,k))
                  fdir(2) = gfdr*yr + 0.5d0*usc5*rr5*
     &                         (sci(4)*uinp(2,i)+scip(4)*uind(2,i)
     &                        + sci(3)*uinp(2,k)+scip(3)*uind(2,k))
                  fdir(3) = gfdr*zr + 0.5d0*usc5*rr5*
     &                         (sci(4)*uinp(3,i)+scip(4)*uind(3,i)
     &                        + sci(3)*uinp(3,k)+scip(3)*uind(3,k))
                  ftm2i(1) = ftm2i(1) + fdir(1) + findmp(1)
                  ftm2i(2) = ftm2i(2) + fdir(2) + findmp(2)
                  ftm2i(3) = ftm2i(3) + fdir(3) + findmp(3)
               end if
c
c     intermediate variables for induced torque terms
c
               gti(2) = 0.5d0 * bn(2) * (sci(4)+scip(4))
               gti(3) = 0.5d0 * bn(2) * (sci(3)+scip(3))
               gti(4) = gfi(4)
               gti(5) = gfi(5)
               gti(6) = gfi(6)
               gtri(2) = 0.5d0 * rr5 * (sci(4)*psc5+scip(4)*dsc5)
               gtri(3) = 0.5d0 * rr5 * (sci(3)*psc5+scip(3)*dsc5)
               gtri(4) = gfri(4)
               gtri(5) = gfri(5)
               gtri(6) = gfri(6)
c
c     get the induced torque with screening
c
               ttm2i(1) = -bn(1)*(dixuk(1)+dixukp(1))*0.5d0
     &           + gti(2)*dixr(1) + gti(4)*(ukxqir(1)+rxqiuk(1)
     &           + ukxqirp(1)+rxqiukp(1))*0.5d0 - gti(5)*rxqir(1)
               ttm2i(2) = -bn(1)*(dixuk(2)+dixukp(2))*0.5d0
     &           + gti(2)*dixr(2) + gti(4)*(ukxqir(2)+rxqiuk(2)
     &           + ukxqirp(2)+rxqiukp(2))*0.5d0 - gti(5)*rxqir(2)
               ttm2i(3) = -bn(1)*(dixuk(3)+dixukp(3))*0.5d0
     &           + gti(2)*dixr(3) + gti(4)*(ukxqir(3)+rxqiuk(3)
     &           + ukxqirp(3)+rxqiukp(3))*0.5d0 - gti(5)*rxqir(3)
               ttm3i(1) = -bn(1)*(dkxui(1)+dkxuip(1))*0.5d0
     &           + gti(3)*dkxr(1) - gti(4)*(uixqkr(1)+rxqkui(1)
     &           + uixqkrp(1)+rxqkuip(1))*0.5d0 - gti(6)*rxqkr(1)
               ttm3i(2) = -bn(1)*(dkxui(2)+dkxuip(2))*0.5d0
     &           + gti(3)*dkxr(2) - gti(4)*(uixqkr(2)+rxqkui(2)
     &           + uixqkrp(2)+rxqkuip(2))*0.5d0 - gti(6)*rxqkr(2)
               ttm3i(3) = -bn(1)*(dkxui(3)+dkxuip(3))*0.5d0
     &           + gti(3)*dkxr(3) - gti(4)*(uixqkr(3)+rxqkui(3)
     &           + uixqkrp(3)+rxqkuip(3))*0.5d0 - gti(6)*rxqkr(3)
c
c     get the induced torque without screening
c
               ttm2ri(1) = -rr3*(dixuk(1)*psc3+dixukp(1)*dsc3)*0.5d0
     &           + gtri(2)*dixr(1) + gtri(4)*((ukxqir(1)+rxqiuk(1))*psc5
     &           +(ukxqirp(1)+rxqiukp(1))*dsc5)*0.5d0 - gtri(5)*rxqir(1)
               ttm2ri(2) = -rr3*(dixuk(2)*psc3+dixukp(2)*dsc3)*0.5d0
     &           + gtri(2)*dixr(2) + gtri(4)*((ukxqir(2)+rxqiuk(2))*psc5
     &           +(ukxqirp(2)+rxqiukp(2))*dsc5)*0.5d0 - gtri(5)*rxqir(2)
               ttm2ri(3) = -rr3*(dixuk(3)*psc3+dixukp(3)*dsc3)*0.5d0
     &           + gtri(2)*dixr(3) + gtri(4)*((ukxqir(3)+rxqiuk(3))*psc5
     &           +(ukxqirp(3)+rxqiukp(3))*dsc5)*0.5d0 - gtri(5)*rxqir(3)
               ttm3ri(1) = -rr3*(dkxui(1)*psc3+dkxuip(1)*dsc3)*0.5d0
     &           + gtri(3)*dkxr(1) - gtri(4)*((uixqkr(1)+rxqkui(1))*psc5
     &           +(uixqkrp(1)+rxqkuip(1))*dsc5)*0.5d0 - gtri(6)*rxqkr(1)
               ttm3ri(2) = -rr3*(dkxui(2)*psc3+dkxuip(2)*dsc3)*0.5d0
     &           + gtri(3)*dkxr(2) - gtri(4)*((uixqkr(2)+rxqkui(2))*psc5
     &           +(uixqkrp(2)+rxqkuip(2))*dsc5)*0.5d0 - gtri(6)*rxqkr(2)
               ttm3ri(3) = -rr3*(dkxui(3)*psc3+dkxuip(3)*dsc3)*0.5d0
     &           + gtri(3)*dkxr(3) - gtri(4)*((uixqkr(3)+rxqkui(3))*psc5
     &           +(uixqkrp(3)+rxqkuip(3))*dsc5)*0.5d0 - gtri(6)*rxqkr(3)
c
c     handle the case where scaling is used
c
               do j = 1, 3
                  ftm2i(j) = f * (ftm2i(j)-ftm2ri(j))
                  ttm2i(j) = f * (ttm2i(j)-ttm2ri(j))
                  ttm3i(j) = f * (ttm3i(j)-ttm3ri(j))
               end do
c
c     increment gradient due to force and torque on first site
c
               dep(1,ii) = dep(1,ii) + ftm2i(1)
               dep(2,ii) = dep(2,ii) + ftm2i(2)
               dep(3,ii) = dep(3,ii) + ftm2i(3)
               call torque_3b (i,ttm2,ttm2i,frcxi,frcyi,frczi)

c
c     increment gradient due to force and torque on second site
c
               dep(1,kk) = dep(1,kk) - ftm2i(1)
               dep(2,kk) = dep(2,kk) - ftm2i(2)
               dep(3,kk) = dep(3,kk) - ftm2i(3)
               call torque_3b (k,ttm3,ttm3i,frcxk,frcyk,frczk)
c
c     increment the internal virial tensor components
c
               iaz = zaxis(i)
               iax = xaxis(i)
               iay = yaxis(i)
               kaz = zaxis(k)
               kax = xaxis(k)
               kay = yaxis(k)
               if (iaz .eq. 0)  iaz = ii
               if (iax .eq. 0)  iax = ii
               if (iay .eq. 0)  iay = ii
               if (kaz .eq. 0)  kaz = kk
               if (kax .eq. 0)  kax = kk
               if (kay .eq. 0)  kay = kk
               xiz = x(iaz) - x(ii)
               yiz = y(iaz) - y(ii)
               ziz = z(iaz) - z(ii)
               xix = x(iax) - x(ii)
               yix = y(iax) - y(ii)
               zix = z(iax) - z(ii)
               xiy = x(iay) - x(ii)
               yiy = y(iay) - y(ii)
               ziy = z(iay) - z(ii)
               xkz = x(kaz) - x(kk)
               ykz = y(kaz) - y(kk)
               zkz = z(kaz) - z(kk)
               xkx = x(kax) - x(kk)
               ykx = y(kax) - y(kk)
               zkx = z(kax) - z(kk)
               xky = x(kay) - x(kk)
               yky = y(kay) - y(kk)
               zky = z(kay) - z(kk)
               vxx = -xr*ftm2i(1)+xix*frcxi(1)+xiy*frcyi(1)+xiz*frczi(1)
               vyx = -yr*ftm2i(1)+yix*frcxi(1)+yiy*frcyi(1)+yiz*frczi(1)
               vzx = -zr*ftm2i(1)+zix*frcxi(1)+ziy*frcyi(1)+ziz*frczi(1)
               vyy = -yr*ftm2i(2)+yix*frcxi(2)+yiy*frcyi(2)+yiz*frczi(2)
               vzy = -zr*ftm2i(2)+zix*frcxi(2)+ziy*frcyi(2)+ziz*frczi(2)
               vzz = -zr*ftm2i(3)+zix*frcxi(3)+ziy*frcyi(3)+ziz*frczi(3)
               vir(1,1) = vir(1,1) + vxx
               vir(2,1) = vir(2,1) + vyx
               vir(3,1) = vir(3,1) + vzx
               vir(1,2) = vir(1,2) + vyx
               vir(2,2) = vir(2,2) + vyy
               vir(3,2) = vir(3,2) + vzy
               vir(1,3) = vir(1,3) + vzx
               vir(2,3) = vir(2,3) + vzy
               vir(3,3) = vir(3,3) + vzz
            end if
         end do
c
c     reset interaction scaling coefficients for connected atoms
c
         do j = 1, n12(ii)
            mscale(i12(j,ii)) = 1.0d0
            pscale(i12(j,ii)) = 1.0d0
         end do
         do j = 1, n13(ii)
            mscale(i13(j,ii)) = 1.0d0
            pscale(i13(j,ii)) = 1.0d0
         end do
         do j = 1, n14(ii)
            mscale(i14(j,ii)) = 1.0d0
            pscale(i14(j,ii)) = 1.0d0
         end do
         do j = 1, n15(ii)
            mscale(i15(j,ii)) = 1.0d0
            pscale(i15(j,ii)) = 1.0d0
         end do
         do j = 1, np11(ii)
            dscale(ip11(j,ii)) = 1.0d0
            uscale(ip11(j,ii)) = 1.0d0
         end do
         do j = 1, np12(ii)
            dscale(ip12(j,ii)) = 1.0d0
            uscale(ip12(j,ii)) = 1.0d0
         end do
         do j = 1, np13(ii)
            dscale(ip13(j,ii)) = 1.0d0
            uscale(ip13(j,ii)) = 1.0d0
         end do
         do j = 1, np14(ii)
            dscale(ip14(j,ii)) = 1.0d0
            uscale(ip14(j,ii)) = 1.0d0
         end do
      end do
c
c     dor periodic boundary conditions with large cutoffs
c
      if (use_replica) then
c
c     calculate interactions with other unit cells
c
      do l1 = 1, npole3b
         i = pnum(l1)
         ii = ipole(i)
         pdi = pdamp(i)
         pti = thole(i)
         ci = rpole(1,i)
         di(1) = rpole(2,i)
         di(2) = rpole(3,i)
         di(3) = rpole(4,i)
         qi(1) = rpole(5,i)
         qi(2) = rpole(6,i)
         qi(3) = rpole(7,i)
         qi(4) = rpole(8,i)
         qi(5) = rpole(9,i)
         qi(6) = rpole(10,i)
         qi(7) = rpole(11,i)
         qi(8) = rpole(12,i)
         qi(9) = rpole(13,i)
c
c     set interaction scaling coefficients for connected atoms
c
         do j = 1, n12(ii)
            mscale(i12(j,ii)) = m2scale
            pscale(i12(j,ii)) = p2scale
         end do
         do j = 1, n13(ii)
            mscale(i13(j,ii)) = m3scale
            pscale(i13(j,ii)) = p3scale
         end do
         do j = 1, n14(ii)
            mscale(i14(j,ii)) = m4scale
            pscale(i14(j,ii)) = p4scale
            do k = 1, np11(ii)
                if (i14(j,ii) .eq. ip11(k,ii))
     &            pscale(i14(j,ii)) = p4scale * p41scale
            end do
         end do
         do j = 1, n15(ii)
            mscale(i15(j,ii)) = m5scale
            pscale(i15(j,ii)) = p5scale
         end do
         do j = 1, np11(ii)
            dscale(ip11(j,ii)) = d1scale
            uscale(ip11(j,ii)) = u1scale
         end do
         do j = 1, np12(ii)
            dscale(ip12(j,ii)) = d2scale
            uscale(ip12(j,ii)) = u2scale
         end do
         do j = 1, np13(ii)
            dscale(ip13(j,ii)) = d3scale
            uscale(ip13(j,ii)) = u3scale
         end do
         do j = 1, np14(ii)
            dscale(ip14(j,ii)) = d4scale
            uscale(ip14(j,ii)) = u4scale
         end do
         do l3 = l1, npole3b
            k = pnum(l3)
            kk = ipole(k)
            do jcell=1,ncell
            xr = x(kk) - x(ii)
            yr = y(kk) - y(ii)
            zr = z(kk) - z(ii)
            call imager (xr,yr,zr,jcell)
            r2 = xr*xr + yr*yr + zr*zr
            if (.not.(use_polymer .and. r2.le.polycut2)) then
              mscale(kk)=1.0d0
              pscale(kk)=1.0d0
              dscale(kk)=1.0d0
              uscale(kk)=1.0d0
            end if
            if (r2 .le. off2) then
               r = sqrt(r2)
               ck = rpole(1,k)
               dk(1) = rpole(2,k)
               dk(2) = rpole(3,k)
               dk(3) = rpole(4,k)
               qk(1) = rpole(5,k)
               qk(2) = rpole(6,k)
               qk(3) = rpole(7,k)
               qk(4) = rpole(8,k)
               qk(5) = rpole(9,k)
               qk(6) = rpole(10,k)
               qk(7) = rpole(11,k)
               qk(8) = rpole(12,k)
               qk(9) = rpole(13,k)
c
c     calculate the real space error function terms
c
               ralpha = aewald * r
               bn(0) = erfc(ralpha) / r
               alsq2 = 2.0d0 * aewald**2
               alsq2n = 0.0d0
               if (aewald .gt. 0.0d0)  alsq2n = 1.0d0 / (sqrtpi*aewald)
               exp2a = exp(-ralpha**2)
               do j = 1, 5
                  bfac = dble(2*j-1)
                  alsq2n = alsq2 * alsq2n
                  bn(j) = (bfac*bn(j-1)+alsq2n*exp2a) / r2
               end do
c
c     apply Thole polarization damping to scale factors
c
               rr1 = 1.0d0 / r
               rr3 = rr1 / r2
               rr5 = 3.0d0 * rr3 / r2
               rr7 = 5.0d0 * rr5 / r2
               rr9 = 7.0d0 * rr7 / r2
               rr11 = 9.0d0 * rr9 / r2
               scale3 = 1.0d0
               scale5 = 1.0d0
               scale7 = 1.0d0
               do j = 1, 3
                  ddsc3(j) = 0.0d0
                  ddsc5(j) = 0.0d0
                  ddsc7(j) = 0.0d0
               end do
               damp = pdi * pdamp(k)
               if (damp .ne. 0.0d0) then
                  pgamma = min(pti,thole(k))
                  damp = -pgamma * (r/damp)**3
                  if (damp .gt. -50.0d0) then
                     expdamp = exp(damp)
                     scale3 = 1.0d0 - expdamp
                     scale5 = 1.0d0 - (1.0d0-damp)*expdamp
                     scale7 = 1.0d0 - (1.0d0-damp+0.6d0*damp**2)
     &                                       *expdamp
                     temp3 = -3.0d0 * damp * expdamp / r2
                     temp5 = -damp
                     temp7 = -0.2d0 - 0.6d0*damp
                     ddsc3(1) = temp3 * xr
                     ddsc3(2) = temp3 * yr
                     ddsc3(3) = temp3 * zr
                     ddsc5(1) = temp5 * ddsc3(1)
                     ddsc5(2) = temp5 * ddsc3(2)
                     ddsc5(3) = temp5 * ddsc3(3)
                     ddsc7(1) = temp7 * ddsc5(1)
                     ddsc7(2) = temp7 * ddsc5(2)
                     ddsc7(3) = temp7 * ddsc5(3)
                  end if
               end if
               dsc3 = 1.0d0 - scale3*dscale(kk)
               dsc5 = 1.0d0 - scale5*dscale(kk)
               dsc7 = 1.0d0 - scale7*dscale(kk)
               psc3 = 1.0d0 - scale3*pscale(kk)
               psc5 = 1.0d0 - scale5*pscale(kk)
               psc7 = 1.0d0 - scale7*pscale(kk)
               usc3 = 1.0d0 - scale3*uscale(kk)
               usc5 = 1.0d0 - scale5*uscale(kk)
c
c     construct necessary auxiliary vectors
c
               dixdk(1) = di(2)*dk(3) - di(3)*dk(2)
               dixdk(2) = di(3)*dk(1) - di(1)*dk(3)
               dixdk(3) = di(1)*dk(2) - di(2)*dk(1)
               dixuk(1) = di(2)*uind(3,k) - di(3)*uind(2,k)
               dixuk(2) = di(3)*uind(1,k) - di(1)*uind(3,k)
               dixuk(3) = di(1)*uind(2,k) - di(2)*uind(1,k)
               dkxui(1) = dk(2)*uind(3,i) - dk(3)*uind(2,i)
               dkxui(2) = dk(3)*uind(1,i) - dk(1)*uind(3,i)
               dkxui(3) = dk(1)*uind(2,i) - dk(2)*uind(1,i)
               dixukp(1) = di(2)*uinp(3,k) - di(3)*uinp(2,k)
               dixukp(2) = di(3)*uinp(1,k) - di(1)*uinp(3,k)
               dixukp(3) = di(1)*uinp(2,k) - di(2)*uinp(1,k)
               dkxuip(1) = dk(2)*uinp(3,i) - dk(3)*uinp(2,i)
               dkxuip(2) = dk(3)*uinp(1,i) - dk(1)*uinp(3,i)
               dkxuip(3) = dk(1)*uinp(2,i) - dk(2)*uinp(1,i)
               dixr(1) = di(2)*zr - di(3)*yr
               dixr(2) = di(3)*xr - di(1)*zr
               dixr(3) = di(1)*yr - di(2)*xr
               dkxr(1) = dk(2)*zr - dk(3)*yr
               dkxr(2) = dk(3)*xr - dk(1)*zr
               dkxr(3) = dk(1)*yr - dk(2)*xr
               qir(1) = qi(1)*xr + qi(4)*yr + qi(7)*zr
               qir(2) = qi(2)*xr + qi(5)*yr + qi(8)*zr
               qir(3) = qi(3)*xr + qi(6)*yr + qi(9)*zr
               qkr(1) = qk(1)*xr + qk(4)*yr + qk(7)*zr
               qkr(2) = qk(2)*xr + qk(5)*yr + qk(8)*zr
               qkr(3) = qk(3)*xr + qk(6)*yr + qk(9)*zr
               qiqkr(1) = qi(1)*qkr(1) + qi(4)*qkr(2) + qi(7)*qkr(3)
               qiqkr(2) = qi(2)*qkr(1) + qi(5)*qkr(2) + qi(8)*qkr(3)
               qiqkr(3) = qi(3)*qkr(1) + qi(6)*qkr(2) + qi(9)*qkr(3)
               qkqir(1) = qk(1)*qir(1) + qk(4)*qir(2) + qk(7)*qir(3)
               qkqir(2) = qk(2)*qir(1) + qk(5)*qir(2) + qk(8)*qir(3)
               qkqir(3) = qk(3)*qir(1) + qk(6)*qir(2) + qk(9)*qir(3)
               qixqk(1) = qi(2)*qk(3) + qi(5)*qk(6) + qi(8)*qk(9)
     &                       - qi(3)*qk(2) - qi(6)*qk(5) - qi(9)*qk(8)
               qixqk(2) = qi(3)*qk(1) + qi(6)*qk(4) + qi(9)*qk(7)
     &                       - qi(1)*qk(3) - qi(4)*qk(6) - qi(7)*qk(9)
               qixqk(3) = qi(1)*qk(2) + qi(4)*qk(5) + qi(7)*qk(8)
     &                       - qi(2)*qk(1) - qi(5)*qk(4) - qi(8)*qk(7)
               rxqir(1) = yr*qir(3) - zr*qir(2)
               rxqir(2) = zr*qir(1) - xr*qir(3)
               rxqir(3) = xr*qir(2) - yr*qir(1)
               rxqkr(1) = yr*qkr(3) - zr*qkr(2)
               rxqkr(2) = zr*qkr(1) - xr*qkr(3)
               rxqkr(3) = xr*qkr(2) - yr*qkr(1)
               rxqikr(1) = yr*qiqkr(3) - zr*qiqkr(2)
               rxqikr(2) = zr*qiqkr(1) - xr*qiqkr(3)
               rxqikr(3) = xr*qiqkr(2) - yr*qiqkr(1)
               rxqkir(1) = yr*qkqir(3) - zr*qkqir(2)
               rxqkir(2) = zr*qkqir(1) - xr*qkqir(3)
               rxqkir(3) = xr*qkqir(2) - yr*qkqir(1)
               qkrxqir(1) = qkr(2)*qir(3) - qkr(3)*qir(2)
               qkrxqir(2) = qkr(3)*qir(1) - qkr(1)*qir(3)
               qkrxqir(3) = qkr(1)*qir(2) - qkr(2)*qir(1)
               qidk(1) = qi(1)*dk(1) + qi(4)*dk(2) + qi(7)*dk(3)
               qidk(2) = qi(2)*dk(1) + qi(5)*dk(2) + qi(8)*dk(3)
               qidk(3) = qi(3)*dk(1) + qi(6)*dk(2) + qi(9)*dk(3)
               qkdi(1) = qk(1)*di(1) + qk(4)*di(2) + qk(7)*di(3)
               qkdi(2) = qk(2)*di(1) + qk(5)*di(2) + qk(8)*di(3)
               qkdi(3) = qk(3)*di(1) + qk(6)*di(2) + qk(9)*di(3)
               qiuk(1) = qi(1)*uind(1,k) + qi(4)*uind(2,k)
     &                      + qi(7)*uind(3,k)
               qiuk(2) = qi(2)*uind(1,k) + qi(5)*uind(2,k)
     &                      + qi(8)*uind(3,k)
               qiuk(3) = qi(3)*uind(1,k) + qi(6)*uind(2,k)
     &                      + qi(9)*uind(3,k)
               qkui(1) = qk(1)*uind(1,i) + qk(4)*uind(2,i)
     &                      + qk(7)*uind(3,i)
               qkui(2) = qk(2)*uind(1,i) + qk(5)*uind(2,i)
     &                      + qk(8)*uind(3,i)
               qkui(3) = qk(3)*uind(1,i) + qk(6)*uind(2,i)
     &                      + qk(9)*uind(3,i)
               qiukp(1) = qi(1)*uinp(1,k) + qi(4)*uinp(2,k)
     &                       + qi(7)*uinp(3,k)
               qiukp(2) = qi(2)*uinp(1,k) + qi(5)*uinp(2,k)
     &                       + qi(8)*uinp(3,k)
               qiukp(3) = qi(3)*uinp(1,k) + qi(6)*uinp(2,k)
     &                       + qi(9)*uinp(3,k)
               qkuip(1) = qk(1)*uinp(1,i) + qk(4)*uinp(2,i)
     &                       + qk(7)*uinp(3,i)
               qkuip(2) = qk(2)*uinp(1,i) + qk(5)*uinp(2,i)
     &                       + qk(8)*uinp(3,i)
               qkuip(3) = qk(3)*uinp(1,i) + qk(6)*uinp(2,i)
     &                       + qk(9)*uinp(3,i)
               dixqkr(1) = di(2)*qkr(3) - di(3)*qkr(2)
               dixqkr(2) = di(3)*qkr(1) - di(1)*qkr(3)
               dixqkr(3) = di(1)*qkr(2) - di(2)*qkr(1)
               dkxqir(1) = dk(2)*qir(3) - dk(3)*qir(2)
               dkxqir(2) = dk(3)*qir(1) - dk(1)*qir(3)
               dkxqir(3) = dk(1)*qir(2) - dk(2)*qir(1)
               uixqkr(1) = uind(2,i)*qkr(3) - uind(3,i)*qkr(2)
               uixqkr(2) = uind(3,i)*qkr(1) - uind(1,i)*qkr(3)
               uixqkr(3) = uind(1,i)*qkr(2) - uind(2,i)*qkr(1)
               ukxqir(1) = uind(2,k)*qir(3) - uind(3,k)*qir(2)
               ukxqir(2) = uind(3,k)*qir(1) - uind(1,k)*qir(3)
               ukxqir(3) = uind(1,k)*qir(2) - uind(2,k)*qir(1)
               uixqkrp(1) = uinp(2,i)*qkr(3) - uinp(3,i)*qkr(2)
               uixqkrp(2) = uinp(3,i)*qkr(1) - uinp(1,i)*qkr(3)
               uixqkrp(3) = uinp(1,i)*qkr(2) - uinp(2,i)*qkr(1)
               ukxqirp(1) = uinp(2,k)*qir(3) - uinp(3,k)*qir(2)
               ukxqirp(2) = uinp(3,k)*qir(1) - uinp(1,k)*qir(3)
               ukxqirp(3) = uinp(1,k)*qir(2) - uinp(2,k)*qir(1)
               rxqidk(1) = yr*qidk(3) - zr*qidk(2)
               rxqidk(2) = zr*qidk(1) - xr*qidk(3)
               rxqidk(3) = xr*qidk(2) - yr*qidk(1)
               rxqkdi(1) = yr*qkdi(3) - zr*qkdi(2)
               rxqkdi(2) = zr*qkdi(1) - xr*qkdi(3)
               rxqkdi(3) = xr*qkdi(2) - yr*qkdi(1)
               rxqiuk(1) = yr*qiuk(3) - zr*qiuk(2)
               rxqiuk(2) = zr*qiuk(1) - xr*qiuk(3)
               rxqiuk(3) = xr*qiuk(2) - yr*qiuk(1)
               rxqkui(1) = yr*qkui(3) - zr*qkui(2)
               rxqkui(2) = zr*qkui(1) - xr*qkui(3)
               rxqkui(3) = xr*qkui(2) - yr*qkui(1)
               rxqiukp(1) = yr*qiukp(3) - zr*qiukp(2)
               rxqiukp(2) = zr*qiukp(1) - xr*qiukp(3)
               rxqiukp(3) = xr*qiukp(2) - yr*qiukp(1)
               rxqkuip(1) = yr*qkuip(3) - zr*qkuip(2)
               rxqkuip(2) = zr*qkuip(1) - xr*qkuip(3)
               rxqkuip(3) = xr*qkuip(2) - yr*qkuip(1)
c
c     calculate the scalar products for permanent components
c
               sc(2) = di(1)*dk(1) + di(2)*dk(2) + di(3)*dk(3)
               sc(3) = di(1)*xr + di(2)*yr + di(3)*zr
               sc(4) = dk(1)*xr + dk(2)*yr + dk(3)*zr
               sc(5) = qir(1)*xr + qir(2)*yr + qir(3)*zr
               sc(6) = qkr(1)*xr + qkr(2)*yr + qkr(3)*zr
               sc(7) = qir(1)*dk(1) + qir(2)*dk(2) + qir(3)*dk(3)
               sc(8) = qkr(1)*di(1) + qkr(2)*di(2) + qkr(3)*di(3)
               sc(9) = qir(1)*qkr(1) + qir(2)*qkr(2) + qir(3)*qkr(3)
               sc(10) = qi(1)*qk(1) + qi(2)*qk(2) + qi(3)*qk(3)
     &                     + qi(4)*qk(4) + qi(5)*qk(5) + qi(6)*qk(6)
     &                     + qi(7)*qk(7) + qi(8)*qk(8) + qi(9)*qk(9)
c
c     calculate the scalar products for induced components
c
               sci(1) = uind(1,i)*dk(1) + uind(2,i)*dk(2)
     &                     + uind(3,i)*dk(3) + di(1)*uind(1,k)
     &                     + di(2)*uind(2,k) + di(3)*uind(3,k)
               sci(2) = uind(1,i)*uind(1,k) + uind(2,i)*uind(2,k)
     &                     + uind(3,i)*uind(3,k)
               sci(3) = uind(1,i)*xr + uind(2,i)*yr + uind(3,i)*zr
               sci(4) = uind(1,k)*xr + uind(2,k)*yr + uind(3,k)*zr
               sci(7) = qir(1)*uind(1,k) + qir(2)*uind(2,k)
     &                     + qir(3)*uind(3,k)
               sci(8) = qkr(1)*uind(1,i) + qkr(2)*uind(2,i)
     &                     + qkr(3)*uind(3,i)
               scip(1) = uinp(1,i)*dk(1) + uinp(2,i)*dk(2)
     &                      + uinp(3,i)*dk(3) + di(1)*uinp(1,k)
     &                      + di(2)*uinp(2,k) + di(3)*uinp(3,k)
               scip(2) = uind(1,i)*uinp(1,k)+uind(2,i)*uinp(2,k)
     &                      + uind(3,i)*uinp(3,k)+uinp(1,i)*uind(1,k)
     &                      + uinp(2,i)*uind(2,k)+uinp(3,i)*uind(3,k)
               scip(3) = uinp(1,i)*xr + uinp(2,i)*yr + uinp(3,i)*zr
               scip(4) = uinp(1,k)*xr + uinp(2,k)*yr + uinp(3,k)*zr
               scip(7) = qir(1)*uinp(1,k) + qir(2)*uinp(2,k)
     &                      + qir(3)*uinp(3,k)
               scip(8) = qkr(1)*uinp(1,i) + qkr(2)*uinp(2,i)
     &                      + qkr(3)*uinp(3,i)
c
c     calculate the gl functions for permanent components
c
               gl(0) = ci*ck
               gl(1) = ck*sc(3) - ci*sc(4)
               gl(2) = ci*sc(6) + ck*sc(5) - sc(3)*sc(4)
               gl(3) = sc(3)*sc(6) - sc(4)*sc(5)
               gl(4) = sc(5)*sc(6)
               gl(5) = -4.0d0 * sc(9)
               gl(6) = sc(2)
               gl(7) = 2.0d0 * (sc(7)-sc(8))
               gl(8) = 2.0d0 * sc(10)
c
c     calculate the gl functions for induced components
c
               gli(1) = ck*sci(3) - ci*sci(4)
               gli(2) = -sc(3)*sci(4) - sci(3)*sc(4)
               gli(3) = sci(3)*sc(6) - sci(4)*sc(5)
               gli(6) = sci(1)
               gli(7) = 2.0d0 * (sci(7)-sci(8))
               glip(1) = ck*scip(3) - ci*scip(4)
               glip(2) = -sc(3)*scip(4) - scip(3)*sc(4)
               glip(3) = scip(3)*sc(6) - scip(4)*sc(5)
               glip(6) = scip(1)
               glip(7) = 2.0d0 * (scip(7)-scip(8))
c
c     compute the energy contributions for this interaction
c
               ei = 0.5d0 * (bn(1)*(gli(1)+gli(6))
     &                      + bn(2)*(gli(2)+gli(7)) + bn(3)*gli(3))
c
c     get the real energy without any screening function
c
               erli = 0.5d0*(rr3*(gli(1)+gli(6))*psc3
     &                   + rr5*(gli(2)+gli(7))*psc5
     &                   + rr7*gli(3)*psc7)
               ei = ei - erli
               ei = f * ei
               if (ii .eq. kk) then
                  ei = 0.5d0 * ei
               end if
               ep = ep + ei
c
c     increment the total intramolecular energy; assumes
c     intramolecular distances are less than half of cell
c     length and less than the ewald cutoff
c
               if (molcule(ii) .eq. molcule(kk)) then
                  eintra = eintra + 0.5d0*pscale(kk)
     &                        * (rr3*(gli(1)+gli(6))*scale3
     &                              + rr5*(gli(2)+gli(7))*scale5
     &                              + rr7*gli(3)*scale7)
               end if
c
c     intermediate variables for induced force terms
c
               gfi(1) = 0.5d0*bn(2)*(gli(1)+glip(1)+gli(6)+glip(6))
     &                     + 0.5d0*bn(2)*scip(2)
     &                     + 0.5d0*bn(3)*(gli(2)+glip(2)+gli(7)+glip(7))
     &                     - 0.5d0*bn(3)*(sci(3)*scip(4)+scip(3)*sci(4))
     &                     + 0.5d0*bn(4)*(gli(3)+glip(3))
               gfi(2) = -ck*bn(1) + sc(4)*bn(2) - sc(6)*bn(3)
               gfi(3) = ci*bn(1) + sc(3)*bn(2) + sc(5)*bn(3)
               gfi(4) = 2.0d0 * bn(2)
               gfi(5) = bn(3) * (sci(4)+scip(4))
               gfi(6) = -bn(3) * (sci(3)+scip(3))
               gfri(1) = 0.5d0*rr5*((gli(1)+gli(6))*psc3
     &                            + (glip(1)+glip(6))*dsc3
     &                            + scip(2)*usc3)
     &                 + 0.5d0*rr7*((gli(7)+gli(2))*psc5
     &                            + (glip(7)+glip(2))*dsc5
     &                     - (sci(3)*scip(4)+scip(3)*sci(4))*usc5)
     &                 + 0.5d0*rr9*(gli(3)*psc7+glip(3)*dsc7)
               gfri(2) = -rr3*ck + rr5*sc(4) - rr7*sc(6)
               gfri(3) = rr3*ci + rr5*sc(3) + rr7*sc(5)
               gfri(4) = 2.0d0 * rr5
               gfri(5) = rr7 * (sci(4)*psc7+scip(4)*dsc7)
               gfri(6) = -rr7 * (sci(3)*psc7+scip(3)*dsc7)
c
c     get the induced force with screening
c
               ftm2i(1) = gfi(1)*xr + 0.5d0*
     &             (gfi(2)*(uind(1,i)+uinp(1,i))
     &            + bn(2)*(sci(4)*uinp(1,i)+scip(4)*uind(1,i))
     &            + gfi(3)*(uind(1,k)+uinp(1,k))
     &            + bn(2)*(sci(3)*uinp(1,k)+scip(3)*uind(1,k))
     &            + (sci(4)+scip(4))*bn(2)*di(1)
     &            + (sci(3)+scip(3))*bn(2)*dk(1)
     &            + gfi(4)*(qkui(1)+qkuip(1)-qiuk(1)-qiukp(1)))
     &            + gfi(5)*qir(1) + gfi(6)*qkr(1)
               ftm2i(2) = gfi(1)*yr + 0.5d0*
     &             (gfi(2)*(uind(2,i)+uinp(2,i))
     &            + bn(2)*(sci(4)*uinp(2,i)+scip(4)*uind(2,i))
     &            + gfi(3)*(uind(2,k)+uinp(2,k))
     &            + bn(2)*(sci(3)*uinp(2,k)+scip(3)*uind(2,k))
     &            + (sci(4)+scip(4))*bn(2)*di(2)
     &            + (sci(3)+scip(3))*bn(2)*dk(2)
     &            + gfi(4)*(qkui(2)+qkuip(2)-qiuk(2)-qiukp(2)))
     &            + gfi(5)*qir(2) + gfi(6)*qkr(2)
               ftm2i(3) = gfi(1)*zr + 0.5d0*
     &             (gfi(2)*(uind(3,i)+uinp(3,i))
     &            + bn(2)*(sci(4)*uinp(3,i)+scip(4)*uind(3,i))
     &            + gfi(3)*(uind(3,k)+uinp(3,k))
     &            + bn(2)*(sci(3)*uinp(3,k)+scip(3)*uind(3,k))
     &            + (sci(4)+scip(4))*bn(2)*di(3)
     &            + (sci(3)+scip(3))*bn(2)*dk(3)
     &            + gfi(4)*(qkui(3)+qkuip(3)-qiuk(3)-qiukp(3)))
     &            + gfi(5)*qir(3) + gfi(6)*qkr(3)
c
c     get the induced force without screening
c
               ftm2ri(1) = gfri(1)*xr + 0.5d0*
     &           (- rr3*ck*(uind(1,i)*psc3+uinp(1,i)*dsc3)
     &            + rr5*sc(4)*(uind(1,i)*psc5+uinp(1,i)*dsc5)
     &            - rr7*sc(6)*(uind(1,i)*psc7+uinp(1,i)*dsc7))
     &            + (rr3*ci*(uind(1,k)*psc3+uinp(1,k)*dsc3)
     &            + rr5*sc(3)*(uind(1,k)*psc5+uinp(1,k)*dsc5)
     &            + rr7*sc(5)*(uind(1,k)*psc7+uinp(1,k)*dsc7))*0.5d0
     &            + rr5*usc5*(sci(4)*uinp(1,i)+scip(4)*uind(1,i)
     &            + sci(3)*uinp(1,k)+scip(3)*uind(1,k))*0.5d0
     &            + 0.5d0*(sci(4)*psc5+scip(4)*dsc5)*rr5*di(1)
     &            + 0.5d0*(sci(3)*psc5+scip(3)*dsc5)*rr5*dk(1)
     &            + 0.5d0*gfri(4)*((qkui(1)-qiuk(1))*psc5
     &            + (qkuip(1)-qiukp(1))*dsc5)
     &            + gfri(5)*qir(1) + gfri(6)*qkr(1)
               ftm2ri(2) = gfri(1)*yr + 0.5d0*
     &           (- rr3*ck*(uind(2,i)*psc3+uinp(2,i)*dsc3)
     &            + rr5*sc(4)*(uind(2,i)*psc5+uinp(2,i)*dsc5)
     &            - rr7*sc(6)*(uind(2,i)*psc7+uinp(2,i)*dsc7))
     &            + (rr3*ci*(uind(2,k)*psc3+uinp(2,k)*dsc3)
     &            + rr5*sc(3)*(uind(2,k)*psc5+uinp(2,k)*dsc5)
     &            + rr7*sc(5)*(uind(2,k)*psc7+uinp(2,k)*dsc7))*0.5d0
     &            + rr5*usc5*(sci(4)*uinp(2,i)+scip(4)*uind(2,i)
     &            + sci(3)*uinp(2,k)+scip(3)*uind(2,k))*0.5d0
     &            + 0.5d0*(sci(4)*psc5+scip(4)*dsc5)*rr5*di(2)
     &            + 0.5d0*(sci(3)*psc5+scip(3)*dsc5)*rr5*dk(2)
     &            + 0.5d0*gfri(4)*((qkui(2)-qiuk(2))*psc5
     &            + (qkuip(2)-qiukp(2))*dsc5)
     &            + gfri(5)*qir(2) + gfri(6)*qkr(2)
               ftm2ri(3) = gfri(1)*zr + 0.5d0*
     &           (- rr3*ck*(uind(3,i)*psc3+uinp(3,i)*dsc3)
     &            + rr5*sc(4)*(uind(3,i)*psc5+uinp(3,i)*dsc5)
     &            - rr7*sc(6)*(uind(3,i)*psc7+uinp(3,i)*dsc7))
     &            + (rr3*ci*(uind(3,k)*psc3+uinp(3,k)*dsc3)
     &            + rr5*sc(3)*(uind(3,k)*psc5+uinp(3,k)*dsc5)
     &            + rr7*sc(5)*(uind(3,k)*psc7+uinp(3,k)*dsc7))*0.5d0
     &            + rr5*usc5*(sci(4)*uinp(3,i)+scip(4)*uind(3,i)
     &            + sci(3)*uinp(3,k)+scip(3)*uind(3,k))*0.5d0
     &            + 0.5d0*(sci(4)*psc5+scip(4)*dsc5)*rr5*di(3)
     &            + 0.5d0*(sci(3)*psc5+scip(3)*dsc5)*rr5*dk(3)
     &            + 0.5d0*gfri(4)*((qkui(3)-qiuk(3))*psc5
     &            + (qkuip(3)-qiukp(3))*dsc5)
     &            + gfri(5)*qir(3) + gfri(6)*qkr(3)
c
c     account for partially excluded induced interactions
c
               temp3 = 0.5d0 * rr3 * ((gli(1)+gli(6))*pscale(kk)
     &                                  +(glip(1)+glip(6))*dscale(kk))
               temp5 = 0.5d0 * rr5 * ((gli(2)+gli(7))*pscale(kk)
     &                                  +(glip(2)+glip(7))*dscale(kk))
               temp7 = 0.5d0 * rr7 * (gli(3)*pscale(kk)
     &                                  +glip(3)*dscale(kk))
               fridmp(1) = temp3*ddsc3(1) + temp5*ddsc5(1)
     &                        + temp7*ddsc7(1)
               fridmp(2) = temp3*ddsc3(2) + temp5*ddsc5(2)
     &                        + temp7*ddsc7(2)
               fridmp(3) = temp3*ddsc3(3) + temp5*ddsc5(3)
     &                        + temp7*ddsc7(3)
c
c     find some scaling terms for induced-induced force
c
               temp3 = 0.5d0 * rr3 * uscale(kk) * scip(2)
               temp5 = -0.5d0 * rr5 * uscale(kk)
     &                    * (sci(3)*scip(4)+scip(3)*sci(4))
               findmp(1) = temp3*ddsc3(1) + temp5*ddsc5(1)
               findmp(2) = temp3*ddsc3(2) + temp5*ddsc5(2)
               findmp(3) = temp3*ddsc3(3) + temp5*ddsc5(3)
c
c     modify the forces for partially excluded interactions
c
               ftm2i(1) = ftm2i(1) - fridmp(1) - findmp(1)
               ftm2i(2) = ftm2i(2) - fridmp(2) - findmp(2)
               ftm2i(3) = ftm2i(3) - fridmp(3) - findmp(3)
c
c     correction to convert mutual to direct polarization force
c
               if (poltyp .eq. 'DIRECT') then
                  gfd = 0.5d0 * (bn(2)*scip(2)
     &                     - bn(3)*(scip(3)*sci(4)+sci(3)*scip(4)))
                  gfdr = 0.5d0 * (rr5*scip(2)*usc3
     &                     - rr7*(scip(3)*sci(4)
     &                           +sci(3)*scip(4))*usc5)

                  ftm2i(1) = ftm2i(1) - gfd*xr - 0.5d0*bn(2)*
     &                          (sci(4)*uinp(1,i)+scip(4)*uind(1,i)
     &                          +sci(3)*uinp(1,k)+scip(3)*uind(1,k))
                  ftm2i(2) = ftm2i(2) - gfd*yr - 0.5d0*bn(2)*
     &                          (sci(4)*uinp(2,i)+scip(4)*uind(2,i)
     &                          +sci(3)*uinp(2,k)+scip(3)*uind(2,k))
                  ftm2i(3) = ftm2i(3) - gfd*zr - 0.5d0*bn(2)*
     &                          (sci(4)*uinp(3,i)+scip(4)*uind(3,i)
     &                          +sci(3)*uinp(3,k)+scip(3)*uind(3,k))
                  fdir(1) = gfdr*xr + 0.5d0*usc5*rr5*
     &                         (sci(4)*uinp(1,i)+scip(4)*uind(1,i)
     &                        + sci(3)*uinp(1,k)+scip(3)*uind(1,k))
                  fdir(2) = gfdr*yr + 0.5d0*usc5*rr5*
     &                         (sci(4)*uinp(2,i)+scip(4)*uind(2,i)
     &                        + sci(3)*uinp(2,k)+scip(3)*uind(2,k))
                  fdir(3) = gfdr*zr + 0.5d0*usc5*rr5*
     &                         (sci(4)*uinp(3,i)+scip(4)*uind(3,i)
     &                        + sci(3)*uinp(3,k)+scip(3)*uind(3,k))
                  ftm2i(1) = ftm2i(1) + fdir(1) + findmp(1)
                  ftm2i(2) = ftm2i(2) + fdir(2) + findmp(2)
                  ftm2i(3) = ftm2i(3) + fdir(3) + findmp(3)
               end if
c
c     intermediate variables for induced torque terms
c
               gti(2) = 0.5d0 * bn(2) * (sci(4)+scip(4))
               gti(3) = 0.5d0 * bn(2) * (sci(3)+scip(3))
               gti(4) = gfi(4)
               gti(5) = gfi(5)
               gti(6) = gfi(6)
               gtri(2) = 0.5d0 * rr5 * (sci(4)*psc5+scip(4)*dsc5)
               gtri(3) = 0.5d0 * rr5 * (sci(3)*psc5+scip(3)*dsc5)
               gtri(4) = gfri(4)
               gtri(5) = gfri(5)
               gtri(6) = gfri(6)
c
c     get the induced torque with screening
c
               ttm2i(1) = -bn(1)*(dixuk(1)+dixukp(1))*0.5d0
     &           + gti(2)*dixr(1) + gti(4)*(ukxqir(1)+rxqiuk(1)
     &           + ukxqirp(1)+rxqiukp(1))*0.5d0 - gti(5)*rxqir(1)
               ttm2i(2) = -bn(1)*(dixuk(2)+dixukp(2))*0.5d0
     &           + gti(2)*dixr(2) + gti(4)*(ukxqir(2)+rxqiuk(2)
     &           + ukxqirp(2)+rxqiukp(2))*0.5d0 - gti(5)*rxqir(2)
               ttm2i(3) = -bn(1)*(dixuk(3)+dixukp(3))*0.5d0
     &           + gti(2)*dixr(3) + gti(4)*(ukxqir(3)+rxqiuk(3)
     &           + ukxqirp(3)+rxqiukp(3))*0.5d0 - gti(5)*rxqir(3)
               ttm3i(1) = -bn(1)*(dkxui(1)+dkxuip(1))*0.5d0
     &           + gti(3)*dkxr(1) - gti(4)*(uixqkr(1)+rxqkui(1)
     &           + uixqkrp(1)+rxqkuip(1))*0.5d0 - gti(6)*rxqkr(1)
               ttm3i(2) = -bn(1)*(dkxui(2)+dkxuip(2))*0.5d0
     &           + gti(3)*dkxr(2) - gti(4)*(uixqkr(2)+rxqkui(2)
     &           + uixqkrp(2)+rxqkuip(2))*0.5d0 - gti(6)*rxqkr(2)
               ttm3i(3) = -bn(1)*(dkxui(3)+dkxuip(3))*0.5d0
     &           + gti(3)*dkxr(3) - gti(4)*(uixqkr(3)+rxqkui(3)
     &           + uixqkrp(3)+rxqkuip(3))*0.5d0 - gti(6)*rxqkr(3)
c
c     get the induced torque without screening
c
               ttm2ri(1) = -rr3*(dixuk(1)*psc3+dixukp(1)*dsc3)*0.5d0
     &           + gtri(2)*dixr(1) + gtri(4)*((ukxqir(1)+rxqiuk(1))*psc5
     &           +(ukxqirp(1)+rxqiukp(1))*dsc5)*0.5d0 - gtri(5)*rxqir(1)
               ttm2ri(2) = -rr3*(dixuk(2)*psc3+dixukp(2)*dsc3)*0.5d0
     &           + gtri(2)*dixr(2) + gtri(4)*((ukxqir(2)+rxqiuk(2))*psc5
     &           +(ukxqirp(2)+rxqiukp(2))*dsc5)*0.5d0 - gtri(5)*rxqir(2)
               ttm2ri(3) = -rr3*(dixuk(3)*psc3+dixukp(3)*dsc3)*0.5d0
     &           + gtri(2)*dixr(3) + gtri(4)*((ukxqir(3)+rxqiuk(3))*psc5
     &           +(ukxqirp(3)+rxqiukp(3))*dsc5)*0.5d0 - gtri(5)*rxqir(3)
               ttm3ri(1) = -rr3*(dkxui(1)*psc3+dkxuip(1)*dsc3)*0.5d0
     &           + gtri(3)*dkxr(1) - gtri(4)*((uixqkr(1)+rxqkui(1))*psc5
     &           +(uixqkrp(1)+rxqkuip(1))*dsc5)*0.5d0 - gtri(6)*rxqkr(1)
               ttm3ri(2) = -rr3*(dkxui(2)*psc3+dkxuip(2)*dsc3)*0.5d0
     &           + gtri(3)*dkxr(2) - gtri(4)*((uixqkr(2)+rxqkui(2))*psc5
     &           +(uixqkrp(2)+rxqkuip(2))*dsc5)*0.5d0 - gtri(6)*rxqkr(2)
               ttm3ri(3) = -rr3*(dkxui(3)*psc3+dkxuip(3)*dsc3)*0.5d0
     &           + gtri(3)*dkxr(3) - gtri(4)*((uixqkr(3)+rxqkui(3))*psc5
     &           +(uixqkrp(3)+rxqkuip(3))*dsc5)*0.5d0 - gtri(6)*rxqkr(3)
c
c     handle the case where scaling is used
c
              if (use_polymer .and. r2.le.polycut2) then
               do j = 1, 3
                  ftm2i(j) = f * (ftm2i(j)-ftm2ri(j))
                  ttm2i(j) = f * (ttm2i(j)-ttm2ri(j))
                  ttm3i(j) = f * (ttm3i(j)-ttm3ri(j))
               end do
              else
                  do j = 1, 3
                     ftm2i(j) = f * (ftm2i(j)-ftm2ri(j))
                     ttm2i(j) = f * (ttm2i(j)-ttm2ri(j))
                     ttm3i(j) = f * (ttm3i(j)-ttm3ri(j))
                  end do
               end if
               if (ii .eq. kk) then
                  do j = 1, 3
                     ftm2i(j) = 0.5d0 * ftm2i(j)
                     ttm2i(j) = 0.5d0 * ttm2i(j)
                     ttm3i(j) = 0.5d0 * ttm3i(j)
                  end do
               end if
c
c     increment gradient due to force and torque on first site
c
               dep(1,ii) = dep(1,ii) + ftm2i(1)
               dep(2,ii) = dep(2,ii) + ftm2i(2)
               dep(3,ii) = dep(3,ii) + ftm2i(3)
               call torque_3b (i,ttm2,ttm2i,frcxi,frcyi,frczi)

c
c     increment gradient due to force and torque on second site
c
               dep(1,kk) = dep(1,kk) - ftm2i(1)
               dep(2,kk) = dep(2,kk) - ftm2i(2)
               dep(3,kk) = dep(3,kk) - ftm2i(3)
               call torque_3b (k,ttm3,ttm3i,frcxk,frcyk,frczk)
c
c     increment the internal virial tensor components
c
               iaz = zaxis(i)
               iax = xaxis(i)
               iay = yaxis(i)
               kaz = zaxis(k)
               kax = xaxis(k)
               kay = yaxis(k)
               if (iaz .eq. 0)  iaz = ii
               if (iax .eq. 0)  iax = ii
               if (iay .eq. 0)  iay = ii
               if (kaz .eq. 0)  kaz = kk
               if (kax .eq. 0)  kax = kk
               if (kay .eq. 0)  kay = kk
               xiz = x(iaz) - x(ii)
               yiz = y(iaz) - y(ii)
               ziz = z(iaz) - z(ii)
               xix = x(iax) - x(ii)
               yix = y(iax) - y(ii)
               zix = z(iax) - z(ii)
               xiy = x(iay) - x(ii)
               yiy = y(iay) - y(ii)
               ziy = z(iay) - z(ii)
               xkz = x(kaz) - x(kk)
               ykz = y(kaz) - y(kk)
               zkz = z(kaz) - z(kk)
               xkx = x(kax) - x(kk)
               ykx = y(kax) - y(kk)
               zkx = z(kax) - z(kk)
               xky = x(kay) - x(kk)
               yky = y(kay) - y(kk)
               zky = z(kay) - z(kk)
               vxx = -xr*ftm2i(1)+xix*frcxi(1)+xiy*frcyi(1)+xiz*frczi(1)
               vyx = -yr*ftm2i(1)+yix*frcxi(1)+yiy*frcyi(1)+yiz*frczi(1)
               vzx = -zr*ftm2i(1)+zix*frcxi(1)+ziy*frcyi(1)+ziz*frczi(1)
               vyy = -yr*ftm2i(2)+yix*frcxi(2)+yiy*frcyi(2)+yiz*frczi(2)
               vzy = -zr*ftm2i(2)+zix*frcxi(2)+ziy*frcyi(2)+ziz*frczi(2)
               vzz = -zr*ftm2i(3)+zix*frcxi(3)+ziy*frcyi(3)+ziz*frczi(3)
               vir(1,1) = vir(1,1) + vxx
               vir(2,1) = vir(2,1) + vyx
               vir(3,1) = vir(3,1) + vzx
               vir(1,2) = vir(1,2) + vyx
               vir(2,2) = vir(2,2) + vyy
               vir(3,2) = vir(3,2) + vzy
               vir(1,3) = vir(1,3) + vzx
               vir(2,3) = vir(2,3) + vzy
               vir(3,3) = vir(3,3) + vzz
            end if
         end do
         end do
c
c     reset interaction scaling coefficients for connected atoms
c
         do j = 1, n12(ii)
            mscale(i12(j,ii)) = 1.0d0
            pscale(i12(j,ii)) = 1.0d0
         end do
         do j = 1, n13(ii)
            mscale(i13(j,ii)) = 1.0d0
            pscale(i13(j,ii)) = 1.0d0
         end do
         do j = 1, n14(ii)
            mscale(i14(j,ii)) = 1.0d0
            pscale(i14(j,ii)) = 1.0d0
         end do
         do j = 1, n15(ii)
            mscale(i15(j,ii)) = 1.0d0
            pscale(i15(j,ii)) = 1.0d0
         end do
         do j = 1, np11(ii)
            dscale(ip11(j,ii)) = 1.0d0
            uscale(ip11(j,ii)) = 1.0d0
         end do
         do j = 1, np12(ii)
            dscale(ip12(j,ii)) = 1.0d0
            uscale(ip12(j,ii)) = 1.0d0
         end do
         do j = 1, np13(ii)
            dscale(ip13(j,ii)) = 1.0d0
            uscale(ip13(j,ii)) = 1.0d0
         end do
         do j = 1, np14(ii)
            dscale(ip14(j,ii)) = 1.0d0
            uscale(ip14(j,ii)) = 1.0d0
         end do
      end do
      end if
c
c     for periodic boundary conditions with large cutoffs
c     neighbors must be found by the replicates method
c
c
c     perform deallocation of some local arrays
c
      deallocate (mscale)
      deallocate (pscale)
      deallocate (dscale)
      deallocate (uscale)
      return
      end

c
c     ##################################################################
c     ##                                                              ##
c     ##  subroutine udirect2a_3b   Ewald real direct field via loop  ##
c     ##  Liam O-Suilleabhain, Omar Demerdash, Teresa Head-Gordon     ##
c     ##                 Spring 2013                                  ##
c     ##                                                              ##
c     ##################################################################
c
c
c     "udirect2a_3b" computes the real space contribution of the permanent
c     atomic multipole moments to the field via a double loop
c
c
      subroutine udirect2a_3b(field,fieldp,npole3b,pnum)
      implicit none
      include 'sizes.i'
      include 'atoms.i'
      include 'boxes.i'
      include 'bound.i'
      include 'cell.i'
      include 'couple.i'
c      include 'combo.i'
      include 'ewald.i'
      include 'math.i'
      include 'mpole.i'
c      include 'polar.i'
      include 'polar2.i'
      include 'polgrp.i'
      include 'polpot.i'
      include 'shunt.i'
      include 'units.i'
      integer i,j,k,m
      integer ii,kk,l1,l2
      real*8 xr,yr,zr,r,r2
      real*8 erfc,bfac,exp2a
      real*8 drr3,drr5,drr7
      real*8 prr3,prr5,prr7
      real*8 ci,dix,diy,diz
      real*8 qixx,qiyy,qizz
      real*8 qixy,qixz,qiyz
      real*8 ck,dkx,dky,dkz
      real*8 qkxx,qkyy,qkzz
      real*8 qkxy,qkxz,qkyz
      real*8 dir,dkr
      real*8 qix,qiy,qiz,qir
      real*8 qkx,qky,qkz,qkr
      real*8 ralpha
      real*8 alsq2,alsq2n
      real*8 pdi,pti,pgamma
      real*8 damp,expdamp
      real*8 scale3,scale5
      real*8 scale7
      real*8 dsc3,dsc5,dsc7
      real*8 psc3,psc5,psc7
      real*8 bn(0:3)
      real*8 fim(3),fkm(3)
      real*8 fid(3),fkd(3)
      real*8 fip(3),fkp(3)
      real*8, allocatable :: dscale(:)
      real*8, allocatable :: pscale(:)
      real*8 field(3,*),aewald,aewald3b
      real*8 fieldp(3,*)
      integer npole3b,pnum(*)
      character*6 mode
      external erfc
      aewald=aewald3b
c
c
c     check for multipoles and set cutoff coefficients
c
      if (npole3b .eq. 0)  return
      mode = 'EWALD'
      call switch (mode)
c
c     perform dynamic allocation of some local arrays
c
      allocate (dscale(n))
      allocate (pscale(n))
c
c     set arrays needed to scale connected atom interactions
c
      do i = 1, n
         pscale(i) = 1.0d0
         dscale(i) = 1.0d0
      end do
c
c     compute the real space portion of the Ewald summation
c
      do l1 = 1,npole3b-1
         i=pnum(l1)
         ii = ipole(i)
         pdi = pdamp(i)
         pti = thole(i)
         ci = rpole(1,i)
         dix = rpole(2,i)
         diy = rpole(3,i)
         diz = rpole(4,i)
         qixx = rpole(5,i)
         qixy = rpole(6,i)
         qixz = rpole(7,i)
         qiyy = rpole(9,i)
         qiyz = rpole(10,i)
         qizz = rpole(13,i)
         do j = 1, n12(ii)
            pscale(i12(j,ii)) = p2scale
         end do
         do j = 1, n13(ii)
            pscale(i13(j,ii)) = p3scale
         end do
         do j = 1, n14(ii)
            pscale(i14(j,ii)) = p4scale
            do k = 1, np11(ii)
               if (i14(j,ii) .eq. ip11(k,ii))
     &            pscale(i14(j,ii)) = p4scale * p41scale
            end do
         end do
         do j = 1, n15(ii)
            pscale(i15(j,ii)) = p5scale
         end do
         do j = 1, np11(ii)
            dscale(ip11(j,ii)) = d1scale
         end do
         do j = 1, np12(ii)
            dscale(ip12(j,ii)) = d2scale
         end do
         do j = 1, np13(ii)
            dscale(ip13(j,ii)) = d3scale
         end do
         do j = 1, np14(ii)
            dscale(ip14(j,ii)) = d4scale
         end do
         do l2 = l1+1,npole3b
            k=pnum(l2)
            kk = ipole(k)
            xr = x(kk) - x(ii)
            yr = y(kk) - y(ii)
            zr = z(kk) - z(ii)
            call image (xr,yr,zr)
            r2 = xr*xr + yr* yr + zr*zr
            if (r2 .le. cut2) then
               r = sqrt(r2)
               ck = rpole(1,k)
               dkx = rpole(2,k)
               dky = rpole(3,k)
               dkz = rpole(4,k)
               qkxx = rpole(5,k)
               qkxy = rpole(6,k)
               qkxz = rpole(7,k)
               qkyy = rpole(9,k)
               qkyz = rpole(10,k)
               qkzz = rpole(13,k)
c
c     calculate the error function damping terms
c
               ralpha = aewald * r
               bn(0) = erfc(ralpha) / r
               alsq2 = 2.0d0 * aewald**2
               alsq2n = 0.0d0
               if (aewald .gt. 0.0d0)
     &            alsq2n = 1.0d0 / (sqrtpi*aewald)
               exp2a = exp(-ralpha**2)
               do j = 1, 3
                  bfac = dble(j+j-1)
                  alsq2n = alsq2 * alsq2n
                  bn(j) = (bfac*bn(j-1)+alsq2n*exp2a) / r2
               end do
c
c     compute the error function scaled and unscaled terms
c
               scale3 = 1.0d0
               scale5 = 1.0d0
               scale7 = 1.0d0
               damp = pdi * pdamp(k)
               if (damp .ne. 0.0d0) then
                  pgamma = min(pti,thole(k))
                  damp = -pgamma * (r/damp)**3
                  if (damp .gt. -50.0d0) then
                     expdamp = exp(damp)
                     scale3 = 1.0d0 - expdamp
                     scale5 = 1.0d0 - expdamp*(1.0d0-damp)
                     scale7 = 1.0d0 - expdamp
     &                           *(1.0d0-damp+0.6d0*damp**2)
                  end if
               end if
               dsc3 = scale3 * dscale(kk)
               dsc5 = scale5 * dscale(kk)
               dsc7 = scale7 * dscale(kk)
               psc3 = scale3 * pscale(kk)
               psc5 = scale5 * pscale(kk)
               psc7 = scale7 * pscale(kk)
               drr3 = (1.0d0-dsc3) / (r*r2)
               drr5 = 3.0d0 * (1.0d0-dsc5) / (r*r2*r2)
               drr7 = 15.0d0 * (1.0d0-dsc7) / (r*r2*r2*r2)
               prr3 = (1.0d0-psc3) / (r*r2)
               prr5 = 3.0d0 * (1.0d0-psc5) / (r*r2*r2)
               prr7 = 15.0d0 * (1.0d0-psc7) / (r*r2*r2*r2)
               dir = dix*xr + diy*yr + diz*zr
               qix = qixx*xr + qixy*yr + qixz*zr
               qiy = qixy*xr + qiyy*yr + qiyz*zr
               qiz = qixz*xr + qiyz*yr + qizz*zr
               qir = qix*xr + qiy*yr + qiz*zr
               dkr = dkx*xr + dky*yr + dkz*zr
               qkx = qkxx*xr + qkxy*yr + qkxz*zr
               qky = qkxy*xr + qkyy*yr + qkyz*zr
               qkz = qkxz*xr + qkyz*yr + qkzz*zr
               qkr = qkx*xr + qky*yr + qkz*zr
               fim(1) = -xr*(bn(1)*ck-bn(2)*dkr+bn(3)*qkr)
     &                     - bn(1)*dkx + 2.0d0*bn(2)*qkx
               fim(2) = -yr*(bn(1)*ck-bn(2)*dkr+bn(3)*qkr)
     &                     - bn(1)*dky + 2.0d0*bn(2)*qky
               fim(3) = -zr*(bn(1)*ck-bn(2)*dkr+bn(3)*qkr)
     &                     - bn(1)*dkz + 2.0d0*bn(2)*qkz
               fkm(1) = xr*(bn(1)*ci+bn(2)*dir+bn(3)*qir)
     &                     - bn(1)*dix - 2.0d0*bn(2)*qix
               fkm(2) = yr*(bn(1)*ci+bn(2)*dir+bn(3)*qir)
     &                     - bn(1)*diy - 2.0d0*bn(2)*qiy
               fkm(3) = zr*(bn(1)*ci+bn(2)*dir+bn(3)*qir)
     &                     - bn(1)*diz - 2.0d0*bn(2)*qiz
               fid(1) = -xr*(drr3*ck-drr5*dkr+drr7*qkr)
     &                     - drr3*dkx + 2.0d0*drr5*qkx
               fid(2) = -yr*(drr3*ck-drr5*dkr+drr7*qkr)
     &                     - drr3*dky + 2.0d0*drr5*qky
               fid(3) = -zr*(drr3*ck-drr5*dkr+drr7*qkr)
     &                     - drr3*dkz + 2.0d0*drr5*qkz
               fkd(1) = xr*(drr3*ci+drr5*dir+drr7*qir)
     &                     - drr3*dix - 2.0d0*drr5*qix
               fkd(2) = yr*(drr3*ci+drr5*dir+drr7*qir)
     &                     - drr3*diy - 2.0d0*drr5*qiy
               fkd(3) = zr*(drr3*ci+drr5*dir+drr7*qir)
     &                     - drr3*diz - 2.0d0*drr5*qiz
               fip(1) = -xr*(prr3*ck-prr5*dkr+prr7*qkr)
     &                     - prr3*dkx + 2.0d0*prr5*qkx
               fip(2) = -yr*(prr3*ck-prr5*dkr+prr7*qkr)
     &                     - prr3*dky + 2.0d0*prr5*qky
               fip(3) = -zr*(prr3*ck-prr5*dkr+prr7*qkr)
     &                     - prr3*dkz + 2.0d0*prr5*qkz
               fkp(1) = xr*(prr3*ci+prr5*dir+prr7*qir)
     &                     - prr3*dix - 2.0d0*prr5*qix
               fkp(2) = yr*(prr3*ci+prr5*dir+prr7*qir)
     &                     - prr3*diy - 2.0d0*prr5*qiy
               fkp(3) = zr*(prr3*ci+prr5*dir+prr7*qir)
     &                     - prr3*diz - 2.0d0*prr5*qiz
c
c     increment the field at each site due to this interaction
c
               do j = 1, 3
c                  field(j,i) = field(j,i) + fim(j) - fid(j)
c                  field(j,k) = field(j,k) + fkm(j) - fkd(j)
c                  fieldp(j,i) = fieldp(j,i) + fim(j) - fip(j)
c                  fieldp(j,k) = fieldp(j,k) + fkm(j) - fkp(j)
                  field(j,l1) = field(j,l1) + fim(j) - fid(j)
                  field(j,l2) = field(j,l2) + fkm(j) - fkd(j)
                  fieldp(j,l1) = fieldp(j,l1) + fim(j) - fip(j)
                  fieldp(j,l2) = fieldp(j,l2) + fkm(j) - fkp(j)
               end do
            end if
         end do
c
c     reset interaction scaling coefficients for connected atoms
c
         do j = 1, n12(ii)
            pscale(i12(j,ii)) = 1.0d0
         end do
         do j = 1, n13(ii)
            pscale(i13(j,ii)) = 1.0d0
         end do
         do j = 1, n14(ii)
            pscale(i14(j,ii)) = 1.0d0
         end do
         do j = 1, n15(ii)
            pscale(i15(j,ii)) = 1.0d0
         end do
         do j = 1, np11(ii)
            dscale(ip11(j,ii)) = 1.0d0
         end do
         do j = 1, np12(ii)
            dscale(ip12(j,ii)) = 1.0d0
         end do
         do j = 1, np13(ii)
            dscale(ip13(j,ii)) = 1.0d0
         end do
         do j = 1, np14(ii)
            dscale(ip14(j,ii)) = 1.0d0
         end do
      end do
c
c     periodic boundary for large cutoffs via replicates method
c
      if (use_replica) then
         do l1 = 1, npole3b
            i=pnum(l1)
            ii = ipole(i)
            pdi = pdamp(i)
            pti = thole(i)
            ci = rpole(1,i)
            dix = rpole(2,i)
            diy = rpole(3,i)
            diz = rpole(4,i)
            qixx = rpole(5,i)
            qixy = rpole(6,i)
            qixz = rpole(7,i)
            qiyy = rpole(9,i)
            qiyz = rpole(10,i)
            qizz = rpole(13,i)
            do j = 1, n12(ii)
               pscale(i12(j,ii)) = p2scale
            end do
            do j = 1, n13(ii)
               pscale(i13(j,ii)) = p3scale
            end do
            do j = 1, n14(ii)
               pscale(i14(j,ii)) = p4scale
            end do
            do j = 1, n15(ii)
               pscale(i15(j,ii)) = p5scale
            end do
            do j = 1, np11(ii)
               dscale(ip11(j,ii)) = d1scale
            end do
            do j = 1, np12(ii)
               dscale(ip12(j,ii)) = d2scale
            end do
            do j = 1, np13(ii)
               dscale(ip13(j,ii)) = d3scale
            end do
            do j = 1, np14(ii)
               dscale(ip14(j,ii)) = d4scale
            end do
            do l2 = l1, npole3b
               k=pnum(l2)
               kk = ipole(k)
               ck = rpole(1,k)
               dkx = rpole(2,k)
               dky = rpole(3,k)
               dkz = rpole(4,k)
               qkxx = rpole(5,k)
               qkxy = rpole(6,k)
               qkxz = rpole(7,k)
               qkyy = rpole(9,k)
               qkyz = rpole(10,k)
               qkzz = rpole(13,k)
               do m = 1, ncell
                  xr = x(kk) - x(ii)
                  yr = y(kk) - y(ii)
                  zr = z(kk) - z(ii)
                  call imager (xr,yr,zr,m)
                  r2 = xr*xr + yr* yr + zr*zr
c
c     calculate the error function damping terms
c
                  if (r2 .le. cut2) then
                     r = sqrt(r2)
                     ralpha = aewald * r
                     bn(0) = erfc(ralpha) / r
                     alsq2 = 2.0d0 * aewald**2
                     alsq2n = 0.0d0
                     if (aewald .gt. 0.0d0)
     &                  alsq2n = 1.0d0 / (sqrtpi*aewald)
                     exp2a = exp(-ralpha**2)
                     do j = 1, 3
                        bfac = dble(j+j-1)
                        alsq2n = alsq2 * alsq2n
                        bn(j) = (bfac*bn(j-1)+alsq2n*exp2a) / r2
                     end do
c
c     compute the error function scaled and unscaled terms
c
                     scale3 = 1.0d0
                     scale5 = 1.0d0
                     scale7 = 1.0d0
                     damp = pdi * pdamp(k)
                     if (damp .ne. 0.0d0) then
                        pgamma = min(pti,thole(k))
                        damp = -pgamma * (r/damp)**3
                        if (damp .gt. -50.0d0) then
                           expdamp = exp(damp)
                           scale3 = 1.0d0 - expdamp
                           scale5 = 1.0d0 - expdamp*(1.0d0-damp)
                           scale7 = 1.0d0 - expdamp
     &                                 *(1.0d0-damp+0.6d0*damp**2)
                        end if
                     end if
                     dsc3 = scale3
                     dsc5 = scale5
                     dsc7 = scale7
                     psc3 = scale3
                     psc5 = scale5
                     psc7 = scale7
                     if (use_polymer) then
                        if (r2 .le. polycut2) then
                           dsc3 = scale3 * dscale(kk)
                           dsc5 = scale5 * dscale(kk)
                           dsc7 = scale7 * dscale(kk)
                           psc3 = scale3 * pscale(kk)
                           psc5 = scale5 * pscale(kk)
                           psc7 = scale7 * pscale(kk)
                        end if
                     end if
                     drr3 = (1.0d0-dsc3) / (r*r2)
                     drr5 = 3.0d0 * (1.0d0-dsc5) / (r*r2*r2)
                     drr7 = 15.0d0 * (1.0d0-dsc7) / (r*r2*r2*r2)
                     prr3 = (1.0d0-psc3) / (r*r2)
                     prr5 = 3.0d0 * (1.0d0-psc5) / (r*r2*r2)
                     prr7 = 15.0d0 * (1.0d0-psc7) / (r*r2*r2*r2)
                     dir = dix*xr + diy*yr + diz*zr
                     qix = qixx*xr + qixy*yr + qixz*zr
                     qiy = qixy*xr + qiyy*yr + qiyz*zr
                     qiz = qixz*xr + qiyz*yr + qizz*zr
                     qir = qix*xr + qiy*yr + qiz*zr
                     dkr = dkx*xr + dky*yr + dkz*zr
                     qkx = qkxx*xr + qkxy*yr + qkxz*zr
                     qky = qkxy*xr + qkyy*yr + qkyz*zr
                     qkz = qkxz*xr + qkyz*yr + qkzz*zr
                     qkr = qkx*xr + qky*yr + qkz*zr
                     fim(1) = -xr*(bn(1)*ck-bn(2)*dkr+bn(3)*qkr)
     &                           - bn(1)*dkx + 2.0d0*bn(2)*qkx
                     fim(2) = -yr*(bn(1)*ck-bn(2)*dkr+bn(3)*qkr)
     &                           - bn(1)*dky + 2.0d0*bn(2)*qky
                     fim(3) = -zr*(bn(1)*ck-bn(2)*dkr+bn(3)*qkr)
     &                           - bn(1)*dkz + 2.0d0*bn(2)*qkz
                     fkm(1) = xr*(bn(1)*ci+bn(2)*dir+bn(3)*qir)
     &                           - bn(1)*dix - 2.0d0*bn(2)*qix
                     fkm(2) = yr*(bn(1)*ci+bn(2)*dir+bn(3)*qir)
     &                           - bn(1)*diy - 2.0d0*bn(2)*qiy
                     fkm(3) = zr*(bn(1)*ci+bn(2)*dir+bn(3)*qir)
     &                           - bn(1)*diz - 2.0d0*bn(2)*qiz
                     fid(1) = -xr*(drr3*ck-drr5*dkr+drr7*qkr)
     &                           - drr3*dkx + 2.0d0*drr5*qkx
                     fid(2) = -yr*(drr3*ck-drr5*dkr+drr7*qkr)
     &                           - drr3*dky + 2.0d0*drr5*qky
                     fid(3) = -zr*(drr3*ck-drr5*dkr+drr7*qkr)
     &                           - drr3*dkz + 2.0d0*drr5*qkz
                     fkd(1) = xr*(drr3*ci+drr5*dir+drr7*qir)
     &                           - drr3*dix - 2.0d0*drr5*qix
                     fkd(2) = yr*(drr3*ci+drr5*dir+drr7*qir)
     &                           - drr3*diy - 2.0d0*drr5*qiy
                     fkd(3) = zr*(drr3*ci+drr5*dir+drr7*qir)
     &                           - drr3*diz - 2.0d0*drr5*qiz
                     fip(1) = -xr*(prr3*ck-prr5*dkr+prr7*qkr)
     &                           - prr3*dkx + 2.0d0*prr5*qkx
                     fip(2) = -yr*(prr3*ck-prr5*dkr+prr7*qkr)
     &                           - prr3*dky + 2.0d0*prr5*qky
                     fip(3) = -zr*(prr3*ck-prr5*dkr+prr7*qkr)
     &                           - prr3*dkz + 2.0d0*prr5*qkz
                     fkp(1) = xr*(prr3*ci+prr5*dir+prr7*qir)
     &                           - prr3*dix - 2.0d0*prr5*qix
                     fkp(2) = yr*(prr3*ci+prr5*dir+prr7*qir)
     &                           - prr3*diy - 2.0d0*prr5*qiy
                     fkp(3) = zr*(prr3*ci+prr5*dir+prr7*qir)
     &                           - prr3*diz - 2.0d0*prr5*qiz
c
c     increment the field at each site due to this interaction
c
                     do j = 1, 3
c                        field(j,i) = field(j,i) + fim(j) - fid(j)
c                        fieldp(j,i) = fieldp(j,i) + fim(j) - fip(j)
                        field(j,l1) = field(j,l1) + fim(j) - fid(j)
                        fieldp(j,l1) = fieldp(j,l1) + fim(j) - fip(j)

                        if (ii .ne. kk) then
c                           field(j,k) = field(j,k) + fkm(j) - fkd(j)
c                           fieldp(j,k) = fieldp(j,k) + fkm(j) - fkp(j)
                           field(j,l2) = field(j,l2) + fkm(j) - fkd(j)
                           fieldp(j,l2) = fieldp(j,l2) + fkm(j) - fkp(j)
                        end if
                     end do
                  end if
               end do
            end do
c
c     reset interaction scaling coefficients for connected atoms
c
            do j = 1, n12(ii)
               pscale(i12(j,ii)) = 1.0d0
            end do
            do j = 1, n13(ii)
               pscale(i13(j,ii)) = 1.0d0
            end do
            do j = 1, n14(ii)
               pscale(i14(j,ii)) = 1.0d0
            end do
            do j = 1, n15(ii)
               pscale(i15(j,ii)) = 1.0d0
            end do
            do j = 1, np11(ii)
               dscale(ip11(j,ii)) = 1.0d0
            end do
            do j = 1, np12(ii)
               dscale(ip12(j,ii)) = 1.0d0
            end do
            do j = 1, np13(ii)
               dscale(ip13(j,ii)) = 1.0d0
            end do
            do j = 1, np14(ii)
               dscale(ip14(j,ii)) = 1.0d0
            end do
         end do
      end if
c
c     perform deallocation of some local arrays
c
      deallocate (dscale)
      deallocate (pscale)
      return
      end
c
c     ################################################################
c     ##                                                            ##
c     ##  subroutine torque_reg3b  --  convert multipole torque to force  ##
c     ##                                                            ##
c     ################################################################
c
c
c     "torque_reg" takes the torque values on sites defined by local
c     coordinate frames and distributes them to convert to forces
c     on the original sites and sites specifying the local frames
c
c     literature reference:
c
c     A. J. Stone and M. Alderton, "Distributed Multipole Analysis -
c     Methods and Applications", Molecular Physics, 56, 1047-1064 (1985)
c
c
c      subroutine torque_reg3b (trq,deptemp,pnum,pnum2,npole3b)
      subroutine torque_reg3b (trq,deptemp,pnum,npole3b)
      implicit none
      include 'sizes.i'
      include 'atoms.i'
      include 'mpole.i'
c      include 'combo.i'
      integer i,j,ia,ib,ic,l1,l3
      real*8 usiz,vsiz,wsiz
      real*8 upsiz,vpsiz
      real*8 dotdu,dotdv
      real*8 dphidu,dphidv,dphidw
      real*8 c,s,uvdis,vudis,du,dv
      real*8 u(3),v(3),w(3)
      real*8 up(3),vp(3),diff(3)
      real*8 trq(3,*)
      real*8 deptemp(3,*)
c      integer pnum(*),npole3b,pnum2(*)
      integer pnum(*),npole3b
      integer ia_l1,ib_l1,ic_l1
c
c     coordinate frame motion described by rotation about u, v and w
c
      do l1 = 1, npole3b
         i = pnum(l1)
         ia = zaxis(i)
         ib = ipole(i)
         ic = xaxis(i)
c         pnum2(ia)=ia_l1
c         pnum2(ib)=ib_l1
c         pnum2(ic)=ic_l1
         u(1) = x(ia) - x(ib)
         u(2) = y(ia) - y(ib)
         u(3) = z(ia) - z(ib)
         usiz = sqrt(u(1)*u(1) + u(2)*u(2) + u(3)*u(3))
         v(1) = x(ic) - x(ib)
         v(2) = y(ic) - y(ib)
         v(3) = z(ic) - z(ib)
         vsiz = sqrt(v(1)*v(1) + v(2)*v(2) + v(3)*v(3))
         w(1) = u(2)*v(3) - u(3)*v(2)
         w(2) = u(3)*v(1) - u(1)*v(3)
         w(3) = u(1)*v(2) - u(2)*v(1)
         wsiz = sqrt(w(1)*w(1) + w(2)*w(2) + w(3)*w(3))
         dotdu = 0.0d0
         dotdv = 0.0d0
         do j = 1, 3
            u(j) = u(j) / usiz
            v(j) = v(j) / vsiz
            w(j) = w(j) / wsiz
            diff(j) = v(j) - u(j)
            dotdu = dotdu + u(j)*diff(j)
            dotdv = dotdv + v(j)*diff(j)
         end do
c
c     get perpendiculars to u,v to get direction of motion of
c     u or v due to rotation about the cross product vector w
c
         upsiz = 0.0d0
         vpsiz = 0.0d0
         do j = 1, 3
            up(j) = diff(j) - dotdu*u(j)
            vp(j) = diff(j) - dotdv*v(j)
            upsiz = upsiz + up(j)*up(j)
            vpsiz = vpsiz + vp(j)*vp(j)
         end do
         upsiz = sqrt(upsiz)
         vpsiz = sqrt(vpsiz)
         do j = 1, 3
            up(j) = up(j) / upsiz
            vp(j) = vp(j) / vpsiz
         end do
c
c     negative of dot product of torque with unit vectors along u, v
c     and w give result of infinitesmal rotation along these vectors
c     i.e. dphi/dtheta = dot product, where dphi is work, dtheta is
c     angle; then dphi/dtheta is torque and the dot product is torque
c     component along unit vector
c
c         dphidu = -trq(1,ib_l1)*u(1) - trq(2,ib_l1)*u(2)
c     &           - trq(3,ib_l1)*u(3)
c         dphidv = -trq(1,ib_l1)*v(1) - trq(2,ib_l1)*v(2) 
c     &           - trq(3,ib_l1)*v(3)
c         dphidw = -trq(1,ib_l1)*w(1) - trq(2,ib_l1)*w(2) 
c     &            - trq(3,ib_l1)*w(3)
         dphidu = -trq(1,ib)*u(1) - trq(2,ib)*u(2) - trq(3,ib)*u(3)
         dphidv = -trq(1,ib)*v(1) - trq(2,ib)*v(2) - trq(3,ib)*v(3)
         dphidw = -trq(1,ib)*w(1) - trq(2,ib)*w(2) - trq(3,ib)*w(3)

c
c     get the projected distances between the vectors
c
         c = u(1)*v(1) + u(2)*v(2) + u(3)*v(3)
         s = sqrt(1.0d0 - c*c)
         uvdis = usiz * s
         vudis = vsiz * s
c
c     force distribution for the bisector local coordinate method
c
         if (polaxe(i) .eq. 'Bisector') then
            do j = 1, 3
               du = -w(j)*dphidv/uvdis + up(j)*dphidw/(2.0d0*usiz)
               dv = w(j)*dphidu/vudis + vp(j)*dphidw/(2.0d0*vsiz)
c               deptemp(j,ia_l1) = deptemp(j,ia_l1) + du
c               deptemp(j,ic_l1) = deptemp(j,ic_l1) + dv
c               deptemp(j,ib_l1) = deptemp(j,ib_l1) - dv - du
               deptemp(j,ia) = deptemp(j,ia) + du
               deptemp(j,ic) = deptemp(j,ic) + dv
               deptemp(j,ib) = deptemp(j,ib) - dv - du
                
            end do
c
c     force distribution for the Z-then-X local coordinate method
c
         else if (polaxe(i) .eq. 'Z-then-X') then
            do j = 1, 3
               du = -w(j)*dphidv/uvdis + up(j)*dphidw/usiz
               dv = w(j)*dphidu/vudis
c               deptemp(j,ia_l1) = deptemp(j,ia_l1) + du
c               deptemp(j,ic_l1) = deptemp(j,ic_l1) + dv
c               deptemp(j,ib_l1) = deptemp(j,ib_l1) - dv - du
               deptemp(j,ia) = deptemp(j,ia) + du
               deptemp(j,ic) = deptemp(j,ic) + dv
               deptemp(j,ib) = deptemp(j,ib) - dv - du
            end do
         end if
      end do
      return
      end


c
c
c     #############################################################
c     ##  COPYRIGHT (C) 2007 by Pengyu Ren & Jay William Ponder  ##
c     ##                   All Rights Reserved                   ##
c     #############################################################
c
c     ##################################################################
c     ##                                                              ##
c     ##  subroutine torque_3b  --  convert single site torque to force  ##
c     ##                                                              ##
c     ##################################################################
c
c
c     "torque" takes the torque values on a single site defined by
c     a local coordinate frame and converts to Cartesian forces on
c     the original site and sites specifying the local frame
c
c     literature reference:
c
c     P. L. Popelier and A. J. Stone, "Formulae for the First and
c     Second Derivatives of Anisotropic Potentials with Respect to
c     Geometrical Parameters", Molecular Physics, 82, 411-425 (1994)
c
c
c               call torque_3b (pnum2,deptemp,i,
c     &           ttm2,ttm2i,frcxi,frcyi,frczi)

c      subroutine torque_3b(pnum2,deptemp,
c     &    i,trq1,trq2,frcx,frcy,frcz)
      subroutine torque_3b(deptemp,
     &    i,trq1,trq2,frcx,frcy,frcz)
      implicit none
      include 'sizes.i'
      include 'atoms.i'
      include 'deriv.i'
      include 'mpole.i'
      integer i,j,l1
      integer ia,ib,ic,id
      real*8 du,dv,dw,random
      real*8 usiz,vsiz,wsiz
      real*8 rsiz,ssiz
      real*8 t1siz,t2siz
      real*8 uvsiz,uwsiz,vwsiz
      real*8 ursiz,ussiz
      real*8 vssiz,wssiz
      real*8 uvcos,uwcos,vwcos
      real*8 urcos,uscos
      real*8 vscos,wscos
      real*8 ut1cos,ut2cos
      real*8 uvsin,uwsin,vwsin
      real*8 ursin,ussin
      real*8 vssin,wssin
      real*8 ut1sin,ut2sin
      real*8 dphidu,dphidv,dphidw
      real*8 dphidr,dphids
      real*8 trq1(3),trq2(3)
      real*8 frcx(3),frcy(3),frcz(3)
      real*8 u(3),v(3),w(3)
      real*8 r(3),s(3)
      real*8 t1(3),t2(3)
      real*8 uv(3),uw(3),vw(3)
      real*8 ur(3),us(3)
      real*8 vs(3),ws(3)
      character*8 axetyp
      real*8 deptemp(3,*)
c      integer pnum2(*)
      integer ia_l1,ib_l1,ic_l1
      integer id_l1
c
c     zero out force components on local frame-defining atoms
c
      do j = 1, 3
         frcz(j) = 0.0d0
         frcx(j) = 0.0d0
         frcy(j) = 0.0d0
      end do
c
c     get the local frame type and the frame-defining atoms
c
c      i=pnum(l1)
      
      ia = zaxis(i)
      ib = ipole(i)
      ic = xaxis(i)
      id = yaxis(i)
c      pnum2(ia)=ia_l1
c      pnum2(ib)=ib_l1
c      pnum2(ic)=ic_l1
c      pnum2(id)=id_l1

      axetyp = polaxe(i)
      if (axetyp .eq. 'None')  return
c
c     construct the three rotation axes for the local frame
c
      u(1) = x(ia) - x(ib)
      u(2) = y(ia) - y(ib)
      u(3) = z(ia) - z(ib)
      if (axetyp .ne. 'Z-Only') then
         v(1) = x(ic) - x(ib)
         v(2) = y(ic) - y(ib)
         v(3) = z(ic) - z(ib)
      else
         v(1) = random ()
         v(2) = random ()
         v(3) = random ()
      end if
      if (axetyp.eq.'Z-Bisect' .or. axetyp.eq.'3-Fold') then
         w(1) = x(id) - x(ib)
         w(2) = y(id) - y(ib)
         w(3) = z(id) - z(ib)
      else
         w(1) = u(2)*v(3) - u(3)*v(2)
         w(2) = u(3)*v(1) - u(1)*v(3)
         w(3) = u(1)*v(2) - u(2)*v(1)
      end if
      usiz = sqrt(u(1)*u(1) + u(2)*u(2) + u(3)*u(3))
      vsiz = sqrt(v(1)*v(1) + v(2)*v(2) + v(3)*v(3))
      wsiz = sqrt(w(1)*w(1) + w(2)*w(2) + w(3)*w(3))
      do j = 1, 3
         u(j) = u(j) / usiz
         v(j) = v(j) / vsiz
         w(j) = w(j) / wsiz
      end do
c
c     build some additional axes needed for the Z-Bisect method
c
      if (axetyp .eq. 'Z-Bisect') then
         r(1) = v(1) + w(1)
         r(2) = v(2) + w(2)
         r(3) = v(3) + w(3)
         s(1) = u(2)*r(3) - u(3)*r(2)
         s(2) = u(3)*r(1) - u(1)*r(3)
         s(3) = u(1)*r(2) - u(2)*r(1)
         rsiz = sqrt(r(1)*r(1) + r(2)*r(2) + r(3)*r(3))
         ssiz = sqrt(s(1)*s(1) + s(2)*s(2) + s(3)*s(3))
         do j = 1, 3
            r(j) = r(j) / rsiz
            s(j) = s(j) / ssiz
         end do
      end if
c
c     find the perpendicular and angle for each pair of axes
c
      uv(1) = v(2)*u(3) - v(3)*u(2)
      uv(2) = v(3)*u(1) - v(1)*u(3)
      uv(3) = v(1)*u(2) - v(2)*u(1)
      uw(1) = w(2)*u(3) - w(3)*u(2)
      uw(2) = w(3)*u(1) - w(1)*u(3)
      uw(3) = w(1)*u(2) - w(2)*u(1)
      vw(1) = w(2)*v(3) - w(3)*v(2)
      vw(2) = w(3)*v(1) - w(1)*v(3)
      vw(3) = w(1)*v(2) - w(2)*v(1)
      uvsiz = sqrt(uv(1)*uv(1) + uv(2)*uv(2) + uv(3)*uv(3))
      uwsiz = sqrt(uw(1)*uw(1) + uw(2)*uw(2) + uw(3)*uw(3))
      vwsiz = sqrt(vw(1)*vw(1) + vw(2)*vw(2) + vw(3)*vw(3))
      do j = 1, 3
         uv(j) = uv(j) / uvsiz
         uw(j) = uw(j) / uwsiz
         vw(j) = vw(j) / vwsiz
      end do
      if (axetyp .eq. 'Z-Bisect') then
         ur(1) = r(2)*u(3) - r(3)*u(2)
         ur(2) = r(3)*u(1) - r(1)*u(3)
         ur(3) = r(1)*u(2) - r(2)*u(1)
         us(1) = s(2)*u(3) - s(3)*u(2)
         us(2) = s(3)*u(1) - s(1)*u(3)
         us(3) = s(1)*u(2) - s(2)*u(1)
         vs(1) = s(2)*v(3) - s(3)*v(2)
         vs(2) = s(3)*v(1) - s(1)*v(3)
         vs(3) = s(1)*v(2) - s(2)*v(1)
         ws(1) = s(2)*w(3) - s(3)*w(2)
         ws(2) = s(3)*w(1) - s(1)*w(3)
         ws(3) = s(1)*w(2) - s(2)*w(1)
         ursiz = sqrt(ur(1)*ur(1) + ur(2)*ur(2) + ur(3)*ur(3))
         ussiz = sqrt(us(1)*us(1) + us(2)*us(2) + us(3)*us(3))
         vssiz = sqrt(vs(1)*vs(1) + vs(2)*vs(2) + vs(3)*vs(3))
         wssiz = sqrt(ws(1)*ws(1) + ws(2)*ws(2) + ws(3)*ws(3))
         do j = 1, 3
            ur(j) = ur(j) / ursiz
            us(j) = us(j) / ussiz
            vs(j) = vs(j) / vssiz
            ws(j) = ws(j) / wssiz
         end do
      end if
c
c     get sine and cosine of angles between the rotation axes
c
      uvcos = u(1)*v(1) + u(2)*v(2) + u(3)*v(3)
      uvsin = sqrt(1.0d0 - uvcos*uvcos)
      uwcos = u(1)*w(1) + u(2)*w(2) + u(3)*w(3)
      uwsin = sqrt(1.0d0 - uwcos*uwcos)
      vwcos = v(1)*w(1) + v(2)*w(2) + v(3)*w(3)
      vwsin = sqrt(1.0d0 - vwcos*vwcos)
      if (axetyp .eq. 'Z-Bisect') then
         urcos = u(1)*r(1) + u(2)*r(2) + u(3)*r(3)
         ursin = sqrt(1.0d0 - urcos*urcos)
         uscos = u(1)*s(1) + u(2)*s(2) + u(3)*s(3)
         ussin = sqrt(1.0d0 - uscos*uscos)
         vscos = v(1)*s(1) + v(2)*s(2) + v(3)*s(3)
         vssin = sqrt(1.0d0 - vscos*vscos)
         wscos = w(1)*s(1) + w(2)*s(2) + w(3)*s(3)
         wssin = sqrt(1.0d0 - wscos*wscos)
      end if
c
c     compute the projection of v and w onto the ru-plane
c
      if (axetyp .eq. 'Z-Bisect') then
         do j = 1, 3
           t1(j) = v(j) - s(j)*vscos
           t2(j) = w(j) - s(j)*wscos
         end do
         t1siz = sqrt(t1(1)*t1(1)+t1(2)*t1(2)+t1(3)*t1(3))
         t2siz = sqrt(t2(1)*t2(1)+t2(2)*t2(2)+t2(3)*t2(3))
         do j = 1, 3
           t1(j) = t1(j) / t1siz
           t2(j) = t2(j) / t2siz
         end do
         ut1cos = u(1)*t1(1) + u(2)*t1(2) + u(3)*t1(3)
         ut1sin = sqrt(1.0d0 - ut1cos*ut1cos)
         ut2cos = u(1)*t2(1) + u(2)*t2(2) + u(3)*t2(3)
         ut2sin = sqrt(1.0d0 - ut2cos*ut2cos)
      end if
c
c     negative of dot product of torque with unit vectors gives
c     result of infinitesimal rotation along these vectors
c
      dphidu = -trq2(1)*u(1) - trq2(2)*u(2) - trq2(3)*u(3)
      dphidv = -trq2(1)*v(1) - trq2(2)*v(2) - trq2(3)*v(3)
      dphidw = -trq2(1)*w(1) - trq2(2)*w(2) - trq2(3)*w(3)
      if (axetyp .eq. 'Z-Bisect') then
         dphidr = -trq2(1)*r(1) - trq2(2)*r(2) - trq2(3)*r(3)
         dphids = -trq2(1)*s(1) - trq2(2)*s(2) - trq2(3)*s(3)
      end if
c
c     force distribution for the Z-Only local coordinate method
c
      if (axetyp .eq. 'Z-Only') then
         do j = 1, 3
            du = uv(j)*dphidv/(usiz*uvsin) + uw(j)*dphidw/usiz
            deptemp(j,ia) = deptemp(j,ia) + du
            deptemp(j,ib) = deptemp(j,ib) - du
            frcz(j) = frcz(j) + du
         end do
c
c     force distribution for the Z-then-X local coordinate method
c
      else if (axetyp .eq. 'Z-then-X') then
         do j = 1, 3
            du = uv(j)*dphidv/(usiz*uvsin) + uw(j)*dphidw/usiz
            dv = -uv(j)*dphidu/(vsiz*uvsin)
            deptemp(j,ia) = deptemp(j,ia) + du
            deptemp(j,ic) = deptemp(j,ic) + dv
            deptemp(j,ib) = deptemp(j,ib) - du - dv
            frcz(j) = frcz(j) + du
            frcx(j) = frcx(j) + dv
         end do
c
c     force distribution for the Bisector local coordinate method
c
      else if (axetyp .eq. 'Bisector') then
         do j = 1, 3
            du = uv(j)*dphidv/(usiz*uvsin) + 0.5d0*uw(j)*dphidw/usiz
            dv = -uv(j)*dphidu/(vsiz*uvsin) + 0.5d0*vw(j)*dphidw/vsiz
            deptemp(j,ia) = deptemp(j,ia) + du
            deptemp(j,ic) = deptemp(j,ic) + dv
            deptemp(j,ib) = deptemp(j,ib) - du - dv
            frcz(j) = frcz(j) + du
            frcx(j) = frcx(j) + dv
         end do
c
c     force distribution for the Z-Bisect local coordinate method
c
      else if (axetyp .eq. 'Z-Bisect') then
         do j = 1, 3
            du = ur(j)*dphidr/(usiz*ursin) + us(j)*dphids/usiz
            dv = (vssin*s(j)-vscos*t1(j))*dphidu
     &              / (vsiz*(ut1sin+ut2sin))
            dw = (wssin*s(j)-wscos*t2(j))*dphidu
     &              / (wsiz*(ut1sin+ut2sin))
            deptemp(j,ia) = deptemp(j,ia) + du
            deptemp(j,ic) = deptemp(j,ic) + dv
            deptemp(j,id) = deptemp(j,id) + dw
            deptemp(j,ib) = deptemp(j,ib) - du - dv - dw
            frcz(j) = frcz(j) + du
            frcx(j) = frcx(j) + dv
            frcy(j) = frcy(j) + dw
         end do
c
c     force distribution for the 3-Fold local coordinate method
c        (correct for uv, uw and vw angles all equal to 90)
c
      else if (axetyp .eq. '3-Fold') then
         do j = 1, 3
            du = uw(j)*dphidw/(usiz*uwsin)
     &              + uv(j)*dphidv/(usiz*uvsin)
     &              - uw(j)*dphidu/(usiz*uwsin)
     &              - uv(j)*dphidu/(usiz*uvsin)
            dv = vw(j)*dphidw/(vsiz*vwsin)
     &              - uv(j)*dphidu/(vsiz*uvsin)
     &              - vw(j)*dphidv/(vsiz*vwsin)
     &              + uv(j)*dphidv/(vsiz*uvsin)
            dw = -uw(j)*dphidu/(wsiz*uwsin)
     &              - vw(j)*dphidv/(wsiz*vwsin)
     &              + uw(j)*dphidw/(wsiz*uwsin)
     &              + vw(j)*dphidw/(wsiz*vwsin)
            du = du / 3.0d0
            dv = dv / 3.0d0
            dw = dw / 3.0d0
            deptemp(j,ia) = deptemp(j,ia) + du
            deptemp(j,ic) = deptemp(j,ic) + dv
            deptemp(j,id) = deptemp(j,id) + dw
            deptemp(j,ib) = deptemp(j,ib) - du - dv - dw
            frcz(j) = frcz(j) + du
            frcx(j) = frcx(j) + dv
            frcy(j) = frcy(j) + dw
         end do
      end if
      return
      end

c
c
c     #############################################################
c     ##  COPYRIGHT (C) 2007 by Pengyu Ren & Jay William Ponder  ##
c     ##                   All Rights Reserved                   ##
c     #############################################################
c
c     ##################################################################
c     ##                                                              ##
c     ##  subroutine torque  --  convert single site torque to force  ##
c     ##                                                              ##
c     ##################################################################
c
c
c     "torque" takes the torque values on a single site defined by
c     a local coordinate frame and converts to Cartesian forces on
c     the original site and sites specifying the local frame
c
c     literature reference:
c
c     P. L. Popelier and A. J. Stone, "Formulae for the First and
c     Second Derivatives of Anisotropic Potentials with Respect to
c     Geometrical Parameters", Molecular Physics, 82, 411-425 (1994)
c
c
      subroutine torque_3b_Perm (i,trq1,trq2,frcx,frcy,frcz)
      implicit none
      include 'sizes.i'
      include 'atoms.i'
      include 'deriv.i'
      include 'mpole.i'
      integer i,j
      integer ia,ib,ic,id
      real*8 du,dv,dw,random
      real*8 usiz,vsiz,wsiz
      real*8 rsiz,ssiz
      real*8 t1siz,t2siz
      real*8 uvsiz,uwsiz,vwsiz
      real*8 ursiz,ussiz
      real*8 vssiz,wssiz
      real*8 uvcos,uwcos,vwcos
      real*8 urcos,uscos
      real*8 vscos,wscos
      real*8 ut1cos,ut2cos
      real*8 uvsin,uwsin,vwsin
      real*8 ursin,ussin
      real*8 vssin,wssin
      real*8 ut1sin,ut2sin
      real*8 dphidu,dphidv,dphidw
      real*8 dphidr,dphids
      real*8 trq1(3),trq2(3)
      real*8 frcx(3),frcy(3),frcz(3)
      real*8 u(3),v(3),w(3)
      real*8 r(3),s(3)
      real*8 t1(3),t2(3)
      real*8 uv(3),uw(3),vw(3)
      real*8 ur(3),us(3)
      real*8 vs(3),ws(3)
      character*8 axetyp
c
c
c     zero out force components on local frame-defining atoms
c
      do j = 1, 3
         frcz(j) = 0.0d0
         frcx(j) = 0.0d0
         frcy(j) = 0.0d0
      end do
c
c     get the local frame type and the frame-defining atoms
c
      ia = zaxis(i)
      ib = ipole(i)
      ic = xaxis(i)
      id = yaxis(i)
      axetyp = polaxe(i)
      if (axetyp .eq. 'None')  return
c
c     construct the three rotation axes for the local frame
c
      u(1) = x(ia) - x(ib)
      u(2) = y(ia) - y(ib)
      u(3) = z(ia) - z(ib)
      if (axetyp .ne. 'Z-Only') then
         v(1) = x(ic) - x(ib)
         v(2) = y(ic) - y(ib)
         v(3) = z(ic) - z(ib)
      else
         v(1) = random ()
         v(2) = random ()
         v(3) = random ()
      end if
      if (axetyp.eq.'Z-Bisect' .or. axetyp.eq.'3-Fold') then
         w(1) = x(id) - x(ib)
         w(2) = y(id) - y(ib)
         w(3) = z(id) - z(ib)
      else
         w(1) = u(2)*v(3) - u(3)*v(2)
         w(2) = u(3)*v(1) - u(1)*v(3)
         w(3) = u(1)*v(2) - u(2)*v(1)
      end if
      usiz = sqrt(u(1)*u(1) + u(2)*u(2) + u(3)*u(3))
      vsiz = sqrt(v(1)*v(1) + v(2)*v(2) + v(3)*v(3))
      wsiz = sqrt(w(1)*w(1) + w(2)*w(2) + w(3)*w(3))
      do j = 1, 3
         u(j) = u(j) / usiz
         v(j) = v(j) / vsiz
         w(j) = w(j) / wsiz
      end do
c
c     build some additional axes needed for the Z-Bisect method
c
      if (axetyp .eq. 'Z-Bisect') then
         r(1) = v(1) + w(1)
         r(2) = v(2) + w(2)
         r(3) = v(3) + w(3)
         s(1) = u(2)*r(3) - u(3)*r(2)
         s(2) = u(3)*r(1) - u(1)*r(3)
         s(3) = u(1)*r(2) - u(2)*r(1)
         rsiz = sqrt(r(1)*r(1) + r(2)*r(2) + r(3)*r(3))
         ssiz = sqrt(s(1)*s(1) + s(2)*s(2) + s(3)*s(3))
         do j = 1, 3
            r(j) = r(j) / rsiz
            s(j) = s(j) / ssiz
         end do
      end if
c
c     find the perpendicular and angle for each pair of axes
c
      uv(1) = v(2)*u(3) - v(3)*u(2)
      uv(2) = v(3)*u(1) - v(1)*u(3)
      uv(3) = v(1)*u(2) - v(2)*u(1)
      uw(1) = w(2)*u(3) - w(3)*u(2)
      uw(2) = w(3)*u(1) - w(1)*u(3)
      uw(3) = w(1)*u(2) - w(2)*u(1)
      vw(1) = w(2)*v(3) - w(3)*v(2)
      vw(2) = w(3)*v(1) - w(1)*v(3)
      vw(3) = w(1)*v(2) - w(2)*v(1)
      uvsiz = sqrt(uv(1)*uv(1) + uv(2)*uv(2) + uv(3)*uv(3))
      uwsiz = sqrt(uw(1)*uw(1) + uw(2)*uw(2) + uw(3)*uw(3))
      vwsiz = sqrt(vw(1)*vw(1) + vw(2)*vw(2) + vw(3)*vw(3))
      do j = 1, 3
         uv(j) = uv(j) / uvsiz
         uw(j) = uw(j) / uwsiz
         vw(j) = vw(j) / vwsiz
      end do
      if (axetyp .eq. 'Z-Bisect') then
         ur(1) = r(2)*u(3) - r(3)*u(2)
         ur(2) = r(3)*u(1) - r(1)*u(3)
         ur(3) = r(1)*u(2) - r(2)*u(1)
         us(1) = s(2)*u(3) - s(3)*u(2)
         us(2) = s(3)*u(1) - s(1)*u(3)
         us(3) = s(1)*u(2) - s(2)*u(1)
         vs(1) = s(2)*v(3) - s(3)*v(2)
         vs(2) = s(3)*v(1) - s(1)*v(3)
         vs(3) = s(1)*v(2) - s(2)*v(1)
         ws(1) = s(2)*w(3) - s(3)*w(2)
         ws(2) = s(3)*w(1) - s(1)*w(3)
         ws(3) = s(1)*w(2) - s(2)*w(1)
         ursiz = sqrt(ur(1)*ur(1) + ur(2)*ur(2) + ur(3)*ur(3))
         ussiz = sqrt(us(1)*us(1) + us(2)*us(2) + us(3)*us(3))
         vssiz = sqrt(vs(1)*vs(1) + vs(2)*vs(2) + vs(3)*vs(3))
         wssiz = sqrt(ws(1)*ws(1) + ws(2)*ws(2) + ws(3)*ws(3))
         do j = 1, 3
            ur(j) = ur(j) / ursiz
            us(j) = us(j) / ussiz
            vs(j) = vs(j) / vssiz
            ws(j) = ws(j) / wssiz
         end do
      end if
c
c     get sine and cosine of angles between the rotation axes
c
      uvcos = u(1)*v(1) + u(2)*v(2) + u(3)*v(3)
      uvsin = sqrt(1.0d0 - uvcos*uvcos)
      uwcos = u(1)*w(1) + u(2)*w(2) + u(3)*w(3)
      uwsin = sqrt(1.0d0 - uwcos*uwcos)
      vwcos = v(1)*w(1) + v(2)*w(2) + v(3)*w(3)
      vwsin = sqrt(1.0d0 - vwcos*vwcos)
      if (axetyp .eq. 'Z-Bisect') then
         urcos = u(1)*r(1) + u(2)*r(2) + u(3)*r(3)
         ursin = sqrt(1.0d0 - urcos*urcos)
         uscos = u(1)*s(1) + u(2)*s(2) + u(3)*s(3)
         ussin = sqrt(1.0d0 - uscos*uscos)
         vscos = v(1)*s(1) + v(2)*s(2) + v(3)*s(3)
         vssin = sqrt(1.0d0 - vscos*vscos)
         wscos = w(1)*s(1) + w(2)*s(2) + w(3)*s(3)
         wssin = sqrt(1.0d0 - wscos*wscos)
      end if
c
c     compute the projection of v and w onto the ru-plane
c
      if (axetyp .eq. 'Z-Bisect') then
         do j = 1, 3
           t1(j) = v(j) - s(j)*vscos
           t2(j) = w(j) - s(j)*wscos
         end do
         t1siz = sqrt(t1(1)*t1(1)+t1(2)*t1(2)+t1(3)*t1(3))
         t2siz = sqrt(t2(1)*t2(1)+t2(2)*t2(2)+t2(3)*t2(3))
         do j = 1, 3
           t1(j) = t1(j) / t1siz
           t2(j) = t2(j) / t2siz
         end do
         ut1cos = u(1)*t1(1) + u(2)*t1(2) + u(3)*t1(3)
         ut1sin = sqrt(1.0d0 - ut1cos*ut1cos)
         ut2cos = u(1)*t2(1) + u(2)*t2(2) + u(3)*t2(3)
         ut2sin = sqrt(1.0d0 - ut2cos*ut2cos)
      end if
c
c     negative of dot product of torque with unit vectors gives
c     result of infinitesimal rotation along these vectors
c
      dphidu = -trq1(1)*u(1) - trq1(2)*u(2) - trq1(3)*u(3)
      dphidv = -trq1(1)*v(1) - trq1(2)*v(2) - trq1(3)*v(3)
      dphidw = -trq1(1)*w(1) - trq1(2)*w(2) - trq1(3)*w(3)
      if (axetyp .eq. 'Z-Bisect') then
         dphidr = -trq1(1)*r(1) - trq1(2)*r(2) - trq1(3)*r(3)
         dphids = -trq1(1)*s(1) - trq1(2)*s(2) - trq1(3)*s(3)
      end if
c
c     force distribution for the Z-Only local coordinate method
c
      if (axetyp .eq. 'Z-Only') then
         do j = 1, 3
            du = uv(j)*dphidv/(usiz*uvsin) + uw(j)*dphidw/usiz
            dem(j,ia) = dem(j,ia) + du
            dem(j,ib) = dem(j,ib) - du
c           frcz(j) = frcz(j) + du
         end do
c
c     force distribution for the Z-then-X local coordinate method
c
      else if (axetyp .eq. 'Z-then-X') then
         do j = 1, 3
            du = uv(j)*dphidv/(usiz*uvsin) + uw(j)*dphidw/usiz
            dv = -uv(j)*dphidu/(vsiz*uvsin)
            dem(j,ia) = dem(j,ia) + du
            dem(j,ic) = dem(j,ic) + dv
            dem(j,ib) = dem(j,ib) - du - dv
            frcz(j) = frcz(j) + du
            frcx(j) = frcx(j) + dv
         end do
c
c     force distribution for the Bisector local coordinate method
c
      else if (axetyp .eq. 'Bisector') then
         do j = 1, 3
            du = uv(j)*dphidv/(usiz*uvsin) + 0.5d0*uw(j)*dphidw/usiz
            dv = -uv(j)*dphidu/(vsiz*uvsin) + 0.5d0*vw(j)*dphidw/vsiz
            dem(j,ia) = dem(j,ia) + du
            dem(j,ic) = dem(j,ic) + dv
            dem(j,ib) = dem(j,ib) - du - dv
            frcz(j) = frcz(j) + du
            frcx(j) = frcx(j) + dv
         end do
c
c     force distribution for the Z-Bisect local coordinate method
c
      else if (axetyp .eq. 'Z-Bisect') then
         do j = 1, 3
            du = ur(j)*dphidr/(usiz*ursin) + us(j)*dphids/usiz
            dv = (vssin*s(j)-vscos*t1(j))*dphidu
     &              / (vsiz*(ut1sin+ut2sin))
            dw = (wssin*s(j)-wscos*t2(j))*dphidu
     &              / (wsiz*(ut1sin+ut2sin))
            dem(j,ia) = dem(j,ia) + du
            dem(j,ic) = dem(j,ic) + dv
            dem(j,id) = dem(j,id) + dw
            dem(j,ib) = dem(j,ib) - du - dv - dw
            frcz(j) = frcz(j) + du
            frcx(j) = frcx(j) + dv
            frcy(j) = frcy(j) + dw
         end do
c
c     force distribution for the 3-Fold local coordinate method
c        (correct for uv, uw and vw angles all equal to 90)
c
      else if (axetyp .eq. '3-Fold') then
         do j = 1, 3
            du = uw(j)*dphidw/(usiz*uwsin)
     &              + uv(j)*dphidv/(usiz*uvsin)
     &              - uw(j)*dphidu/(usiz*uwsin)
     &              - uv(j)*dphidu/(usiz*uvsin)
            dv = vw(j)*dphidw/(vsiz*vwsin)
     &              - uv(j)*dphidu/(vsiz*uvsin)
     &              - vw(j)*dphidv/(vsiz*vwsin)
     &              + uv(j)*dphidv/(vsiz*uvsin)
            dw = -uw(j)*dphidu/(wsiz*uwsin)
     &              - vw(j)*dphidv/(wsiz*vwsin)
     &              + uw(j)*dphidw/(wsiz*uwsin)
     &              + vw(j)*dphidw/(wsiz*vwsin)
            du = du / 3.0d0
            dv = dv / 3.0d0
            dw = dw / 3.0d0
            dem(j,ia) = dem(j,ia) + du
            dem(j,ic) = dem(j,ic) + dv
            dem(j,id) = dem(j,id) + dw
            dem(j,ib) = dem(j,ib) - du - dv - dw
            frcz(j) = frcz(j) + du
            frcx(j) = frcx(j) + dv
            frcy(j) = frcy(j) + dw
         end do
      end if
      return
      end
