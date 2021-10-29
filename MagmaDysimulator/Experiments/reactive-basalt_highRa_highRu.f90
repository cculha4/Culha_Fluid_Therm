! Navier-Stokes solution
!
! density difference
! interfacial tension
!
! symmetry or no slip bc at y=0,1
! inlet velocity x=0
! outlet Pressure x=1
!
! needs curvture smoothing
!
! bubbles,planer interface,solids
!
!-----------------------------------------------------------------------------

     program FVNS
        implicit real*8(a-h,o-z)
        common/param/g,sp,ubc,uout,mint
        common/grid/dx,dy,xmx,xmn,ymx,ymn,nx,ny,tt,tb,crmx,ns,nsp,nsl
        allocatable :: x(:),y(:),uin(:),u(:,:),v(:,:),p(:,:), &
                       rho1x(:,:),Tem(:,:),rho(:,:),&
                       rho1y(:,:),rho2x(:,:),rho2y(:,:),intf(:,:),&
                       isize(:),control(:),xout(:),xout0(:),xmu(:,:),&
                       rpx(:,:),rpy(:,:),x1(:),Tem1(:,:),phisolid(:,:),&
                       phiw(:,:,:),xs(:,:),ys(:,:),vf(:,:),phigas(:,:),phi_solid(:,:,:),&
                       xo(:),yo(:),ubr(:),vbr(:),wbr(:),nmcluster(:,:),nsize(:),&
                       nxbound(:),nybound(:),nwall_l(:),nwall_r(:),rho1xs(:,:),&
                       rho2xs(:,:),rho1ys(:,:),rho2ys(:,:),ixmn(:),ixmx(:),iymn(:),iymx(:),&
                       w(:,:),sig(:,:,:),u_TVD(:,:),v_TVD(:,:),u_pre(:,:),&
                       rhof(:,:),xmuf(:,:),rhog(:,:),xmug(:,:),phi_eq(:,:),&
                       xo_new(:),yo_new(:),phi_eq_before(:,:),p_vol(:),&
                       p_rho_f(:),p_mu_f(:),p_rho_m(:),p_mu_m(:),itracker(:),&
                       x_tr(:),y_tr(:),itracker_tr(:),u_tr(:),v_tr(:),x_loc1(:),y_loc1(:),&
                       x_tr2(:),y_tr2(:),itracker_tr2(:),u_tr2(:),v_tr2(:),x_loc2(:),y_loc2(:),&
                       x_tr3(:),y_tr3(:),itracker_tr3(:),u_tr3(:),v_tr3(:),x_loc3(:),y_loc3(:),&
                       x_tr4(:),y_tr4(:),itracker_tr4(:),u_tr4(:),v_tr4(:),x_loc4(:),y_loc4(:),&
                       x_tr5(:),y_tr5(:),itracker_tr5(:),u_tr5(:),v_tr5(:),x_loc5(:),y_loc5(:),&
                       cnctr(:,:),xlct(:,:),ylct(:,:),u_lct(:,:),v_lct(:,:)


        integer, dimension(:,:), allocatable :: ns_needed_before, ns_needed, ns_in_box
        character (len=3) num,num0,nss,nsc
        integer*8 inumeric,isymbolic,istatus,isys

        call cpu_time(stime)
        !---------------------------------------------------------
        nx = 210 + 1
        ny = 210 + 1
        nsp = 1 !percent of crystals at bottom boundary
        nsl = 800!16*nsp/5  !percent of crystal at the top boundary
        ns = 0 !actual number of crystals from set up
        rdx = 0.003!0.002
        !-----------------------------
        ninlet = ny                  ! inlet grid
        nwrite = 200
        !-----------------------------
        g = -9.8d0                 ! gravitaional constant
        !-----------------------------
        ifelsic = 0                  !magma type
        if (ifelsic.eq.0) then
             rhos = 3000.d0                ! density of crystal
             capg = 1367.                ! heat capibility of the initially lower flow
             capf = 1367.                ! heat capibility of the initially higher flow
             tkthg = 1.53               ! heat conductivity of the initially lower flow
             tkthf = 1.53                ! heat conductivity of the initially higher flow
             xit = 10.d0/100.d0          ! interfacial tension
            !-----------------------------
             tt = 800.                   !top temperature
             tb = 1040.            !bottom temperature
            !--------------------------------------------------------
             cfl1 = 1.d0/5.d0          ! cfl number for advection
             cfl2 = 1.d0/1000.d0
            !--------------------------------------------------------
             tst = 0.0                  ! starting time
             tmax = 1000.                 ! max simulation time
         elseif (ifelsic.eq.1) then
            rhos = 2600.d0                ! density of crystal
            capg = 1367.                ! heat capibility of the initially lower flow
            capf = 1367.                ! heat capibility of the initially higher flow
            tkthg = 1.53               ! heat conductivity of the initially lower flow
            tkthf = 1.53                ! heat conductivity of the initially higher flow
            xit = 10.d0/100.d0          ! interfacial tension
            !-----------------------------
            tt = 700.                   !top temperature
            tb = 940.            !bottom temperature
            !-----------------------------
            cfl1 = 4.d0/5.d0          ! cfl number for advection
            cfl2 = 1.d0/1000.d0
            tst = 0.0                  ! starting time
            tmax = 5000.                 ! max simulation time
        end if
        !--------------------------------------------------------
        ubc =  1.d0                ! free slip ubc=1, no slip ubc=-1
        uinn = 0.d0                ! inlet velocity
        uout = 0.d0                ! outlet velocity(i=nx), for P=0, uout=1
        !--------------------------------------------------------
        iout_max = 1E9             ! maximum iterations
        tol = 1E-8                 ! steady-state convergence
        isw = 1                    ! write to screen
        ivrd = 0                   ! read initial value data
        iwp_icr = 10               ! delta itr for param.dat output
        !--------------------------------------------------------
        xmx =   0.3                 ! maximum x coordinate value
        xmn =   0.                  ! min x
        ymx =   0.3                 ! max y
        ymn =   0.                  ! min y
        dx = (xmx-xmn)/dble(nx-1)   ! delta x
        dy = (ymx-ymn)/dble(ny-1)   ! delta y
        ar = dx/dy                  ! aspect ratio
        umx = 0.00001
        vmx = 0.00001
        iout = 1
        iout0 = 1
        sp = 1.5d0*dx               ! spread around interface
        nav = 20
        nxnav = int(dble(nx)/dble(nav))+1
        nynav = int(dble(ny)/dble(nav))+1
        !--------------------------------------------------------
        !levelset redistance controls
        itrmax_ls = 10  !maximum itrations
        tol_ls = 1E-10 !redistancing tolerance
        tmx_ls = 1.    !max ir time
        cfl_ls = 0.01  !redist cfl
        err_ls = 0.    !redist error
        itr_ls = 1     !itr counter
        mls_incr = 1  !delta itr for redist
        mls = 1       !delta itr for redist start
        !Turn Cansu's subroutines on or off
        ccon = 1.
        !lagrangian trackers
        nolgr = 20
        navxlct = nx*3
        navylct = ny*3


        write(nss,'(i3)') nsp+100
        write(nsc,'(i3)') nsl+100
        open(30,file="param_magma_chamber_ev_hghRa_"// nss //"_"// nsc //".dat")
        write(30,*) 'variables=t,urms,umx,vmx,pavg,Tem'

        allocate(x(nx),y(ny),uin(ny+1),u(nx,ny),v(nx,ny),Tem(nx,ny), &
                 p(nx,ny),rho(nx,ny), &
                 rho1x(nx,ny),rho1y(nx,ny),rho2x(nx,ny),rho2y(nx,ny), &
                 intf(nx,ny),xmu(nx,ny),control(20),&
                 rpx(nx,2),rpy(ny,2),xout(nwrite),x1(nx+1),Tem1(nx+1,ny),&
                 xout0(31),phiw(nx,ny,3),xs(nx,3),ys(ny,3),phisolid(nx,ny),&
                 vf(nx,ny),phigas(nx,ny),phi_solid(nx,ny,3),&
                 xo(10000),yo(10000),ubr(10000),vbr(10000),wbr(10000),&
                 nmcluster(10000,10000),nsize(10000),&
                 nxbound(10000),nybound(10000),nwall_l(10000),nwall_r(10000),&
                 rho1xs(nx,ny),rho1ys(nx,ny),rho2xs(nx,ny),rho2ys(nx,ny),&
                 ixmn(10000),ixmx(10000),iymn(10000),iymx(10000),w(nx,ny),sig(nx,ny,2),&
                 u_TVD(nx,ny),v_TVD(nx,ny),u_pre(nx,ny),rhof(nx,ny),xmuf(nx,ny),&
                rhog(nx,ny),xmug(nx,ny),phi_eq(nx,ny),xo_new(10000),yo_new(10000),&
                 phi_eq_before(nx,ny),p_vol(4),p_rho_f(2),&
                 p_mu_f(2),p_rho_m(3),p_mu_m(2),itracker(10000),x_tr(nolgr),y_tr(nolgr),&
                 itracker_tr(nolgr),u_tr(nolgr),v_tr(nolgr),x_loc1(2),y_loc1(2),&
                 x_tr2(nolgr),y_tr2(nolgr),itracker_tr2(nolgr),u_tr2(nolgr),v_tr2(nolgr),&
                 x_loc2(2),y_loc2(2),x_tr3(nolgr),y_tr3(nolgr),itracker_tr3(nolgr),u_tr3(nolgr),&
                 v_tr3(nolgr),x_loc3(2),y_loc3(2),x_tr4(nolgr),y_tr4(nolgr),itracker_tr4(nolgr),&
                 u_tr4(nolgr),v_tr4(nolgr),x_loc4(2),y_loc4(2),x_tr5(nolgr),y_tr5(nolgr),&
                 itracker_tr5(nolgr),u_tr5(nolgr),v_tr5(nolgr),x_loc5(2),y_loc5(2),cnctr(navxlct,navylct),&
                 xlct(navxlct,navylct),ylct(navxlct,navylct),u_lct(navxlct,navylct),v_lct(navxlct,navylct),&
                 ns_needed_before(nxnav,nynav),ns_needed(nxnav,nynav),ns_in_box(nxnav,nynav))



!
!!--------------------------------------------------------
!! rheology constants
!!!!!Lassen
!       if (ccon.eq.1.) then
!          p_rho_f = 0.
!          p_rho_m = 0.
!          p_mu_f = 0.
!          p_mu_f = 0.
!          p_rho_f(1) = 6.486720266388353e-01
!          p_rho_f(2) = 1.647388263558799e3
!          p_rho_m(1) = 8.436758856006292e-01
!          p_rho_m(2) = 1.460876276533728e3
!          p_mu_f(1) = -0.004703370118443
!          p_mu_f(2) = 7.824398014164816
!          p_mu_m(1) = -0.007596738709235
!          p_mu_m(2) = 9.606732807516890
!         p_vol = 0.
!         p_vol(1) = 2.118278271446838e-8
!         p_vol(2) = -6.528044709758813e-05
!         p_vol(3) = 0.063483689955344
!         p_vol(4) = -19.240081905719467
!      end if
!!!!Unzen
       if (ifelsic.eq.0) then
          if (ccon.eq.1.) then
             p_rho_f = 0.
             p_rho_m = 0.
             p_mu_f = 0.
             p_mu_f = 0.

             p_rho_f(1) = 6.486720266388353e-01
             p_rho_f(2) = 1.647388263558799e3
             p_rho_m(1) = 0.000005361328234
             p_rho_m(2) = -0.009120372457138
             p_rho_m(3) = 6.083106105210905
             p_mu_f(1) = -0.004703370118443
             p_mu_f(2) = 7.824398014164816
             p_mu_m(1) = -0.0083
             p_mu_m(2) = 9.8978
             p_vol = 0.
             p_vol(1) = 1.195349843624094e-08
             p_vol(2) = -4.577245303322344e-05
             p_vol(3) = 0.051005252900785
             p_vol(4) = -16.954103757044030
          end if
      elseif (ifelsic.eq.1) then
         if (ccon.eq.1.) then
            p_rho_f = 0.
            p_rho_m = 0.
            p_mu_f = 0.
            p_mu_m = 0.

            p_rho_f(1) = 0.000003965780087
            p_rho_f(2) = 1.647388263558799e3
            p_rho_m(1) = 0.000003965780087
            p_rho_m(2) = -0.006181681951894
            p_rho_m(3) = 4.580188157467195
            p_mu_f(1) = -0.004703370118443
            p_mu_f(2) = 7.824398014164816
            p_mu_m(1) = -0.006880209123201
            p_mu_m(2) = 9.754526113628000
            p_vol = 0.
            p_vol(1) = -0.000000097564990
            p_vol(2) = 0.000246154051750
            p_vol(3) = -0.208262866735695
            p_vol(4) = 59.365280578985910
         end if
      end if
           x_tr = 0.
           y_tr = 0.
           u_tr = 0.
           v_tr = 0.
           itracker_tr = 0.





!----------initial location and velocity of crystals------------
        ixcntr = 0
        itracker(1) = ixcntr

        axmx = 0.1                         ! maximum x coordinate value
        axmn = 0.08                       ! min x
        aymx = 0.1                         !max y
        aymn = 0.0                         ! min y
!
        pii = 4.*atan(1.d0)
        vcrystal = pii*rdx**2.
        alayervol = (axmx-axmn)*(aymx-aymn)
        ns = 0

        phigas = -100.

       ns_needed1 = 1
!       call smart_crystal_adding(rdx,rdx,ns_needed1,xo,yo,axmx,&
!            axmn,aymx,aymn,phigas,itracker,xcntr)
       call smart_crystal_adding(rdx,ns_needed1,xo,yo,axmx,&
                  axmn,aymx,aymn,phigas,itracker,ixcntr)

        print *,'Added Crystals'


        ubr = 0.
        vbr = 0.
        wbr = 0.
!--------------------------------------------------------------

        pi = 4.*atan(1.d0)
        t = tst
        umd = 1.
        vmd = 1.
        urd = 1.
        uin = uinn

!----- initialize -------------

        lm = 0  !main loop counter init
        ns_needed = 0

    call init(p,u,v,rho1x,rho1y,rho2x,rho2y,xmu,rpx,rpy,xout,xout0,x,y,uin,tmax,intf,Tem,lm,rho,&
              nwrite,vf,phigas,rhog,xmug,rhof,xmuf,xit,w,sig,rdx,&
              p_rho_f,p_mu_f,p_rho_m,p_mu_m,cnctr,phisolid,xo,yo,axmn,navxlct,navylct,xlct,&
              ylct,u_lct,v_lct)




!----------------main loop------------------------------------

        iwp = 1 !param.dat output counter init

        do while(t.le.tmax.and.lm.lt.iout_max)
           lm = lm + 1
           ngas = 0
         
           !------------------------------------------------------------------

           if (iout.eq.1) then
              write(num,'(i3)') iout+100
           end if
!------------Adding or subtracting crystals--------------------------
            if (ccon.eq.1.) then
            call determining_volume(Tem,phigas,phi_eq,ns_needed,iout,lm,iout_max,t,&
                                     tmax,xout,nwrite,nav,rdx,xo,yo,num,p_vol,ns_in_box,nxnav,nynav)

              do k = 1,int(dble(ny)/dble(nav))+1
                do l =1,int(dble(nx)/dble(nav))+1
                   phi_eq_before(l,k) = phi_eq(l,k)
                   ns_needed_before(l,k) = ns_needed(l,k)
                end do
              end do

              do j = 1,int(dble(ny)/dble(nav))+1
               do i = 1,int(dble(nx)/dble(nav))+1


                if (ns_needed(i,j).ne.0) then
                 write(*,*)'negative needed',j,i,ny/nav,nx/nav,ns_needed(i,j)

                 print *, ns_needed(i,j),xmx/nx,dx,dble(i)*dx*nav,dble(i-1)*dx*nav
                 axmn = dble(i-1)*dx*dble(nav)! + dx
                 axmx = dble(i)*dx*dble(nav)! - dx
                 aymn = dble(j-1)*dy*dble(nav)! +dy
                 aymx = dble(j)*dy*dble(nav)! - dy

!                 call smart_crystal_adding(rdx,rdx,ns_needed(i,j),xo,yo,axmx,&
!                      axmn,aymx,aymn,phigas,itracker,xcntr)
                  call smart_crystal_adding(rdx,ns_needed(i,j),xo,yo,axmx,&
                       axmn,aymx,aymn,phigas,itracker,ixcntr)

                 WRITE(*,*)'Input complete. Number of records: ',ns

                end if
               end do
              end do
           end if
!----------------------------------------------------------------------------
           do k = 1, ns

              !--------------------------------limit the grid for crystals on level 2------------------------------------
              ixmn(k) = floor(((xo(k)-rdx)-x(1))/dx)-1-1
              ixmx(k) = floor(((xo(k)+rdx)-x(1))/dx)+2+1
              iymn(k) = floor(((yo(k)-rdx)-y(1))/dx)-1-1
              iymx(k) = floor(((yo(k)+rdx)-y(1))/dx)+2+1

           end do

           rho1xs = rho1x
           rho2xs = rho2x
           rho1ys = rho1y
           rho2ys = rho2y
           if (ns.gt.0) then
              call solid_property_level(x,y,xo,yo,rdx,rho1xs,rho2xs,rho1ys,rho2ys,ixmn,ixmx,iymn,iymx,rhos)
           end if
           !----------------------solve Stoke equation by iteration method-----------------------
           TOL1 = 1E-3
           TOL2 = 1E-3
           k_Stoke = 0
           error = 1E5
           error1 = 1E0
           error2 = 1E5
           rhomn = 1.0E6
           xmumx  = 10000.

           do j = 1,ny
              do i = 1,nx
                 rhomn = min(rhomn,rhog(i,j))
                 rhomn = min(rhomn,rhof(i,j))
              end do
           end do
           do j = 1,ny
              do i = 1,nx
                 xmumx = max(xmumx,xmu(i,j))
                 xmumx = max(xmumx,xmu(i,j))
              end do
           end do
           C_cfl=umx/dx+vmx/dy
           V_cfl=0.
           G_cfl=sqrt(abs(g)/dx)
           S_cfl=sqrt(xit*crmx/rhomn/dx**2.)
           dt=cfl2/((C_cfl+V_cfl)+sqrt((C_cfl+V_cfl)**2.+4.*G_cfl**2.))!+4.*S_cfl**2.))
           do while (error.gt.TOL2)
                    u_pre = u
                    k_Stoke = k_Stoke + 1
                    !---------------------solve Stokes equation------------------------
                    call solve1(u,v,p,uin,rho1x,rho1y,rho2x,rho2y,xmu,rpx,rpy, &
                                dt,dx,dy,intf,control,div,err,res,errm,nx,ny,mu, &
                                inumeric,isymbolic,itr,t,x,y,rho1xs,rho1ys,&
                                rho2xs,rho2ys,w,sig,k_Stoke,phigas)
                    umx = 0.0000001
                    vmx = 0.0000001
                    do j = 1,ny
                       do i = 1,nx
                          umx = max(umx,abs(u(i,j)))  !max velocity
                          vmx = max(vmx,abs(v(i,j)))  !min velocity
                       end do
                    end do
                    C_cfl=umx/dx+vmx/dy
                    V_cfl=0.
                    G_cfl=sqrt(abs(g)/dx)
                    S_cfl=sqrt(xit*crmx/rhomn/dx**2.)
                    dt=cfl2/((C_cfl+V_cfl)+sqrt((C_cfl+V_cfl)**2.+4.*G_cfl**2.))!+4.*S_cfl**2.))
                    !-------------calculate the error at every iteration---------------
                    !error = sqrt(abs(sum(u**2-u_pre**2)/dble(nx*ny)))/umx
                    error = 0.
                    do j = 1,ny
                       do i = 1,nx
                          error = max(error,abs((u(i,j)-u_pre(i,j))/dt))
                       end do
                    end do
                    Re_test = max(umx,vmx)*rdx*rhomn/xmumx
                    error2 = abs(error-error1)/error1
                    error1 = error
                    write(*,*) k_Stoke, error, error1, error2, t, Re_test
           end do

           call umf4fnum(inumeric)
           call umf4fsym(isymbolic)

           dt=min(cfl1*dx/(umx+vmx), cfl1*dx**2./(umx*dx+vmx*dy+4.*tkthg/(rhomn*capg)))
            if (ns.gt.0) then
           !------------------crystal computation---------------------
           call collid_find(xo,yo,ncluster,nsize,nmcluster,rdx,ns)
           !---------------------------------update the velocity for crystal in level1-----------------------------------
           call solve_solid(u,v,dt,rdx,xo,yo,ubr,vbr,wbr,ixmn,ixmx,iymn,iymx,&
                            nmcluster,nsize,ncluster,ninlet,phi_solid)
!            call determining_volume(Tem,phigas,phi_eq,ns_needed,iout,lm,iout_max,t,tmax,xout,nwrite,nav,rdx,&
!                                    xo,yo,1,p_vol,ns_in_box)

            call solid_move(xo,yo,dt,ubr,vbr,ns)
           end if
           !----------------------------------------------------------

           ! solve the energy equation
           call energy(phigas,u,v,uin,Tem,capg,capf,tkthg,tkthf,rhof,xmuf,rhog,xmug,dt,ninlet,&
                      lm,x,y,intf,&
                      p_rho_f,p_mu_f,p_rho_m,p_mu_m)

!           call advect_WENO31(phigas,u,v,x,y,cfl,dt,phimx,dx,dy,xmn,ymn,nx,ny,uin,phigas(2,:),intf)

           !------------------------------------------------------------------------------

           call curvature(x,y,rho1x,rho2x,rho1y,rho2y,xmu,u,v,uin,intf,Tem,dt,lm,rho,&
                          vf,phigas,rhog,xmug,rhof,xmuf,xit,w,sig)

           !------------------------------------------------------------------------------
           t = t + dt
           
           xmerr = abs(xim-xcm)/xim !mass error

           if(iwp.eq.lm) then
              write(30,10) t,umx,vmx,amass_L,amass_L_int,amass_L-amass_L_int
              iwp = iwp + iwp_icr
           end if

           if(isw.eq.1) then
              write(6,'(i8,i5,3f10.5,3e10.2,2f10.5)') lm,isw,t,dt,umx,amass_L,amass_L_int,amass_L-amass_L_int
              write(6,'(i8,i5,3f10.5,3e10.2,2f10.5)') lm,isw,t,dt,umx,vmx
           end if

           !2D object
           phisolid = 1E9
           do j = 1,ny
              do i = 1,nx
                 do k = 1,ns
                    !bubble
                    phisolid(i,j) =  min(phisolid(i,j), - rdx + sqrt((x(i)-xo(k))**2. + (y(j)-yo(k))**2.))
                 end do
              end do
           end do


!          call concentration_tracking(navxlct,navylct,xlct,ylct,u,v,u_lct,v_lct,dt)

           if (lm.eq.iout_max.or.t.ge.tmax.or.t.ge.xout(iout)) then

              write(num,'(i3)') iout+100
              iout = iout + 1



              open(10,file="data_magma_chamber_trc_ev_hghRa_" // nss // "_" // nsc //"_" // num // ".dat")
              write(10,*) 'TITLE = "velocity field"'
              write(10,*) 'VARIABLES = x,y,p,u,v,rho,xmu,Tem,solid,gas,reaction,cnctr'
              write(10,*) "ZONE T = ",'"Rectanguler zone"'
              write(10,*) "time=",t,"I=",nx,",J=",ny,",F=POINT"
              write(10,*) "DT=(SINGLE SINGLE SINGLE SINGLE)"
              do j=1,ny
                 do i=1,nx
                    write(10,*) x(i),y(j),p(i,j),u(i,j),v(i,j), &
                                rho(i,j),xmu(i,j),Tem(i,j),&
                                phisolid(i,j),phigas(i,j)!,phivf(i,j,1),phivf(i,j,2)
                 end do
              end do
              close(10)
!              open(10,file="data_magma_chamber_phi_trc_" // nss // "_" // nsc //"_" // num // ".dat")
!               do k = 1,int(ny/nav)+1
!                   do l =1,int(nx/nav)+1
!                      write(10,*) phi_eq_before(l,k),ns_needed_before(l,k),&
!                                  phi_eq(l,k),ns_needed(l,k)
!                   end do
!               end do
!              close(10)
              open(10,file="data_magma_chamber_xtl_trc_ev_hghRa_" // nss // "_" // nsc //"_" // num // ".dat")
               write(10,*) 'TITLE = "velocity field"'
               write(10,*) 'VARIABLES = x,y,u,v,itracker'
               write(10,*) "ZONE T = ",'"Rectanguler zone"'
               write(10,*) "ns=",ns,",J=",ny,",F=POINT"
               write(10,*) "DT=(SINGLE SINGLE SINGLE SINGLE)"
               do k = 1,ns
                      write(10,*) xo(k),yo(k),ubr(k),vbr(k),itracker(k)
               end do
              close(10)

                  open(10,file="data_magma_chamber_gctrc_ev_hghRa_" // nss // "_" // nsc //"_" // num // ".dat")
                   write(10,*) 'TITLE = "Lagragian Grid Tracers"'
                   write(10,*) 'VARIABLES = xlct,ylct,u_lct,v_lct,cnctr'
                   write(10,*) "ZONE T = ",'"Rectanguler zone"'
                   write(10,*) "nx=",navxlct,",ny=",navylct,",F=POINT"
                   write(10,*) "DT=(SINGLE SINGLE SINGLE SINGLE SINGLE)"
                   do j = 1,navylct
                      do i = 1,navxlct
                          write(10,*) xlct(i,j),ylct(i,j),u_lct(i,j),v_lct(i,j),cnctr(i,j)
                      end do
                   end do
                   close(10)
!               if (iout.le.2) then
!                  open(10,file="data_magma_chamber_xtltrc_ev_" // nss // "_" // nsc //"_" // num // ".dat")
!                   write(10,*) 'TITLE = "Lagragian Tracers Around Crystals"'
!                   write(10,*) 'VARIABLES = ns,navxlct,navylct,xo,yo,xlct,ylct,u_lct,v_lct,cnctr'
!                   write(10,*) "ZONE T = ",'"Rectanguler zone"'
!                   write(10,*) "nx=",navxlct,",ny=",navylct,",F=POINT"
!                   write(10,*) "DT=(SINGLE SINGLE SINGLE SINGLE SINGLE)"
!                   do j = 1,navylct
!                      do i = 1,navxlct
!                       do l = 1,ns
!                          if (sqrt((xlct(i,j)-xo(l))**2.+(ylct(i,j)-yo(l))**2.).le.6.*rdx) then
!                            write(10,*) l,i,j,xo(l),yo(l),xlct(i,j),ylct(i,j),u_lct(i,j),v_lct(i,j),cnctr(i,j)
!                          end if
!                        end do
!                      end do
!                   end do
!                   close(10)
!              end if


           end if

        end do

        call cpu_time(etime)

        if(isw.eq.1) then
           write(6,'(a20,f20.10,i5)') "elapsed time in mins = ",(etime-stime)/60.
        end if

10      format(12e25.15)
20      format(i6,2f10.5,4e15.5)


      end program FVNS

!---------------------------------------------------------------------------------
!*********************************************************************************
!
!*********************************************************************************
!---------------------------------------------------------------------------------




           subroutine concentration_tracking(navxlct,navylct,xlct,ylct,u,v,u_lct,v_lct,dt)
             implicit real*8(a-h,o-z)
             common/grid/dx,dy,xmx,xmn,ymx,ymn,nx,ny,tt,tb,crmx,ns,nsp,nsl
             common/param/g,sp,ubc,uout,mint
             dimension :: xlct(navxlct,navylct),ylct(navxlct,navylct),&
                          u_lct(navxlct,navylct),v_lct(navxlct,navylct),&
                          u(nx,ny),v(nx,ny)
             allocatable :: x_tr(:),y_tr(:),u_tr(:),v_tr(:)
             allocate(x_tr(navxlct*navylct),y_tr(navxlct*navylct),&
                      u_tr(navxlct*navylct),v_tr(navxlct*navylct))
           l = 0
           do j = 1,navylct
              do i = 1,navxlct
                 l = l+1
                 x_tr(l) = xlct(i,j)
                 y_tr(l) = ylct(i,j)
                 u_tr(l) = u_lct(i,j)
                 v_tr(l) = v_lct(i,j)
              end do
           end do



           call tracermover(x_tr,y_tr,u_tr,v_tr,navxlct*navylct,u,v,dt)
           do j = 1,navylct
              do i = 1,navxlct
                 xlct(i,j) = x_tr((j-1)*navxlct+i)
                 ylct(i,j) = y_tr((j-1)*navxlct+i)
                 u_lct(i,j) = u_tr((j-1)*navxlct+i)
                 v_lct(i,j) = v_tr((j-1)*navxlct+i)
              end do
           end do
!           print *, 'Made it through concentration tracking'
       end subroutine concentration_tracking

!---------------------------------------------------------------------------------
!*********************************************************************************
!
!*********************************************************************************
!---------------------------------------------------------------------------------


     subroutine concentration_tracking_continuum(phisolid,u,v,axmx,axmn,aymx,aymn,cnctr,ninlet,intf,&
                                       uin,dt)
       implicit real*8(a-h,o-z)
       common/grid/dx,dy,xmx,xmn,ymx,ymn,nx,ny,tt,tb,crmx,ns,nsp,nsl
       common/param/g,sp,ubc,uout,mint
       dimension :: phisolid(nx,ny),u(nx,ny),v(nx,ny),cnctr(nx,ny),intf(nx,ny),uin(ny+1)
        allocatable :: cnctrT(:,:),tkth(:,:)
        allocate(cnctrT(nx,ny),tkth(nx,ny))

        pii = 4.*atan(1.d0)
!
!        !-------------------update density&thermal diffusivity------------------------

        cnctrT = cnctr
        tkth = 0.

        do j = 1,ny
           do i = 1,nx
              if (intf(i,j).eq.0) then
                 if (i.eq.1) then
                       sxl = cnctrT(i,j)
                       ssxl= tkth(i,j)*(cnctrT(i,j)-cnctrT(i,j))
                       ul = uin(j)
                 else
                       sxl = cnctrT(i-1,j)
                       ssxl= (tkth(i,j)+tkth(i-1,j))/2.*(cnctrT(i,j)-cnctrT(i-1,j))
                       ul = u(i-1,j)
                 end if
                 if(i.eq.nx) then
                    sxr = cnctrT(i,j)
                    ssxr= tkth(i,j)*((cnctrT(i,j))-cnctrT(i,j))
                 else
                       sxr = cnctrT(i+1,j)
                       ssxr= (tkth(i+1,j)+tkth(i,j))/2.*(cnctrT(i+1,j)-cnctrT(i,j))
                 end if
                 ur = u(i,j)

                 if(j.eq.1) then
                    syl = cnctrT(i,j)
                    ssyl= tkth(i,j)*(cnctrT(i,j)-cnctrT(i,j))
                    vl = 0.
                 else
                       syl = cnctrT(i,j-1)
                       ssyl= (tkth(i,j)+tkth(i,j-1))/2.*(cnctrT(i,j)-cnctrT(i,j-1))
                       vl = v(i,j-1)
                 end if

                 if(j.eq.ny) then
                    syr = cnctrT(i,j)
                    ssyr= tkth(i,j)*(cnctrT(i,j)-cnctrT(i,j))
                 else
                       syr = cnctrT(i,j+1)
                       ssyr= (tkth(i,j)+tkth(i,j+1))/2.*(cnctrT(i,j+1)-cnctrT(i,j))
                 end if
                 vr = v(i,j)

                 sm = cnctrT(i,j)



                 if((ur+ul).ge.0.) then
                    fxr = sm*(ur+ul)/2.
                 else
                    fxr = sxr*(ur+ul)/2.
                 end if

                 if((ur+ul).ge.0.) then
                    fxl = sxl*(ur+ul)/2.
                 else
                    fxl = sm*(ur+ul)/2.
                 end if


                 if((vr+vl).ge.0.) then
                    fyr = sm*(vr+vl)/2.
                 else
                    fyr = syr*(vr+vl)/2.
                 end if

                 if((vr+vl).ge.0.) then
                    fyl = syl*(vr+vl)/2.
                 else
                    fyl = sm*(vr+vl)/2.
                 end if
                 cnctr(i,j) = cnctrT(i,j) - dt * ( (fxr-fxl)/dx + (fyr-fyl)/dy )
!                 write(10,*) ur, ul, sm,sxl,sxr,syr,syl, &
!                             fyl,dt * ( (fxr-fxl)/dx + (fyr-fyl)/dy ),dt

             end if
           end do
        end do
!       close(10)
!        do j = 1,ny
!          do i = 1,nx
!             if (phisolid(i,j).ge.0.) then
!                cnctr(i,j) = 1.0
!              end if
!          end do
!        end do

     end subroutine concentration_tracking_continuum
!---------------------------------------------------------------------------------
!*********************************************************************************
!
!*********************************************************************************
!---------------------------------------------------------------------------------
!Lagrangian Trackers
     subroutine lgntrc(xo,yo,rdx,u,v,x_loc,y_loc,ime_dp,t,dt,xout,nwrite,&
                       x_tr,y_tr,u_tr,v_tr,itracker_tr,nolgr)
        !xo and yo:  crystal locations (in-real)
        !rdx: crystal size (in-real)
        !u and v: liquid velocities (in-real)
        !ime_dp: when the tracers are added (in-integer)
        !t: the current time (in-integer)
        !dt: time step
        !xout: times of output (in-real)
        !nwrite: number of outputs (in-integer)
        !x_tr and y_tr: lagrangian tracer locations (out-real)
        !x_loc and y_loc: initiation location of the tracers (in-real)
        !itracker_tr: lagrangian tracer counter (out-integer)
        !u_tr and v_tr: velocities of the trackers (out-real)
        implicit real*8(a-h,o-z)
        common/grid/dx,dy,xmx,xmn,ymx,ymn,nx,ny,tt,tb,crmx,ns,nsp,nsl
        common/param/g,sp,ubc,uout,mint
        dimension :: xo(ns),yo(ns),u(nx,ny),v(nx,ny),x_loc(2),y_loc(2),xout(nwrite),&
                     x_tr(nolgr),y_tr(nolgr),itracker_tr(nolgr),u_tr(nolgr),v_tr(nolgr)
!        allocatable :: x_tr(:),y_tr(:),itracker_tr(:)
!       allocate(x_tr(10),y_tr(10),itracker_tr(10))


       if (t.ge.xout(ime_dp)) then
!          print *, t,xout(ime_dp),itracker_tr
          if (itracker_tr(1).eq.0.) then
          axmn = x_loc(1)
          axmx = x_loc(2)
          aymn = y_loc(1)
          aymx = y_loc(2)
          rlagrdx = dx
            call smart_tracker_adding(rdx,rlagrdx,nolgr,xo,yo,x_tr,y_tr,&
                                      axmx,axmn,aymx,aymn,itracker_tr)

            u_tr = 0.
            v_tr = 0.
            !print and save the initial x_tr and y_tr
            open(10,file="Lag_Tracker_Init.dat")
               write(10,*) 'TITLE = "Lagragian Tracers Initial"'
               write(10,*) 'VARIABLES = x_trc,y_trc,itracker_trc'

               do k = 1,nolgr
                  write(10,*) x_tr(k),y_tr(k),u_tr(k),v_tr(k),itracker_tr(k)
               end do
            close(10)
           end if
       end if


       if (t.ge.xout(ime_dp)) then


         call tracermover(x_tr,y_tr,u_tr,v_tr,nolgr,u,v,dt)
       end if

!      print *, 'Made it through lgntrc'


     end subroutine lgntrc

     subroutine tracermover(x_tr,y_tr,u_tr,v_tr,nolgr,u,v,dt)
         implicit real*8(a-h,o-z)



          common/grid/dx,dy,xmx,xmn,ymx,ymn,nx,ny,tt,tb,crmx,ns,nsp,nsl
          common/param/g,sp,ubc,uout,mint
          dimension :: x_tr(nolgr),y_tr(nolgr),u(nx,ny),v(nx,ny),&
                       u_tr(nolgr),v_tr(nolgr)
!          allocatable ::
!          allocate()

           do k = 1,nolgr
                !identifying the grid points closest to the tracker
                ixm = floor(x_tr(k)/dx)
                ixn = ceiling(x_tr(k)/dx)
                iym = floor(y_tr(k)/dy)
                iyn = ceiling(y_tr(k)/dy)

                if (ixm.eq.ixn) then
                     if (ixm.le.1) then
                       ixn = ixn + 1
                     elseif (ixm.ge.nx) then
                       ixm = nx - 1
                     end if
                end if
                if (iym.eq.iyn) then
                     if (iym.le.1) then
                       iyn = iyn + 1
                     elseif (iym.ge.ny) then
                       iym = ny - 1
                     end if
                end if

!                print *, ixm,ixn,iym,iyn

!                if (ixn.lt.1) then
!                   ixn = 1
!                end
!              if (ixn.gt.nx) then
!                  ixn = nx
!              end if
!              if (iyn.lt.1) then
!                  iyn = 1
!              end if
!              if (iyn.lt.ny) then
!                  iyn = ny
!              end if
                !identify how close the grids are to the tracker
                axf = abs(x_tr(k) - dx*ixm)/dx
                bxc = abs(x_tr(k) - dx*ixn)/dx
                ayf = abs(y_tr(k) - dy*iym)/dy
                byc = abs(y_tr(k) - dy*iyn)/dy
                if (axf.eq.0.) then
                   axf = 0.0000000001
                end if
                if (bxc.eq.0.) then
                   bxc = 0.0000000001
                end if
                if (ayf.eq.0.) then
                   ayf = 0.0000000001
                end if
                if (byc.eq.0.) then
                   byc = 0.0000000001
                end if

!                print *, x_tr(k), dx*ixm, dx*ixn, y_tr(k), dy*iym, dy*iyn
!                print *, axf, bxc, ayf, byc
                ii = int(floor( min(axf,bxc)/axf) * dble(ixm) + &
                        floor( min(axf,bxc)/bxc) * dble(ixn))
!                print *, axf,bxc, ayf,byc
                jj = int(floor( min(ayf,byc)/ayf) * dble(iym) + &
                        floor( min(ayf,byc)/byc) * dble(iyn))
!                print *, 'i,j', ii, jj

                !identify the velocity at the tracker
                u_tr(k) = axf*u(ixm,jj) + bxc*u(ixn,jj)
                v_tr(k) = ayf*v(ii,iym) + byc*v(ii,iyn)
!                print *,u_tr(k), u(ii,jj), v_tr(k), v(ii,jj)
                !move the tracer
                x_tr(k) = x_tr(k) + dt*u_tr(k)
                y_tr(k) = y_tr(k) + dt*v_tr(k)

                if (x_tr(k).gt.xmx) then
                  x_tr(k) = xmx
                end if
                if (y_tr(k).gt.ymx) then
                  y_tr(k) = ymx
                end if
                if (x_tr(k).lt.xmn) then
                  x_tr(k) = xmn
                end if
                if (y_tr(k).lt.ymn) then
                  y_tr(k) = ymn
                end if
           end do
!          print *, 'Made it through tracermover'
      end subroutine tracermover


!find the random value
       subroutine true_random(r,isize,ilength)
         implicit real*8(a-h,o-z)
         integer :: values(1:8), k
         integer, dimension(:), allocatable :: seed
         dimension :: r(isize,ilength)

         call date_and_time(values=values)

         call random_seed(size=k)
         allocate(seed(1:k))
         seed(:) = values(8)
         call random_seed(put=seed)
         call random_number(r)

      end subroutine

      subroutine determining_volume(Tem,sliq,phi_eq,ns_needed,iout,lm,iout_max,&
                                    t,tmax,xout,nwrite,nav,rdx,xo,yo,mboa,p,ns_in_box,nxnav,nynav)
!          !!!!It takes temperature (Tem), the levelset value (sliq), output conditions (iout,lm,iout_max,t,tmax,xout,nwrite), it takes how many cells in each direction it should be averaged (nav), and the crystal radius (rdx)
!          !!!!!! Out put is the volume of crystals present in the cell (phi_eq) and the number of crystals needed to get to that equilibrium (ns_needed)

          implicit real*8(a-h,o-z)
          integer :: n



          common/grid/dx,dy,xmx,xmn,ymx,ymn,nx,ny,tt,tb,crmx,ns,nsp,nsl
          common/param/g,sp,ubc,uout,mint
          dimension :: Tem(nx,ny),sliq(nx,ny),phi_eq(nx,ny),xout(nwrite),&
                       xo(10000),yo(10000),p(4)
          integer, dimension(nxnav,nynav) :: ns_needed,ns_in_box
          allocatable :: phivol(:,:),phivol10(:,:),ks(:),nxs(:),nys(:),phivf1(:),&
                         phivol10av(:,:),Tem10(:,:),Tem10av(:,:),phi_T(:,:),&
                         sliq10(:,:),sliq10av(:,:),x(:,:,:),y(:,:,:),s_part(:,:),&
                         ixmn(:),ixmx(:),iymn(:),iymx(:),n_s(:),phivf(:,:,:)
          character (len=3) num,num0,nss,nsc,mboa
          character (len=1) mum
!          ! get as many random number as you want by change 15000 below
          allocate(phivol(nx,ny),phivol10(nx,ny),phivol10av(nx,ny),Tem10(nx,ny),&
                  Tem10av(nx,ny),ks(nx*ny),nxs(nx*ny),nys(nx*ny),phivf(nx,ny,3),&
                  phivf1(nx*ny),&
                  phi_T(nx,ny),sliq10av(nx,ny),sliq10(nx,ny),x(nx,3,3),y(ny,3,3),s_part(nx,ny),&
                  ixmn(10000),ixmx(10000),iymn(10000),iymx(10000),n_s(ns))

!
!
!
!
!         !!!!!!!!!!!!!Implement the look up table
!
!         p(1) = 1.488733e-21!92390307e-21
!         p(2) = -1.35073236010659e-17
!         p(3) = 5.509921943400374e-14
!         p(4) = -1.330736688175155e-10
!         p(5) = 2.107295473308341e-07
!         p(6) = -0.0002286242346953226
!         p(7) = 0.172099650405596
!         p(8) = -88.757970744807030
!         p(9) = 30014.67085608857
!         p(10) = -6009642.1910502655
!         p(11) = 5.410180995646929e8
!         p = 0.
!         p(1) = 1.195349843624094e-08
!         p(2) = -4.577245303322344e-05
!         p(3) = 0.051005252900785
!         p(4) = -16.954103757044030
!
!         open(10,file="p.dat")
!          do l = 1,4
!             write(10,*) p(l)
!          end do
!          close(10)
          n = 0
          phivol = 0.
          phivol10 = 0.
          phivf = 0.

        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        !!!!!!Defining crystal area volume!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        pii = 4.*atan(1.d0)
        Area = pii*rdx**2.
        Area1 = pii*rdx**4.

        do i = 1,nx
           x(i,1,1) = dx*float(i-1)+xmn                !coordinate for bubble on level 1
           x(i,2,1) = dx*float(i-1)+xmn                !coordinate for particle on level 1
           x(i,1,2) = dx*float(i-1)+xmn                !coordinate for volume fraction for u on level 1
           x(i,2,2) = dx*float(i-1)+xmn-dx/2.          !coordinate for volume fraction for v on level 1
           x(i,3,2) = dx*float(i-1)+xmn-dx/2.          !coordinate for volume fraction for w on level 1
           x(i,1,3) = dx*float(i-1)+xmn+dx/2.          !coordinate for u update on level 1
           x(i,2,3) = dx*float(i-1)+xmn                !coordinate for v update on level 1
           x(i,3,3) = dx*float(i-1)+xmn
        end do
        do j = 1,ny
           y(j,1,1) = dy*float(j-1)+ymn                !coordinate for bubble on level 1
           y(j,2,1) = dy*float(j-1)+ymn                !coordinate for parkicle on level 1
           y(j,1,2) = dy*float(j-1)+ymn-dy/2.          !coordinate for volume fraction for u on level 1
           y(j,2,2) = dy*float(j-1)+ymn                !coordinate for volume fraction for v on level 1
           y(j,3,2) = dy*float(j-1)+ymn-dy/2.          !coordinate for volume fraction for w on level 1
           y(j,1,3) = dy*float(j-1)+ymn                !coordinate for u update on level 1
           y(j,2,3) = dy*float(j-1)+ymn+dy/2.          !coordinate for v update on level 1
           y(j,3,3) = dy*float(j-1)+ymn
        end do


          do k = 1, ns

!--------------------------------limit the grid for crystals on level 2------------------------------------
             ixmn(k) = floor(((xo(k)-rdx)-xmn)/dx)-1-1
             ixmx(k) = floor(((xo(k)+rdx)-xmn)/dx)+2+1
             iymn(k) = floor(((yo(k)-rdx)-ymn)/dx)-1-1
             iymx(k) = floor(((yo(k)+rdx)-ymn)/dx)+2+1

         end do
         do k = 1,ns
            n_s(k) = k
         end do
        write(nss,'(i3)') nsp+100
!        open(10,file="data_phivf" // nss // ".dat")
        do k1 = 1,ns
           do k2 = 1,3
              do j = iymn(k1),iymx(k1)
                 do i = ixmn(k1),ixmx(k1)
                    s_part(i,j) = - rdx + sqrt((x(i,k2,2)-xo(n_s(k1)))**2. + (y(j,k2,2)-yo(n_s(k1)))**2.)
                    phivf(i,j,k2) = 0.
                 end do
              end do
              !-------------volume fraction for particle------------
              call volume_frac(s_part(:,:),phivf(:,:,k2),ixmn(k1),ixmx(k1),iymn(k1),iymx(k1),&
                               x(:,k2,2),y(:,k2,2),rdx,xo(n_s(k1)),yo(n_s(k1)),k1)
           end do
!       write(mum,'(i3)') mboa
!        open(10,file="data_phivf" // nss // "_" // mum // ".dat")
            do j=iymn(k1),iymx(k1)
               do i=ixmn(k1),ixmx(k1)
                ! write (10,*) k1,i,j,phivf(i,j,3)
                 phivol(i,j) = phivf(i,j,2)
               end do
             end do
         end do
 !        close(10)

           !!! Identify the averaged areas


!           open(10,file="data_phivf2"//mum//".dat")
!         open(10,file="data_phivf" // nss // "_" // mboa // ".dat")
           do k = 1,int(dble(ny)/dble(nav))+1
            do l =1,int(dble(nx)/dble(nav))+1
               phivol10(l,k) = 0.
               phivol10av(l,k) = 0.
               Tem10(l,k) = 0.
               Tem10av(l,k) = 0.
               sliq10(l,k) = 0.
               sliq10av(l,k) = 0.
               phi_eq(l,k) = 0.
               phi_T(l,k) = 0.
               ns_needed(l,k) = 0
               m = k*nav-nav+1
               n = l*nav-nav+1

               axmn = dble(l-1)*dx*dble(nav)
               axmx = dble(l)*dx*dble(nav)
               aymn = dble(k-1)*dy*dble(nav)
               aymx = dble(k)*dy*dble(nav)

               do j = 0,nav-1
                !avj = avj + 1
                do i = 0,nav-1
                  if ((n+i).le.nx.and.(m+j).le.ny) then
                    if (x(n+i,1,1).gt.xmn+6.*dx+rdx.and.x(n+i,1,1).le.xmx-6.*dx-rdx.and.&
                        y(m+j,1,1).gt.ymn+6.*dx+rdx.and.y(m+j,1,1).le.ymx-6.*dx-rdx) then
                       if (sliq(n+i,m+j).lt.0.) then !check that we are in lower liquid

                      !Determine the volume of solids
                      phivol10(l,k) = phivol10(l,k) + phivol(n+i,m+j)
                      !Determine the Temperature
                      Tem10(l,k) = Tem10(l,k) + Tem(n+i,m+j)
                      !Determine the average liquid
                      sliq10(l,k) = sliq10(l,k) + 1.
                     endif
                   endif
                  endif
!                  sliq10(l,k) = sliq10(l,k) + sliq(n+i,m+j)
                end do
               end do
               !Identify the average
!               if (sliq10(l,k).gt.nav) then
!                if (sliq10(l,k).gt.0.) then
                  phivol10av(l,k) = phivol10(l,k)/sliq10(l,k)!/(nav*nav)
                  Tem10av(l,k) = Tem10(l,k)/sliq10(l,k)!(nav*nav)
!                else
!                  phivol10av(l,k) = 0.
!                  Tem10av(l,k) = Tem10(l,k)
!                endif

                do i = 1,ns
                 if (xo(i).gt.(axmn).and.xo(i).lt.(axmx).and.yo(i).gt.(aymn).and.yo(i).lt.(aymx)) then
                     ns_in_box(l,k) = ns_in_box(l,k) + 1
                 endif
                end do
                if (sliq10(l,k).gt.0.) then
                  do ll = 1,4
                     phi_T(l,k) = phi_T(l,k) + p(ll)*Tem10av(l,k)**(4.-dble(ll))
                  end do
                else
                 phi_T(l,k) = 0.
                endif
!                  if (Tem10av(l,k).ge.1015.) then
!                      phi_T(l,k) = 0.
!                  elseif (Tem10av(l,k).lt.1015.) then
!                      phi_T(l,k) = dble(nsp)/100.
!                  endif
!---------------------Top liquid
!                     phi_T(l,k) = (950.-Tem10av(l,k))*0.5/(950.-800.)
!                     if (Tem10av(l,k).ge.950.) then
!                        phi_T(l,k) = 0.
!                     end if

                  phi_eq(l,k) = (phi_T(l,k) - phivol10av(l,k))


               !Identify how many crystals are needed
                  pii = 4.*atan(1.d0)
                  vcrystal = pii*rdx**2.
                  cellvol = sliq10(l,k)*dx*dy
                  cryst_ns =  phi_eq(l,k)*cellvol/vcrystal
                  ns_needed(l,k) = nint(cryst_ns)!sign(abs(cryst_ns),phi_eq(l,k)))
!                  write (10,*) phi_eq(l,k),cryst_ns,ns_needed(l,k)
            end do
          end do
!          close(10)

          !!!!!!!!Look up table for the temperature and if the volume matches


       end subroutine determining_volume





!---------------------------------------------------------------------------------
!*********************************************************************************
!
!*********************************************************************************
!---------------------------------------------------------------------------------
!
!
!      subroutine smart_crystal_adding(rdx1,rdx2,ns_needed,xo,yo,axmx,axmn,aymx,aymn,phigas,itracker,&
!                                      xcntr)
!     !!!It takes in the radius of the crystal (rdx1,rdx2-lagrangian tracker), the number of crystals needed or subtracted to stay at equilibrium (ns_needed), current crystal locations (xo,yo), the bounds of the cell we are filling (xmx,ymn,xmx,ymx), and the size of the domain the crystals are averaged over
!     !!!!Output is crystal locations (xo_new, yo_new)
!
!          implicit real*8(a-h,o-z)
!          common/grid/dx,dy,xmx,xmn,ymx,ymn,nx,ny,tt,tb,crmx,ns,nsp,nsl
!          common/param/g,sp,ubc,uout,mint
!!          integer, intent(in) :: ns_needed,nav
!!          real, intent(in) :: rdx,axmn,aymn,xo,yo
!!          real, intent(out) :: new_crystals
!          dimension :: xo(10000),yo(10000),phigas(nx,ny),itracker(10000)
!          allocatable :: x_level1(:),y_level1(:),r(:,:),r1(:,:),r2(:,:),r3(:,:),&
!                         new_crystals(:,:),xo_new(:),yo_new(:),inew_tracker(:),&
!                         itracker_new(:)
!          ! get as many random number as you want by change 15000 below
!          allocate(r(15000,2),r1(15000,2),r2(15000,2),r3(15000,2),xo_new(15000),&
!                  yo_new(15000),new_crystals(10000,2),inew_tracker(15000),itracker_new(15000))
!          ! creat the random number at first
!          CALL true_random(r,15000,2)
!
!             nsx = 0
!             li = 0.
!             pii = 4.*atan(1.d0)
!!
!!             print *, pii, axmx, aymx, axmn, aymn, axmn + 2.*rdx,aymn + 2.*rdx
!
!             if (ns_needed.gt.0) then
!
!                do i = 1,15000
!                  !!placing the crystals into the domain
!                  r1(i,1) = (axmx-axmn)*r(i,1)+axmn
!                  r1(i,2) = (aymx-aymn)*r(i,2) + aymn
!
!                end do
!
!
!                ! delete the number overlaping with each other
!                nhave = 1
!
!               r2(nhave,:) = r1(nhave,:)
!               do i = 2,15000
!               !! make sure the crystals are not on top of each other (was 2.5)
!                  ncheck = 0
!                  do j = 1,nhave
!                    if (sqrt((r1(i,1)-r2(j,1))**2.+(r1(i,2)-r2(j,2))**2.).le.2.5*rdx2) then !!originally 3
!                       ncheck = 1
!
!                       exit
!                    end if
!                  end do
!                  if (ncheck.eq.0) then
!                    nhave = nhave + 1
!                    r2(nhave,:) = r1(i,:)
!                  end if
!               end do
!               !check if any are overlapping with current crysals within 3rd (was 2.5)
!               nhave1 = 0
!               do i = 1,nhave
!                 ncheck = 0
!                 do j = 1,ns
!                   if (sqrt((r2(i,1)-xo(j))**2.+(r2(i,2)-yo(j))**2.).le.2.5*rdx1) then
!                      ncheck = 1
!                      exit
!                   end if
!                 end do
!                !check if the crystals are in mafic magma, and delete ones in felsic magma
!                 do j = 1,ny
!                    do k = 1,nx
!!                       if (dble(k)*dx.ge.(r2(i,1)-rdx)).and.dble(k)*dx
!                     if ((dble(k)*dx).ge.(r2(i,1)-3.*rdx1).and.(dble(k)*dx).le.(r2(i,1)+3.*rdx1).and.&
!                         (dble(j)*dy).ge.(r2(i,2)-3.*rdx1).and.(dble(j)*dy).le.(r2(i,2)+3.*rdx1)) then
!                         if (phigas(k,j).ge.-0.) then
!                            ncheck = 1
!                         end if
!                       end if
!                     end do
!                 end do
!                 !check if the crystals are by the walls
!                 if ((r2(i,1).le.(xmn+2.*rdx1+1.*dx)).or.(r2(i,1).ge.(xmx-2.*rdx1-1.*dx)).or.&
!                    (r2(i,2).le.(ymn+2.*rdx1+1.*dy)).or.(r2(i,2).ge.(ymx-2.*rdx1-1.*dy))) then
!                     ncheck = 1
!                 end if
!
!                 if (ncheck.eq.0) then
!                    nhave1 = nhave1+1
!                    r3(nhave1,:) =r2(i,:)
!                 end if
!                 !check if the crystals are in mafic magma, and delete ones in felsic magma
!
!               end do
!
!
!
!           ! check whether this is overlaping
!           !allocate(dis(nhave,nhave))
!              dis = 1000
!
!              do k1 = 1,nhave1
!                do k2 = 1,nhave1
!                   if (k2.ne.k1) then
!                     dis = sqrt((r3(k1,1)-r3(k2,1))**2.+(r3(k1,2)-r3(k2,2))**2.) !(k1,k2)
!                   end if
!                end do
!              end do
!
!
!
!               !track the crystals and identify the crystal count in total (xcnt)
!               if (ns_needed.gt.nhave1) then
!                  ns_needed = nhave1
!               end if
!               do i = 1,ns
!                  xo_new(i) = xo(i)
!                  yo_new(i) = yo(i)
!                  itracker_new(i) = itracker(i)
!               end do
!               do i = 1,ns_needed
!                  ns = ns +1
!                  xcntr = xcntr+1.
!                  xo_new(ns) = r3(i,1)
!                  yo_new(ns) = r3(i,2)
!                  itracker_new(ns) = xcntr
!               end do
!           end if
!
!           ns_new = 0
!
!
!           if (ns_needed.lt.0) then
!             ns_check = 0
!             ns_pres = 0
!
!                do i = 1,ns
!
!
!                 if (xo(i).gt.(axmn).and.xo(i).lt.(axmx).and.yo(i).gt.(aymn).and.yo(i).lt.(aymx)) then
!                   ns_pres = ns_pres +1
!                   if (ns_check.lt.abs(ns_needed)) then
!                       ns_check = ns_check + 1
!                    else
!                       ns_new = ns_new+1
!                       xo_new(ns_new) = xo(i)
!                       yo_new(ns_new) = yo(i)
!                       itracker_new(ns_new) = itracker(i)
!                   end if
!                 else
!                    ns_new = ns_new+1
!                    xo_new(ns_new) = xo(i)
!                    yo_new(ns_new) = yo(i)
!                    itracker_new(ns_new) = itracker(i)
!                 end if
!               end do
!!             write(*,*)'made it into subroutine'
!!             print *, ns_check, ns_needed,ns_pres
!
!             ns = ns - ns_check
!!         print *, ns_check, ns
!          end if
!          if (ns_needed.eq.0) then
!            do i = 1,ns
!             xo_new(i) = xo(i)
!             yo_new(i) = yo(i)
!             itracker_new(i) = itracker(i)
!            end do
!          end if
!
!!         open(10,file="data_crystal.dat")
!
!
!           do i = 1,(ns)
!             new_crystals(i,1) = xo_new(i)
!             new_crystals(i,2) = yo_new(i)
!             inew_tracker(i) = itracker_new(i)
!!             write (10,*) xo_new(i),yo_new(i)
!           end do
!!         close(10)
!         do i = 1,ns
!           xo(i) = xo_new(i)
!           yo(i) = yo_new(i)
!           itracker(i) = itracker_new(i)
!         end do
!
!     end subroutine smart_crystal_adding

     subroutine smart_crystal_adding(rdx,ns_needed,xo,yo,axmx,axmn,aymx,aymn,phigas,itracker,&
                                      ixcntr)
     !!!It takes in the radius of the crystal (rdx), the number of crystals needed or subtracted to stay at equilibrium (ns_needed), current crystal locations (xo,yo), the bounds of the cell we are filling (xmx,ymn,xmx,ymx), and the size of the domain the crystals are averaged over
     !!!!Output is crystal locations (xo_new, yo_new)

          implicit real*8(a-h,o-z)

          common/grid/dx,dy,xmx,xmn,ymx,ymn,nx,ny,tt,tb,crmx,ns,nsp,nsl
          common/param/g,sp,ubc,uout,mint

          dimension :: xo(10000),yo(10000),phigas(nx,ny),itracker(10000)
          allocatable :: x_level1(:),y_level1(:),r(:,:),r1(:,:),r2(:,:),r3(:,:),&
                         new_crystals(:,:),xo_new(:),yo_new(:),inew_tracker(:),&
                         itracker_new(:)
          ! get as many random number as you want by change 15000 below
          allocate(r(15000,2),r1(15000,2),r2(15000,2),r3(15000,2),xo_new(15000),&
                  yo_new(15000),new_crystals(10000,2),inew_tracker(15000),itracker_new(15000))
          ! creat the random number at first
          CALL true_random(r,15000,2)
          ! Then define the geometry


             nsx = 0
             li = 0.
             pii = 4.*atan(1.d0)
!

             ns_new = 0
             print *, 'ns_needed', ns_needed
!             print *, pii, axmx, aymx, axmn, aymn, axmn + 2.*rdx,aymn + 2.*rdx

             if (ns_needed.gt.0) then

                do i = 1,15000
                  !!placing the crystals into the domain
                  r1(i,1) = (axmx-axmn)*r(i,1)+axmn!((axmx-4.*rdx-2.*dble(dx))-axmn)*r(i,1) + &
                             !axmn + 2.*rdx+1.*dble(dx)
                  r1(i,2) = (aymx-aymn)*r(i,2) + aymn!((aymx-4.*rdx-2.*dble(dy))-aymn)*r(i,2) + &
                             !aymn + 2.*rdx+1.*dble(dy)

                end do


                ! delete the number overlaping with each other
                nhave = 1

               r2(nhave,:) = r1(nhave,:)
               do i = 2,15000
               !! make sure the crystals are not on top of each other
                  ncheck = 0
                  do j = 1,nhave
                    if (sqrt((r1(i,1)-r2(j,1))**2+(r1(i,2)-r2(j,2))**2).le.2.5*rdx) then !!originally 3
                       ncheck = 1

                       exit
                    end if
                  end do
                  if (ncheck.eq.0) then
                    nhave = nhave + 1
                    r2(nhave,:) = r1(i,:)
                  end if
               end do
               print *, 'nhave', nhave
               !check if any are overlapping with current crysals within 3rd
               nhave1 = 0
               do i = 1,nhave
                 ncheck = 0
                 do j = 1,ns
                   if (sqrt((r2(i,1)-xo(j))**2+(r2(i,2)-yo(j))**2).le.2.5*rdx) then
                      ncheck = 1
                      print *, 'crystals are too close to existing crystals', sqrt((r2(i,1)-xo(j))**2+(r2(i,2)-yo(j))**2)
                      exit
                   end if
                 end do
                !check if the crystals are in mafic magma, and delete ones in felsic magma
                 do j = 1,ny
                    do k = 1,nx
!                       if (dble(k)*dx.ge.(r2(i,1)-rdx)).and.dble(k)*dx
                     if ((dble(k)*dx).ge.(r2(i,1)-2.5*rdx).and.(dble(k)*dx).le.(r2(i,1)+2.5*rdx).and.&
                         (dble(j)*dy).ge.(r2(i,2)-2.5*rdx).and.(dble(j)*dy).le.(r2(i,2)+2.5*rdx)) then
                         if (phigas(k,j).ge.-0.) then
                            ncheck = 1
                            print *, 'phigas is greater than -0', phigas(k,j)
                         end if
                       end if
                     end do
                 end do
                 !check if the crystals are by the walls
                 if (r2(i,1).le.(xmn+rdx+6.*dx)) then
                     ncheck = 1
                     print *, 'crystal by xmn', r2(i,1)
                 end if
                 if (r2(i,1).ge.(xmx-rdx-6.*dx)) then
                     ncheck = 1
                     print *, 'crystal by xmx', r2(i,1)
                 end if
                 if (r2(i,2).le.(ymn+rdx+6.*dy)) then
                     ncheck = 1
                     print *, 'crystal by ymn', r2(i,2)
                end if
                if (r2(i,2).ge.(ymx-rdx-6.*dy)) then
                    ncheck = 1
                     print *, 'crystal by ymx', r2(i,2)
                end if
!if ((r2(i,1).le.(xmn+1.5*rdx+1.*dx)).or.(r2(i,1).ge.(xmx-1.5*rdx-1.*dx)).or.&
!   (r2(i,2).le.(ymn+1.5*rdx+1.*dy)).or.(r2(i,2).ge.(ymx-1.5*rdx-1.*dy))) then
!    ncheck = 1
!end if

                 if (ncheck.eq.0) then
                    nhave1 = nhave1+1
                    r3(nhave1,:) =r2(i,:)
                 end if
               end do



!           ! check whether this is overlaping
!           !allocate(dis(nhave,nhave))
!              dis = 1000
!
!              do k1 = 1,nhave1
!                do k2 = 1,nhave1
!                   if (k2.ne.k1) then
!                     dis = sqrt((r3(k1,1)-r3(k2,1))**2+(r3(k1,2)-r3(k2,2))**2) !(k1,k2)
!                   end if
!                end do
!              end do


             print *, 'nhave1', nhave1
               !track the crystals and identify the crystal count in total (xcnt)
               if (ns_needed.gt.nhave1) then
                  ns_needed = nhave1
               end if
               do i = 1,ns
                  xo_new(i) = xo(i)
                  yo_new(i) = yo(i)
                  itracker_new(i) = itracker(i)
               end do
               do i = 1,ns_needed
                  ns = ns +1
                  ixcntr = ixcntr+1
                  xo_new(ns) = r3(i,1)
                  yo_new(ns) = r3(i,2)
                  itracker_new(ns) = ixcntr
               end do
!           end if



!           if (ns_needed.lt.0) then
            elseif (ns_needed.lt.0) then
             ns_check = 0
             ns_pres = 0

                do i = 1,ns


                 if (xo(i)+rdx.gt.(axmn).and.xo(i)-rdx.lt.(axmx).and.yo(i)+rdx.gt.(aymn).and.yo(i)-rdx.lt.(aymx)) then
                   ns_pres = ns_pres +1
                   if (ns_check.lt.abs(ns_needed)) then
                       ns_check = ns_check + 1
                    else
                       ns_new = ns_new+1
                       xo_new(ns_new) = xo(i)
                       yo_new(ns_new) = yo(i)
                       itracker_new(ns_new) = itracker(i)
                   end if
                 else
                    ns_new = ns_new+1
                    xo_new(ns_new) = xo(i)
                    yo_new(ns_new) = yo(i)
                    itracker_new(ns_new) = itracker(i)
                 end if
               end do
             write(*,*)'made it into negative ns_needed sub'
             print *, ns_check, ns_needed,ns_pres

             ns = ns - ns_check
             print *, ns_check, ns, ns_new
!          end if
!          if (ns_needed.eq.0) then
          elseif (ns_needed.eq.0) then
            do i = 1,ns
             xo_new(i) = xo(i)
             yo_new(i) = yo(i)
             itracker_new(i) = itracker(i)
            end do
          end if


           do i = 1,(ns)
             new_crystals(i,1) = xo_new(i)
             new_crystals(i,2) = yo_new(i)
             inew_tracker(i) = itracker_new(i)
           end do
         print *, ns_needed
         do i = 1,ns
           xo(i) = xo_new(i)
           yo(i) = yo_new(i)
           itracker(i) = itracker_new(i)
         end do

     end subroutine smart_crystal_adding



!---------------------------------------------------------------------------------
!*********************************************************************************
!
!*********************************************************************************
!---------------------------------------------------------------------------------
!
!
      subroutine smart_tracker_adding(rdx1,rdx2,ns_needed,xo,yo,x_tr,y_tr,axmx,axmn,aymx,aymn,&
                                     itracker_new)
     !!!It takes in the radius of the crystal (rdx1,rdx2-lagrangian tracker), the number of crystals needed or subtracted to stay at equilibrium (ns_needed), current crystal locations (xo,yo), the bounds of the cell we are filling (xmx,ymn,xmx,ymx), and the size of the domain the crystals are averaged over
     !!!!Output is crystal locations (xo_new, yo_new)

          implicit real*8(a-h,o-z)
          common/grid/dx,dy,xmx,xmn,ymx,ymn,nx,ny,tt,tb,crmx,ns,nsp,nsl
          common/param/g,sp,ubc,uout,mint
!          integer, intent(in) :: ns_needed,nav
!          real, intent(in) :: rdx,axmn,aymn,xo,yo
!          real, intent(out) :: new_crystals
          dimension :: xo(ns),yo(ns),x_tr(ns_needed),y_tr(ns_needed),itracker_new(ns_needed)
          allocatable :: x_level1(:),y_level1(:),r(:,:),r1(:,:),r2(:,:),r3(:,:)
          ! get as many random number as you want by change 15000 below
          allocate(r(15000,2),r1(15000,2),r2(15000,2),r3(15000,2))
          ! creat the random number at first
          CALL true_random(r,15000,2)
          ! Then define the geometry



             nsx = 0
             li = 0.
             pii = 4.*atan(1.d0)


                do i = 1,15000
                  !!placing the crystals into the domain
                  r1(i,1) = (axmx-axmn)*r(i,1)+axmn
                  r1(i,2) = (aymx-aymn)*r(i,2) + aymn

                end do


                ! delete the number overlaping with each other
                nhave = 1

               r2(nhave,:) = r1(nhave,:)
               do i = 2,15000
               !! make sure the crystals are not on top of each other
                  ncheck = 0
                  do j = 1,nhave
                    if (sqrt((r1(i,1)-r2(j,1))**2.+(r1(i,2)-r2(j,2))**2.).le.3.*rdx2) then !!originally 3
                       ncheck = 1

                       exit
                    end if
                  end do
                  if (ncheck.eq.0) then
                    nhave = nhave + 1
                    r2(nhave,:) = r1(i,:)
                  end if
               end do
               !check if any are overlapping with current crysals within 3rd
               nhave1 = 0

               do i = 1,nhave
                 ncheck = 0
                 do j = 1,ns
                   if (sqrt((r2(i,1)-xo(j))**2.+(r2(i,2)-yo(j))**2.).le.1.1*rdx1) then
                      ncheck = 1
                      exit
                   end if
                 end do

                 !check if the crystals are by the walls
                 if ((r2(i,1).le.(xmn+rdx1+7.*dx)).or.(r2(i,1).ge.(xmx-rdx1-7.*dx)).or.&
                    (r2(i,2).le.(ymn+rdx1+7.*dy)).or.(r2(i,2).ge.(ymx-rdx1-7.*dy))) then
                     ncheck = 1
                 end if

                 if (ncheck.eq.0) then
                    nhave1 = nhave1+1
                    r3(nhave1,:) =r2(i,:)
                 end if
                 !check if the crystals are in mafic magma, and delete ones in felsic magma

               end do



           ! check whether this is overlaping
           !allocate(dis(nhave,nhave))
              dis = 1000

              do k1 = 1,nhave1
                do k2 = 1,nhave1
                   if (k2.ne.k1) then
                     dis = sqrt((r3(k1,1)-r3(k2,1))**2.+(r3(k1,2)-r3(k2,2))**2.) !(k1,k2)
                   end if
                end do
              end do



               !track the crystals and identify the crystal count in total (xcnt)
               if (ns_needed.gt.nhave1) then
                  ns_needed = nhave1
               end if
               do i = 1,ns_needed
                  ixcntr = ixcntr+1
                  x_tr(ixcntr) = r3(i,1)
                  y_tr(ixcntr) = r3(i,2)
                  itracker_new(ixcntr) = ixcntr
               end do

     end subroutine smart_tracker_adding









!---------------------------------------------------------------------------------
!*********************************************************************************
!
!*********************************************************************************
!---------------------------------------------------------------------------------

      subroutine random_interface(rdx,xo,yo,ixcntr,itracker,axmx,axmn,aymx,aymn,iller)

          implicit real*8(a-h,o-z)
          common/grid/dx,dy,xmx,xmn,ymx,ymn,nx,ny,tt,tb,crmx,ns,nsp,nsl
          common/param/g,sp,ubc,uout,mint
          dimension :: xo(ns),yo(ns),phigas(nx,ny),itracker(ns)
          allocatable :: x_level1(:),y_level1(:),r(:,:),r1(:,:),r2(:,:),amod(:)
          ! get as many random number as you want by change 15000 below
          allocate(r(15000,2),r1(15000,2),r2(15000,2),amod(nx*ny))
          ! creat the random number at first
          CALL true_random(r,15000,2)
          ! Then define the geometry

!             axmx = 0.09                         ! maximum x coordinate value
!             axmn = 0.07                       ! min x
!             aymx = 0.1                         !max y
!             aymn = 0.                         ! min y

             nsx = 0
             li = 0.
             pii = 4.*atan(1.d0)
!             vcrystal = pii*rdx**2.
!             alayervol = (axmx-axmn)*(aymx-aymn)
!             ns = int(dble(nsp)/100.*alayervol/vcrystal) !dble(nsp)
!             print *, iller
             !allocate(amod(ny_level1))
!             print *, pii
             ! delete the number overlaping the wall
!             open(10,file="data_particle_amod.dat")
              do j=1,ny!81
                ay=dble(dy)*dble(j)
                   !amod(j)=xmx/5.+xmx/10.*cos(1./2.*pii*ay/ymx)+min(dx,dy)*2.
                amod(j) = axmx !- axmx/2.
!                write (10,*) amod(j), ay
              end do
!             close(10)

!             open(10,file="data_particle_amod_loc1.dat")
!               do i = 1,15000
!                 write (10,*) r(i,1), r(i,2)
!               end do
!             close(10)

!             open(10,file="data_particle_amod_loc2.dat")
                do i = 1,15000
                  !!placing the crystals into the domain
                  r(i,1) = (axmx-axmn-2.*rdx)*r(i,1)+axmn!
                  r(i,2) = (aymx-aymn-2.*rdx)*r(i,2) + 3.*rdx+aymn!

                  if (r(i,1).ge.(axmn + 2.*rdx+1.*dble(dx)).and.r(i,2).ge.(aymn + 2.*rdx+1.*dble(dy))) then
                    li = li +1.
                    doman = r(i,1)+2.*rdx
                    do j = 1,(ny-1)
                      if (doman<amod(j).and.r(i,2).ge.dble(dy)*dble(j).and.r(i,2).lt.dble(dy)*dble(j+1)) then

                         nsx = nsx + 1
                         r1(nsx,1) = r(i,1)
                         r1(nsx,2) = r(i,2)
                      end if
                    end do
                  end if
                end do
!            close(10)
!            print *, nsx
!            print *, li
  !          open(10,file="data_particle_initial1.dat")
  !          do k=1,nsx
  !             write(10,*) r1(k,1),r1(k,2)
  !          end do
  !          close(10)
!            write(*,*) nsx
            ! delete the number overlaping with each other
            nhave = 1
            r2(nhave,:) = r1(nhave,:)
            do i = 2,nsx
            !! make sure the crystals are not on top of each other
               ncheck = 0
               do j = 1,nhave
                 if (sqrt((r1(i,1)-r2(j,1))**2.+(r1(i,2)-r2(j,2))**2.).le.3.*rdx) then !!originally 3
                   ncheck = 1
                   exit
                 end if
               end do
               if (ncheck.eq.0) then
                 nhave = nhave + 1
                 r2(nhave,:) = r1(i,:)
               end if
           end do
!           write(*,*) nhave

           ! check whether this is overlaping
           !allocate(dis(nhave,nhave))
           dis = 1000

           do k1 = 1,nhave
             do k2 = 1,nhave
                if (k2.ne.k1) then
                  dis = sqrt((r2(k1,1)-r2(k2,1))**2.+(r2(k1,2)-r2(k2,2))**2.) !(k1,k2)
                end if
             end do
           end do

!3.1415926535897931
!62384
!131686
!62384
!62384

          do k1 = 1,ns
             ixcntr = ixcntr+1
             itracker(k1) = ixcntr
             xo(k1) = r2(k1,1)
             yo(k1) = r2(k1,2)
          end do

           !write(*,*) dis

           ! out put the locations
!           open(10,file="data_particle_initial_RT.dat")
!            do k=1,nhave
!              write(10,*) r2(k,1),r2(k,2)
!            end do
!           close(10)


     end subroutine random_interface

!---------------------------------------------------------------------------------
!*********************************************************************************
!
!*********************************************************************************
!---------------------------------------------------------------------------------

      subroutine init(p,u,v,rho1x,rho1y,rho2x,rho2y,xmu,rpx,rpy,xout,xout0,x,y,uin,tmax,intf,Tem,lm,rho,&
                nwrite,vf,phigas,rhog,xmug,rhof,xmuf,xit,w,sig,rdx,&
                p_rho_f,p_mu_f,p_rho_m,p_mu_m,cnctr,phisolid,xo,yo,axmn,navxlct,navylct,xlct,&
                ylct,u_lct,v_lct)

        implicit real*8(a-h,o-z)
        common/grid/dx,dy,xmx,xmn,ymx,ymn,nx,ny,tt,tb,crmx,ns,nsp,nsl
        common/param/g,sp,ubc,uout,mint
        dimension :: p(nx,ny),u(nx,ny),v(nx,ny),xmu(nx,ny),&
                     x(nx),y(ny),b(ny),uin(ny+1),Tem(nx,ny),&
                     intf(nx,ny),rho1x(nx,ny),rho1y(nx,ny),&
                     rho2x(nx,ny),rho2y(nx,ny),rho(nx,ny),&
                     rpx(nx,2),rpy(ny,2),xout(nwrite),xout0(31),&
                     vf(nx,ny),phigas(nx,ny),w(nx,ny),sig(nx,ny,2),&
                     xmuf(nx,ny),rhof(nx,ny),xmug(nx,ny),rhog(nx,ny),&
                     axmuf(nx,ny),arhof(nx,ny),axmug(nx,ny),arhog(nx,ny),&
                     p_rho_f(2),p_mu_f(2),p_rho_m(3),p_mu_m(2),cnctr(navxlct,navylct),&
                     phisolid(ns,ny),xo(ns),yo(ns),xlct(navxlct,navylct),ylct(navxlct,navylct),&
                     u_lct(navxlct,navylct),v_lct(navxlct,navylct)

        allocatable::vol(:),xco(:),yco(:),phicorner(:,:,:),s_part(:,:),xx(:,:,:),yy(:,:,:),&
                     ixmn(:),ixmx(:),iymn(:),iymx(:),n_s(:),phivf(:,:,:),phivol(:,:),&
                     randm_perpx(:,:),randm_perpy(:,:)
        allocate(vol(nx),xco(4),yco(4),phicorner(nx,ny,3),xx(nx,3,3),yy(ny,3,3),s_part(nx,ny),&
                 ixmn(10000),ixmx(10000),iymn(10000),iymx(10000),n_s(ns),phivf(nx,ny,3),&
                 phivol(nx,ny),randm_perpx(navxlct,navylct),randm_perpy(navxlct,navylct))

        p = 0.
        u = 0.
        v = 0.


        !output times
        do i = 1,nwrite
           xout(i)=0.+1./dble(nwrite-1)*dble(i-1)
        end do

        xout = xout * tmax

        do j = 1,31
           xout0(j) = 22.5+1.*dble(j-1)
        end do

!----- initialize parabolic velocity profile at the inlet -------------
        do i = 1,nx
           x(i) = dx*float(i-1)+xmn
        end do
        do j = 1,ny
           y(j) = dy*float(j-1)+ymn
        end do

!--------bounadry conditions on pressure---------------------------


        do i = 1,nx
           rpx(i,1) = 0. !bc at j=1
           rpx(i,2) = 0. !bc at j=ny
        end do
        do j = 1,ny
           rpy(j,1) = 0. !bc at i=1
           rpy(j,2) = 0. !bc at i=nx
        end do

!-------2D solid objects-----------------------------------------------------------
        !intf marks cells that are solid (intf(i,j) = 1)
        !mint counts cells that are solid
         mint = 0
         intf = 0

!-------------2D fluid objects-------------------------------------
        ! calculate the volume fraction for three phases
        pii = 4.*atan(1.d0)
        do j = 1,ny
        b(j) = 4.*(xmx-xmn)/4. !+ xmx/10.*cos(3./2.*pii*y(j)/ymx)
        end do
!        if (nsp.eq.5) then
         at = 2250.!reference is 750 !160. !200 at 10 crystal diameter distance, 160 at 30 dx distance
!        else if (nsp.eq.10.) then
!        at = 750.
!        else if (nsp.eq.0.)  then
!        at = 10.
!        else if (nsp.eq.15) then
!        at = 1200.
!        end if
        alp = 1E-7

        do j = 1,ny
           do i = 1,nx
              !phigas(i,j) = x(i)-(xmx+xmn)/2. + dx*sin((y(j)-ymn)/(ymx-ymn)*2.*pii)
              phigas(i,j) = -100.!x(i)-b(j) !if + then top, if - then bottom
              Tem(i,j) = (tt+tb)/2.-(tb-tt)/2.*erf((x(i)-b(j))/(2.*sqrt(alp*at)))!xmx)/(2.*sqrt(alp*at)))!
!              if (x(i).ge.b+rdx*10.) then
!                Tem(i,j) = tt
!              else if (x(i).ge.b.and.x(i).le.b+rdx*10.) then
!                Tem(i,j) = 814.
!              else if (x(i).lt.b.and.x(i).ge.b-rdx*10.) then
!                Tem(i,j) = 1026.
!              else if (x(i).lt.b-rdx*10.) then
!                Tem(i,j) = tb
!             end if
           end do
        end do
!      Idenitify particle volume distribution
         phivf = 0.
!
!        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!        !!!!!!Defining crystal area volume!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        pii = 4.*atan(1.d0)

!
        do i = 1,nx
           xx(i,1,1) = dx*float(i-1)+xmn                !coordinate for bubble on level 1
           xx(i,2,1) = dx*float(i-1)+xmn                !coordinate for particle on level 1
           xx(i,1,2) = dx*float(i-1)+xmn                !coordinate for volume fraction for u on level 1
           xx(i,2,2) = dx*float(i-1)+xmn-dx/2.          !coordinate for volume fraction for v on level 1
           xx(i,3,2) = dx*float(i-1)+xmn-dx/2.          !coordinate for volume fraction for w on level 1
           xx(i,1,3) = dx*float(i-1)+xmn+dx/2.          !coordinate for u update on level 1
           xx(i,2,3) = dx*float(i-1)+xmn                !coordinate for v update on level 1
           xx(i,3,3) = dx*float(i-1)+xmn
        end do
        do j = 1,ny
           yy(j,1,1) = dy*float(j-1)+ymn                !coordinate for bubble on level 1
           yy(j,2,1) = dy*float(j-1)+ymn                !coordinate for parkicle on level 1
           yy(j,1,2) = dy*float(j-1)+ymn-dy/2.          !coordinate for volume fraction for u on level 1
           yy(j,2,2) = dy*float(j-1)+ymn                !coordinate for volume fraction for v on level 1
           yy(j,3,2) = dy*float(j-1)+ymn-dy/2.          !coordinate for volume fraction for w on level 1
           yy(j,1,3) = dy*float(j-1)+ymn                !coordinate for u update on level 1
           yy(j,2,3) = dy*float(j-1)+ymn+dy/2.          !coordinate for v update on level 1
           yy(j,3,3) = dy*float(j-1)+ymn
        end do
!
!
        do k = 1, ns
!
!!--------------------------------limit the grid for crystals on level 2------------------------------------
             ixmn(k) = floor(((xo(k)-rdx)-xmn)/dx)-1-1
             ixmx(k) = floor(((xo(k)+rdx)-xmn)/dx)+2+1
             iymn(k) = floor(((yo(k)-rdx)-ymn)/dx)-1-1
             iymx(k) = floor(((yo(k)+rdx)-ymn)/dx)+2+1
!
         end do
         do k = 1,ns
            n_s(k) = k
         end do

        do k1 = 1,ns
           do k2 = 1,3
              do j = iymn(k1),iymx(k1)
                 do i = ixmn(k1),ixmx(k1)
                    s_part(i,j) = - rdx + sqrt((xx(i,k2,2)-xo(n_s(k1)))**2. + (yy(j,k2,2)-yo(n_s(k1)))**2.)
                    phivf(i,j,k2) = 0.
                 end do
              end do
              !-------------volume fraction for particle------------
              call volume_frac(s_part(:,:),phivf(:,:,k2),ixmn(k1),ixmx(k1),iymn(k1),iymx(k1),&
                               xx(:,k2,2),yy(:,k2,2),rdx,xo(n_s(k1)),yo(n_s(k1)),k1)
           end do
            do j=iymn(k1),iymx(k1)
               do i=ixmn(k1),ixmx(k1)
                 phivol(i,j) = phivf(i,j,2)
               end do
             end do
         end do

!      Concentration Distribution
!       cnctr = 0.
!       cnctcx = 1.
!       cnctcl = 0.03
!       cnctcb = cnctcl*(1-dble(nsp)/100.) + cnctcx*(dble(nsp)/100.)
!       do j = 1,ny
!          do i = 1,nx
!            if (x(i).le.axmn) then
!                cnctr(i,j) = cnctcb
!            else
!                cnctr(i,j) = cnctcx*phivol(i,j) + cnctcl*(1-phivol(i,j))
!            end if
!          end do
!       end do
              dxlct = (xmx-4.*dx)/dble(navxlct)
              dylct = (ymx-4.*dy)/dble(navylct)
              cnctr = 0.
              cnctcx = 1.
              cnctcl = 0.03
              cnctcb = cnctcl*(1-dble(nsp)/100.) + cnctcx*(dble(nsp)/100.)
              u_lct = 0.
              v_lct = 0.
              xlct  = 0.
              ylct  = 0.
!              print *, axmn
              call true_random(randm_perpx,navxlct,navylct)
              call true_random(randm_perpy,navxlct,navylct)


              do j = 1, navylct
                 do i = 1,navxlct
                    xlct(i,j) = 2.*dx + dble(dxlct)*dble(i) + dble(dxlct)*3./4.*randm_perpx(i,j) !!!perturbing the spacing
                    ylct(i,j) = 2.*dy + dble(dylct)*dble(j) + dble(dylct)*3./4.*randm_perpy(i,j)

                    !!! Need to assign a concentration for each point


                   if (xlct(i,j).le.axmn) then
                        cnctr(i,j) = cnctcb
                    else
                        cnctr(i,j) = cnctcl

                       do l = 1,ns

                        if (sqrt((xlct(i,j)-xo(l))**2.+(ylct(i,j)-yo(l))**2.).le.0.8*rdx) then
                           cnctr(i,j) = cnctcx
                        end if

                      end do

                    end if

                end do
             end do



!        open(10,file="data_magma_chamber_init_cnctr.dat")
!         do j=1,navylct
!          do i=1,navxlct
!            write(10,*) i, j, xlct(i,j), ylct(i,j), cnctr(i,j)
!           end do
!          end do
!        close(10)

        ! Temperature distribution

        do j = 1,ny
           do i = 1,nx
             if (phigas(i,j).ge.0.) then
                Temm = Tem(i,j)
              else
                Temm = Tem(i,j)
              end if

            arhof(i,j) = 1000.*(p_rho_m(1)*Temm + p_rho_m(2))
            axmuf(i,j) = 10.**(p_mu_m(1)*Temm + p_mu_m(2))
            arhog(i,j) = 1000.*(p_rho_m(1)*Temm**2 + p_rho_m(2)*Temm + p_rho_m(3))
            axmug(i,j) = 10.**(p_mu_m(1)*Temm + p_mu_m(2))
           end do
        end do
      do j = 1,ny
        do i = 1,nx
          vf(i,j) = (1.d0 + derf(-phigas(i,j)/sp))/2.d0
          xmuf(i,j) =  axmuf(i,j) - vf(i,j)*axmuf(i,j) + vf(i,j)*axmug(i,j)
          rhof(i,j) =  arhof(i,j) - vf(i,j)*arhof(i,j) + vf(i,j)*arhog(i,j)
          xmug(i,j) =  axmuf(i,j) - vf(i,j)*axmuf(i,j) + vf(i,j)*axmug(i,j)
          rhog(i,j) =  arhof(i,j) - vf(i,j)*arhof(i,j) + vf(i,j)*arhog(i,j)
        end do
      end do


        do j = 1,ny
          do i = 1,nx
             vf(i,j) = (1.d0 + derf(-phigas(i,j)/sp))/2.d0
             xmu(i,j) =  xmuf(i,j) - vf(i,j)*xmuf(i,j) + vf(i,j)*xmug(i,j)
             rho(i,j) =  rhof(i,j) - vf(i,j)*rhof(i,j) + vf(i,j)*rhog(i,j)
          end do
        end do

        open(10,file="data_magma_chamber_init.dat")
         write(10,*) 'TITLE = "Initial Condition"'
         write(10,*) 'VARIABLES = x,y,xmuf,xmug,xmu,rhof,rhog,rho'
         write(10,*) "ZONE T = ",'"Rectanguler zone"'
         write(10,*) "I=",nx,",J=",ny,",F=POINT"
         write(10,*) "DT=(SINGLE SINGLE SINGLE SINGLE)"
         do j=1,ny
          do i=1,nx
            write(10,*) x(i),y(j),axmuf(i,j),axmug(i,j),xmu(i,j), &
                        arhof(i,j),arhog(i,j),rho(i,j),vf(i,j),phivol(i,j)
            end do
          end do
        close(10)


!------------apply interfacial jump conditions--------------------------




        call curvature(x,y,rho1x,rho2x,rho1y,rho2y,xmu,u,v,uin,intf,Tem,dt,lm,rho,&
                       vf,phigas,rhog,xmug,rhof,xmuf,xit,w,sig)

      end subroutine init



!------------------------------------------------------------------------
!------------------------------------------------------------------------
!------------------------------------------------------------------------


      subroutine curvature(x,y,rho1x,rho2x,rho1y,rho2y,xmu,u,v,uin,intf,Tem,dt,lm,rho,&
                           vf,sgas,rhog,xmug,rhof,xmuf,xit,w,sig)

        implicit real*8(a-h,o-z)
        common/param/g,sp,ubc,uout,mint
        common/grid/dx,dy,xmx,xmn,ymx,ymn,nx,ny,tt,tb,crmx,ns,nsp,nsl10
        dimension :: x(nx),y(ny),rho1x(nx,ny),uin(ny+1),rho(nx,ny),&
                     rho2x(nx,ny),rho1y(nx,ny),rho2y(nx,ny),u(nx,ny),v(nx,ny),&
                     xmu(nx,ny),Tem(nx,ny),intf(nx,ny),phiw(nx,ny,3),&
                     vf(nx,ny),sgas(nx,ny),ssolid(nx,ny),w(nx,ny),sig(nx,ny,2),&
                     rhof(nx,ny),xmuf(nx,ny),xmug(nx,ny),rhog(nx,ny)
        allocatable :: Temt(:,:),convT(:,:),diffT(:,:),vff(:,:),cap(:,:),&
                       tkth(:,:),pf(:,:),dns(:,:),xip(:,:),adf(:,:),crv(:,:),adfx(:,:),adfy(:,:)

        nipp = max(nx,ny) * min(nx,ny)/10
        allocate(Temt(nx,ny),convT(nx,ny),diffT(nx,ny),vff(nx,ny),cap(nx,ny),tkth(nx,ny),&
                 pf(nx,ny),dns(nx,ny),xip(nipp,6),adf(nx,ny),crv(nx,ny),adfx(nx,ny),adfy(nx,ny))

        do j = 1,ny
          do i = 1,nx
             vf(i,j) = (1.d0 + derf(-sgas(i,j)/sp))/2.d0
             xmu(i,j) =  xmuf(i,j) - vf(i,j)*xmuf(i,j) + vf(i,j)*xmug(i,j)
             rho(i,j) =  rhof(i,j) - vf(i,j)*rhof(i,j) + vf(i,j)*rhog(i,j)
          end do
        end do

        !---------------------smooth density--------------------
        rho1x = 0.
        rho2x = 0.
        do j = 1,ny
           do i = 1,nx-1
              if (intf(i,j).eq.0) then
                 if (intf(i+1,j).eq.0) then
                    rho1x(i,j) = 0.
                    rho2x(i,j) = 2.d0/(rho(i+1,j)+rho(i,j))
                 else
                    rho1x(i,j) = 0.
                    rho2x(i,j) = 1./rho(i,j)
                 end if
              end if
           end do
           rho1x(nx,j) = 0.
           rho2x(nx,j) = 1./rho(nx,j)
        end do


        rho1y = 0.
        rho2y = 0.
        do i = 1,nx
           do j = 1,ny-1
              if (intf(i,j).eq.0) then
                 if (intf(i,j+1).eq.0) then
                    rho1y(i,j) = 0.
                    rho2y(i,j) = 2.d0/(rho(i,j+1)+rho(i,j))
                 else
                    rho1y(i,j) = 0.
                    rho2y(i,j) = 1./rho(i,j)
                 end if
              end if
           end do
           rho1y(i,ny) = 0.
           rho2y(i,ny) = 1./rho(i,ny)
        end do

      end subroutine curvature

!------------------------------------------------------------------------
!------------------------------------------------------------------------
!------------------------------------------------------------------------


      subroutine energy(sgas,u,v,uin,Tem,capg,capf,tkthg,tkthf,rhof,xmuf,rhog,xmug,dt,ninlet,&
                 lm,x,y,intf,&
                 p_rho_f,p_mu_f,p_rho_m,p_mu_m)
        implicit real*8(a-h,o-z)
        common/param/g,sp,ubc,uout,mint
        common/grid/dx,dy,xmx,xmn,ymx,ymn,nx,ny,tt,tb,crmx,ns,nsp,nsl
        dimension :: sgas(nx,ny),Tem(nx,ny),u(nx,ny),&
                     v(nx,ny),uin(ny+1),rhof(nx,ny),xmuf(nx,ny),&
                     x(nx),y(ny),intf(nx,ny),rhog(nx,ny),xmug(nx,ny),&
                     arhog(nx,ny),axmug(nx,ny),arhof(nx,ny),axmuf(nx,ny),&
                     p_rho_f(2),p_mu_f(2),p_rho_m(3),p_mu_m(2)


        allocatable :: vf(:,:),Temt(:,:),rho(:,:),cap(:,:),tkth(:,:)

        nipp = max(nx,ny) * min(nx,ny)/10
        allocate(vf(nx,ny),Temt(nx,ny),rho(nx,ny),cap(nx,ny),tkth(nx,ny))

        pii = 4.*atan(1.d0)

        !-------------------update density&thermal diffusivity------------------------
        do j = 1,ny
           do i = 1,nx
              vf(i,j)   = (1.d0 + derf(-sgas(i,j)/sp))/2.d0
              rho(i,j)  = vf(i,j)*(rhog(i,j)-rhof(i,j)) + rhof(i,j)
              cap(i,j)  = vf(i,j)*(capg-capf) + capf
              tkth(i,j) = vf(i,j)*(tkthg-tkthf)+tkthf
           end do
        end do
        Temt = Tem
        do j = 1,ny
           do i = 1,nx
              if (intf(i,j).eq.0) then
                 if (i.eq.1) then
                    if (j.ge.(ny-ninlet)/2+1.and.j.le.(ny-ninlet)/2+ninlet) then
                       sxl = tb*2.-Temt(i,j)
                       ssxl= tkth(i,j)*(Temt(i,j)-(tb*2.-Temt(i,j)))
                       ul = uin(j)
                    else
                       sxl = Temt(i,j)
                       ssxl= tkth(i,j)*(Temt(i,j)-Temt(i,j))
                       ul = uin(j)
                    end if
                 else
                    if (intf(i-1,j).eq.0) then
                       sxl = Temt(i-1,j)
                       ssxl= (tkth(i,j)+tkth(i-1,j))/2.*(Temt(i,j)-Temt(i-1,j))
                       ul = u(i-1,j)
                    else
                       sxl = Temt(i,j)
                       ssxl= (tkth(i,j)+tkth(i-1,j))/2.*(Temt(i,j)-Temt(i,j))
                       ul = 0.
                    end if
                 end if
                 if(i.eq.nx) then
                    sxr = tt*2.-Temt(i,j)
                    ssxr= tkth(i,j)*((tt*2.-Temt(i,j))-Temt(i,j))
                 else
                    if (intf(i+1,j).eq.0) then
                       sxr = Temt(i+1,j)
                       ssxr= (tkth(i+1,j)+tkth(i,j))/2.*(Temt(i+1,j)-Temt(i,j))
                    else
                       sxr = Temt(i,j)
                       ssxr= (tkth(i+1,j)+tkth(i,j))/2.*(Temt(i,j)-Temt(i,j))
                    end if
                 end if
                 ur = u(i,j)

                 if(j.eq.1) then
                    syl = Temt(i,j)
                    ssyl= tkth(i,j)*(Temt(i,j)-Temt(i,j))
                    vl = 0.
                 else
                    if (intf(i,j-1).eq.0) then
                       syl = Temt(i,j-1)
                       ssyl= (tkth(i,j)+tkth(i,j-1))/2.*(Temt(i,j)-Temt(i,j-1))
                       vl = v(i,j-1)
                    else
                       syl = Temt(i,j)
                       ssyl= (tkth(i,j)+tkth(i,j-1))/2.*(Temt(i,j)-Temt(i,j))
                       vl = 0.
                    end if
                 end if

                 if(j.eq.ny) then
                    syr = Temt(i,j)
                    ssyr= tkth(i,j)*(Temt(i,j)-Temt(i,j))
                 else
                    if (intf(i,j+1).eq.0) then
                       syr = Temt(i,j+1)
                       ssyr= (tkth(i,j)+tkth(i,j+1))/2.*(Temt(i,j+1)-Temt(i,j))
                    else
                       syr = Temt(i,j)
                       ssyr= (tkth(i,j)+tkth(i,j+1))/2.*(Temt(i,j)-Temt(i,j))
                    end if
                 end if
                 vr = v(i,j)

                 sm = Temt(i,j)

                 if((ur+ul).ge.0.) then
                    fxr = sm*(ur+ul)/2.
                 else
                    fxr = sxr*(ur+ul)/2.
                 end if

                 if((ur+ul).ge.0.) then
                    fxl = sxl*(ur+ul)/2.
                 else
                    fxl = sm*(ur+ul)/2.
                 end if


                 if((vr+vl).ge.0.) then
                    fyr = sm*(vr+vl)/2.
                 else
                    fyr = syr*(vr+vl)/2.
                 end if

                 if((vr+vl).ge.0.) then
                    fyl = syl*(vr+vl)/2.
                 else
                    fyl = sm*(vr+vl)/2.
                 end if

                 Tem(i,j) = Temt(i,j) - dt * ( (fxr-fxl)/dx + (fyr-fyl)/dy ) &
                            + 1./(rho(i,j)*cap(i,j)) * dt * ( (ssxr-ssxl)/dx**2. + (ssyr-ssyl)/dy**2. )

              end if
           end do
        end do

        do j = 1,ny
          do i = 1,nx
             if (sgas(i,j).ge.0.) then
                Temm = Tem(i,j)
              else
                Temm = Tem(i,j)
              end if
             arhof(i,j) = 1000.*(p_rho_m(1)*Temm + p_rho_m(2))
             axmuf(i,j) = 10.**(p_mu_m(1)*Temm + p_mu_m(2))
             arhog(i,j) = 1000.*(p_rho_m(1)*Temm**2 + p_rho_m(2)*Temm + p_rho_m(3))
             axmug(i,j) = 10.**(p_mu_m(1)*Temm + p_mu_m(2))
          end do
        end do
        do j = 1,ny
         do i = 1,nx
          vf(i,j) = (1.d0 + derf(-sgas(i,j)/sp))/2.d0
          xmuf(i,j) =  axmuf(i,j) - vf(i,j)*axmuf(i,j) + vf(i,j)*axmug(i,j)
          rhof(i,j) =  arhof(i,j) - vf(i,j)*arhof(i,j) + vf(i,j)*arhog(i,j)
          xmug(i,j) =  axmuf(i,j) - vf(i,j)*axmuf(i,j) + vf(i,j)*axmug(i,j)
          rhog(i,j) =  arhof(i,j) - vf(i,j)*arhof(i,j) + vf(i,j)*arhog(i,j)
         end do
       end do

      end subroutine energy

!---------------------------------------------------------------------------------
!---------------------------------------------------------------------------------
      subroutine solve1(u,v,p,uin,rho1x,rho1y,rho2x,rho2y,xmu,rpx,rpy, &
                       dt,dx,dy,intf,control,div,err,res,errm,nx,ny,mu, &
                       inumeric,isymbolic,itr,t,x,y,rho1xs,&
                       rho1ys,rho2xs,rho2ys,w,sig,k_Stoke,phigas)
        implicit real*8(a-h,o-z)
        common/param/g,sp,ubc,uout,mint
        dimension :: control(20),xinfo(90),u(nx,ny),v(nx,ny),p(nx,ny),&
                     intf(nx,ny),x(nx),y(ny),uin(ny+1),&
                     rho1x(nx,ny),rho1y(nx,ny),rho2x(nx,ny), &
                     rho2y(nx,ny),rpx(nx,2),rpy(ny,2),&
                     xmu(nx,ny),rho1xs(nx,ny),rho1ys(nx,ny),rho2xs(nx,ny), &
                     rho2ys(nx,ny),w(nx,ny),sig(nx,ny,2),phigas(nx,ny)
        allocatable :: cmat(:,:),rhs(:),rhsU(:,:),sol(:),xmuh(:,:),fup(:,:),fvp(:,:),tmp(:,:),sov_medi(:,:,:),&
                       ub(:),ut(:),utt(:),ubb(:),ubbb(:),vbb(:),vtt(:),vb(:),vt(:),ul(:),ur(:),uttt(:),&
                       vl(:),vr(:),urr(:),ull(:),vrr(:),vll(:),vlll(:),coeflx(:,:,:),coefrx(:,:,:),vrrr(:),&
                       coefly(:,:,:),coefry(:,:,:),aux(:),bux(:),cux(:),dux(:),cuxp(:),duxp(:),&
                       um(:),umm(:),ummm(:),utr(:,:),auy(:),buy(:),cuy(:),duy(:),cuyp(:),duyp(:),vttt(:),&
                       ustar(:,:),vstar(:,:)

        integer*8 inumeric,isymbolic,istatus,isys
        
        nm = nx*ny - mint
        nmm = nx*ny
        dx2 = 1./dx**2.
        dy2 = 1./dy**2.
        dtx = dt*dx2
        dty = dt*dy2

        eps = 10E-8
        re = 1.

        mumf = 9*nm
        allocate(cmat(nm,18),rhs(nm),rhsU(nmm,2),sol(nm),xmuh(nx+1,ny+1),fup(nx,ny),fvp(nx,ny),tmp(nx,ny),sov_medi(nx,ny,3),&
                 ub(ny+1),ut(ny+1),utt(ny),ubb(ny),ubbb(ny),vbb(ny),vtt(ny),vb(ny),vt(ny),ul(nx),ur(nx),uttt(ny),&
                 vl(nx),vr(nx),urr(nx),ull(nx),vrr(nx),vll(nx),vlll(nx),coeflx(nx,ny,2),coefrx(nx,ny,2),vrrr(nx),&
                 coefly(nx,ny,2),coefry(nx,ny,2),aux(nx),bux(nx),cux(nx),dux(nx),cuxp(nx),duxp(nx),&
                 um(nx),umm(nmm),ummm(ny),utr(nx,ny),auy(ny),buy(ny),cuy(ny),duy(ny),cuyp(ny),duyp(ny),vttt(ny),&
                 ustar(nx,ny),vstar(nx,ny))

        fup = 0.d0
        fvp = 0.d0
        rhs = 0.d0

        coeflx = 0.d0
        coefrx = 0.d0
        coefly = 0.d0
        coefry = 0.d0

        ii = 0

        do j = 1,ny + 1
           ub(j)  = uin(j)                                              !u@i=1/2    t@n
           ut(j)  = 0.                                                  !u@i=nx+3/2 t@n
        end do
        do j = 1,ny
           vb(j)  = v(1,j)*ubc                                             !v@i=0      t@n
           vt(j)  = v(nx,j)*ubc                                            !v@i=nx+1   t@n
        end do

        do i = 1,nx
           ul(i)  = u(i,1)*ubc                                          !u@j=0      t@n
           ur(i)  = u(i,ny)*ubc                                         !u@j=ny+1   t@n

           vl(i)  = 0.                                                  !v@j=1/2    t@n
           vr(i)  = 0.                                                  !v@j=ny+3/2 t@n
        end do

        vcorn = 0.                                                      !v@i=nx+1,j=1/2

!****************************** j = 1, i = 1 ********************************


        j = 1
        i = 1

           ii = ii + 1

           rho = rho1x(i,j) + rho2x(i,j)

           xmuc = (xmu(i,j)+xmu(i+1,j)*(1.-dble(intf(i+1,j)))+xmu(i,j)*dble(intf(i+1,j))&
                  +xmu(i,j+1)*(1.-dble(intf(i,j+1)))+xmu(i,j)*dble(intf(i,j+1))+(xmu(i+1,j+1)*&
                  (1.-dble(intf(i+1,j+1)))+xmu(i+1,j)*dble(intf(i+1,j+1)))*(1.-dble(intf(i+1,j)))+&
                  (xmu(i,j+1)*(1.-dble(intf(i,j+1)))+xmu(i,j)*dble(intf(i,j+1)))*dble(intf(i+1,j)))/4.d0
           xmub = (xmu(i,j)+xmu(i+1,j)*(1.-dble(intf(i+1,j)))+xmu(i,j)*dble(intf(i+1,j)))/2.d0
           xmul = (xmu(i,j)+xmu(i,j+1)*(1.-dble(intf(i,j+1)))+xmu(i,j)*dble(intf(i,j+1)))/2.d0

           xmuh(i+1,j+1) = xmuc
           xmuh(i+1,j) = xmub
           xmuh(i,j+1) = xmul

           coeflx(i,j,1)=dtx*xmu(i,j)*rho/re
           coefrx(i,j,1)=dtx*xmu(i+1,j)*rho/re
           coefly(i,j,1)=dty*xmuh(i+1,j)*rho/re
           coefry(i,j,1)=dty*xmuh(i+1,j+1)*rho/re

           fup(i,j) = ( u(i,j)/dt + g + &
                   (coefrx(i,j,1)-coeflx(i,j,1))*(u(i+1,j)*(1.-dble(intf(i+1,j)))-ub(j))/2./dt + &
                   (coefry(i,j,1)-coefly(i,j,1))*((v(i+1,j)+vl(i+1))*(1.-dble(intf(i+1,j)))-&
                   (v(i,j)+vl(i)))/2./dt) * dt
           
           rhsU(ii,1) = fup(i,j)

           rho = rho1y(i,j) + rho2y(i,j)

           coeflx(i,j,2)=dtx*xmuh(i,j+1)*rho/re
           coefrx(i,j,2)=dtx*xmuh(i+1,j+1)*rho/re
           coefly(i,j,2)=dty*xmu(i,j)*rho/re
           coefry(i,j,2)=dty*xmu(i,j+1)*rho/re

           fvp(i,j) = ( v(i,j)/dt + &
                   (coefrx(i,j,2)-coeflx(i,j,2))*((u(i,j+1)+ub(j+1))-(u(i,j)+ub(j)))/2./dt + &
                   (coefry(i,j,2)-coefly(i,j,2))*(v(i,j+1)*(1.-dble(intf(i,j+1)))-vl(i))/2./dt ) * dt

           rhsU(ii,2) = fvp(i,j)


!-------------------------- j = 1, i = 2,nx-1 --------------------------------

        do i = 2,nx-1


           ii = ii + 1
          
           rho = rho1x(i,j) + rho2x(i,j)

           xmuc = (xmu(i,j)+xmu(i+1,j)*(1.-dble(intf(i+1,j)))+xmu(i,j)*dble(intf(i+1,j))&
                  +xmu(i,j+1)*(1.-dble(intf(i,j+1)))+xmu(i,j)*dble(intf(i,j+1))+(xmu(i+1,j+1)*&
                  (1.-dble(intf(i+1,j+1)))+xmu(i+1,j)*dble(intf(i+1,j+1)))*(1.-dble(intf(i+1,j)))+&
                  (xmu(i,j+1)*(1.-dble(intf(i,j+1)))+xmu(i,j)*dble(intf(i,j+1)))*dble(intf(i+1,j)))/4.d0
           xmub = (xmu(i,j)+xmu(i+1,j)*(1.-dble(intf(i+1,j)))+xmu(i,j)*dble(intf(i+1,j)))/2.d0
           xmul = (xmu(i,j)+xmu(i-1,j)*(1.-dble(intf(i-1,j)))+xmu(i,j)*dble(intf(i-1,j))&
                  +xmu(i,j+1)*(1.-dble(intf(i,j+1)))+xmu(i,j)*dble(intf(i,j+1))+(xmu(i-1,j+1)*&
                  (1.-dble(intf(i-1,j+1)))+xmu(i-1,j)*dble(intf(i-1,j+1)))*(1.-dble(intf(i-1,j)))+&
                  (xmu(i,j+1)*(1.-dble(intf(i,j+1)))+xmu(i,j)*dble(intf(i,j+1)))*dble(intf(i-1,j)))/4.d0

           xmuh(i+1,j+1) = xmuc
           xmuh(i+1,j) = xmub
           xmuh(i,j+1) = xmul

           coeflx(i,j,1)=dtx*xmu(i,j)*rho/re
           coefrx(i,j,1)=dtx*xmu(i+1,j)*rho/re
           coefly(i,j,1)=dty*xmuh(i+1,j)*rho/re
           coefry(i,j,1)=dty*xmuh(i+1,j+1)*rho/re

           fup(i,j) = ( u(i,j)/dt + g + &
                   (coefrx(i,j,1)-coeflx(i,j,1))*(u(i+1,j)*(1.-dble(intf(i+1,j)))-u(i-1,j)&
                   *(1.-dble(intf(i-1,j))))/2./dt + &
                   (coefry(i,j,1)-coefly(i,j,1))*((v(i+1,j)+vl(i+1))*(1.-dble(intf(i+1,j)))&
                   -(v(i,j)+vl(i)))/2./dt ) * dt


           rhsU(ii,1) = fup(i,j)

           rho = rho1y(i,j) + rho2y(i,j)

           coeflx(i,j,2)=dtx*xmuh(i,j+1)*rho/re
           coefrx(i,j,2)=dtx*xmuh(i+1,j+1)*rho/re
           coefly(i,j,2)=dty*xmu(i,j)*rho/re
           coefry(i,j,2)=dty*xmu(i,j+1)*rho/re

           fvp(i,j) = ( v(i,j)/dt  + &
                   (coefrx(i,j,2)-coeflx(i,j,2))*((u(i,j+1)+u(i-1,j+1))*(1.-dble(intf(i,j+1)))&
                   -(u(i,j)+u(i-1,j)))/2./dt + &
                   (coefry(i,j,2)-coefly(i,j,2))*(v(i,j+1)*(1.-dble(intf(i,j+1)))&
                   -vl(i))/2./dt )* dt

           rhsU(ii,2) = fvp(i,j)



        end do

!------------------j = 1, i = nx ----------------------------------------------

        i = nx

           ii = ii + 1

           rho = rho1x(i,j) + rho2x(i,j)

           xmuc = (xmu(i,j)+xmu(i,j+1)*(1.-dble(intf(i,j+1)))+xmu(i,j)*dble(intf(i,j+1)))/2.d0
           xmub = xmu(i,j)
           xmul = (xmu(i,j)+xmu(i-1,j)*(1.-dble(intf(i-1,j)))+xmu(i,j)*dble(intf(i-1,j))&
                  +xmu(i,j+1)*(1.-dble(intf(i,j+1)))+xmu(i,j)*dble(intf(i,j+1))+(xmu(i-1,j+1)*&
                  (1.-dble(intf(i-1,j+1)))+xmu(i-1,j)*dble(intf(i-1,j+1)))*(1.-dble(intf(i-1,j)))+&
                  (xmu(i,j+1)*(1.-dble(intf(i,j+1)))+xmu(i,j)*dble(intf(i,j+1)))*dble(intf(i-1,j)))/4.d0

           xmuh(i+1,j+1) = xmuc
           xmuh(i+1,j) = xmub
           xmuh(i,j+1) = xmul

           coeflx(i,j,1)=dtx*xmu(i,j)*rho/re
           coefrx(i,j,1)=dtx*xmu(i,j)*rho/re
           coefly(i,j,1)=dty*xmuh(i+1,j)*rho/re
           coefry(i,j,1)=dty*xmuh(i+1,j+1)*rho/re

           fup(i,j) = ( u(i,j)/dt + g + &
                   (coefrx(i,j,1)-coeflx(i,j,1))*(ut(j)-u(i-1,j)*(1.-dble(intf(i-1,j))))/2./dt + &
                   (coefry(i,j,1)-coefly(i,j,1))*((vt(j)+vcorn)-(v(i,j)+vl(i)))/2./dt) * dt


           rhsU(ii,1) = fup(i,j)


           rho = rho1y(i,j) + rho2y(i,j)

           coeflx(i,j,2)=dtx*xmuh(i,j+1)*rho/re
           coefrx(i,j,2)=dtx*xmuh(i+1,j+1)*rho/re
           coefly(i,j,2)=dty*xmu(i,j)*rho/re
           coefry(i,j,2)=dty*xmu(i,j+1)*rho/re

           fvp(i,j) = ( v(i,j)/dt + &
                   (coefrx(i,j,2)-coeflx(i,j,2))*((u(i,j+1)+u(i-1,j+1))*(1.-dble(intf(i,j+1)))-&
                   (u(i,j)+u(i-1,j)))/2./dt + &
                   (coefry(i,j,2)-coefly(i,j,2))*(v(i,j+1)*(1.-dble(intf(i,j+1)))-vl(i))/2./dt ) * dt

           rhsU(ii,2) = fvp(i,j)



!***************************** j = 2,ny-1, i = 1 *********************************

        do j = 2,ny-1

           i = 1


              ii = ii + 1

              rho = rho1x(i,j) + rho2x(i,j)

              xmuc = (xmu(i,j)+xmu(i+1,j)*(1.-dble(intf(i+1,j)))+xmu(i,j)*dble(intf(i+1,j))&
                     +xmu(i,j+1)*(1.-dble(intf(i,j+1)))+xmu(i,j)*dble(intf(i,j+1))+(xmu(i+1,j+1)*&
                     (1.-dble(intf(i+1,j+1)))+xmu(i+1,j)*dble(intf(i+1,j+1)))*(1.-dble(intf(i+1,j)))+&
                     (xmu(i,j+1)*(1.-dble(intf(i,j+1)))+xmu(i,j)*dble(intf(i,j+1)))*dble(intf(i+1,j)))/4.d0
              xmub = (xmu(i,j)+xmu(i+1,j)*(1.-dble(intf(i+1,j)))+xmu(i,j)*dble(intf(i+1,j))&
                     +xmu(i,j-1)*(1.-dble(intf(i,j-1)))+xmu(i,j)*dble(intf(i,j-1))+(xmu(i+1,j-1)*&
                     (1.-dble(intf(i+1,j-1)))+xmu(i+1,j)*dble(intf(i+1,j-1)))*(1.-dble(intf(i+1,j)))+&
                     (xmu(i,j-1)*(1.-dble(intf(i,j-1)))+xmu(i,j)*dble(intf(i,j-1)))*dble(intf(i+1,j)))/4.d0
              xmul = (xmu(i,j)+xmu(i,j+1)*(1.-dble(intf(i,j+1)))+xmu(i,j)*dble(intf(i,j+1)))/2.d0

              xmuh(i+1,j+1) = xmuc
              xmuh(i+1,j) = xmub
              xmuh(i,j+1) = xmul

              coeflx(i,j,1)=dtx*xmu(i,j)*rho/re
              coefrx(i,j,1)=dtx*xmu(i+1,j)*rho/re
              coefly(i,j,1)=dty*xmuh(i+1,j)*rho/re
              coefry(i,j,1)=dty*xmuh(i+1,j+1)*rho/re

              fup(i,j) = ( u(i,j)/dt + g  + &
                   (coefrx(i,j,1)-coeflx(i,j,1))*(u(i+1,j)*(1.-dble(intf(i+1,j)))-ub(j))/2./dt + &
                   (coefry(i,j,1)-coefly(i,j,1))*((v(i+1,j)+v(i+1,j-1))*(1.-dble(intf(i+1,j)))-&
                   (v(i,j)+v(i,j-1)))/2./dt) * dt

              rhsU(ii,1) = fup(i,j)
              
              rho = rho1y(i,j) + rho2y(i,j)

              coeflx(i,j,2)=dtx*xmuh(i,j+1)*rho/re
              coefrx(i,j,2)=dtx*xmuh(i+1,j+1)*rho/re
              coefly(i,j,2)=dty*xmu(i,j)*rho/re
              coefry(i,j,2)=dty*xmu(i,j+1)*rho/re
 
              fvp(i,j) = ( v(i,j)/dt + &
                   (coefrx(i,j,2)-coeflx(i,j,2))*((u(i,j+1)+ub(j+1))-(u(i,j)+ub(j)))/2./dt + &
                   (coefry(i,j,2)-coefly(i,j,2))*(v(i,j+1)*(1.-dble(intf(i,j+1)))-&
                   v(i,j-1)*(1.-dble(intf(i,j-1))))/2./dt) *dt

              rhsU(ii,2) = fvp(i,j)




!------------------------- j = 2,ny-1, i = 2,nx-1 ------------------------------
           
           do i = 2,nx-1


                 ii = ii + 1

                 rho = rho1x(i,j) + rho2x(i,j)

                 xmuc = (xmu(i,j)+xmu(i+1,j)*(1.-dble(intf(i+1,j)))+xmu(i,j)*dble(intf(i+1,j))&
                        +xmu(i,j+1)*(1.-dble(intf(i,j+1)))+xmu(i,j)*dble(intf(i,j+1))+(xmu(i+1,j+1)*&
                        (1.-dble(intf(i+1,j+1)))+xmu(i+1,j)*dble(intf(i+1,j+1)))*(1.-dble(intf(i+1,j)))+&
                        (xmu(i,j+1)*(1.-dble(intf(i,j+1)))+xmu(i,j)*dble(intf(i,j+1)))*dble(intf(i+1,j)))/4.d0
                 xmub = (xmu(i,j)+xmu(i+1,j)*(1.-dble(intf(i+1,j)))+xmu(i,j)*dble(intf(i+1,j))&
                        +xmu(i,j-1)*(1.-dble(intf(i,j-1)))+xmu(i,j)*dble(intf(i,j-1))+(xmu(i+1,j-1)*&
                        (1.-dble(intf(i+1,j-1)))+xmu(i+1,j)*dble(intf(i+1,j-1)))*(1.-dble(intf(i+1,j)))+&
                        (xmu(i,j-1)*(1.-dble(intf(i,j-1)))+xmu(i,j)*dble(intf(i,j-1)))*dble(intf(i+1,j)))/4.d0
                 xmul = (xmu(i,j)+xmu(i-1,j)*(1.-dble(intf(i-1,j)))+xmu(i,j)*dble(intf(i-1,j))&
                        +xmu(i,j+1)*(1.-dble(intf(i,j+1)))+xmu(i,j)*dble(intf(i,j+1))+(xmu(i-1,j+1)*&
                        (1.-dble(intf(i-1,j+1)))+xmu(i-1,j)*dble(intf(i-1,j+1)))*(1.-dble(intf(i-1,j)))+&
                        (xmu(i,j+1)*(1.-dble(intf(i,j+1)))+xmu(i,j)*dble(intf(i,j+1)))*dble(intf(i-1,j)))/4.d0

                 xmuh(i+1,j+1) = xmuc
                 xmuh(i+1,j) = xmub
                 xmuh(i,j+1) = xmul

                 coeflx(i,j,1)=dtx*xmu(i,j)*rho/re
                 coefrx(i,j,1)=dtx*(xmu(i+1,j)*(1.-dble(intf(i+1,j)))+xmu(i,j)*dble(intf(i+1,j)))*rho/re
                 coefly(i,j,1)=dty*xmuh(i+1,j)*rho/re
                 coefry(i,j,1)=dty*xmuh(i+1,j+1)*rho/re

                 fup(i,j) = ( u(i,j)/dt + g + &
                   (coefrx(i,j,1)-coeflx(i,j,1))*(u(i+1,j)*(1.-dble(intf(i+1,j)))-&
                   u(i-1,j)*(1.-dble(intf(i-1,j))))/2./dt + &
                   (coefry(i,j,1)-coefly(i,j,1))*((v(i+1,j)+v(i+1,j-1))*(1.-dble(intf(i+1,j)))-&
                   (v(i,j)+v(i,j-1)))/2./dt) * dt

                 rhsU(ii,1) = fup(i,j)

                 rho = rho1y(i,j) + rho2y(i,j)

                 coeflx(i,j,2)=dtx*xmuh(i,j+1)*rho/re
                 coefrx(i,j,2)=dtx*xmuh(i+1,j+1)*rho/re
                 coefly(i,j,2)=dty*xmu(i,j)*rho/re
                 coefry(i,j,2)=dty*xmu(i,j+1)*rho/re

                 fvp(i,j) = ( v(i,j)/dt + &
                   (coefrx(i,j,2)-coeflx(i,j,2))*((u(i,j+1)+u(i-1,j+1))*(1.-dble(intf(i,j+1)))-&
                   (u(i,j)+u(i-1,j)))/2./dt + &
                   (coefry(i,j,2)-coefly(i,j,2))*(v(i,j+1)*(1.-dble(intf(i,j+1)))-&
                   v(i,j-1)*(1.-dble(intf(i,j-1))))/2./dt) * dt

                 rhsU(ii,2) = fvp(i,j)

           end do

!-----------------------j = 2,ny-1, i = nx ----------------------------------


           i = nx

              ii = ii + 1

              rho = rho1x(i,j) + rho2x(i,j)

              xmuc = (xmu(i,j)+xmu(i,j+1)*(1.-dble(intf(i,j+1)))+xmu(i,j)*dble(intf(i,j+1)))/2.d0
              xmub = (xmu(i,j)+xmu(i,j-1)*(1.-dble(intf(i,j-1)))+xmu(i,j)*dble(intf(i,j-1)))/2.d0
              xmul = (xmu(i,j)+xmu(i-1,j)*(1.-dble(intf(i-1,j)))+xmu(i,j)*dble(intf(i-1,j))&
                     +xmu(i,j+1)*(1.-dble(intf(i,j+1)))+xmu(i,j)*dble(intf(i,j+1))+(xmu(i-1,j+1)*&
                     (1.-dble(intf(i-1,j+1)))+xmu(i-1,j)*dble(intf(i-1,j+1)))*(1.-dble(intf(i-1,j)))+&
                     (xmu(i,j+1)*(1.-dble(intf(i,j+1)))+xmu(i,j)*dble(intf(i,j+1)))*dble(intf(i-1,j)))/4.d0

              xmuh(i+1,j+1) = xmuc
              xmuh(i+1,j) = xmub
              xmuh(i,j+1) = xmul

              coeflx(i,j,1)=dtx*xmu(i,j)*rho/re
              coefrx(i,j,1)=dtx*xmu(i,j)*rho/re
              coefly(i,j,1)=dty*xmuh(i+1,j)*rho/re
              coefry(i,j,1)=dty*xmuh(i+1,j+1)*rho/re

              fup(i,j) = ( u(i,j)/dt + g + &
                   (coefrx(i,j,1)-coeflx(i,j,1))*(ut(j)-u(i-1,j)*(1.-dble(intf(i-1,j))))/2./dt + &
                   (coefry(i,j,1)-coefly(i,j,1))*((vt(j)+vt(j-1))-(v(i,j)+v(i,j-1)))/2./dt) * dt

              rhsU(ii,1) = fup(i,j)

              rho = rho1y(i,j) + rho2y(i,j)

              coeflx(i,j,2)=dtx*xmuh(i,j+1)*rho/re
              coefrx(i,j,2)=dtx*xmuh(i+1,j+1)*rho/re
              coefly(i,j,2)=dty*xmu(i,j)*rho/re
              coefry(i,j,2)=dty*xmu(i,j+1)*rho/re

              fvp(i,j) = ( v(i,j)/dt  + &
                   (coefrx(i,j,2)-coeflx(i,j,2))*((u(i,j+1)+u(i-1,j+1))*(1.-dble(intf(i,j+1)))-&
                   (u(i,j)+u(i-1,j)))/2./dt + &
                   (coefry(i,j,2)-coefly(i,j,2))*(v(i,j+1)*(1.-dble(intf(i,j+1)))-&
                   v(i,j-1)*(1.-dble(intf(i,j-1))))/2./dt) *dt

              rhsU(ii,2) = fvp(i,j)

        end do


!*************************** j = ny, i = 1 ***************************************

        j = ny
        i = 1

           ii = ii + 1

           rho = rho1x(i,j) + rho2x(i,j)

           xmuc = (xmu(i,j)+xmu(i+1,j)*(1.-dble(intf(i+1,j)))+xmu(i,j)*dble(intf(i+1,j)))/2.d0
           xmub = (xmu(i,j)+xmu(i+1,j)*(1.-dble(intf(i+1,j)))+xmu(i,j)*dble(intf(i+1,j))&
                  +xmu(i,j-1)*(1.-dble(intf(i,j-1)))+xmu(i,j)*dble(intf(i,j-1))+(xmu(i+1,j-1)*&
                  (1.-dble(intf(i+1,j-1)))+xmu(i+1,j)*dble(intf(i+1,j-1)))*(1.-dble(intf(i+1,j)))+&
                  (xmu(i,j-1)*(1.-dble(intf(i,j-1)))+xmu(i,j)*dble(intf(i,j-1)))*dble(intf(i+1,j)))/4.d0
           xmul = xmu(i,j)

           xmuh(i+1,j+1) = xmuc
           xmuh(i+1,j) = xmub
           xmuh(i,j+1) = xmul

           coeflx(i,j,1)=dtx*xmu(i,j)*rho/re
           coefrx(i,j,1)=dtx*xmu(i+1,j)*rho/re
           coefly(i,j,1)=dty*xmuh(i+1,j)*rho/re
           coefry(i,j,1)=dty*xmuh(i+1,j+1)*rho/re

           fup(i,j) = ( u(i,j)/dt + g + &
                   (coefrx(i,j,1)-coeflx(i,j,1))*(u(i+1,j)*(1.-dble(intf(i+1,j)))-ub(j))/2./dt + &
                   (coefry(i,j,1)-coefly(i,j,1))*((v(i+1,j)+v(i+1,j-1))*(1.-dble(intf(i+1,j)))-&
                   (v(i,j)+v(i,j-1)))/2./dt) * dt

           rhsU(ii,1) = fup(i,j)


           rho = rho1y(i,j) + rho2y(i,j)

           coeflx(i,j,2)=dtx*xmuh(i,j+1)*rho/re
           coefrx(i,j,2)=dtx*xmuh(i+1,j+1)*rho/re
           coefly(i,j,2)=dty*xmu(i,j)*rho/re
           coefry(i,j,2)=dty*xmu(i,j)*rho/re

           fvp(i,j) = ( v(i,j)/dt + &
                   (coefrx(i,j,2)-coeflx(i,j,2))*((ul(i)+ucorn)-(u(i,j)+ub(j)))/2./dt + &
                   (coefry(i,j,2)-coefly(i,j,2))*(vr(i)-v(i,j-1)*(1.-dble(intf(i,j-1))))/2./dt) * dt

           rhsU(ii,2) = fvp(i,j)



!-------------------------- j = ny, i = 2,nx-1 ------------------------------------

        do i = 2,nx-1

              ii = ii + 1

              rho = rho1x(i,j) + rho2x(i,j)

              xmuc = (xmu(i,j)+xmu(i+1,j)*(1.-dble(intf(i+1,j)))+xmu(i,j)*dble(intf(i+1,j)))/2.d0
              xmub = (xmu(i,j)+xmu(i+1,j)*(1.-dble(intf(i+1,j)))+xmu(i,j)*dble(intf(i+1,j))&
                     +xmu(i,j-1)*(1.-dble(intf(i,j-1)))+xmu(i,j)*dble(intf(i,j-1))+(xmu(i+1,j-1)*&
                     (1.-dble(intf(i+1,j-1)))+xmu(i+1,j)*dble(intf(i+1,j-1)))*(1.-dble(intf(i+1,j)))+&
                     (xmu(i,j-1)*(1.-dble(intf(i,j-1)))+xmu(i,j)*dble(intf(i,j-1)))*dble(intf(i+1,j)))/4.d0
              xmul = (xmu(i,j)+xmu(i-1,j)*(1.-dble(intf(i-1,j)))+xmu(i,j)*dble(intf(i-1,j)))/2.d0

              xmuh(i+1,j+1) = xmuc
              xmuh(i+1,j) = xmub
              xmuh(i,j+1) = xmul

              coeflx(i,j,1)=dtx*xmu(i,j)*rho/re
              coefrx(i,j,1)=dtx*xmu(i+1,j)*rho/re
              coefly(i,j,1)=dty*xmuh(i+1,j)*rho/re
              coefry(i,j,1)=dty*xmuh(i+1,j+1)*rho/re

              fup(i,j) = ( u(i,j)/dt + g + &
                   (coefrx(i,j,1)-coeflx(i,j,1))*(u(i+1,j)*(1.-dble(intf(i+1,j)))-&
                   u(i-1,j)*(1.-dble(intf(i-1,j))))/2./dt + &
                   (coefry(i,j,1)-coefly(i,j,1))*((v(i+1,j)+v(i+1,j-1))*(1.-dble(intf(i+1,j)))-&
                   (v(i,j)+v(i,j-1)))/2./dt) * dt


              rhsU(ii,1) = fup(i,j)

              rho = rho1y(i,j) + rho2y(i,j)

              coeflx(i,j,2)=dtx*xmuh(i,j+1)*rho/re
              coefrx(i,j,2)=dtx*xmuh(i+1,j+1)*rho/re
              coefly(i,j,2)=dty*xmu(i,j)*rho/re
              coefry(i,j,2)=dty*xmu(i,j)*rho/re
 
              fvp(i,j) = ( v(i,j)/dt + &
                   (coefrx(i,j,2)-coeflx(i,j,2))*((ur(i)+ur(i-1))-(u(i,j)+u(i-1,j)))/2./dt + &
                   (coefry(i,j,2)-coefly(i,j,2))*(vr(i)-v(i,j-1)*(1.-dble(intf(i,j-1))))/2./dt) * dt

              rhsU(ii,2) = fvp(i,j)


        end do

!------------------ j = ny, i = nx -------------------------------------------------

        i = nx


           ii = ii + 1

           rho = rho1x(i,j) + rho2x(i,j)

           xmuc = xmu(i,j)
           xmub = (xmu(i,j)+xmu(i,j-1)*(1.-dble(intf(i,j-1)))+xmu(i,j)*dble(intf(i,j-1)))/2.d0
           xmul = (xmu(i,j)+xmu(i-1,j)*(1.-dble(intf(i-1,j)))+xmu(i,j)*dble(intf(i-1,j)))/2.d0

           xmuh(i+1,j+1) = xmuc
           xmuh(i+1,j) = xmub
           xmuh(i,j+1) = xmul

           coeflx(i,j,1)=dtx*xmu(i,j)*rho/re
           coefrx(i,j,1)=dtx*xmu(i,j)*rho/re
           coefly(i,j,1)=dty*xmuh(i+1,j)*rho/re
           coefry(i,j,1)=dty*xmuh(i+1,j+1)*rho/re

           fup(i,j) = ( u(i,j)/dt + g + &
                   (coefrx(i,j,1)-coeflx(i,j,1))*(ut(j)-u(i-1,j)*(1.-dble(intf(i-1,j))))/2./dt + &
                   (coefry(i,j,1)-coefly(i,j,1))*((vt(j)+vt(j-1))-(v(i,j)+v(i,j-1)))/2./dt) * dt

           rhsU(ii,1) = fup(i,j)

           rho = rho1y(i,j) + rho2y(i,j)

           coeflx(i,j,2)=dtx*xmuh(i,j+1)*rho/re
           coefrx(i,j,2)=dtx*xmuh(i+1,j+1)*rho/re
           coefly(i,j,2)=dty*xmu(i,j)*rho/re
           coefry(i,j,2)=dty*xmu(i,j)*rho/re

           fvp(i,j) = ( v(i,j)/dt + &
                   (coefrx(i,j,2)-coeflx(i,j,2))*((ur(i)+ur(i-1))-(u(i,j)+u(i-1,j)))/2./dt + &
                   (coefry(i,j,2)-coefly(i,j,2))*(vr(i)-v(i,j-1)*(1.-dble(intf(i,j-1))))/2./dt) * dt

           rhsU(ii,2) = fvp(i,j)


!*********************************************************************************

        sov_medi = 0.
        uttt = 0.
        ubbb = 0.
        vttt = 0.
        vlll = 0.
        vrrr = 0.

        ! Thomas Algorithm for intermit velocity u and v
        do k = 1,2
           if (k.eq.1) then
              nymax = ny
           else
              nymax = ny-1
           end if
           if (k.eq.1) then
              nxmax = nx-1
           else
              nxmax = nx
           end if
           utr = 0.
           !--------------solve the linear system in x direction--------------
           !build the main diagonal matrix for x axis
           do j = 1,nymax
              iicx = 0
              do i = 1,nxmax
                 if (k.eq.1) then
                    if (intf(i+1,j).eq.0.and.intf(i,j).eq.0) then
                       iicx = iicx + 1
                       if (i.eq.1) then
                          aux(iicx) = 0.
                          bux(iicx) = 1. + coeflx(i,j,k) + coefrx(i,j,k)
                          if (intf(i+2,j).eq.0) then
                             cux(iicx) = - coefrx(i,j,k)
                          else
                             cux(iicx) = 0.
                          end if
                          dux(iicx) = rhsU(nx*(j-1)+i,k)+coeflx(i,j,k)*uin(j)
                       else if (i.eq.nxmax) then
                          if (intf(i-1,j).eq.0) then
                             aux(iicx) = - coeflx(i,j,k)
                          else
                             aux(iicx) = 0.
                          end if
                          bux(iicx) = 1. + coeflx(i,j,k) + coefrx(i,j,k)
                          cux(iicx) = 0.
                          dux(iicx) = rhsU(nx*(j-1)+i,k)+coefrx(i,j,k)*uttt(j)
                       else
                          if (intf(i-1,j).eq.0) then
                             aux(iicx) = - coeflx(i,j,k)
                          else
                             aux(iicx) = 0.
                          end if
                          bux(iicx) = 1. + coeflx(i,j,k) + coefrx(i,j,k)
                          if (intf(i+2,j).eq.0) then
                             cux(iicx) = - coefrx(i,j,k)
                          else
                             cux(iicx) = 0.
                          end if
                          dux(iicx) = rhsU(nx*(j-1)+i,k)
                       end if
                    end if
                 else
                    if (intf(i,j+1).eq.0.and.intf(i,j).eq.0) then
                       iicx = iicx + 1
                       if (i.eq.1) then
                          if (ubc.eq.-1.) then
                             aux(iicx) = 0.
                             bux(iicx) = 1.  + coefrx(i,j,k) + coeflx(i,j,k) + coeflx(i,j,k)
                             dux(iicx) = rhsU(nx*(j-1)+i,k) + coeflx(i,j,k)*dt/dy*(rho1y(i,j)+rho2y(i,j))*(p(i,j+1)-p(i,j))*2.
                          else
                             aux(iicx) = 0.
                             bux(iicx) = 1.  + coefrx(i,j,k) + coeflx(i,j,k) - coeflx(i,j,k)
                             dux(iicx) = rhsU(nx*(j-1)+i,k)
                          end if
                          if (intf(i+1,j+1).eq.0.and.intf(i+1,j).eq.0) then
                             cux(iicx) = - coefrx(i,j,k)
                          else
                             cux(iicx) = 0.
                             bux(iicx) = bux(iicx) + coefrx(i,j,k)
                             dux(iicx) = dux(iicx) + coefrx(i,j,k)*dt/dy*(rho1y(i,j)+rho2y(i,j))*(p(i,j+1)-p(i,j))*2.
                          end if
                       else if (i.eq.nxmax) then
                          if (ubc.eq.-1) then
                             bux(iicx) = 1. + coeflx(i,j,k) + coefrx(i,j,k) + coefrx(i,j,k)
                             cux(iicx) = 0.
                             dux(iicx) = rhsU(nx*(j-1)+i,k) + coefrx(i,j,k)*dt/dy*(rho1y(i,j)+rho2y(i,j))*(p(i,j+1)-p(i,j))*2.
                          else
                             bux(iicx) = 1. + coeflx(i,j,k) + coefrx(i,j,k) - coefrx(i,j,k)
                             cux(iicx) = 0.
                             dux(iicx) = rhsU(nx*(j-1)+i,k)
                          end if
                          if (intf(i-1,j+1).eq.0.and.intf(i-1,j).eq.0) then
                             aux(iicx) = - coeflx(i,j,k)
                          else
                             aux(iicx) = 0.
                             bux(iicx) = bux(iicx) + coeflx(i,j,k)
                             dux(iicx) = dux(iicx) + coeflx(i,j,k)*dt/dy*(rho1y(i,j)+rho2y(i,j))*(p(i,j+1)-p(i,j))*2.
                          end if
                       else
                          bux(iicx) = 1. + coeflx(i,j,k) + coefrx(i,j,k)
                          dux(iicx) = rhsU(nx*(j-1)+i,k)
                          if (intf(i-1,j+1).eq.0.and.intf(i-1,j).eq.0) then
                             aux(iicx) = - coeflx(i,j,k)
                          else
                             aux(iicx) = 0.
                             bux(iicx) = bux(iicx) + coeflx(i,j,k)
                             dux(iicx) = dux(iicx) + coeflx(i,j,k)*dt/dy*(rho1y(i,j)+rho2y(i,j))*(p(i,j+1)-p(i,j))*2.
                          end if
                          if (intf(i+1,j+1).eq.0.and.intf(i+1,j).eq.0) then
                             cux(iicx) = - coefrx(i,j,k)
                          else
                             cux(iicx) = 0.
                             bux(iicx) = bux(iicx) + coefrx(i,j,k)
                             dux(iicx) = dux(iicx) + coefrx(i,j,k)*dt/dy*(rho1y(i,j)+rho2y(i,j))*(p(i,j+1)-p(i,j))*2.
                          end if

                       end if
                    end if
                 end if
              end do
           
             ! forward sweep
              cuxp(1) = cux(1)/(bux(1)+eps)
              duxp(1) = dux(1)/(bux(1)+eps)
              do i = 2,iicx
                 cuxp(i) = cux(i)/(bux(i)-aux(i)*cuxp(i-1)+eps)
                 duxp(i) = (dux(i)-aux(i)*duxp(i-1))/(bux(i)-aux(i)*cuxp(i-1)+eps)
              end do
        
             ! backward sweep
              um(iicx) = duxp(iicx)
              do i = iicx-1,1,-1
                 um(i) = duxp(i)-um(i+1)*cuxp(i)
              end do
  
              iicx = 0
              do i= 1,nxmax
                 if (k.eq.1) then
                    if (intf(i,j).eq.0.and.intf(i+1,j).eq.0) then
                        iicx = iicx + 1
                        utr(i,j) = um(iicx)
                    end if
                 else
                    if (intf(i,j).eq.0.and.intf(i,j+1).eq.0) then
                        iicx = iicx + 1
                        utr(i,j) = um(iicx)
                    end if
                 end if
              end do
           end do
           !------------------------------------------------------------------

           do i = 1,nxmax
              do j = 1,nymax
                 umm(nymax*(i-1)+j) = utr(i,j)
              end do
           end do

           !--------------solve the linear system in y direction--------------
           do i = 1,nxmax
              jjcy = 0
              do j = 1,nymax
                 if (k.eq.1) then
                    if (intf(i,j).eq.0.and.intf(i+1,j).eq.0) then
                       jjcy = jjcy + 1
                       if (j.eq.1) then
                          if (ubc.eq.-1) then
                             auy(jjcy) = 0.
                             buy(jjcy) = 1. + coefly(i,j,k) + coefry(i,j,k) + coefly(i,j,k)
                             duy(jjcy) = umm(nymax*(i-1)+j) + coefly(i,j,k)*dt/dx*&
                                         (rho1x(i,j)+rho2x(i,j))*(p(i+1,j)-p(i,j))*2.
                          else
                             auy(jjcy) = 0.
                             buy(jjcy) = 1. + coefly(i,j,k) + coefry(i,j,k) - coefly(i,j,k)
                             duy(jjcy) = umm(nymax*(i-1)+j)
                          end if
                          if (intf(i,j+1).eq.0.and.intf(i+1,j+1).eq.0) then
                             cuy(jjcy) = - coefry(i,j,k)
                          else
                             cuy(jjcy) = 0.
                             buy(jjcy) = buy(jjcy) + coefry(i,j,k)
                             duy(jjcy) = duy(jjcy) + coefry(i,j,k)*dt/dx*&
                                         (rho1x(i,j)+rho2x(i,j))*(p(i+1,j)-p(i,j))*2.
                          end if
                       elseif (j.eq.nymax) then
                          if (ubc.eq.-1) then
                             buy(jjcy) = 1. + coefly(i,j,k) + coefry(i,j,k) + coefry(i,j,k)
                             cuy(jjcy) = 0.
                             duy(jjcy) = umm(nymax*(i-1)+j) + coefry(i,j,k)*dt/dx*&
                                         (rho1x(i,j)+rho2x(i,j))*(p(i+1,j)-p(i,j))*2.
                          else
                             buy(jjcy) = 1. + coefly(i,j,k) + coefry(i,j,k) - coefry(i,j,k)
                             cuy(jjcy) = 0.
                             duy(jjcy) = umm(nymax*(i-1)+j)
                          end if
                          if (intf(i,j-1).eq.0.and.intf(i+1,j-1).eq.0) then
                             auy(jjcy) = - coefly(i,j,k)
                          else
                             auy(jjcy) = 0.
                             buy(jjcy) = buy(jjcy) + coefly(i,j,k)
                             duy(jjcy) = duy(jjcy) + coefly(i,j,k)*dt/dx*&
                                         (rho1x(i,j)+rho2x(i,j))*(p(i+1,j)-p(i,j))*2.

                          end if
                       else
                          buy(jjcy) = 1. + coefly(i,j,k) + coefry(i,j,k)
                          duy(jjcy) = umm(nymax*(i-1)+j)
                          if (intf(i,j-1).eq.0.and.intf(i+1,j-1).eq.0) then
                             auy(jjcy) = - coefly(i,j,k)
                          else
                             auy(jjcy) = 0.
                             buy(jjcy) = buy(jjcy) + coefly(i,j,k)
                             duy(jjcy) = duy(jjcy) + coefly(i,j,k)*dt/dx*&
                                         (rho1x(i,j)+rho2x(i,j))*(p(i+1,j)-p(i,j))*2.
                          end if
                          if (intf(i,j+1).eq.0.and.intf(i+1,j+1).eq.0) then
                             cuy(jjcy) = -coefry(i,j,k)
                          else
                             cuy(jjcy) = 0.
                             buy(jjcy) = buy(jjcy) + coefry(i,j,k)
                             duy(jjcy) = duy(jjcy) + coefry(i,j,k)*dt/dx*&
                                         (rho1x(i,j)+rho2x(i,j))*(p(i+1,j)-p(i,j))*2.
                          end if
                       end if
                    end if
                 !if (j.eq.11) then
                 !   write(*,*) 20, i, auy(jjcy),buy(jjcy),cuy(jjcy),duy(jjcy)
                 !end if
                 else
                    if (intf(i,j).eq.0.and.intf(i,j+1).eq.0) then
                       jjcy = jjcy + 1
                       if (j.eq.1) then
                          auy(jjcy) = 0.
                          buy(jjcy) = 1. + coefry(i,j,k) + coefly(i,j,k)
                          if (intf(i,j+2).eq.0) then
                             cuy(jjcy) = - coefry(i,j,k)
                          else
                             cuy(jjcy) = 0.
                          end if
                          duy(jjcy) = umm(nymax*(i-1)+j)
                       elseif (j.eq.nymax) then
                          if (intf(i,j-1).eq.0) then
                             auy(jjcy) = -coefly(i,j,k)
                          else
                             auy(jjcy) = 0.
                          end if
                          buy(jjcy) = 1.+ coefry(i,j,k) + coefly(i,j,k)
                          cuy(jjcy) = 0.
                          duy(jjcy) = umm(nymax*(i-1)+j)
                       else
                          if (intf(i,j-1).eq.0) then
                             auy(jjcy) = - coefly(i,j,k)
                          else
                             auy(jjcy) = 0.
                          end if
                          buy(jjcy) = 1. + coefly(i,j,k) + coefry(i,j,k)
                          if (intf(i,j+2).eq.0) then
                             cuy(jjcy) = -coefry(i,j,k)
                          else
                             cuy(jjcy) = 0.
                          end if
                          duy(jjcy) = umm(nymax*(i-1)+j)
                       end if
                    end if
                 end if
              end do
             ! forward sweep
              cuyp(1) = cuy(1)/(buy(1)+eps)
              duyp(1) = duy(1)/(buy(1)+eps)
              do j = 2,jjcy
                 cuyp(j) = cuy(j)/(buy(j)-auy(j)*cuyp(j-1)+eps)
                 duyp(j) = (duy(j)-auy(j)*duyp(j-1))/(buy(j)-auy(j)*cuyp(j-1)+eps)
              end do

             ! backward sweep
              ummm(jjcy) = duyp(jjcy)
              do j = jjcy-1,1,-1
                 ummm(j) = duyp(j)-ummm(j+1)*cuyp(j)
              end do

              jjcy = 0
              do j = 1,nymax
                 if (k.eq.1) then
                    if (intf(i,j).eq.0.and.intf(i+1,j).eq.0) then
                       jjcy = jjcy + 1
                       sov_medi(i,j,k) = ummm(jjcy)
                    end if
                 else
                    if (intf(i,j).eq.0.and.intf(i,j+1).eq.0) then
                       jjcy = jjcy + 1
                       sov_medi(i,j,k) = ummm(jjcy)
                    end if
                 end if
              end do
    
           end do
           !------------------------------------------------------------------
     
        end do

        !solve for pressure using LU factorization computed initially
        ii = 0
        j = 1
        i = 1
        if(intf(i,j).eq.0) then
           ii = ii + 1
           rhs(ii) = (sov_medi(i,j,1)* (1.d0 - dble(intf(i+1,j)))-ubbb(j))/dt/dx&
                    +(sov_medi(i,j,2)* (1.d0 - dble(intf(i,j+1)))-vlll(i))/dt/dy + w(i,j)
        end if
        do i = 2,nx-1
           if(intf(i,j).eq.0) then
              ii = ii + 1
              rhs(ii) = (sov_medi(i,j,1)* (1.d0 - dble(intf(i+1,j)))-sov_medi(i-1,j,1)* (1.d0 - dble(intf(i-1,j))))/dt/dx&
                        +(sov_medi(i,j,2)* (1.d0 - dble(intf(i,j+1)))-vlll(i))/dt/dy + w(i,j)
           end if
        end do
        i = nx
        if(intf(i,j).eq.0) then
           ii = ii + 1
           rhs(ii) = (uttt(j)-sov_medi(i-1,j,1)* (1.d0 - dble(intf(i-1,j))))/dt/dx&
                       +(sov_medi(i,j,2)* (1.d0 - dble(intf(i,j+1)))-vlll(i))/dt/dy + w(i,j)
        end if
        do j = 2,ny-1
           i = 1
           if(intf(i,j).eq.0) then
             ii = ii + 1
             rhs(ii) = (sov_medi(i,j,1)* (1.d0 - dble(intf(i+1,j)))-ubbb(j))/dt/dx&
                         +(sov_medi(i,j,2)* (1.d0 - dble(intf(i,j+1)))-sov_medi(i,j-1,2)* (1.d0 - dble(intf(i,j-1))))/dt/dy &
                         + w(i,j)
           end if
           do i = 2,nx-1
              if(intf(i,j).eq.0) then
                ii = ii + 1
                rhs(ii) = (sov_medi(i,j,1)* (1.d0 - dble(intf(i+1,j)))-sov_medi(i-1,j,1)* (1.d0 - dble(intf(i-1,j))))/dt/dx&
                           +(sov_medi(i,j,2)* (1.d0 - dble(intf(i,j+1)))-sov_medi(i,j-1,2)* (1.d0 - dble(intf(i,j-1))))/dt/dy&
                           + w(i,j)
              end if
           end do
           i = nx
           if(intf(i,j).eq.0) then
              ii = ii + 1
              rhs(ii) = (uttt(j)-sov_medi(i-1,j,1)* (1.d0 - dble(intf(i-1,j))))/dt/dx+&
                        (sov_medi(i,j,2)* (1.d0 - dble(intf(i,j+1)))-sov_medi(i,j-1,2)* (1.d0 - dble(intf(i,j-1))))/dt/dy&
                        + w(i,j)
           end if
        end do
        j = ny
        i = 1
        if(intf(i,j).eq.0) then
           ii = ii + 1
           rhs(ii) = (sov_medi(i,j,1)* (1.d0 - dble(intf(i+1,j)))-ubbb(j))/dt/dx&
                     + (vrrr(i)-sov_medi(i,j-1,2)* (1.d0 - dble(intf(i,j-1))))/dt/dy &
                     + w(i,j)
        end if
        do i = 2,nx-1
           if(intf(i,j).eq.0) then
              ii = ii + 1
              rhs(ii) = (sov_medi(i,j,1)* (1.d0 - dble(intf(i+1,j)))-sov_medi(i-1,j,1)* (1.d0 - dble(intf(i-1,j))))/dt/dx&
                      +(vrrr(i)-sov_medi(i,j-1,2)* (1.d0 - dble(intf(i,j-1))))/dt/dy&
                      + w(i,j)
           end if
        end do
        i = nx
        if(intf(i,j).eq.0) then
           ii = ii + 1
           rhs(ii) = (uttt(j)-sov_medi(i-1,j,1)* (1.d0 - dble(intf(i-1,j))))/dt/dx+&
                     (vrrr(i)-sov_medi(i,j-1,2)* (1.d0 - dble(intf(i,j-1))))/dt/dy &
                     + w(i,j)
        end if

        if (k_Stoke.eq.1) then
           call umf4def(control)
           control(1) = 1
           call umf4pcon(control)
           call LU(rho1xs,rho1ys,rho2xs,rho2ys,dx,dy,intf,control,nx,ny,mu,inumeric,isymbolic)
        end if

        isys = 0
        sol = 0.d0
        call umf4sol(isys, sol, rhs, inumeric, control, xinfo)
        if (xinfo (1) .lt. 0) then
           print *, 'Error occurred in umf4sol: ', xinfo (1)
           stop
        end if

        ii = 0
        do j = 1,ny
           do i = 1,nx
              if(intf(i,j).eq.0) then
                 ii = ii + 1
                 sov_medi(i,j,3) = sol(ii)
              end if
           end do
        end do

!-------------check residual and compute velocity---------------------------

        eps = 1E-15
        err = 0.
        div = 0.
        res = 0.
        u = 0.
        v = 0.
        up = 0.
        vp = 0.
        ii = 0
        xm1 = 0.
        xm2 = 0.
        do j = 1,ny
           do i = 1,nx
              if(intf(i,j).eq.0) then

                 ii = ii + 1

                 dpxr = 0.
                 dpxl = 0.
                 dpyr = 0.
                 dpyl = 0.
                 up1r = 0.
                 up1l = 0.
                 vp1r = 0.
                 vp1l = 0.
                 up2r = 0.
                 up2l = 0.
                 vp2r = 0.
                 vp2l = 0.
                 ur1 = 0.
                 ul1 = 0.
                 vr1 = 0.
                 vl1 = 0.

                 !--veolcity field--------------------------------

                 if(i.eq.1) then
                    dpxl = (rho1xs(i,j) + rho2xs(i,j)) * 0.
                 else
                    dpxl = (rho1xs(i-1,j) + rho2xs(i-1,j)) &
                           * (sov_medi(i,j,3) - sov_medi(i-1,j,3)) * &
                           (1.d0 - dble(intf(i-1,j)))
                 end if
                 if(i.eq.nx) then
                    dpxr = - (rho1xs(i,j) + rho2xs(i,j)) * sov_medi(i,j,3) * uout
                 else
                    dpxr = (rho1xs(i,j) + rho2xs(i,j)) &
                         * (sov_medi(i+1,j,3) - sov_medi(i,j,3)) * &
                           (1.d0 - dble(intf(i+1,j)))
                 end if

                 if(j.eq.1) then
                    dpyl = (rho1ys(i,j) + rho2ys(i,j)) * sov_medi(i,j,3) * 0.
                 else
                    dpyl = (rho1ys(i,j-1) + rho2ys(i,j-1)) &
                         * (sov_medi(i,j,3) - sov_medi(i,j-1,3)) * &
                           (1.d0 - dble(intf(i,j-1)))
                 end if
                 if(j.eq.ny) then
                    dpyr = 0.
                 else
                    dpyr = (rho1ys(i,j) + rho2ys(i,j)) &
                         * (sov_medi(i,j+1,3) - sov_medi(i,j,3)) * &
                           (1.d0 - dble(intf(i,j+1)))
                 end if

                 xx = (dpxr - dpxl)/dx**2. + (dpyr - dpyl)/dy**2.

                 if (i.eq.nx) then
                    u(i,j) = uttt(j) - dpxr/dx * dt + sig(i,j,1) * dt
                 else
                    u(i,j) = sov_medi(i,j,1) - dpxr/dx * dt + sig(i,j,1) * dt
                 end if
                 if (j.eq.ny) then
                    v(i,j) = 0. + sig(i,j,2) * dt
                 else
                    v(i,j) = sov_medi(i,j,2) - dpyr/dy * dt + sig(i,j,2) * dt
                 end if

                 p(i,j) = sov_medi(i,j,3)

                 !-phase 1 velocities----------------------------

                 if (i.lt.nx) then
                    ur1 = u(i,j)
                 else
                    ur1 = u(i,j)
                 end if

                 if(i.gt.1) then
                    ul1 = u(i-1,j)
                 else
                    ul1 = ubbb(j)
                 end if

                 if (j.lt.ny) then
                    vr1 = v(i,j)
                 else
                    vr1 = vrrr(i)
                 end if

                 if(j.gt.1) then
                    vl1 = v(i,j-1)
                 else
                    vl1 = vlll(i)
                 end if


                 !----------residual, divergence --------------------

                 div1 = abs((ur1-ul1)/dx + (vr1-vl1)/dy)
                 div = div + div1 / dble(nx*ny)
                 res = max(res,abs(xx-rhs(ii)))

                 if(res.gt.1E-2) then
                    write(6,'(a20,3i5,e20.10)') "presure residual not zero ",i,j,ii,res
                    write(*,*) xx, rhs(ii), dpyr,dpxr,dpyl,dpxl
                    stop
                 end if
                 if(div1.gt.1E-5) then
                    write(6,'(a20,3i5,e20.10)') "divergence not zero ",i,j,ii,div1
                    write(*,*) u(i,j)
                    stop
                 end if

              end if
           end do
          
        end do

        errm = abs(xm1-xm2)/dble(nx)/dble(ny)


10      format(2i5,10e20.10)
20      format(3i5,10e20.10)
30      format(3i5,15f13.5)


      end subroutine solve1

!-----------------------------------------------------------------------------------
!-----------------------------------------------------------------------------------
!-----------------------------------------------------------------------------------

                     
      subroutine LU(rho1x,rho1y,rho2x,rho2y,dx,dy,intf,control,nx,ny,mu,inumeric,isymbolic)
        implicit real*8(a-h,o-z)
        common/param/g,sp,ubc,uout,mint
        dimension :: intf(nx,ny),control(20),xinfo(90), &
                     rho1x(nx,ny),rho1y(nx,ny),rho2x(nx,ny),rho2y(nx,ny)
        allocatable :: cmat(:,:),iAp(:),iAp1(:),iAi(:,:),iAi1(:),Ax(:), &
                       sol(:),cmat1(:,:)
        integer*8 inumeric,isymbolic,istatus,isys
        
        nm = nx*ny - mint
        dx2 = 1./dx**2.
        dy2 = 1./dy**2.
        eps = 10E-8

        mumf = 9*nm
!        allocate(cmat(nm,18),iAp(nm+1),iAi(nm,18),iAi1(mumf),sol(nm),Ax(mumf), &
!                 cmat1(nm,nm))
        allocate(cmat(nm,18),iAp(nm+1),iAi(nm,18),iAi1(mumf),sol(nm),Ax(mumf))

        iAp = 0
        iAi1 = 0
        iAi = 0
        cmat = 0.
        ii = 0

!*********************************************************************************

        j = 1
        i = 1

        if(intf(i,j).eq.0) then

           ii = ii + 1

           fl = (rho1x(i,j) + rho2x(i,j))*dx2
           fr = (rho1x(i,j) + rho2x(i,j))*dx2

           call cmat_cfx(cmat,fr,fl,iAp,iAi,intf,dx,nm,nx,ny,ii,i,j, &
                0,0,1,1,1)

           fl = (rho1y(i,j) + rho2y(i,j))*dy2
           fr = (rho1y(i,j) + rho2y(i,j))*dy2

           call cmat_cfy(cmat,fr,fl,iAp,iAi,intf,dy,nm,nx,ny,ii,i,j, &
                0,0,1,1,1)

        end if
!---------------------------------j=1---------------------------------------------

        do i = 2,nx-1

           if(intf(i,j).eq.0) then

              ii = ii + 1

              fl = (rho1x(i-1,j) + rho2x(i-1,j))*dx2
              fr = (rho1x(i,j) + rho2x(i,j))*dx2

              k1=3+(1-3/i)+(1-4/i); k2=2+(1-2/i)+(1-3/i)
              k3=2+(1-2/i); k4=2; k5=1
              call cmat_cfx(cmat,fr,fl,iAp,iAi,intf,dx,nm,nx,ny,ii,i,j, &
                   k1,k2,k3,k4,k5)!,nchoose)

              fl = (rho1y(i,j) + rho2y(i,j))*dy2
              fr = (rho1y(i,j) + rho2y(i,j))*dy2

              k3=2+(1-2/i)
              call cmat_cfy(cmat,fr,fl,iAp,iAi,intf,dy,nm,nx,ny,ii,i,j, &
                   0,0,k3,1,1)!,nchoose)

           end if

        end do

!---------------------------------j=1---------------------------------------------

        i = nx

        if(intf(i,j).eq.0) then

           ii = ii + 1

           fl = (rho1x(i-1,j) + rho2x(i-1,j))*dx2
           fr = (rho1x(i,j) + rho2x(i,j))*dx2

           call cmat_cfx(cmat,fr,fl,iAp,iAi,intf,dx,nm,nx,ny,ii,i,j, &
                5,4,3,0,0)


           fl = (rho1y(i,j) + rho2y(i,j))*dy2
           fr = (rho1y(i,j) + rho2y(i,j))*dy2

           call cmat_cfy(cmat,fr,fl,iAp,iAi,intf,dy,nm,nx,ny,ii,i,j, &
                0,0,3,1,1)

        end if

!*********************************************************************************

        do j = 2,ny-1

           jj = 1 - 2/j
           jj1 = (1 - 2/j) + (1 - 3/j)
           jj2 = (1 - 3/j) + (1 - 4/j)

           i = 1

           if(intf(i,j).eq.0) then

              ii = ii + 1

              fl = (rho1x(i,j) + rho2x(i,j))*dx2
              fr = (rho1x(i,j) + rho2x(i,j))*dx2

              k1=0; k2=0; k3=2+jj; k4=2+jj; k5=2+jj
              call cmat_cfx(cmat,fr,fl,iAp,iAi,intf,dx,nm,nx,ny,ii,i,j, &
                            k1,k2,k3,k4,k5)
   
              fl = (rho1y(i,j-1) + rho2y(i,j-1))*dy2
              fr = (rho1y(i,j) + rho2y(i,j))*dy2

              k1=5+jj2; k2=4+jj1; k3=2+jj; k4=2; k5=1
              call cmat_cfy(cmat,fr,fl,iAp,iAi,intf,dy,nm,nx,ny,ii,i,j, &
                            k1,k2,k3,k4,k5)

           end if

!---------------------------------j=2,ny-1----------------------------------------

           do i = 2,nx-1
              
              if(intf(i,j).eq.0) then

                 ii = ii + 1

                 if (intf(i-1,j).eq.0) then
                    fl = (rho1x(i-1,j) + rho2x(i-1,j))*dx2
                 else
                    fl = (rho1x(i,j) + rho2x(i,j))*dx2
                 end if
                 fr = (rho1x(i,j) + rho2x(i,j))*dx2
                 
                 k1=4+(1-3/i)+(1-4/i)+jj; k2=3+(1-2/i)+(1-3/i)+jj
                 k3=3+(1-2/i)+jj; k4=3+jj; k5=2+jj
                 call cmat_cfx(cmat,fr,fl,iAp,iAi,intf,dx,nm,nx,ny,ii,i,j, &
                      k1,k2,k3,k4,k5)
   
                 if (intf(i,j-1).eq.0) then
                    fl = (rho1y(i,j-1) + rho2y(i,j-1))*dy2
                 else
                    fl = (rho1y(i,j) + rho2y(i,j))*dy2
                 end if
                 fr = (rho1y(i,j) + rho2y(i,j))*dy2

                 k1=7-2/i-i/(nx-1)+jj2; k2=6-2/i-i/(nx-1)+jj1
                 k3=3+(1-2/i)+jj; k4=2; k5=1
                 call cmat_cfy(cmat,fr,fl,iAp,iAi,intf,dy,nm,nx,ny,ii,i,j, &
                               k1,k2,k3,k4,k5)

              end if

           end do

!---------------------------------j=2,ny-1----------------------------------------

           i = nx

           if(intf(i,j).eq.0) then

              ii = ii + 1

              fl = (rho1x(i-1,j) + rho2x(i-1,j))*dx2
              fr = (rho1x(i,j) + rho2x(i,j))*dx2

              k1=6+jj; k2=5+jj; k3=4+jj; k4=0; k5=0
              call cmat_cfx(cmat,fr,fl,iAp,iAi,intf,dx,nm,nx,ny,ii,i,j, &
                   k1,k2,k3,k4,k5)
   
              fl = (rho1y(i,j-1) + rho2y(i,j-1))*dy2
              fr = (rho1y(i,j) + rho2y(i,j))*dy2

              k1=5+jj2; k2=4+jj1; k3=4+jj; k4=2; k5=1
              call cmat_cfy(cmat,fr,fl,iAp,iAi,intf,dy,nm,nx,ny,ii,i,j, &
                            k1,k2,k3,k4,k5)

           end if

        end do



!*********************************************************************************

        j = ny
        i = 1

        if(intf(i,j).eq.0) then

           ii = ii + 1

           fl = (rho1x(i,j) + rho2x(i,j))*dx2
           fr = (rho1x(i,j) + rho2x(i,j))*dx2

           k1=0; k2=0; k3=3; k4=3; k5=3
           call cmat_cfx(cmat,fr,fl,iAp,iAi,intf,dx,nm,nx,ny,ii,i,j, &
                k1,k2,k3,k4,k5)
   
           fl = (rho1y(i,j-1) + rho2y(i,j-1))*dy2
           fr = (rho1y(i,j) + rho2y(i,j))*dy2

           k1=7; k2=6; k3=3; k4=0; k5=0
           call cmat_cfy(cmat,fr,fl,iAp,iAi,intf,dy,nm,nx,ny,ii,i,j, &
                k1,k2,k3,k4,k5)

        end if

!---------------------------------j=ny------------------------------------------ &
                                        

        do i = 2,nx-1

           if(intf(i,j).eq.0) then

              ii = ii + 1

              fl = (rho1x(i-1,j) + rho2x(i-1,j))*dx2
              fr = (rho1x(i,j) + rho2x(i,j))*dx2
           
              k1=5+(1-3/i)+(1-4/i); k2=4+(1-2/i)+(1-3/i); k3=4+(1-2/i);
              k4=4; k5=3
              call cmat_cfx(cmat,fr,fl,iAp,iAi,intf,dx,nm,nx,ny,ii,i,j, &
                   k1,k2,k3,k4,k5)
    
              fl = (rho1y(i,j-1) + rho2y(i,j-1))*dy2
              fr = (rho1y(i,j) + rho2y(i,j))*dy2

              k1=9-2/i-i/(nx-1); k2=8-2/i-i/(nx-1); k3=4+(1-2/i); k4=0; k5=0
              call cmat_cfy(cmat,fr,fl,iAp,iAi,intf,dy,nm,nx,ny,ii,i,j, &
                            k1,k2,k3,k4,k5)

           end if

        end do

!---------------------------------j=ny--------------------------------------------


        i = nx

        if(intf(i,j).eq.0) then

           ii = ii + 1

           fl = (rho1x(i-1,j) + rho2x(i-1,j))*dx2
           fr = (rho1x(i,j) + rho2x(i,j))*dx2

           k1=7; k2=6; k3=5; k4=0; k5=0
           call cmat_cfx(cmat,fr,fl,iAp,iAi,intf,dx,nm,nx,ny,ii,i,j, &
                k1,k2,k3,k4,k5)
  
           fl = (rho1y(i,j-1) + rho2y(i,j-1))*dy2
           fr = (rho1y(i,j) + rho2y(i,j))*dy2

           k1=7; k2=6; k3=5; k4=0; k5=0
           call cmat_cfy(cmat,fr,fl,iAp,iAi,intf,dy,nm,nx,ny,ii,i,j, &
                k1,k2,k3,k4,k5)

        end if

!---------------------------------------------------------------------------------


        iAp = 0
        mu = 0
        do i = 1,nm
           do j = 1,18
              if(cmat(i,j).ne.0.) then
                 mu = mu + 1
                 Ax(mu) = cmat(i,j)
                 iAi1(mu) = iAi(i,j) - 1
              end if
           end do
           iAp(i+1) = mu
        end do

!!$        open(10,file="tmp.dat")
!!$        do i = 1,mu
!!$           write(10,'(i5,5E15.5)') i,cmat(i,1:5)
!!$        end do
!!$        close(10)



        call umf4sym(nm, nm, iAp, iAi1, Ax, isymbolic, control, xinfo)
        if (xinfo(1) .lt. 0) then
           print *, 'Error occurred in umf4sys: ', xinfo (1)
           stop
        endif
        call umf4num(iAp, iAi1, Ax, isymbolic, inumeric, control, xinfo)
        if (xinfo(1) .lt. 0) then
           print *, 'Error occurred in umf4num: ', xinfo (1)
           stop
        endif

10      format(2i5,10e20.10)

      end subroutine LU



!------------------------------------------------------------------------
!------------------------------------------------------------------------
!------------------------------------------------------------------------

      subroutine cmat_cfx(cmat,fr,fl,iAp,iAi,intf,dx,nm,nx,ny,ii,i,j, &
                          k1,k2,k3,k4,k5)
        implicit real*8(a-h,o-z)
        common/param/g,sp,ubc,uout,mint
        dimension :: cmat(nm,18),iAp(nm+1),iAi(nm,18),intf(nx,ny),xnj(nx,ny)

        if(i.eq.1) then
           x1 = 0.d0
        else
           x1 = 1.d0 - dble(intf(i-1,j))
        end if

        if(i.eq.nx) then
           xn = uout
        else
           xn = 1.d0 - dble(intf(i+1,j))
        end if

        if(i.gt.1) then

           if(intf(i-1,j).eq.0) then
              iAp(ii-1) = iAp(ii-1) + 1
              iAi(ii-1,k2) = ii
              cmat(ii-1,k2)  =  x1*fl
           end if

        end if

        iAp(ii) = iAp(ii) + 1
        iAi(ii,k3) = ii
        cmat(ii,k3)    = - (xn*fr + x1*fl)

        if(i.lt.nx) then

           if(intf(i+1,j).eq.0) then
              iAp(ii+1) = iAp(ii+1) + 1
              iAi(ii+1,k4) = ii
              cmat(ii+1,k4)  = xn*fr
           end if

        end if

      end subroutine cmat_cfx


!------------------------------------------------------------------------
!------------------------------------------------------------------------
!------------------------------------------------------------------------

      subroutine cmat_cfy(cmat,fr,fl,iAp,iAi,intf,dx,nm,nx,ny,ii,i,j, &
                          k1,k2,k3,k4,k5)
        implicit real*8(a-h,o-z)
        common/param/g,sp,ubc,uout,mint
        dimension :: cmat(nm,10),iAp(nm+1),iAi(nm,10),intf(nx,ny)

        if(j.eq.1) then
           x1 = 0.d0
        else
           x1 = 1.d0 - dble(intf(i,j-1))
        end if

        if(j.eq.ny) then
           xn = 0.d0
        else
           xn = 1.d0 - dble(intf(i,j+1))
        end if

        if(j.gt.1) then

           kk = 0
           do k = i-1,1,-1
              if(intf(k,j).eq.1) then
                 kk = kk + 1
              end if
           end do
           do k = nx,i+1,-1
              if(intf(k,j-1).eq.1) then
                 kk = kk + 1
              end if
           end do

           if(intf(i,j-1).eq.0) then
              iAp(ii-nx+kk) = iAp(ii-nx+kk) + 1
              iAi(ii-nx+kk,k2) = ii
              cmat(ii-nx+kk,k2)  = x1*fl
           end if

        end if

        iAp(ii) = iAp(ii) + 1

        cmat(ii,k3) = cmat(ii,k3) - (xn*fr + x1*fl)

        if(j.lt.ny) then

           kk = 0
           do k = i+1,nx
              if(intf(k,j).eq.1) then
                 kk = kk + 1
              end if
           end do
           do k = 1,i-1
              if(intf(k,j+1).eq.1) then
                 kk = kk + 1
              end if
           end do

           if(intf(i,j+1).eq.0) then
              iAp(ii+nx-kk) = iAp(ii+nx-kk) + 1
              iAi(ii+nx-kk,k4) = ii
              cmat(ii+nx-kk,k4)  =  xn*fr
           end if

        end if

 
      end subroutine cmat_cfy

!------------------------------------------------------------------------
!------------------------------------------------------------------------
!------------------------------------------------------------------------

   
      subroutine advect_WENO3(s,u,v,x,y,cfl,dt,smx,dx,dy,xmn,ymn,nx,ny,uin,s1,intf,phiw)
        implicit real*8(a-h,o-z)
        dimension :: s(nx,ny),u(nx,ny),v(nx,ny),x(nx),y(ny),uin(ny+1),s1(ny),intf(nx,ny),phiw(nx,ny)
        allocatable :: so(:,:),cm(:,:)
        allocate(so(nx,ny),cm(3,3))

        eps = 1E-15
        smx = 0.

        so = s

        cm = reshape( &
             (/ 1./3. , -1./6. ,   1./3.  , &
                5./6. ,  5./6. ,  -7./6.  , &
               -1./6. ,  1./3. ,  11./6.  /) &
             , (/3,3/) )

        d2 = 1./10.
        d1 = 3./5.
        d0 = 3./10.

        do j = 1,3
           do i = 1,nx
              if(i.eq.1) then
                 sxl = s1(j)
                 ul = uin(j)
              else
                 sxl = so(i-1,j)
                 ul = u(i-1,j)
              end if
              sxr = so(i+1,j)
              ur = u(i,j)

              if(j.eq.1) then
                 syl = so(i,j+1)
                 vl = 0.
              else
                 syl = so(i,j-1)
                 vl = v(i,j-1)
              end if
              if(j.eq.ny) then
                 syr = so(i,ny-1)
              else
                 syr = so(i,j+1)
              end if
              vr = v(i,j)

              sm = so(i,j)

              if(ur.ge.0.) then
                 fxr = sm*ur
              else
                 fxr = sxr*ur
              end if

              if(ul.ge.0.) then
                 fxl = sxl*ul
              else
                 fxl = sm*ul
              end if

              if(vr.ge.0.) then
                 fyr = sm*vr
              else
                 fyr = syr*vr
              end if

              if(vl.ge.0.) then
                 fyl = syl*vl
              else
                 fyl = sm*vl
              end if

              s(i,j) = so(i,j) - dt * ( (fxr-fxl)/dx + (fyr-fyl)/dy )
           end do
        end do

        do j = ny-2,ny
           do i = 1,nx
              if(i.eq.1) then
                 sxl = s1(j)
                 ul = uin(j)
              else
                 sxl = so(i-1,j)
                 ul = u(i-1,j)
              end if
              sxr = so(i+1,j)
              ur = u(i,j)

              if(j.eq.1) then
                 syl = so(i,j+1)
                 vl = 0.
              else
                 syl = so(i,j-1)
                 vl = v(i,j-1)
              end if
              if(j.eq.ny) then
                 syr = so(i,ny-1)
              else
                 syr = so(i,j+1)
              end if
              vr = v(i,j)

              sm = so(i,j)

              if(ur.ge.0.) then
                 fxr = sm*ur
              else
                 fxr = sxr*ur
              end if

              if(ul.ge.0.) then
                 fxl = sxl*ul
              else
                 fxl = sm*ul
              end if

              if(vr.ge.0.) then
                 fyr = sm*vr
              else
                 fyr = syr*vr
              end if

              if(vl.ge.0.) then
                 fyl = syl*vl
              else
                 fyl = sm*vl
              end if

              s(i,j) = so(i,j) - dt * ( (fxr-fxl)/dx + (fyr-fyl)/dy )
           end do
        end do

        do j = 4,ny-3

           if (phiw(i,j-3).lt.0.or.phiw(i,j+3).lt.0) then
              if(i.eq.1) then
                 sxl = s1(j)
                 ul = uin(j)
              else
                 sxl = so(i-1,j)
                 ul = u(i-1,j)
              end if
              sxr = so(i+1,j)
              ur = u(i,j)

              if(j.eq.1) then
                 syl = so(i,j+1)
                 vl = 0.
              else
                 syl = so(i,j-1)
                 vl = v(i,j-1)
              end if
              if(j.eq.ny) then
                 syr = so(i,ny-1)
              else
                 syr = so(i,j+1)
              end if
              vr = v(i,j)

              sm = so(i,j)

              if(ur.ge.0.) then
                 fxr = sm*ur
              else
                 fxr = sxr*ur
              end if

              if(ul.ge.0.) then
                 fxl = sxl*ul
              else
                 fxl = sm*ul
              end if

              if(vr.ge.0.) then
                 fyr = sm*vr
              else
                 fyr = syr*vr
              end if

              if(vl.ge.0.) then
                 fyl = syl*vl
              else
                 fyl = sm*vl
              end if

              s(i,j) = so(i,j) - dt * ( (fxr-fxl)/dx + (fyr-fyl)/dy )

           else

           do i = 1,3

              if(i.eq.1) then
                 sxl = s1(j)
                 ul = uin(j)
              else
                 sxl = so(i-1,j)
                 ul = u(i-1,j)
              end if
              sxr = so(i+1,j)
              ur = u(i,j)

              if(j.eq.1) then
                 syl = so(i,j+1)
                 vl = 0.
              else
                 syl = so(i,j-1)
                 vl = v(i,j-1)
              end if
              if(j.eq.ny) then
                 syr = so(i,ny-1)
              else
                 syr = so(i,j+1)
              end if
              vr = v(i,j)

              sm = so(i,j)

              if(ur.ge.0.) then
                 fxr = sm*ur
              else
                 fxr = sxr*ur
              end if

              if(ul.ge.0.) then
                 fxl = sxl*ul
              else
                 fxl = sm*ul
              end if

              if(vr.ge.0.) then
                 fyr = sm*vr
              else
                 fyr = syr*vr
              end if

              if(vl.ge.0.) then
                 fyl = syl*vl
              else
                 fyl = sm*vl
              end if

              s(i,j) = so(i,j) - dt * ( (fxr-fxl)/dx + (fyr-fyl)/dy )

           end do

           do i = 4,nx-3

           if (phiw(i-3,j).lt.0.or.phiw(i+3,j).lt.0) then
              if(i.eq.1) then
                 sxl = s1(j)
                 ul = uin(j)
              else
                 sxl = so(i-1,j)
                 ul = u(i-1,j)
              end if
              sxr = so(i+1,j)
              ur = u(i,j)

              if(j.eq.1) then
                 syl = so(i,j+1)
                 vl = 0.
              else
                 syl = so(i,j-1)
                 vl = v(i,j-1)
              end if
              if(j.eq.ny) then
                 syr = so(i,ny-1)
              else
                 syr = so(i,j+1)
              end if
              vr = v(i,j)

              sm = so(i,j)

              if(ur.ge.0.) then
                 fxr = sm*ur
              else
                 fxr = sxr*ur
              end if

              if(ul.ge.0.) then
                 fxl = sxl*ul
              else
                 fxl = sm*ul
              end if

              if(vr.ge.0.) then
                 fyr = sm*vr
              else
                 fyr = syr*vr
              end if

              if(vl.ge.0.) then
                 fyl = syl*vl
              else
                 fyl = sm*vl
              end if

              s(i,j) = so(i,j) - dt * ( (fxr-fxl)/dx + (fyr-fyl)/dy )

           else

              !if (intf(i,j).eq.0) then

              if(i.eq.1) then
                 ul = uin(j)
              else
                 ul = u(i-1,j)
              end if
              if(i.eq.nx) then
                 ur = u(i,j)
              else
                 ur = u(i,j)
              end if

              if(j.eq.1) then
                 vl = 0
              else
                 vl = v(i,j-1)
              end if
              if(j.eq.ny) then
                 vr = 0.
              else
                 vr = v(i,j)
              end if

              !x flux----------------------


              n = nx

              isgn = int((u(i,j)+eps)/abs(u(i,j)+eps))
           

              m1r = (1+isgn)/2 * ((1/i)*n + (2/i)*(i/2)*(n+1) + &
                   i*(i/max(i,3))) + &
                   (1-isgn)/2 * (i/n - ((n-2)/i)*(i/(n-2)) + &
                   i*((n-3)/max(i,n-3)))
           
              m2r = (1+isgn)/2 * ((1/i)*n + i*(i/max(i,2))) + &
                   (1-isgn)/2*(i/n + i*((n-2)/max(i,n-2)))
              
              m3r = i*(1+isgn)/2 + (1-isgn)/2*(i/n + i*((n-1)/max(i,n-1)))
              
              m4r = i*(1-isgn)/2 + (1+isgn)/2*(i/n + i*((n-1)/max(i,n-1)))

              m5r = (1-isgn)/2 * ((1/i)*n + i*(i/max(i,2))) + &
                   (1+isgn)/2 * (i/n + i*((n-2)/max(i,n-2)))

              m1l = (1+isgn)/2 * ((1/i)*n + (2/i)*(i/2)*(n+1) + &
                   (3/i)*(i/3)*(n+2) + i*(i/max(i,4))) + &
                    (1-isgn)/2 * (i/n + i*((n-2)/max(i,n-2)))

              m2l = (1+isgn)/2 * ((1/i)*n + (2/i)*(i/2)*(n+1) + &
                   i*(i/max(i,3))) + (1-isgn)/2 * (i/n + i*((n-1)/max(i,n-1)))

              m3l = (1+isgn)/2 * ((1/i)*n + i*(i/max(i,2))) + i*(1-isgn)/2

              m4l = (1-isgn)/2 * ((1/i)*n + i*(i/max(i,2))) + i*(1+isgn)/2

              m5l = (1-isgn)/2 * ((1/i)*n + (2/i)*(i/2)*(n+1) + &
                   i*(i/max(i,3))) + (1+isgn)/2 * (i/n + i*((n-1)/max(i,n-1)))

              i1r = m1r - 2 + 5*(1-isgn)/2
              i2r = m2r - 1 + 3*(1-isgn)/2
              i3r = m3r     + 1*(1-isgn)/2
              i4r = m4r + 1 - 1*(1-isgn)/2
              i5r = m5r + 2 - 3*(1-isgn)/2
              
              i1l = m1l - 3 + 5*(1-isgn)/2
              i2l = m2l - 2 + 3*(1-isgn)/2
              i3l = m3l - 1 + 1*(1-isgn)/2
              i4l = m4l     - 1*(1-isgn)/2
              i5l = m5l + 1 - 3*(1-isgn)/2

              exr = 1.
              exl = 1.

              f1r = so(i1r,j)**exr
              f2r = so(i2r,j)**exr
              f3r = so(i3r,j)**exr
              f4r = so(i4r,j)**exr
              f5r = so(i5r,j)**exr
              
              f1l = so(i1l,j)**exl
              f2l = so(i2l,j)**exl
              f3l = so(i3l,j)**exl
              f4l = so(i4l,j)**exl
              f5l = so(i5l,j)**exl

              !WENO3 smoothness functions in x-coordinate

              b0r = 13./12.*(f1r - 2.*f2r + f3r)**2. + &
                   1./4.*(f1r - 4.*f2r + 3.*f3r)**2.
              b1r = 13./12.*(f2r - 2.*f3r + f4r)**2. + 1./4.*(f2r - f4r)**2.
              b2r = 13./12.*(f3r - 2.*f4r + f5r)**2. + &
                   1./4.*(3.*f3r - 4.*f4r + f5r)**2.

              b0l = 13./12.*(f1l - 2.*f2l + f3l)**2. + &
                   1./4.*(f1l - 4.*f2l + 3.*f3l)**2.
              b1l = 13./12.*(f2l - 2.*f3l + f4l)**2. + 1./4.*(f2l - f4l)**2.
              b2l = 13./12.*(f3l - 2.*f4l + f5l)**2. + &
                   1./4.*(3.*f3l - 4.*f4l + f5l)**2.
              
              aar = d0/(eps+b0r)**2. + d1/(eps+b1r)**2. + d2/(eps+b2r)**2.
              aal = d0/(eps+b0l)**2. + d1/(eps+b1l)**2. + d2/(eps+b2l)**2.

              w0r = d0 / (eps + b0r)**2. / aar
              w1r = d1 / (eps + b1r)**2. / aar
              w2r = d2 / (eps + b2r)**2. / aar

              w0l = d0 / (eps + b0l)**2. / aal
              w1l = d1 / (eps + b1l)**2. / aal
              w2l = d2 / (eps + b2l)**2. / aal

              c1r =   1./3.*w0r
              c2r = - 1./6.*w1r - 7./6.*w0r
              c3r =   1./3.*w2r + 5./6.*w1r + 11./6.*w0r
              c4r =   5./6.*w2r + 1./3.*w1r
              c5r = - 1./6.*w2r

              c1l =   1./3.*w0l
              c2l = - 1./6.*w1l - 7./6.*w0l
              c3l =   1./3.*w2l + 5./6.*w1l + 11./6.*w0l
              c4l =   5./6.*w2l + 1./3.*w1l
              c5l = - 1./6.*w2l

              frx = c5r*f5r + c4r*f4r + c3r*f3r + c2r*f2r + c1r*f1r
              flx = c5l*f5l + c4l*f4l + c3l*f3l + c2l*f2l + c1l*f1l

              if(i.eq.1) then
                 flx = 0.
              elseif(i.eq.nx) then
                 frx = 0.
              end if

              !y flux----------------------

              n = ny
          
              isgn = int((v(i,j)+eps)/abs(v(i,j)+eps))

              m1r = (1+isgn)/2 * ((1/j)*n + (2/j)*(j/2)*(n+1) + &
                   j*(j/max(j,3))) + &
                    (1-isgn)/2 * (j/n - ((n-2)/j)*(j/(n-2)) + &
                    j*((n-3)/max(j,n-3)))
           
              m2r = (1+isgn)/2 * ((1/j)*n + j*(j/max(j,2))) + &
                     (1-isgn)/2*(j/n + j*((n-2)/max(j,n-2)))
              
              m3r = j*(1+isgn)/2 + (1-isgn)/2*(j/n + j*((n-1)/max(j,n-1)))
              
              m4r = j*(1-isgn)/2 + (1+isgn)/2*(j/n + j*((n-1)/max(j,n-1)))

              m5r = (1-isgn)/2 * ((1/j)*n + j*(j/max(j,2))) + &
                    (1+isgn)/2 * (j/n + j*((n-2)/max(j,n-2)))

              m1l = (1+isgn)/2 * ((1/j)*n + (2/j)*(j/2)*(n+1) + &
                   (3/j)*(j/3)*(n+2) + j*(j/max(j,4))) + &
                    (1-isgn)/2 * (j/n + j*((n-2)/max(j,n-2)))

              m2l = (1+isgn)/2 * ((1/j)*n + (2/j)*(j/2)*(n+1) + &
                   j*(j/max(j,3))) + &
                    (1-isgn)/2 * (j/n + j*((n-1)/max(j,n-1)))

              m3l = (1+isgn)/2 * ((1/j)*n + j*(j/max(j,2))) + j*(1-isgn)/2

              m4l = (1-isgn)/2 * ((1/j)*n + j*(j/max(j,2))) + j*(1+isgn)/2

              m5l = (1-isgn)/2 * ((1/j)*n + (2/j)*(j/2)*(n+1) + &
                   j*(j/max(j,3))) + (1+isgn)/2 * (j/n + j*((n-1)/max(j,n-1)))


              i1r = m1r - 2 + 5*(1-isgn)/2
              i2r = m2r - 1 + 3*(1-isgn)/2
              i3r = m3r     + 1*(1-isgn)/2
              i4r = m4r + 1 - 1*(1-isgn)/2
              i5r = m5r + 2 - 3*(1-isgn)/2
              
              i1l = m1l - 3 + 5*(1-isgn)/2
              i2l = m2l - 2 + 3*(1-isgn)/2
              i3l = m3l - 1 + 1*(1-isgn)/2
              i4l = m4l     - 1*(1-isgn)/2
              i5l = m5l + 1 - 3*(1-isgn)/2

              eyr = 1.
              eyl = 1.

              f1r = so(i,i1r)**eyr
              f2r = so(i,i2r)**eyr
              f3r = so(i,i3r)**eyr
              f4r = so(i,i4r)**eyr
              f5r = so(i,i5r)**eyr
              
              f1l = so(i,i1l)**eyl
              f2l = so(i,i2l)**eyl
              f3l = so(i,i3l)**eyl
              f4l = so(i,i4l)**eyl
              f5l = so(i,i5l)**eyl

              !WENO3 smoothness functions in y-coordinate

              b0r = 13./12.*(f1r - 2.*f2r + f3r)**2. + &
                   1./4.*(f1r - 4.*f2r + 3.*f3r)**2.
              b1r = 13./12.*(f2r - 2.*f3r + f4r)**2. + 1./4.*(f2r - f4r)**2.
              b2r = 13./12.*(f3r - 2.*f4r + f5r)**2. + &
                   1./4.*(3.*f3r - 4.*f4r + f5r)**2.

              b0l = 13./12.*(f1l - 2.*f2l + f3l)**2. + &
                   1./4.*(f1l - 4.*f2l + 3.*f3l)**2.
              b1l = 13./12.*(f2l - 2.*f3l + f4l)**2. + 1./4.*(f2l - f4l)**2.
              b2l = 13./12.*(f3l - 2.*f4l + f5l)**2. + &
                   1./4.*(3.*f3l - 4.*f4l + f5l)**2.
              
              aar = d0/(eps+b0r)**2. + d1/(eps+b1r)**2. + d2/(eps+b2r)**2.
              aal = d0/(eps+b0l)**2. + d1/(eps+b1l)**2. + d2/(eps+b2l)**2.

              w0r = d0 / (eps + b0r)**2. / aar
              w1r = d1 / (eps + b1r)**2. / aar
              w2r = d2 / (eps + b2r)**2. / aar

              w0l = d0 / (eps + b0l)**2. / aal
              w1l = d1 / (eps + b1l)**2. / aal
              w2l = d2 / (eps + b2l)**2. / aal

              c1r =   1./3.*w0r
              c2r = - 1./6.*w1r - 7./6.*w0r
              c3r =   1./3.*w2r + 5./6.*w1r + 11./6.*w0r
              c4r =   5./6.*w2r + 1./3.*w1r
              c5r = - 1./6.*w2r

              c1l =   1./3.*w0l
              c2l = - 1./6.*w1l - 7./6.*w0l
              c3l =   1./3.*w2l + 5./6.*w1l + 11./6.*w0l
              c4l =   5./6.*w2l + 1./3.*w1l
              c5l = - 1./6.*w2l

              fry = c5r*f5r + c4r*f4r + c3r*f3r + c2r*f2r + c1r*f1r
              fly = c5l*f5l + c4l*f4l + c3l*f3l + c2l*f2l + c1l*f1l

              if(j.eq.1) then
                 fly = 0.
              elseif(j.eq.ny) then
                 fry = 0.
              end if


              s(i,j) = so(i,j) - dt*(frx*ur - flx*ul)/dx &
                               - dt*(fry*vr - fly*vl)/dy

              smx = max(smx,abs(s(i,j)))

              end if

           end do

           do i = nx-2,nx

              sxl = so(i-1,j)
              ul = u(i-1,j)

              if(i.eq.nx) then
                 sxr = so(i-1,j)
              else
                 sxr = so(i+1,j)
              end if
              ur = u(i,j)

              if(j.eq.1) then
                 syl = so(i,j+1)
                 vl = 0.
              else
                 syl = so(i,j-1)
                 vl = v(i,j-1)
              end if
              if(j.eq.ny) then
                 syr = so(i,ny-1)
              else
                 syr = so(i,j+1)
              end if
              vr = v(i,j)

              sm = so(i,j)

              if(ur.ge.0.) then
                 fxr = sm*ur
              else
                 fxr = sxr*ur
              end if

              if(ul.ge.0.) then
                 fxl = sxl*ul
              else
                 fxl = sm*ul
              end if

              if(vr.ge.0.) then
                 fyr = sm*vr
              else
                 fyr = syr*vr
              end if

              if(vl.ge.0.) then
                 fyl = syl*vl
              else
                 fyl = sm*vl
              end if

              s(i,j) = so(i,j) - dt * ( (fxr-fxl)/dx + (fyr-fyl)/dy )

           end do

           end if

        end do


        err = sqrt(sum((s - so)**2.)/dble(nx-1)/dble(ny-1))


      end subroutine advect_WENO3

!------------------------------------------------------------------------
!------------------------------------------------------------------------
!------------------------------------------------------------------------

   
      subroutine advect_WENO31(s,u,v,x,y,cfl,dt,smx,dx,dy,xmn,ymn,nx,ny,uin,s1,intf)
        implicit real*8(a-h,o-z)
        dimension :: s(nx,ny),u(nx,ny),v(nx,ny),x(nx),y(ny),uin(ny+1),s1(ny),intf(nx,ny)
        allocatable :: so(:,:),cm(:,:),nxstart(:,:),nxend(:,:),nystart(:,:),nyend(:,:),nxgap(:),nygap(:)
        allocate(so(nx,ny),cm(3,3),nxstart(ny,10),nxend(ny,10),nystart(nx,10),nyend(nx,10),nxgap(ny),nygap(nx))

        eps = 1E-15
        smx = 0.

        so = s

        cm = reshape( &
             (/ 1./3. , -1./6. ,   1./3.  , &
                5./6. ,  5./6. ,  -7./6.  , &
               -1./6. ,  1./3. ,  11./6.  /) &
             , (/3,3/) )

        d2 = 1./10.
        d1 = 3./5.
        d0 = 3./10.
      
        do j = 1,ny
           nxgap(j) = 1
           nxstart(j,nxgap(j)) = 0
           do i = 2,nx
              if (intf(i,j).eq.0.and.intf(i-1,j).ne.0) then
                 nxgap(j) = nxgap(j) + 1
                 nxstart(j,nxgap(j)) = i-1
              end if
           end do
           nxgap(j) = 0
           do i = 1,nx-1
              if (intf(i,j).eq.0.and.intf(i+1,j).ne.0) then
                 nxgap(j) = nxgap(j) + 1
                 nxend(j,nxgap(j)) = i
              end if
           end do
           nxgap(j) = nxgap(j) + 1
           nxend(j,nxgap(j)) = nx
        end do

        do i = 1,nx
           nygap(i) = 1
           nystart(i,nygap(i)) = 0
           do j = 2,ny
              if (intf(i,j).eq.0.and.intf(i,j-1).ne.0) then
                 nygap(i) = nygap(i) + 1
                 nystart(i,nygap(i)) = j-1
              end if
           end do
           nygap(i) = 1
           nyend(i,nygap(i)) = ny
           do j = 1,ny-1
              if (intf(i,j).eq.0.and.intf(i,j+1).ne.0) then
                 nygap(i) = nygap(i) + 1
                 nyend(i,nygap(i)) = j
              end if
           end do
        end do

        !write(*,*) nygap(121),nystart(121,1),nystart(121,2),nyend(121,1),nyend(121,2)

        !stop
           
        do j = 1,ny

           do i = 1,nx

              if (intf(i,j).eq.0) then


              !              if(j.eq.1.or.intf(i,j-1).ne.0) then
              !                 vl = 0.
              !              else
              !                 vl = v(i,j-1)
              !              end if
                             if(j.eq.1) then
                                vl = 0.
                             elseif (intf(i,j-1).ne.0) then
                                vl = 0.
                             else
                                vl = v(i,j-1)
                             end if
              
              
              !              if(j.eq.ny.or.intf(i,j+1).ne.0) then
              !                 vr = 0.
              !              else
              !                 vr = v(i,j)
              !              end if
              
                             if(j.eq.ny) then
                                vr = 0.
                             elseif (intf(i,j+1).ne.0) then
                                vr = 0.
                             else
                                vr = v(i,j)
                             end if

              !x flux----------------------

              do k = 1,nxgap(j)
                 if (i.ge.nxstart(j,k).and.i.le.nxend(j,k)) then
                    ik = k
                 end if
              end do

              n = nxend(j,ik)-nxstart(j,ik)

              isgn = int((u(i,j)+eps)/abs(u(i,j)+eps))

              m1r = (1+isgn)/2 * ((1/(i-nxstart(j,ik)))*n &
                    + (2/(i-nxstart(j,ik)))*((i-nxstart(j,ik))/2)*(n+1) + &
                   (i-nxstart(j,ik))*((i-nxstart(j,ik))/max((i-nxstart(j,ik)),3))) + &
                   (1-isgn)/2 * ((i-nxstart(j,ik))/n - ((n-2)/(i-nxstart(j,ik)))&
                   *((i-nxstart(j,ik))/(n-2)) + &
                   (i-nxstart(j,ik))*((n-3)/max((i-nxstart(j,ik)),n-3))) + nxstart(j,ik)
           
              m2r = (1+isgn)/2 * ((1/(i-nxstart(j,ik)))*n &
                   + (i-nxstart(j,ik))*((i-nxstart(j,ik))/max((i-nxstart(j,ik)),2))) + &
                   (1-isgn)/2*((i-nxstart(j,ik))/n &
                   + (i-nxstart(j,ik))*((n-2)/max((i-nxstart(j,ik)),n-2))) + nxstart(j,ik)
              
              m3r = (i-nxstart(j,ik))*(1+isgn)/2 + (1-isgn)/2*((i-nxstart(j,ik))/n&
                    + (i-nxstart(j,ik))*((n-1)/max((i-nxstart(j,ik)),n-1))) + nxstart(j,ik)
              
              m4r = (i-nxstart(j,ik))*(1-isgn)/2 + (1+isgn)/2*((i-nxstart(j,ik))/n&
                    + (i-nxstart(j,ik))*((n-1)/max((i-nxstart(j,ik)),n-1))) + nxstart(j,ik)

              m5r = (1-isgn)/2 * ((1/(i-nxstart(j,ik)))*n &
                   + (i-nxstart(j,ik))*((i-nxstart(j,ik))/max((i-nxstart(j,ik)),2))) + &
                   (1+isgn)/2 * ((i-nxstart(j,ik))/n &
                   + (i-nxstart(j,ik))*((n-2)/max((i-nxstart(j,ik)),n-2))) + nxstart(j,ik)

              m1l = (1+isgn)/2 * ((1/(i-nxstart(j,ik)))*n &
                   + (2/(i-nxstart(j,ik)))*((i-nxstart(j,ik))/2)*(n+1) + &
                   (3/(i-nxstart(j,ik)))*((i-nxstart(j,ik))/3)*(n+2)&
                   + (i-nxstart(j,ik))*((i-nxstart(j,ik))/max((i-nxstart(j,ik)),4))) + &
                   (1-isgn)/2 * ((i-nxstart(j,ik))/n &
                   + (i-nxstart(j,ik))*((n-2)/max((i-nxstart(j,ik)),n-2))) + nxstart(j,ik)

              m2l = (1+isgn)/2 * ((1/(i-nxstart(j,ik)))*n &
                   + (2/(i-nxstart(j,ik)))*((i-nxstart(j,ik))/2)*(n+1) + &
                   (i-nxstart(j,ik))*((i-nxstart(j,ik))/max((i-nxstart(j,ik)),3)))&
                   + (1-isgn)/2 * ((i-nxstart(j,ik))/n &
                   + (i-nxstart(j,ik))*((n-1)/max((i-nxstart(j,ik)),n-1))) + nxstart(j,ik)

              m3l = (1+isgn)/2 * ((1/(i-nxstart(j,ik)))*n &
                    + (i-nxstart(j,ik))*((i-nxstart(j,ik))/max((i-nxstart(j,ik)),2)))&
                    + (i-nxstart(j,ik))*(1-isgn)/2 + nxstart(j,ik)

              m4l = (1-isgn)/2 * ((1/(i-nxstart(j,ik)))*n &
                    + (i-nxstart(j,ik))*((i-nxstart(j,ik))/max((i-nxstart(j,ik)),2)))&
                    + (i-nxstart(j,ik))*(1+isgn)/2 + nxstart(j,ik)

              m5l = (1-isgn)/2 * ((1/(i-nxstart(j,ik)))*n &
                   + (2/(i-nxstart(j,ik)))*((i-nxstart(j,ik))/2)*(n+1) + &
                   (i-nxstart(j,ik))*((i-nxstart(j,ik))/max((i-nxstart(j,ik)),3))) &
                   + (1+isgn)/2 * ((i-nxstart(j,ik))/n &
                   + (i-nxstart(j,ik))*((n-1)/max((i-nxstart(j,ik)),n-1))) + nxstart(j,ik)

              i1r = m1r - 2 + 5*(1-isgn)/2
              i2r = m2r - 1 + 3*(1-isgn)/2
              i3r = m3r     + 1*(1-isgn)/2
              i4r = m4r + 1 - 1*(1-isgn)/2
              i5r = m5r + 2 - 3*(1-isgn)/2
              
              i1l = m1l - 3 + 5*(1-isgn)/2
              i2l = m2l - 2 + 3*(1-isgn)/2
              i3l = m3l - 1 + 1*(1-isgn)/2
              i4l = m4l     - 1*(1-isgn)/2
              i5l = m5l + 1 - 3*(1-isgn)/2

              exr = 1.
              exl = 1.

              f1r = so(i1r,j)**exr
              f2r = so(i2r,j)**exr
              f3r = so(i3r,j)**exr
              f4r = so(i4r,j)**exr
              f5r = so(i5r,j)**exr
              
              f1l = so(i1l,j)**exl
              f2l = so(i2l,j)**exl
              f3l = so(i3l,j)**exl
              f4l = so(i4l,j)**exl
              f5l = so(i5l,j)**exl

              !WENO3 smoothness functions in x-coordinate

              b0r = 13./12.*(f1r - 2.*f2r + f3r)**2. + &
                   1./4.*(f1r - 4.*f2r + 3.*f3r)**2.
              b1r = 13./12.*(f2r - 2.*f3r + f4r)**2. + 1./4.*(f2r - f4r)**2.
              b2r = 13./12.*(f3r - 2.*f4r + f5r)**2. + &
                   1./4.*(3.*f3r - 4.*f4r + f5r)**2.

              b0l = 13./12.*(f1l - 2.*f2l + f3l)**2. + &
                   1./4.*(f1l - 4.*f2l + 3.*f3l)**2.
              b1l = 13./12.*(f2l - 2.*f3l + f4l)**2. + 1./4.*(f2l - f4l)**2.
              b2l = 13./12.*(f3l - 2.*f4l + f5l)**2. + &
                   1./4.*(3.*f3l - 4.*f4l + f5l)**2.
              
              aar = d0/(eps+b0r)**2. + d1/(eps+b1r)**2. + d2/(eps+b2r)**2.
              aal = d0/(eps+b0l)**2. + d1/(eps+b1l)**2. + d2/(eps+b2l)**2.

              w0r = d0 / (eps + b0r)**2. / aar
              w1r = d1 / (eps + b1r)**2. / aar
              w2r = d2 / (eps + b2r)**2. / aar

              w0l = d0 / (eps + b0l)**2. / aal
              w1l = d1 / (eps + b1l)**2. / aal
              w2l = d2 / (eps + b2l)**2. / aal

              c1r =   1./3.*w0r
              c2r = - 1./6.*w1r - 7./6.*w0r
              c3r =   1./3.*w2r + 5./6.*w1r + 11./6.*w0r
              c4r =   5./6.*w2r + 1./3.*w1r
              c5r = - 1./6.*w2r

              c1l =   1./3.*w0l
              c2l = - 1./6.*w1l - 7./6.*w0l
              c3l =   1./3.*w2l + 5./6.*w1l + 11./6.*w0l
              c4l =   5./6.*w2l + 1./3.*w1l
              c5l = - 1./6.*w2l

              frx = c5r*f5r + c4r*f4r + c3r*f3r + c2r*f2r + c1r*f1r
              flx = c5l*f5l + c4l*f4l + c3l*f3l + c2l*f2l + c1l*f1l

              !              if(i.eq.1.or.intf(i-1,j).ne.0) then
              !                 flx = 0.
              !              elseif(i.eq.nx.or.intf(i+1,j).ne.0) then
              !                 frx = 0.
              !              end if
                            if(i.eq.1) then
                               flx = 0.
                            elseif (intf(i-1,j).ne.0) then
                               flx = 0.
                            elseif(i.eq.nx) then
                               frx = 0.
                            elseif (intf(i+1,j).ne.0) then
                               frx = 0.
                            end if

              !y flux----------------------
              do k = 1,nygap(i)
                 if (j.ge.nystart(i,k).and.j.le.nyend(i,k)) then
                    jk = k
                 end if
              end do

              n = nyend(i,jk)-nystart(i,jk)

              !write(*,*) i, j, n, jk, nystart(i,jk)
          
              isgn = int((v(i,j)+eps)/abs(v(i,j)+eps))

              m1r = (1+isgn)/2 * ((1/(j-nystart(i,jk)))*n &
                    + (2/(j-nystart(i,jk)))*((j-nystart(i,jk))/2)*(n+1) + &
                   (j-nystart(i,jk))*((j-nystart(i,jk))/max((j-nystart(i,jk)),3))) + &
                    (1-isgn)/2 * ((j-nystart(i,jk))/n - ((n-2)/(j-nystart(i,jk)))&
                    *((j-nystart(i,jk))/(n-2)) + &
                    (j-nystart(i,jk))*((n-3)/max((j-nystart(i,jk)),n-3))) + nystart(i,jk)
           
              m2r = (1+isgn)/2 * ((1/(j-nystart(i,jk)))*n &
                     + (j-nystart(i,jk))*((j-nystart(i,jk))/max((j-nystart(i,jk)),2))) + &
                     (1-isgn)/2*((j-nystart(i,jk))/n &
                     + (j-nystart(i,jk))*((n-2)/max((j-nystart(i,jk)),n-2))) + nystart(i,jk)
              
              m3r = (j-nystart(i,jk))*(1+isgn)/2 + (1-isgn)/2*((j-nystart(i,jk))/n &
                    + (j-nystart(i,jk))*((n-1)/max((j-nystart(i,jk)),n-1))) + nystart(i,jk)
              
              m4r = (j-nystart(i,jk))*(1-isgn)/2 + (1+isgn)/2*((j-nystart(i,jk))/n &
                    + (j-nystart(i,jk))*((n-1)/max((j-nystart(i,jk)),n-1))) + nystart(i,jk)

              m5r = (1-isgn)/2 * ((1/(j-nystart(i,jk)))*n &
                    + (j-nystart(i,jk))*((j-nystart(i,jk))/max((j-nystart(i,jk)),2))) + &
                    (1+isgn)/2 * ((j-nystart(i,jk))/n &
                    + (j-nystart(i,jk))*((n-2)/max((j-nystart(i,jk)),n-2))) + nystart(i,jk)

              m1l = (1+isgn)/2 * ((1/(j-nystart(i,jk)))*n &
                    + (2/(j-nystart(i,jk)))*((j-nystart(i,jk))/2)*(n+1) + &
                    (3/(j-nystart(i,jk)))*((j-nystart(i,jk))/3)*(n+2) &
                    + (j-nystart(i,jk))*((j-nystart(i,jk))/max((j-nystart(i,jk)),4))) + &
                    (1-isgn)/2 * ((j-nystart(i,jk))/n &
                    + (j-nystart(i,jk))*((n-2)/max((j-nystart(i,jk)),n-2))) + nystart(i,jk)

              m2l = (1+isgn)/2 * ((1/(j-nystart(i,jk)))*n &
                    + (2/(j-nystart(i,jk)))*((j-nystart(i,jk))/2)*(n+1) + &
                    (j-nystart(i,jk))*((j-nystart(i,jk))/max((j-nystart(i,jk)),3))) + &
                    (1-isgn)/2 * ((j-nystart(i,jk))/n &
                    + (j-nystart(i,jk))*((n-1)/max((j-nystart(i,jk)),n-1))) + nystart(i,jk)

              m3l = (1+isgn)/2 * ((1/(j-nystart(i,jk)))*n &
                    + (j-nystart(i,jk))*((j-nystart(i,jk))/max((j-nystart(i,jk)),2))) &
                    + (j-nystart(i,jk))*(1-isgn)/2 + nystart(i,jk)

              m4l = (1-isgn)/2 * ((1/(j-nystart(i,jk)))*n &
                    + (j-nystart(i,jk))*((j-nystart(i,jk))/max((j-nystart(i,jk)),2))) &
                    + (j-nystart(i,jk))*(1+isgn)/2 + nystart(i,jk)

              m5l = (1-isgn)/2 * ((1/(j-nystart(i,jk)))*n &
                    + (2/(j-nystart(i,jk)))*((j-nystart(i,jk))/2)*(n+1) + &
                    (j-nystart(i,jk))*((j-nystart(i,jk))/max((j-nystart(i,jk)),3))) &
                    + (1+isgn)/2 * ((j-nystart(i,jk))/n &
                    + (j-nystart(i,jk))*((n-1)/max((j-nystart(i,jk)),n-1))) + nystart(i,jk)


              i1r = m1r - 2 + 5*(1-isgn)/2
              i2r = m2r - 1 + 3*(1-isgn)/2
              i3r = m3r     + 1*(1-isgn)/2
              i4r = m4r + 1 - 1*(1-isgn)/2
              i5r = m5r + 2 - 3*(1-isgn)/2
              
              i1l = m1l - 3 + 5*(1-isgn)/2
              i2l = m2l - 2 + 3*(1-isgn)/2
              i3l = m3l - 1 + 1*(1-isgn)/2
              i4l = m4l     - 1*(1-isgn)/2
              i5l = m5l + 1 - 3*(1-isgn)/2

              eyr = 1.
              eyl = 1.

              f1r = so(i,i1r)**eyr
              f2r = so(i,i2r)**eyr
              f3r = so(i,i3r)**eyr
              f4r = so(i,i4r)**eyr
              f5r = so(i,i5r)**eyr
              
              f1l = so(i,i1l)**eyl
              f2l = so(i,i2l)**eyl
              f3l = so(i,i3l)**eyl
              f4l = so(i,i4l)**eyl
              f5l = so(i,i5l)**eyl

              !WENO3 smoothness functions in y-coordinate

              b0r = 13./12.*(f1r - 2.*f2r + f3r)**2. + &
                   1./4.*(f1r - 4.*f2r + 3.*f3r)**2.
              b1r = 13./12.*(f2r - 2.*f3r + f4r)**2. + 1./4.*(f2r - f4r)**2.
              b2r = 13./12.*(f3r - 2.*f4r + f5r)**2. + &
                   1./4.*(3.*f3r - 4.*f4r + f5r)**2.

              b0l = 13./12.*(f1l - 2.*f2l + f3l)**2. + &
                   1./4.*(f1l - 4.*f2l + 3.*f3l)**2.
              b1l = 13./12.*(f2l - 2.*f3l + f4l)**2. + 1./4.*(f2l - f4l)**2.
              b2l = 13./12.*(f3l - 2.*f4l + f5l)**2. + &
                   1./4.*(3.*f3l - 4.*f4l + f5l)**2.
              
              aar = d0/(eps+b0r)**2. + d1/(eps+b1r)**2. + d2/(eps+b2r)**2.
              aal = d0/(eps+b0l)**2. + d1/(eps+b1l)**2. + d2/(eps+b2l)**2.

              w0r = d0 / (eps + b0r)**2. / aar
              w1r = d1 / (eps + b1r)**2. / aar
              w2r = d2 / (eps + b2r)**2. / aar

              w0l = d0 / (eps + b0l)**2. / aal
              w1l = d1 / (eps + b1l)**2. / aal
              w2l = d2 / (eps + b2l)**2. / aal

              c1r =   1./3.*w0r
              c2r = - 1./6.*w1r - 7./6.*w0r
              c3r =   1./3.*w2r + 5./6.*w1r + 11./6.*w0r
              c4r =   5./6.*w2r + 1./3.*w1r
              c5r = - 1./6.*w2r

              c1l =   1./3.*w0l
              c2l = - 1./6.*w1l - 7./6.*w0l
              c3l =   1./3.*w2l + 5./6.*w1l + 11./6.*w0l
              c4l =   5./6.*w2l + 1./3.*w1l
              c5l = - 1./6.*w2l

              fry = c5r*f5r + c4r*f4r + c3r*f3r + c2r*f2r + c1r*f1r
              fly = c5l*f5l + c4l*f4l + c3l*f3l + c2l*f2l + c1l*f1l

              !              if(j.eq.1.or.intf(i,j-1).ne.0) then
              !                 fly = 0.
              !              elseif(j.eq.ny.or.intf(i,j+1).ne.0) then
              !                 fry = 0.
              !              end if
                            if(j.eq.1) then
                               fly = 0.
                            elseif (intf(i,j-1).ne.0) then
                               fly = 0.
                            elseif(j.eq.ny) then
                               fry = 0.
                            elseif (intf(i,j+1).ne.0) then
                               fry = 0.
                            end if


              s(i,j) = so(i,j) - dt*(frx*ur - flx*ul)/dx &
                               - dt*(fry*vr - fly*vl)/dy

              smx = max(smx,abs(s(i,j)))

           end if

           end do

        end do


        err = sqrt(sum((s - so)**2.)/dble(nx-1)/dble(ny-1))


      end subroutine advect_WENO31


!------------------------------------------------------------------------
!------------------------------------------------------------------------
!------------------------------------------------------------------------

      subroutine advect_UW(phi,u,v,uin,x,y,dx,dy,cfl,dt,phimx,nx,ny,intf)
        implicit real*8(a-h,o-z)
        dimension :: phi(nx,ny),u(nx,ny),v(nx,ny),x(nx),y(ny),uin(ny+1),intf(nx,ny)
        allocatable :: phi1(:,:)
        allocate(phi1(nx,ny))

        phi1 = phi
        phi = 0.

        phimx = -1E10

        do j = 1,ny
           do i = 1,nx

              if(i.eq.1.or.intf(i-1,j).ne.0) then
                 sxl = phi1(i+1,j)
                 ul = uin(j)
              else
                 sxl = phi1(i-1,j)
                 ul = u(i-1,j)
              end if
              if(i.eq.nx.or.intf(i+1,j).ne.0) then
                 sxr = phi1(i-1,j)
              else
                 sxr = phi1(i+1,j)
              end if
              ur = u(i,j)

              if(j.eq.1.or.intf(i,j-1).ne.0) then
                 syl = phi1(i,j+1)
                 vl = 0.
              else
                 syl = phi1(i,j-1)
                 vl = v(i,j-1)
              end if
              if(j.eq.ny.or.intf(i,j+1).ne.0) then
                 syr = phi1(i,j-1)
              else
                 syr = phi1(i,j+1)
              end if
              vr = v(i,j)

              sm = phi1(i,j)

              if(ur.ge.0.) then
                 fxr = sm*ur
              else
                 fxr = sxr*ur
              end if

              if(ul.ge.0.) then
                 fxl = sxl*ul
              else
                 fxl = sm*ul
              end if


              if(vr.ge.0.) then
                 fyr = sm*vr
              else
                 fyr = syr*vr
              end if

              if(vl.ge.0.) then
                 fyl = syl*vl
              else
                 fyl = sm*vl
              end if


              phi(i,j) = phi1(i,j) - dt * ( (fxr-fxl)/dx + (fyr-fyl)/dy )

              phimx = max(abs(phi(i,j)),phimx)

           end do
        end do


      end subroutine advect_UW

!------------------------------------------------------------------------
!------------------------------------------------------------------------
!------------------------------------------------------------------------


      subroutine solid_property_level(x,y,xo,yo,rdx,rho1xs,rho2xs,rho1ys,rho2ys,ixmn,ixmx,iymn,iymx,rhos)

        implicit real*8(a-h,o-z)
        common/param/g,sp,ubc,uout,mint
        common/grid/dx,dy,xmx,xmn,ymx,ymn,nx,ny,tt,tb,crmx,ns,nsp,nsl
        dimension :: x(nx),y(ny),xo(ns),yo(ns),rho1xs(nx,ny),rho2xs(nx,ny),rho1ys(nx,ny),&
                     rho2ys(nx,ny),n_s(ns),ixmn(ns),ixmx(ns),iymn(ns),iymx(ns)

        do k = 1,ns
          !particle
           do j = iymn(k),iymx(k)-1
              do i = ixmn(k),ixmx(k)-1
                 a =  - rdx + sqrt((x(i)-xo(k))**2. + (y(j)-yo(k))**2.)         ! i,  j
                 b =  - rdx + sqrt((x(i+1)-xo(k))**2. + (y(j)-yo(k))**2.)       ! i+1,j
                 c =  - rdx + sqrt((x(i)-xo(k))**2. + (y(j+1)-yo(k))**2.)       ! i,  j+1
                 if (a.le.0..and.b.le.0.) then
                    rho1xs(i,j) = 1./rhos
                    rho2xs(i,j) = 0.
                 end if
                 if (a.le.0..and.c.le.0.) then
                    rho1ys(i,j) = 1./rhos
                    rho2ys(i,j) = 0.
                 end if
              end do
           end do
        end do

      end subroutine solid_property_level

!---------------------------------------------------------------------------------
!---------------------------------------------------------------------------------
      subroutine collid_find(xo,yo,ncluster,nsize,nmcluster,rdx,ns)

        implicit real*8(a-h,o-z)
        dimension :: xo(ns),yo(ns),nmcluster(ns,ns),nsize(ns)
        allocatable :: nmfind(:),nmclusternew(:,:),nsizenew(:)
        allocate(nmfind(ns),nmclusternew(ns,ns),nsizenew(ns))

        !-------------------------------------------find the colliding the crystals-----------------------------------------
        ncluster = 0
        do k1 = 1,ns
           nfind1 = 0
           do k2 = 1,ncluster
              nfind2 = 0
              do k3 = 1,nsize(k2)
                 dis = sqrt((xo(k1)-xo(nmcluster(k3,k2)))**2.+(yo(k1)-yo(nmcluster(k3,k2)))**2.)
                 if (dis.le.2.*rdx) then
                    nfind2 = 1
                    exit
                 end if
              end do
              if (nfind2.eq.1) then
                 nfind1 = nfind1 + 1
                 nmfind(nfind1) = k2
              end if
           end do
           !write(*,*) k1,nfind1
           if (nfind1.eq.0) then
              ncluster = ncluster + 1
              nsize(ncluster) = 1
              nmcluster(nsize(ncluster),ncluster) = k1
           else if (nfind1.eq.1) then
              nsize(nmfind(nfind1)) = nsize(nmfind(nfind1)) + 1
              nmcluster(nsize(nmfind(nfind1)),nmfind(nfind1)) = k1
           else
              nclusternew = 0
              do k4 = 1,ncluster
                 nfind3 = 0
                 do k5 = 1,nfind1
                    if (k4.eq.nmfind(k5)) then
                       nfind3 = 1
                       exit
                    end if
                 end do
                 if (nfind3.eq.0) then
                    nclusternew = nclusternew + 1
                    nsizenew(nclusternew) = nsize(k4)
                    do k51 = 1,nsizenew(nclusternew)
                       nmclusternew(k51,nclusternew) = nmcluster(k51,k4)
                    end do
                 end if
              end do
              nclusternew = nclusternew + 1
              nsizelast = 0
              do k6 = 1,nfind1
                 nmclusternew(nsizelast+1:nsizelast + nsize(nmfind(k6)),nclusternew) = nmcluster(1:nsize(nmfind(k6)),nmfind(k6))
                 nsizelast = nsizelast + nsize(nmfind(k6))
              end do
              nsizelast = nsizelast + 1
              nmclusternew(nsizelast,nclusternew) = k1
              nsizenew(nclusternew) = nsizelast
              ncluster = nclusternew
              do k7 = 1,ncluster
                 nsize(k7) = nsizenew(k7)
                 do k8 = 1, nsize(k7)
                    nmcluster(k8,k7) = nmclusternew(k8,k7)
                 end do
              end do
           end if
        end do
      end subroutine collid_find

!---------------------------------------------------------------------------------
!---------------------------------------------------------------------------------
      subroutine collid_boundary(xo,yo,nxbound,nybound,nwall_l,nwall_r,dx,dy,rdx,nx,ny,xmn,ymn,ns,&
                                 a1,a2,a3,b1,b2,b3)

        implicit real*8(a-h,o-z)

        dimension :: xo(ns),yo(ns),nxbound(ns),nybound(ns),nwall_l(ns),nwall_r(ns),x(nx),y(ny)

        nxbound = 0
        nybound = 0
        do k = 1,ns
           if (floor((xo(k)+2.*rdx-xmn)/dx)+1.ge.nx) then
              nxbound(k) = 1
           end if
           if (floor((xo(k)-2.*rdx-0.)/dx).le.1) then
              nxbound(k) = 2
           end if
           if (floor((yo(k)+2.*rdx-ymn)/dy)+1.ge.ny) then
              nybound(k) = 1
           end if
           if (floor((yo(k)-2.*rdx-ymn)/dy).le.1) then
              nybound(k) = 2
           end if
        end do

        nwall_l = 0
        nwall_r = 0
        !---------------------------collison strategy for wall-crystal case----------------------------
        do k = 1,ns
           x1 = (a2**2.*xo(k)-a1*a2*yo(k)-a1*a3)/(a1**2.+a2**2.)
           y1 = (-a1*a2*xo(k)+a1**2.*yo(k)-a2*a3)/(a1**2.+a2**2.)
           d = sqrt((x1-xo(k))**2.+(y1-yo(k))**2.)
           if (d.le.2.*rdx) then
              nwall_l(k) = 1
           end if
        end do

        do k = 1,ns
           x1 = (b2**2.*xo(k)-b1*b2*yo(k)-b1*b3)/(b1**2.+b2**2.)
           y1 = (-b1*b2*xo(k)+b1**2.*yo(k)-b2*b3)/(b1**2.+b2**2.)
           d = sqrt((x1-xo(k))**2.+(y1-yo(k))**2.)
           if (d.le.2.*rdx) then
              nwall_r(k) = 1
           end if
        end do

    
      end subroutine collid_boundary

!------------------------------------------------------------------------
!------------------------------------------------------------------------
!------------------------------------------------------------------------
      subroutine collid_strategy(rdx,xo,yo,ubr,vbr,dx,dy,nmcluster,nsize,ncluster,nxbound,nybound,&
                                 nwall_l,nwall_r,ns,xmx,xmn,ymx,ymn,a1,a2,a3,b1,b2,b3)

        implicit real*8(a-h,o-z)
        dimension :: xo(ns),yo(ns),ubr(ns),vbr(ns),nmcluster(ns,ns),nsize(ns),&
                     nxbound(ns),nybound(ns),nwall_l(ns),nwall_r(ns)

        allocatable :: nmcoll(:),angmutual1(:),angmutual2(:),umutual(:),angfree1(:),angfree2(:),ufree(:),&
                       nbound1(:),free1(:),free2(:)

        eps = 10E-15

        mumf = 9*nm
        allocate(nmcoll(ns),angmutual1(ns),angmutual2(ns),umutual(ns),angfree1(ns),angfree2(ns),ufree(ns),&
                 nbound1(ns),free1(ns),free2(ns))

        do k = 1,ncluster
           if (nsize(k).eq.1) then

              !-----------------------for the single-crystal cluster, change the direction of the velocity when it hits the wall------------------------
              kk1 = nmcluster(1,k)
              if (nxbound(kk1).eq.1) then
                 ubr(kk1) = -abs(ubr(kk1))
              elseif (nxbound(kk1).eq.2) then
                 ubr(kk1) =  abs(ubr(kk1))
              end if
              if (nybound(kk1).eq.1) then
                 vbr(kk1) = -abs(vbr(kk1))
              elseif (nybound(kk1).eq.2) then
                 vbr(kk1) =  abs(vbr(kk1))
              end if
              if (nwall_l(kk1).eq.1) then
                 x1 = (a2**2.*xo(k)-a1*a2*yo(k)-a1*a3)/(a1**2.+a2**2.)
                 y1 = (-a1*a2*xo(k)+a1**2.*yo(k)-a2*a3)/(a1**2.+a2**2.)
                 dis = sqrt((x1-xo(kk1))**2.+(y1-yo(kk1))**2.)
                 unn = -abs(ubr(kk1)*(x1-xo(kk1))/dis + vbr(kk1)*(y1-yo(kk1))/dis)
                 utt = ubr(kk1)*(y1-yo(kk1))/dis-vbr(kk1)*(x1-xo(kk1))/dis
                 ubr(kk1)=unn*(x1-xo(kk1))/dis + utt*(y1-yo(kk1))/dis
                 vbr(kk1)=unn*(y1-yo(kk1))/dis - utt*(x1-xo(kk1))/dis
              end if
              if (nwall_r(kk1).eq.1) then
                 x1 = (b2**2.*xo(k)-b1*b2*yo(k)-b1*b3)/(b1**2.+b2**2.)
                 y1 = (-b1*b2*xo(k)+b1**2.*yo(k)-b2*b3)/(b1**2.+b2**2.)
                 dis = sqrt((x1-xo(kk1))**2.+(y1-yo(kk1))**2.)
                 unn = -abs(ubr(kk1)*(x1-xo(kk1))/dis + vbr(kk1)*(y1-yo(kk1))/dis)
                 utt = ubr(kk1)*(y1-yo(kk1))/dis-vbr(kk1)*(x1-xo(kk1))/dis
                 ubr(kk1)=unn*(x1-xo(kk1))/dis + utt*(y1-yo(kk1))/dis
                 vbr(kk1)=unn*(y1-yo(kk1))/dis - utt*(x1-xo(kk1))/dis
              end if
              !-----------------------------------------------------------------------------------------------------------------------------------------

           else
              if (nsize(k).eq.2) then
 
                 !----------------for the two-crystals cluster, consider the effect from both wall and the interaction between each other---------------
                 nbound1 = 0
                 do kk = 1,nsize(k)
                    kk1 = nmcluster(kk,k)
                    do kkk = 1,nsize(k)
                       kkkk = nmcluster(kkk,k)
                       if (kkkk.ne.kk1) then
                          kk2 = kkkk
                       end if
                    end do
                    !-----------------free direction for each crystal----------
                    if (nxbound(kk1).eq.1) then
                       nbound1(kk1) = 1
                       free1(kk1) = xo(kk1)-xo(kk2) + xo(kk1)-xmx
                       free2(kk1) = yo(kk1)-yo(kk2) + 0.
                    end if
                    if (nxbound(kk1).eq.2) then
                       nbound1(kk1) = 1
                       free1(kk1) = xo(kk1)-xo(kk2) + xo(kk1)-0.
                       free2(kk1) = yo(kk1)-yo(kk2) + 0.
                    end if
                    if (nybound(kk1).eq.1) then
                       nbound1(kk1) = 1
                       free1(kk1) = xo(kk1)-xo(kk2) + 0.
                       free2(kk1) = yo(kk1)-yo(kk2) + yo(kk1)-ymx
                    end if
                    if (nybound(kk1).eq.2) then
                       nbound1(kk1) = 1
                       free1(kk1) = xo(kk1)-xo(kk2) + 0.
                       free2(kk1) = yo(kk1)-yo(kk2) + yo(kk1)-ymn
                    end if
                    if (nwall_l(kk1).eq.1) then
                       nbound1(kk1) = 1
                       x1 = (a2**2.*xo(k)-a1*a2*yo(k)-a1*a3)/(a1**2.+a2**2.)
                       y1 = (-a1*a2*xo(k)+a1**2.*yo(k)-a2*a3)/(a1**2.+a2**2.)
                       free1(kk1) = xo(kk1)-xo(kk2) + xo(kk1)-x1
                       free2(kk1) = yo(kk1)-yo(kk2) + yo(kk1)-y1
                    end if
                    if (nwall_r(kk1).eq.1) then
                       nbound1(kk1) = 1
                       x1 = (b2**2.*xo(k)-b1*b2*yo(k)-b1*b3)/(b1**2.+b2**2.)
                       y1 = (-b1*b2*xo(k)+b1**2.*yo(k)-b2*b3)/(b1**2.+b2**2.)
                       free1(kk1) = xo(kk1)-xo(kk2) + xo(kk1)-x1
                       free2(kk1) = yo(kk1)-yo(kk2) + yo(kk1)-y1
                    end if
                    !-------------------------------------------------------
                 end do
                 !----------------------check whether there is any crystal in the cluster hit the wall---------------
                 nbound2 = 0
                 do kk = 1,nsize(k)
                    kk1 = nmcluster(kk,k)
                    if (nbound1(kk1).eq.1) then
                       nbound2 = 1
                    end if
                 end do
                 !---------------------------------------------------------------------------------------------------
                 do kk = 1,nsize(k)
                    kk1 = nmcluster(kk,k)
                    do kkk = 1,nsize(k)
                       kkkk = nmcluster(kkk,k)
                       if (kkkk.ne.kk1) then
                          kk2 = kkkk
                       end if
                    end do
                    if (nbound1(kk1).eq.1) then
                       nfree = 1
                       dis = sqrt(free1(kk1)**2.+free2(kk1)**2.)
                       free1(kk1) = free1(kk1)/dis
                       free2(kk1) = free2(kk1)/dis
                       ufreek1 = ubr(kk1)*free1(kk1)+vbr(kk1)*free2(kk1)
                       if (ufreek1.lt.0.) then
                          nfree = 0
                       else
                          if (free1(kk1)*(xo(kk1)-xo(kk2))+free2(kk1)*(yo(kk1)-yo(kk2)).le.0.) then
                             ufreek2 = ubr(kk2)*free1(kk1)+vbr(kk2)*free2(kk1)
                             if ((ufreek1.ge.0.and.ufreek2.le.0).or.(ufreek1.ge.0.and.ufreek2.ge.0.and.ufreek2.le.ufreek1).or.&
                                (ufreek1.le.0.and.ufreek2.le.0.and.ufreek2.le.ufreek1)) then
                                nfree = 0
                             end if
                          else
                             ufreek2 = ubr(kk2)*free1(kk1)+vbr(kk2)*free2(kk1)
                             if (ufreek1.ge.0.and.ufreek2.ge.0.and.ufreek2.ge.ufreek1) then
                                nfree = 0
                             end if
                          end if
                       end if
                       if (nfree.eq.0) then
                          angfree1(kk1) = 0.
                          angfree2(kk1) = 0.
                          ufree(kk1) = 0.
                          dis = sqrt(ubr(kk1)**2.+vbr(kk1)**2.)
                          angmutual1(kk1) = ubr(kk1)/dis
                          angmutual2(kk1) = vbr(kk1)/dis
                          umutual(kk1) = 0.
                       else
                          angfree1(kk1) = free1(kk1)
                          angfree2(kk1) = free2(kk1)
                          ufree(kk1) = ubr(kk1)*angfree1(kk1)+vbr(kk1)*angfree2(kk1)
                          angmutual1(kk1) = -free2(kk1)
                          angmutual2(kk1) =  free1(kk1)
                          umutual(kk1) = 0.
                       end if
                    else
                       if (nbound2.eq.1) then
                          dis = sqrt((xo(kk1)-xo(kk2))**2.+(yo(kk1)-yo(kk2))**2.)
                          angfree1(kk1) = -(yo(kk1)-yo(kk2))/dis
                          angfree2(kk1) = (xo(kk1)-xo(kk2))/dis
                          ufree(kk1) = ubr(kk1)*angfree1(kk1)+vbr(kk1)*angfree2(kk1)
                          angmutual1(kk1) = (xo(kk1)-xo(kk2))/dis
                          angmutual2(kk1) = (yo(kk1)-yo(kk2))/dis
                          umutual(kk1) = 0.
                       else
                          dis = sqrt((xo(kk1)-xo(kk2))**2.+(yo(kk1)-yo(kk2))**2.)
                          angfree1(kk1) = -(yo(kk1)-yo(kk2))/dis
                          angfree2(kk1) = (xo(kk1)-xo(kk2))/dis
                          ufree(kk1) = ubr(kk1)*angfree1(kk1)+vbr(kk1)*angfree2(kk1)
                          angmutual1(kk1) = (xo(kk1)-xo(kk2))/dis
                          angmutual2(kk1) = (yo(kk1)-yo(kk2))/dis
                          umutual(kk1) = (ubr(kk1)*angmutual1(kk1)+vbr(kk1)*angmutual2(kk1))/dble(nsize(k))
                       end if
                    end if
                 end do
                 do kk = 1,nsize(k)
                    kk1 = nmcluster(kk,k)
                    ubr(kk1) = ufree(kk1)*angfree1(kk1)
                    vbr(kk1) = ufree(kk1)*angfree2(kk1)
                    do kkk = 1,nsize(k)
                       kk2 = nmcluster(kkk,k)
                       ubr(kk1) = ubr(kk1) + umutual(kk2)*angmutual1(kk2)
                       vbr(kk1) = vbr(kk1) + umutual(kk2)*angmutual2(kk2)
                    end do
                 end do
                 !---------------------------------------------------------------------------------------------------------------------------
              else
                 !--------------------------------------------for the multi_crystals cluster-------------------------------------------------
                 !-----------------------------------look for the crystals hiting the wall-------------------------------------
                 nbound1 = 0
                 do kk = 1,nsize(k)
                    kk1 = nmcluster(kk,k)
                    if (nxbound(kk1).eq.1) then
                       nbound1(kk1) = 1
                    end if
                    if (nxbound(kk1).eq.2) then
                       nbound1(kk1) = 1
                    end if
                    if (nybound(kk1).eq.1) then
                       nbound1(kk1) = 1
                    end if
                    if (nybound(kk1).eq.2) then
                       nbound1(kk1) = 1
                    end if
                    if (nwall_l(kk1).eq.1) then
                       nbound1(kk1) = 1
                    end if
                    if (nwall_r(kk1).eq.1) then
                       nbound1(kk1) = 1
                    end if
                 end do
                 !-----------------------------------------------------------------------------------------------------------
                 !-------------------------------check whether there is any crystal hiting the wall--------------------------
                 nbound2 = 0
                 do kk = 1,nsize(k)
                    kk1 = nmcluster(kk,k)
                    if (nbound1(kk1).eq.1) then
                       nbound2 = 1
                    end if
                 end do
                 !-----------------------------------------------------------------------------------------------------------

                 do kk = 1,nsize(k)
                    kk1 = nmcluster(kk,k)
                    ncoll = 0
                    do kkk = 1,nsize(k)
                       kk2 = nmcluster(kkk,k)
                       dis = sqrt((xo(kk1)-xo(kk2))**2.+(yo(kk1)-yo(kk2))**2.)
                       if (kk2.ne.kk1.and.dis.le.2.*(rdx)) then
                          ncoll = ncoll + 1
                          nmcoll(ncoll) = kk2
                       end if
                    end do

                    !----------------------------compute the free direction-----------------------------------
                    free1(kk1) = 0.
                    free2(kk1) = 0.
                    do kkk = 1,ncoll
                       kk2 = nmcoll(kkk)
                       free1(kk1) = free1(kk1) + xo(kk1)-xo(kk2)
                       free2(kk1) = free2(kk1) + yo(kk1)-yo(kk2)
                    end do
                    if (nxbound(kk1).eq.1) then
                       free1(kk1) = free1(kk1) + xo(kk1)-xmx
                       free2(kk1) = free2(kk1) + 0.
                    end if
                    if (nxbound(kk1).eq.2) then
                       free1(kk1) = free1(kk1) + xo(kk1)-0.
                       free2(kk1) = free2(kk1) + 0.
                    end if
                    if (nybound(kk1).eq.1) then
                       free1(kk1) = free1(kk1) + 0.
                       free2(kk1) = free2(kk1) + yo(kk1)-ymx
                    end if
                    if (nybound(kk1).eq.2) then
                       free1(kk1) = free1(kk1) + 0.
                       free2(kk1) = free2(kk1) + yo(kk1)-ymn
                    end if
                    if (nwall_l(kk1).eq.1) then
                       x1 = (a2**2.*xo(k)-a1*a2*yo(k)-a1*a3)/(a1**2.+a2**2.)
                       y1 = (-a1*a2*xo(k)+a1**2.*yo(k)-a2*a3)/(a1**2.+a2**2.)
                       free1(kk1) = free1(kk1) + xo(kk1)-x1
                       free2(kk1) = free2(kk1) + yo(kk1)-y1
                    end if
                    if (nwall_r(kk1).eq.1) then
                       x1 = (b2**2.*xo(k)-b1*b2*yo(k)-b1*b3)/(b1**2.+b2**2.)
                       y1 = (-b1*b2*xo(k)+b1**2.*yo(k)-b2*b3)/(b1**2.+b2**2.)
                       free1(kk1) = free1(kk1) + xo(kk1)-x1
                       free2(kk1) = free2(kk1) + yo(kk1)-y1
                    end if
                    !----------------------------------------------------------------------------------------

                    !---------------------------for the crystal who only hit one other crystal---------------
                    if (ncoll.eq.1) then
                       kk2 = nmcoll(1)
                       if (nbound1(kk1).eq.1) then
                          nfree = 1
                          dis = sqrt(free1(kk1)**2.+free2(kk1)**2.)
                          free1(kk1) = free1(kk1)/dis
                          free2(kk1) = free2(kk1)/dis
                          ufreek1 = ubr(kk1)*free1(kk1)+vbr(kk1)*free2(kk1)
                          if (ufreek1.lt.0.) then
                             nfree = 0
                          else
                             if (free1(kk1)*(xo(kk1)-xo(kk2))+free2(kk1)*(yo(kk1)-yo(kk2)).le.0.) then
                                ufreek2 = ubr(kk2)*free1(kk1)+vbr(kk2)*free2(kk1)
                                if ((ufreek1.ge.0.and.ufreek2.le.0).or.(ufreek1.ge.0.and.ufreek2.ge.0.and.ufreek2.le.ufreek1).or.&
                                   (ufreek1.le.0.and.ufreek2.le.0.and.ufreek2.le.ufreek1)) then
                                   nfree = 0
                                end if
                             else
                                ufreek2 = ubr(kk2)*free1(kk1)+vbr(kk2)*free2(kk1)
                                if (ufreek1.ge.0.and.ufreek2.ge.0.and.ufreek2.ge.ufreek1) then
                                   nfree = 0
                                end if
                             end if
                          end if
                          if (nfree.eq.0) then
                             angfree1(kk1) = 0.
                             angfree2(kk1) = 0.
                             ufree(kk1) = 0.
                             dis = sqrt(ubr(kk1)**2.+vbr(kk1)**2.)
                             angmutual1(kk1) = ubr(kk1)/dis
                             angmutual2(kk1) = vbr(kk1)/dis
                             umutual(kk1) = 0.
                          else
                             angfree1(kk1) = free1(kk1)
                             angfree2(kk1) = free2(kk1)
                             ufree(kk1) = ubr(kk1)*angfree1(kk1)+vbr(kk1)*angfree2(kk1)
                             angmutual1(kk1) = -free2(kk1)
                             angmutual2(kk1) =  free1(kk1)
                             umutual(kk1) = 0.
                          end if
                       else
                          dis = sqrt((xo(kk1)-xo(kk2))**2.+(yo(kk1)-yo(kk2))**2.)
                          angfree1(kk1) = -(yo(kk1)-yo(kk2))/dis
                          angfree2(kk1) = (xo(kk1)-xo(kk2))/dis
                          ufree(kk1) = ubr(kk1)*angfree1(kk1)+vbr(kk1)*angfree2(kk1)
                          angmutual1(kk1) = (xo(kk1)-xo(kk2))/dis
                          angmutual2(kk1) = (yo(kk1)-yo(kk2))/dis
                          if (nbound2.eq.1) then
                             umutual(kk1) = 0.
                          else
                             umutual(kk1) = (ubr(kk1)*angmutual1(kk1)+vbr(kk1)*angmutual2(kk1))/dble(nsize(k))
                          end if
                       end if
                    !----------------------------------------------------------------------------------------

                    elseif (ncoll.eq.2) then
                       !-------------------------for the crystl who hit two other crystal-----------------------
                       nfree = 1
                       dis = sqrt(free1(kk1)**2.+free2(kk1)**2.)
                       free1(kk1) = free1(kk1)/dis
                       free2(kk1) = free2(kk1)/dis
                       ufreek1 = ubr(kk1)*free1(kk1)+vbr(kk1)*free2(kk1)
                       if (ufreek1.lt.0.) then
                          nfree = 0
                       else
                          do kkk = 1,ncoll
                             kk2 = nmcoll(kkk)
                             if (free1(kk1)*(xo(kk1)-xo(kk2))+free2(kk1)*(yo(kk1)-yo(kk2)).le.0.) then
                                ufreek2 = ubr(kk2)*free1(kk1)+vbr(kk2)*free2(kk1)
                                if ((ufreek1.ge.0.and.ufreek2.le.0).or.(ufreek1.ge.0.and.ufreek2.ge.0.and.ufreek2.le.ufreek1).or.&
                                   (ufreek1.le.0.and.ufreek2.le.0.and.ufreek2.le.ufreek1)) then
                                   nfree = 0
                                   exit
                                end if
                             else
                                ufreek2 = ubr(kk2)*free1(kk1)+vbr(kk2)*free2(kk1)
                                if (ufreek1.ge.0.and.ufreek2.ge.0.and.ufreek2.ge.ufreek1) then
                                   nfree = 0
                                   exit
                                end if
                             end if
                          end do
                       end if
                       if (nfree.eq.0) then
                          angfree1(kk1) = 0.
                          angfree2(kk1) = 0.
                          ufree(kk1) = 0.
                          dis = sqrt(ubr(kk1)**2.+vbr(kk1)**2.)
                          angmutual1(kk1) = ubr(kk1)/dis
                          angmutual2(kk1) = vbr(kk1)/dis
                          if (nbound2.eq.1) then
                             umutual(kk1) = 0.
                          else
                             umutual(kk1) = (ubr(kk1)*angmutual1(kk1)+vbr(kk1)*angmutual2(kk1))/dble(nsize(k))
                          end if
                       else
                          angfree1(kk1) = free1(kk1)
                          angfree2(kk1) = free2(kk1)
                          ufree(kk1) = ubr(kk1)*angfree1(kk1)+vbr(kk1)*angfree2(kk1)
                          angmutual1(kk1) = -free2(kk1)
                          angmutual2(kk1) =  free1(kk1)
                          if (nbound2.eq.1) then
                             umutual(kk1) = 0.
                          else
                             umutual(kk1) = (ubr(kk1)*angmutual1(kk1)+vbr(kk1)*angmutual2(kk1))/dble(nsize(k))
                          end if
                       end if
                       !-------------------------------------------------------------------------------------

                    else
                       !-------------------------for the crystl who hit more than two other crystal-----------------------
                       dmove = free1(kk1)*(xo(kk1)-xo(nmcoll(1)))+free2(kk1)*(yo(kk1)-yo(nmcoll(1)))
                       if (dmove.lt.0.) then
                          angfree1(kk1) = 0.
                          angfree2(kk1) = 0.
                          ufree(kk1) = 0.
                          dis = sqrt(ubr(kk1)**2.+vbr(kk1)**2.)
                          angmutual1(kk1) = ubr(kk1)/dis
                          angmutual2(kk1) = vbr(kk1)/dis
                          if (nbound2.eq.1) then
                             umutual(kk1) = 0.
                          else
                             umutual(kk1) = (ubr(kk1)*angmutual1(kk1)+vbr(kk1)*angmutual2(kk1))/dble(nsize(k))
                          end if
                       else
                          dis = sqrt(free1(kk1)**2.+free2(kk1)**2.)
                          free1(kk1) = free1(kk1)/dis
                          free2(kk1) = free2(kk1)/dis
                          nfree = 1
                          do kkk = 1,ncoll
                             kk2 = nmcoll(kkk)
                             ufreek1 = ubr(kk1)*free1(kk1)+vbr(kk1)*free2(kk1)
                             if (ufreek1.lt.0.) then
                                nfree = 0
                                exit
                             else
                                if (free1(kk1)*(xo(kk1)-xo(kk2))+free2(kk1)*(yo(kk1)-yo(kk2)).le.0.) then
                                   ufreek2 = ubr(kk2)*free1(kk1)+vbr(kk2)*free2(kk1)
                                   if ((ufreek1.ge.0.and.ufreek2.le.0).or.&
                                       (ufreek1.ge.0.and.ufreek2.ge.0.and.ufreek2.le.ufreek1).or.&
                                       (ufreek1.le.0.and.ufreek2.le.0.and.ufreek2.le.ufreek1)) then
                                      nfree = 0
                                      exit
                                   end if
                                else
                                   ufreek2 = ubr(kk2)*free1(kk1)+vbr(kk2)*free2(kk1)
                                   if (ufreek1.ge.0.and.ufreek2.ge.0.and.ufreek2.ge.ufreek1) then
                                      nfree = 0
                                      exit
                                   end if
                                end if
                             end if
                          end do
                          if (nfree.eq.0) then
                             angfree1(kk1) = 0.
                             angfree2(kk1) = 0.
                             ufree(kk1) = 0.
                             dis = sqrt(ubr(kk1)**2.+vbr(kk1)**2.)
                             angmutual1(kk1) = ubr(kk1)/dis
                             angmutual2(kk1) = vbr(kk1)/dis
                             if (nbound2.eq.1) then
                                umutual(kk1) = 0.
                             else
                                umutual(kk1) = (ubr(kk1)*angmutual1(kk1)+vbr(kk1)*angmutual2(kk1))/dble(nsize(k))
                             end if
                          else
                             angfree1(kk1) = free1(kk1)
                             angfree2(kk1) = free2(kk1)
                             ufree(kk1) = ubr(kk1)*angfree1(kk1)+vbr(kk1)*angfree2(kk1)
                             angmutual1(kk1) = -free2(kk1)
                             angmutual2(kk1) =  free1(kk1)
                             if (nbound2.eq.1) then
                                umutual(kk1) = 0.
                             else
                                umutual(kk1) = (ubr(kk1)*angmutual1(kk1)+vbr(kk1)*angmutual2(kk1))/dble(nsize(k))
                             end if
                          end if
                       end if
                       !-----------------------------------------------------------------------------------------
                    end if

                 end do

                 do kk = 1,nsize(k)
                    kk1 = nmcluster(kk,k)
                    ubr(kk1) = ufree(kk1)*angfree1(kk1)
                    vbr(kk1) = ufree(kk1)*angfree2(kk1)
                    do kkk = 1,nsize(k)
                       kk2 = nmcluster(kkk,k)
                       ubr(kk1) = ubr(kk1) + umutual(kk2)*angmutual1(kk2)
                       vbr(kk1) = vbr(kk1) + umutual(kk2)*angmutual2(kk2)
                    end do
                 end do
              end if
           end if
        end do
    
        !write(*,*) 3000,n2n,n3n,n4n,n5n,nfree

      end subroutine collid_strategy

!--------------------------------------------------------------------------------------------------------------------
!--------------------------------------------------------------------------------------------------------------------

     subroutine solid_move(xo,yo,dt,ubr,vbr,ns)
        implicit real*8(a-h,o-z)

        dimension :: xo(ns),yo(ns),ubr(ns),vbr(ns)

        do k = 1,ns
           xo(k) = xo(k) + ubr(k)*dt
           yo(k) = yo(k) + vbr(k)*dt
        end do

     end subroutine solid_move

!---------------------------------------------------------------------------------
!---------------------------------------------------------------------------------
      subroutine solve_solid(u,v,dt,rdx,xo,yo,ubr,vbr,wbr,ixmn,ixmx,iymn,iymx,&
                            nmcluster,nsize,ncluster,ninlet,phi_solid)

        implicit real*8(a-h,o-z)
        common/param/g,sp,ubc,uout,mint
        common/grid/dx,dy,xmx,xmn,ymx,ymn,nx,ny,tt,tb,crmx,ns,nsp,nsl
        dimension :: u(nx,ny),v(nx,ny),xo(ns),yo(ns),&
                     ubr(ns),vbr(ns),wbr(ns),&
                     ixmn(ns),ixmx(ns),iymn(ns),iymx(ns),&
                     nmcluster(ns,ns),nsize(ns),phi_solid(nx,ny,3)

        allocatable :: tmp(:,:),s_part(:,:),phivf(:,:,:),ndone(:),angn1(:),angn2(:),un(:),&
                       nxbound(:),nybound(:),nwall_l(:),nwall_r(:),nmcoll(:),angmutual1(:),angmutual2(:),&
                       umutual(:),angfree1(:),angfree2(:),ufree(:),n_s(:),&
                       nbound1(:),free1(:),free2(:),x(:,:,:),y(:,:,:)
        character (len=3) num,num0,nss,nsc
        eps = 10E-15

        mumf = 9*nm
        allocate(tmp(nx,ny),s_part(nx,ny),phivf(nx,ny,3),ndone(ns),angn1(ns),angn2(ns),un(ns),&
                 nxbound(ns),nybound(ns),nwall_l(ns),nwall_r(ns),nmcoll(ns),angmutual1(ns),angmutual2(ns),&
                 umutual(ns),angfree1(ns),angfree2(ns),ufree(ns),n_s(ns),&
                 nbound1(ns),free1(ns),free2(ns),x(nx,3,3),y(ny,3,3))

        pii = 4.*atan(1.d0)
        Area = pii*rdx**2.
        Area1 = pii*rdx**4.

        do i = 1,nx
           x(i,1,1) = dx*float(i-1)+xmn                !coordinate for bubble on level 1
           x(i,2,1) = dx*float(i-1)+xmn                !coordinate for particle on level 1
           x(i,1,2) = dx*float(i-1)+xmn                !coordinate for volume fraction for u on level 1
           x(i,2,2) = dx*float(i-1)+xmn-dx/2.          !coordinate for volume fraction for v on level 1
           x(i,3,2) = dx*float(i-1)+xmn-dx/2.          !coordinate for volume fraction for w on level 1
           x(i,1,3) = dx*float(i-1)+xmn+dx/2.          !coordinate for u update on level 1
           x(i,2,3) = dx*float(i-1)+xmn                !coordinate for v update on level 1
           x(i,3,3) = dx*float(i-1)+xmn
        end do
        do j = 1,ny
           y(j,1,1) = dy*float(j-1)+ymn                !coordinate for bubble on level 1
           y(j,2,1) = dy*float(j-1)+ymn                !coordinate for parkicle on level 1
           y(j,1,2) = dy*float(j-1)+ymn-dy/2.          !coordinate for volume fraction for u on level 1
           y(j,2,2) = dy*float(j-1)+ymn                !coordinate for volume fraction for v on level 1
           y(j,3,2) = dy*float(j-1)+ymn-dy/2.          !coordinate for volume fraction for w on level 1
           y(j,1,3) = dy*float(j-1)+ymn                !coordinate for u update on level 1
           y(j,2,3) = dy*float(j-1)+ymn+dy/2.          !coordinate for v update on level 1
           y(j,3,3) = dy*float(j-1)+ymn
        end do

        n_s_n = ns
        do k = 1,ns
           n_s(k) = k
        end do
!        write(nss,'(i3)') nsp+100
!        open(10,file="data_phivf" // nss // ".dat")
        do k1 = 1,n_s_n
           do k2 = 1,3
              do j = iymn(k1),iymx(k1)
                 do i = ixmn(k1),ixmx(k1)
                    s_part(i,j) = - rdx + sqrt((x(i,k2,2)-xo(n_s(k1)))**2. + (y(j,k2,2)-yo(n_s(k1)))**2.)
                    phivf(i,j,k2) = 0.
                 end do
              end do
              !-------------volume fraction for particle------------
              call volume_frac(s_part(:,:),phivf(:,:,k2),ixmn(k1),ixmx(k1),iymn(k1),iymx(k1),&
                               x(:,k2,2),y(:,k2,2),rdx,xo(n_s(k1)),yo(n_s(k1)),k1)
           end do


           ubr(n_s(k1)) = 0.
           vbr(n_s(k1)) = 0.
           wbr(n_s(k1)) = 0.
           do j = iymn(k1)+1,iymx(k1)
              do i = ixmn(k1)+1,ixmx(k1)
                 ubr(n_s(k1)) = ubr(n_s(k1)) + phivf(i,j,1)*u(i,j)*dx*dy/Area
                 vbr(n_s(k1)) = vbr(n_s(k1)) + phivf(i,j,2)*v(i,j)*dx*dy/Area
              end do
           end do





        end do


        !--------------------------collision strategy for boundary-crystal case---------------------------
        nxbound = 0
        nybound = 0
        do k = 1,ns
           if (floor((xo(k)+1.*rdx-xmn)/dx)+2+2.ge.nx) then
              nxbound(k) = 2
           end if
           if (floor((xo(k)-1.*rdx-xmn)/dx)-1-2.le.1) then
              nxbound(k) = 1
           end if
           if (floor((yo(k)+1.*rdx-ymn)/dy)+2+2.ge.ny) then
              nybound(k) = 2
           end if
           if (floor((yo(k)-1.*rdx-ymn)/dy)-1-2.le.1) then
              nybound(k) = 1
           end if
        end do

        nwall_l = 0
        nwall_r = 0

!        write(*,*) 2000

        do k = 1,ns

           !write(*,*) 2000, nxbound(k),nybound(k)

        end do


        do k = 1,ncluster
           if (nsize(k).eq.1) then

              !-----------------------for the single-crystal cluster, change the direction of the velocity when it hits the wall------------------------
              kk1 = nmcluster(1,k)
              if (nxbound(kk1).eq.2) then
                 ubr(kk1) = -abs(ubr(kk1))
              elseif (nxbound(kk1).eq.1) then
                 ubr(kk1) =  abs(ubr(kk1))
              end if
              if (nybound(kk1).eq.2) then
                 vbr(kk1) = -abs(vbr(kk1))
              elseif (nybound(kk1).eq.1) then
                 vbr(kk1) =  abs(vbr(kk1))
              end if
              !-----------------------------------------------------------------------------------------------------------------------------------------

           else
              if (nsize(k).eq.2) then
 
                 !----------------for the two-crystals cluster, consider the effect from both wall and the interaction between each other---------------
                 nbound1 = 0
                 do kk = 1,nsize(k)
                    kk1 = nmcluster(kk,k)
                    do kkk = 1,nsize(k)
                       kkkk = nmcluster(kkk,k)
                       if (kkkk.ne.kk1) then
                          kk2 = kkkk
                       end if
                    end do
                    dis1 = sqrt((xo(kk1)-xo(kk2))**2.+(yo(kk1)-yo(kk2))**2.)+eps
                    !-----------------free direction for each crystal----------
                    if (nxbound(kk1).eq.1) then
                       nbound1(kk1) = 1
                    end if
                    if (nxbound(kk1).eq.2) then
                       nbound1(kk1) = 2
                    end if
                    if (nybound(kk1).eq.1) then
                       nbound1(kk1) = 3
                    end if
                    if (nybound(kk1).eq.2) then
                       nbound1(kk1) = 4
                    end if
                    !-------------------------------------------------------
                 end do
                 !----------------------check whether there is any crystal in the cluster hit the wall---------------
                 nbound2_b = 0
                 nbound2_t = 0
                 do kk = 1,nsize(k)
                    kk1 = nmcluster(kk,k)
                    if (nbound1(kk1).eq.1) then
                       nbound2_b = 1
                    end if
                    if (nbound1(kk1).eq.2) then
                       nbound2_t = 1
                    end if
                 end do
                 if (nbound2_b.eq.1.and.nbound2_t.eq.0) then
                    do kk = 1,nsize(k)
                       kk1 = nmcluster(kk,k)
                       ubr(kk1) = abs(ubr(kk1))
                    end do
                 elseif (nbound2_b.eq.0.and.nbound2_t.eq.1) then
                    do kk = 1,nsize(k)
                       kk1 = nmcluster(kk,k)
                       ubr(kk1) = -abs(ubr(kk1))
                    end do
                 elseif (nbound2_b.eq.1.and.nbound2_t.eq.1) then
                    do kk = 1,nsize(k)
                       kk1 = nmcluster(kk,k)
                       ubr(kk1) = 0.
                    end do
                 end if

                 nbound2_l = 0
                 nbound2_r = 0
                 do kk = 1,nsize(k)
                    kk1 = nmcluster(kk,k)
                    if (nbound1(kk1).eq.3) then
                       nbound2_l = 1
                    end if
                    if (nbound1(kk1).eq.4) then
                       nbound2_r = 1
                    end if
                 end do
                 if (nbound2_l.eq.1.and.nbound2_r.eq.0) then
                    do kk = 1,nsize(k)
                       kk1 = nmcluster(kk,k)
                       vbr(kk1) = abs(vbr(kk1))
                    end do
                 elseif (nbound2_l.eq.0.and.nbound2_r.eq.1) then
                    do kk = 1,nsize(k)
                       kk1 = nmcluster(kk,k)
                       vbr(kk1) = -abs(vbr(kk1))
                    end do
                 elseif (nbound2_l.eq.1.and.nbound2_r.eq.1) then
                    do kk = 1,nsize(k)
                       kk1 = nmcluster(kk,k)
                       vbr(kk1) = 0.
                    end do
                 end if
                 !---------------------------------------------------------------------------------------------------
                 do kk = 1,nsize(k)
                    kk1 = nmcluster(kk,k)
                    do kkk = 1,nsize(k)
                       kkkk = nmcluster(kkk,k)
                       if (kkkk.ne.kk1) then
                          kk2 = kkkk
                       end if
                    end do
                    dis = sqrt((xo(kk1)-xo(kk2))**2.+(yo(kk1)-yo(kk2))**2.)+eps
                    angfree1(kk1) = -(yo(kk1)-yo(kk2))/dis
                    angfree2(kk1) = (xo(kk1)-xo(kk2))/dis
                    ufree(kk1) = ubr(kk1)*angfree1(kk1)+vbr(kk1)*angfree2(kk1)
                    angmutual1(kk1) = (xo(kk1)-xo(kk2))/dis
                    angmutual2(kk1) = (yo(kk1)-yo(kk2))/dis
                    umutual(kk1) = (ubr(kk1)*angmutual1(kk1)+vbr(kk1)*angmutual2(kk1))/dble(nsize(k))
                 end do
                 do kk = 1,nsize(k)
                    kk1 = nmcluster(kk,k)
                    ubr(kk1) = ufree(kk1)*angfree1(kk1)
                    vbr(kk1) = ufree(kk1)*angfree2(kk1)
                    do kkk = 1,nsize(k)
                       kk2 = nmcluster(kkk,k)
                       ubr(kk1) = ubr(kk1) + umutual(kk2)*angmutual1(kk2)
                       vbr(kk1) = vbr(kk1) + umutual(kk2)*angmutual2(kk2)
                    end do
                 end do
                 !---------------------------------------------------------------------------------------------------------------------------
              else
                 !--------------------------------------------for the multi_crystals cluster-------------------------------------------------
                 !-----------------------------------look for the crystals hiting the wall-------------------------------------
                 nbound1 = 0
                 do kk = 1,nsize(k)
                    kk1 = nmcluster(kk,k)
                    if (nxbound(kk1).eq.1) then
                       nbound1(kk1) = 1
                    end if
                    if (nxbound(kk1).eq.2) then
                       nbound1(kk1) = 2
                    end if
                    if (nybound(kk1).eq.1) then
                       nbound1(kk1) = 3
                    end if
                    if (nybound(kk1).eq.2) then
                       nbound1(kk1) = 4
                    end if
                 end do
                 !-----------------------------------------------------------------------------------------------------------
                 !-------------------------------check whether there is any crystal hiting the wall--------------------------
                 !----------------------check whether there is any crystal in the cluster hit the wall---------------
                 nbound2_b = 0
                 nbound2_t = 0
                 do kk = 1,nsize(k)
                    kk1 = nmcluster(kk,k)
                    if (nbound1(kk1).eq.1) then
                       nbound2_b = 1
                    end if
                    if (nbound1(kk1).eq.2) then
                       nbound2_t = 1
                    end if
                 end do
                 if (nbound2_b.eq.1.and.nbound2_t.eq.0) then
                    do kk = 1,nsize(k)
                       kk1 = nmcluster(kk,k)
                       ubr(kk1) = abs(ubr(kk1))
                    end do
                 elseif (nbound2_b.eq.0.and.nbound2_t.eq.1) then
                    do kk = 1,nsize(k)
                       kk1 = nmcluster(kk,k)
                       ubr(kk1) = -abs(ubr(kk1))
                    end do
                 elseif (nbound2_b.eq.1.and.nbound2_t.eq.1) then
                    do kk = 1,nsize(k)
                       kk1 = nmcluster(kk,k)
                       ubr(kk1) = 0.
                    end do
                 end if

                 nbound2_l = 0
                 nbound2_r = 0
                 do kk = 1,nsize(k)
                    kk1 = nmcluster(kk,k)
                    if (nbound1(kk1).eq.3) then
                       nbound2_l = 1
                    end if
                    if (nbound1(kk1).eq.4) then
                       nbound2_r = 1
                    end if
                 end do
                 if (nbound2_l.eq.1.and.nbound2_r.eq.0) then
                    do kk = 1,nsize(k)
                       kk1 = nmcluster(kk,k)
                       vbr(kk1) = abs(vbr(kk1))
                    end do
                 elseif (nbound2_l.eq.0.and.nbound2_r.eq.1) then
                    do kk = 1,nsize(k)
                       kk1 = nmcluster(kk,k)
                       vbr(kk1) = -abs(vbr(kk1))
                    end do
                 elseif (nbound2_l.eq.1.and.nbound2_r.eq.1) then
                    do kk = 1,nsize(k)
                       kk1 = nmcluster(kk,k)
                       vbr(kk1) = 0.
                    end do
                 end if
                 !-----------------------------------------------------------------------------------------------------------

                 do kk = 1,nsize(k)
                    kk1 = nmcluster(kk,k)
                    ncoll = 0
                    do kkk = 1,nsize(k)
                       kk2 = nmcluster(kkk,k)
                       dis = sqrt((xo(kk1)-xo(kk2))**2.+(yo(kk1)-yo(kk2))**2.)+eps
                       if (kk2.ne.kk1.and.dis.le.2.*(rdx)) then
                          ncoll = ncoll + 1
                          nmcoll(ncoll) = kk2
                       end if
                    end do

                    !----------------------------compute the free direction-----------------------------------
                    free1(kk1) = 0.
                    free2(kk1) = 0.
                    do kkk = 1,ncoll
                       kk2 = nmcoll(kkk)
                       free1(kk1) = free1(kk1) + xo(kk1)-xo(kk2)
                       free2(kk1) = free2(kk1) + yo(kk1)-yo(kk2)
                    end do
                    dis1 = sqrt(free1(kk1)**2.+free2(kk1)**2.)+eps
                    !----------------------------------------------------------------------------------------

                    !---------------------------for the crystal who only hit one other crystal---------------
                    if (ncoll.eq.1) then
                       kk2 = nmcoll(1)
                       dis = sqrt((xo(kk1)-xo(kk2))**2.+(yo(kk1)-yo(kk2))**2.)+eps
                       angfree1(kk1) = -(yo(kk1)-yo(kk2))/dis
                       angfree2(kk1) = (xo(kk1)-xo(kk2))/dis
                       ufree(kk1) = ubr(kk1)*angfree1(kk1)+vbr(kk1)*angfree2(kk1)
                       angmutual1(kk1) = (xo(kk1)-xo(kk2))/dis
                       angmutual2(kk1) = (yo(kk1)-yo(kk2))/dis
                       umutual(kk1) = (ubr(kk1)*angmutual1(kk1)+vbr(kk1)*angmutual2(kk1))/dble(nsize(k))
                    !----------------------------------------------------------------------------------------

                    elseif (ncoll.eq.2) then
                       !-------------------------for the crystl who hit two other crystal-----------------------
                       nfree = 1
                       dis = sqrt(free1(kk1)**2.+free2(kk1)**2.)+eps
                       free1(kk1) = free1(kk1)/dis
                       free2(kk1) = free2(kk1)/dis
                       ufreek1 = ubr(kk1)*free1(kk1)+vbr(kk1)*free2(kk1)
                       if (ufreek1.lt.0.) then
                          nfree = 0
                       else
                          do kkk = 1,ncoll
                             kk2 = nmcoll(kkk)
                             if (free1(kk1)*(xo(kk1)-xo(kk2))+free2(kk1)*(yo(kk1)-yo(kk2)).le.0.) then
                                ufreek2 = ubr(kk2)*free1(kk1)+vbr(kk2)*free2(kk1)
                                if ((ufreek1.ge.0.and.ufreek2.le.0).or.(ufreek1.ge.0.and.ufreek2.ge.0.and.ufreek2.le.ufreek1).or.&
                                   (ufreek1.le.0.and.ufreek2.le.0.and.ufreek2.le.ufreek1)) then
                                   nfree = 0
                                   exit
                                end if
                             else
                                ufreek2 = ubr(kk2)*free1(kk1)+vbr(kk2)*free2(kk1)
                                if (ufreek1.ge.0.and.ufreek2.ge.0.and.ufreek2.ge.ufreek1) then
                                   nfree = 0
                                   exit
                                end if
                             end if
                          end do
                       end if
                       if (nfree.eq.0) then
                          angfree1(kk1) = 0.
                          angfree2(kk1) = 0.
                          ufree(kk1) = 0.
                          dis = sqrt(ubr(kk1)**2.+vbr(kk1)**2.)+eps
                          angmutual1(kk1) = ubr(kk1)/dis
                          angmutual2(kk1) = vbr(kk1)/dis
                          umutual(kk1) = (ubr(kk1)*angmutual1(kk1)+vbr(kk1)*angmutual2(kk1))/dble(nsize(k))
                       else
                          angfree1(kk1) = free1(kk1)
                          angfree2(kk1) = free2(kk1)
                          ufree(kk1) = ubr(kk1)*angfree1(kk1)+vbr(kk1)*angfree2(kk1)
                          angmutual1(kk1) = -free2(kk1)
                          angmutual2(kk1) =  free1(kk1)
                          umutual(kk1) = (ubr(kk1)*angmutual1(kk1)+vbr(kk1)*angmutual2(kk1))/dble(nsize(k))
                       end if
                       !-------------------------------------------------------------------------------------
                    else
                       !-------------------------for the crystl who hit more than two other crystal-----------------------
                       ndmove = 0
                       do kkk = 1,ncoll
                          kk2 = nmcoll(kkk)
                          if (free1(kk1)*(xo(kk1)-xo(kk2))+free2(kk1)*(yo(kk1)-yo(kk2)).lt.0) then
                             ndmove = 1
                          end if
                       end do
                       dmove = free1(kk1)*(xo(kk1)-xo(nmcoll(1)))+free2(kk1)*(yo(kk1)-yo(nmcoll(1)))
                       if (ndmove.eq.1) then
                          angfree1(kk1) = 0.
                          angfree2(kk1) = 0.
                          ufree(kk1) = 0.
                          dis = sqrt(ubr(kk1)**2.+vbr(kk1)**2.)+eps
                          angmutual1(kk1) = ubr(kk1)/dis
                          angmutual2(kk1) = vbr(kk1)/dis
                          umutual(kk1) = (ubr(kk1)*angmutual1(kk1)+vbr(kk1)*angmutual2(kk1))/dble(nsize(k))
                       else
                          dis = sqrt(free1(kk1)**2.+free2(kk1)**2.)+eps
                          free1(kk1) = free1(kk1)/dis
                          free2(kk1) = free2(kk1)/dis
                          nfree = 1
                          do kkk = 1,ncoll
                             kk2 = nmcoll(kkk)
                             ufreek1 = ubr(kk1)*free1(kk1)+vbr(kk1)*free2(kk1)
                             if (ufreek1.lt.0.) then
                                nfree = 0
                                exit
                             else
                                if (free1(kk1)*(xo(kk1)-xo(kk2))+free2(kk1)*(yo(kk1)-yo(kk2)).le.0.) then
                                   ufreek2 = ubr(kk2)*free1(kk1)+vbr(kk2)*free2(kk1)
                                   if ((ufreek1.ge.0.and.ufreek2.le.0).or.&
                                       (ufreek1.ge.0.and.ufreek2.ge.0.and.ufreek2.le.ufreek1).or.&
                                       (ufreek1.le.0.and.ufreek2.le.0.and.ufreek2.le.ufreek1)) then
                                      nfree = 0
                                      exit
                                   end if
                                else
                                   ufreek2 = ubr(kk2)*free1(kk1)+vbr(kk2)*free2(kk1)
                                   if (ufreek1.ge.0.and.ufreek2.ge.0.and.ufreek2.ge.ufreek1) then
                                      nfree = 0
                                      exit
                                   end if
                                end if
                             end if
                          end do
                          if (nfree.eq.0) then
                             angfree1(kk1) = 0.
                             angfree2(kk1) = 0.
                             ufree(kk1) = 0.
                             dis = sqrt(ubr(kk1)**2.+vbr(kk1)**2.)+eps
                             angmutual1(kk1) = ubr(kk1)/dis
                             angmutual2(kk1) = vbr(kk1)/dis
                             umutual(kk1) = (ubr(kk1)*angmutual1(kk1)+vbr(kk1)*angmutual2(kk1))/dble(nsize(k))
                          else
                             angfree1(kk1) = free1(kk1)
                             angfree2(kk1) = free2(kk1)
                             ufree(kk1) = ubr(kk1)*angfree1(kk1)+vbr(kk1)*angfree2(kk1)
                             angmutual1(kk1) = -free2(kk1)
                             angmutual2(kk1) =  free1(kk1)
                             umutual(kk1) = (ubr(kk1)*angmutual1(kk1)+vbr(kk1)*angmutual2(kk1))/dble(nsize(k))
                          end if
                       end if
                       !-----------------------------------------------------------------------------------------
                    end if

                 end do

                 do kk = 1,nsize(k)
                    kk1 = nmcluster(kk,k)
                    ubr(kk1) = ufree(kk1)*angfree1(kk1)
                    vbr(kk1) = ufree(kk1)*angfree2(kk1)
                    do kkk = 1,nsize(k)
                       kk2 = nmcluster(kkk,k)
                       ubr(kk1) = ubr(kk1) + umutual(kk2)*angmutual1(kk2)
                       vbr(kk1) = vbr(kk1) + umutual(kk2)*angmutual2(kk2)
                    end do
                 end do
              end if
           end if
        end do

        !write(*,*) 1000,ubr(1),vbr(1)
  
        !write(*,*) 3000,n2n,n3n,n4n,n5n,nfree
  
        phi_solid = 100.
!        write(*,*) 3000
        do kk = 1,3
           do k = 1,n_s_n
              do j = iymn(k),iymx(k)
                 do i = ixmn(k),ixmx(k)
                    !particle
                     phi_solid(i,j,kk) =  min( phi_solid(i,j,kk), &
                                - rdx + sqrt((x(i,kk,3)-xo(n_s(k)))**2. + (y(j,kk,3)-yo(n_s(k)))**2.))    !particle

                 end do
              end do
           end do
        end do



        do k = 1,n_s_n

           !write(*,*) 10000

           do j = iymn(k),iymx(k)
              do i = ixmn(k),ixmx(k)
                 s =  - rdx + sqrt((x(i,1,3)-xo(n_s(k)))**2. + (y(j,1,3)-yo(n_s(k)))**2.)       !kth particle
                 if (s.le.0.) then
                     na = 1
                 else
                     na = 0
                 end if
                 s =  - rdx + sqrt((x(i+1,1,3)-xo(n_s(k)))**2. + (y(j,1,3)-yo(n_s(k)))**2.)     !kth particle
                 if (s.le.0.) then
                     nb = 1
                 else
                     nb = 0
                 end if
                 if (na.eq.1) then
                    u(i,j) = ubr(n_s(k)) - wbr(n_s(k))*(y(j,1,3)-yo(n_s(k)))
                    if (nb.eq.0.and.phi_solid(i+1,j,1).gt.0) then
                       if ((sqrt(rdx**2.-(y(j,1,3)-yo(n_s(k)))**2.)+xo(n_s(k))).gt.x(i,1,3).and.&
                           (sqrt(rdx**2.-(y(j,1,3)-yo(n_s(k)))**2.)+xo(n_s(k))).lt.x(i+1,1,3)) then
                           sign = 1.d0
                       else
                           sign = -1.d0
                       end if
                       xs = sign*sqrt(rdx**2.-(y(j,1,3)-yo(n_s(k)))**2.)+xo(n_s(k))
                       u(i+1,j) = (x(i+1,1,3)-xs)/(x(i+2,1,3)-xs)*(u(i+2,j)-(ubr(n_s(k)) &
                                - wbr(n_s(k))*(y(j,1,3)-yo(n_s(k))))) &
                                + ubr(n_s(k)) - wbr(n_s(k))*(y(j,1,3)-yo(n_s(k)))
                    end if
                 else
                    if (nb.eq.1.and.phi_solid(i,j,1).gt.0) then
                       if ((sqrt(rdx**2.-(y(j,1,3)-yo(n_s(k)))**2.)+xo(n_s(k))).gt.x(i,1,3).and.&
                           (sqrt(rdx**2.-(y(j,1,3)-yo(n_s(k)))**2.)+xo(n_s(k))).lt.x(i+1,1,3)) then
                           sign = 1.d0
                       else
                           sign = -1.d0
                       end if
                       xs = sign*sqrt(rdx**2.-(y(j,1,3)-yo(n_s(k)))**2.)+xo(n_s(k))
                       u(i,j) = (x(i,1,3)-xs)/(x(i-1,1,3)-xs)*(u(i-1,j)-(ubr(n_s(k)) &
                              - wbr(n_s(k))*(y(j,1,3)-yo(n_s(k))))) &
                              + ubr(n_s(k)) - wbr(n_s(k))*(y(j,1,3)-yo(n_s(k)))
                    end if
                 end if

                 s =  - rdx + sqrt((x(i,2,3)-xo(n_s(k)))**2. + (y(j,2,3)-yo(n_s(k)))**2.)       !kth particle
                 if (s.le.0.) then
                    nc = 1
                 else
                    nc = 0
                 end if
                 s =  - rdx + sqrt((x(i,2,3)-xo(n_s(k)))**2. + (y(j+1,2,3)-yo(n_s(k)))**2.)       !kth particle
                 if (s.le.0.) then
                    nd = 1
                 else
                    nd = 0
                 end if
                 if (nc.eq.1) then
                    v(i,j) = vbr(n_s(k)) + wbr(n_s(k))*(x(i,2,3)-xo(n_s(k)))
                    if (nd.eq.0.and.phi_solid(i,j+1,2).gt.0) then
                       if ((sqrt(rdx**2.-(x(i,2,3)-xo(n_s(k)))**2.)+yo(n_s(k))).gt.y(j,2,3).and.&
                           (sqrt(rdx**2.-(x(i,2,3)-xo(n_s(k)))**2.)+yo(n_s(k))).lt.y(j+1,2,3)) then
                           sign = 1.d0
                       else
                           sign = -1.d0
                       end if
                       ys = sign*sqrt(rdx**2.-(x(i,2,3)-xo(n_s(k)))**2.)+yo(n_s(k))
                       v(i,j+1) = (y(j+1,2,3)-ys)/(y(j+2,2,3)-ys)*(v(i,j+2)&
                                - (vbr(n_s(k)) + wbr(n_s(k))*(x(i,2,3)-xo(n_s(k))))) &
                                + vbr(n_s(k)) + wbr(n_s(k))*(x(i,2,3)-xo(n_s(k)))
                    end if
                 else
                    if (nd.eq.1.and.phi_solid(i,j,2).gt.0) then
                       if ((sqrt(rdx**2.-(x(i,2,3)-xo(n_s(k)))**2.)+yo(n_s(k))).gt.y(j,2,3).and.&
                           (sqrt(rdx**2.-(x(i,2,3)-xo(n_s(k)))**2.)+yo(n_s(k))).lt.y(j+1,2,3)) then
                            sign = 1.d0
                       else
                            sign = -1.d0
                       end if
                       ys = sign*sqrt(rdx**2.-(x(i,2,3)-xo(n_s(k)))**2.)+yo(n_s(k))
                       v(i,j) = (y(j,2,3)-ys)/(y(j-1,2,3)-ys)*(v(i,j-1)&
                              - (vbr(n_s(k)) + wbr(n_s(k))*(x(i,2,3)-xo(n_s(k))))) &
                              + vbr(n_s(k)) + wbr(n_s(k))*(x(i,2,3)-xo(n_s(k)))
                    end if
                 end if
              end do
           end do

        end do

!        write(*,*) 4000



 
10      format(2i5,10e20.10)
20      format(3i5,10e20.10)
30      format(3i5,15f13.5)


      end subroutine solve_solid

!---------------------------------------------------------------------------------
!*********************************************************************************
!
!*********************************************************************************
!---------------------------------------------------------------------------------

      subroutine volume_frac(phi,phivf,ixmn,ixmx,iymn,iymx,x,y,rdx,xo,yo,ks)
        implicit real*8(a-h,o-z)
        common/param/g,sp,ubc,uout,mint
        common/grid/dx,dy,xmx,xmn,ymx,ymn,nx,ny,tt,tb,crmx,ns,nsp,nsl
        dimension:: phi(nx,ny),phivf(nx,ny),x(nx),y(ny)
        allocatable:: phicon(:)
        allocate(phicon(4))

        eps = 1E-10
        do j=iymn+1,iymx-1
           do i=ixmn+1,ixmx-1

                 phicon(1)=phi(i  ,j)
                 phicon(2)=phi(i+1,j)
                 phicon(3)=phi(i+1,j+1)
                 phicon(4)=phi(i  ,j+1)

                 if (phicon(1).lt.0.and.phicon(2).lt.0.and.phicon(3).lt.0.and.phicon(4).lt.0) then
                    phivf(i,j)=1.
                 else if (phicon(1).gt.0.and.phicon(2).gt.0.and.phicon(3).gt.0.and.phicon(4).gt.0) then
                    phivf(i,j)=0.
                 else
                    if ((phicon(1)*phicon(2)*phicon(3)*phicon(4)).lt.0.) then
                        do k=1,4
                           n=0
                           do l=1,4
                              if (l.ne.k) then
                                 dm=phicon(k)*phicon(l)
                                 if (dm.lt.0) then
                                    n=n+1
                                 else
                                    exit
                                 end if
                              end if
                           end do
                           if (n.eq.3) then
                              exit
                           end if
                        end do

                        if (k.eq.1) then
                           if (sqrt(rdx**2.-(y(j)-yo)**2.)+xo.gt.x(i).and.sqrt(rdx**2.-(y(j)-yo)**2.)+xo.lt.x(i+1)) then
                              sign1 = 1.
                           else
                              sign1 =-1.
                           end if
                           if (sqrt(rdx**2.-(x(i)-xo)**2.)+yo.gt.y(j).and.sqrt(rdx**2.-(x(i)-xo)**2.)+yo.lt.y(j+1)) then
                              sign2 = 1.
                           else
                              sign2 =-1.
                           end if
                           a=(yo-y(j))*(sign1*sqrt(rdx**2.-(y(j)-yo)**2.)+xo-x(i))&
                             +sign2*(sign1*(abs(y(j)-yo)*sqrt(rdx**2.-(y(j)-yo)**2.)/2.+rdx**2./2.&
                             *asin(sqrt(rdx**2.-(y(j)-yo)**2.)/rdx))&
                             -((x(i)-xo)*sqrt(rdx**2.-(x(i)-xo)**2.)/2.+rdx**2./2.*asin((x(i)-xo)/rdx)))

                        else if (k.eq.2) then
                           if (sqrt(rdx**2.-(y(j)-yo)**2.)+xo.gt.x(i).and.sqrt(rdx**2.-(y(j)-yo)**2.)+xo.lt.x(i+1)) then
                              sign1 = 1.
                           else
                              sign1 =-1.
                           end if
                           if (sqrt(rdx**2.-(x(i+1)-xo)**2.)+yo.gt.y(j).and.sqrt(rdx**2.-(x(i+1)-xo)**2.)+yo.lt.y(j+1)) then
                              sign2 = 1.
                           else
                              sign2 =-1.
                           end if
                           a=(yo-y(j))*(x(i+1)-sign1*sqrt(rdx**2.-(y(j)-yo)**2.)-xo)&
                             +sign2*(((x(i+1)-xo)*sqrt(rdx**2.-(x(i+1)-xo)**2.)/2.+rdx**2./2.*asin((x(i+1)-xo)/rdx))&
                             -sign1*(abs(y(j)-yo)*sqrt(rdx**2.-(y(j)-yo)**2.)/2.+rdx**2./2.&
                             *asin(sqrt(rdx**2.-(y(j)-yo)**2.)/rdx)))

                        else if (k.eq.3) then
                           if (sqrt(rdx**2.-(y(j+1)-yo)**2.)+xo.gt.x(i).and.sqrt(rdx**2.-(y(j+1)-yo)**2.)+xo.lt.x(i+1)) then
                              sign1 = 1.
                           else
                              sign1 =-1.
                           end if
                           if (sqrt(rdx**2.-(x(i+1)-xo)**2.)+yo.gt.y(j).and.sqrt(rdx**2.-(x(i+1)-xo)**2.)+yo.lt.y(j+1)) then
                              sign2 = 1.
                           else
                              sign2 =-1.
                           end if
                           a=(y(j+1)-yo)*(x(i+1)-sign1*sqrt(rdx**2.-(y(j+1)-yo)**2.)-xo)&
                             -sign2*(((x(i+1)-xo)*sqrt(rdx**2.-(x(i+1)-xo)**2.)/2.+rdx**2./2.*asin((x(i+1)-xo)/rdx))&
                             -sign1*(abs(y(j+1)-yo)*sqrt(rdx**2.-(y(j+1)-yo)**2.)/2.+rdx**2./2.&
                             *asin(sqrt(rdx**2.-(y(j+1)-yo)**2.)/rdx)))
       
                        else if (k.eq.4) then
                           if (sqrt(rdx**2.-(y(j+1)-yo)**2.)+xo.gt.x(i).and.sqrt(rdx**2.-(y(j+1)-yo)**2.)+xo.lt.x(i+1)) then
                              sign1 = 1.
                           else
                              sign1 =-1.
                           end if
                           if (sqrt(rdx**2.-(x(i)-xo)**2.)+yo.gt.y(j).and.sqrt(rdx**2.-(x(i)-xo)**2.)+yo.lt.y(j+1)) then
                              sign2 = 1.
                           else
                              sign2 =-1.
                           end if
                           a=(y(j+1)-yo)*(sign1*sqrt(rdx**2.-(y(j+1)-yo)**2.)+xo-x(i))&
                             -sign2*(sign1*(abs(y(j+1)-yo)*sqrt(rdx**2.-(y(j+1)-yo)**2.)/2.+rdx**2./2.&
                             *asin(sqrt(rdx**2.-(y(j+1)-yo)**2.)/rdx))&
                             -((x(i)-xo)*sqrt(rdx**2.-(x(i)-xo)**2.)/2.+rdx**2./2.*asin((x(i)-xo)/rdx)))

                        end if
                        if (phicon(k).lt.0) then
                           phivf(i,j)=a/(dx*dy)
                        else
                           phivf(i,j)=1.-a/(dx*dy)
                        end if
                    else
                        if (phicon(1)*phicon(2).lt.0) then
                           if (sqrt(rdx**2.-(y(j+1)-yo)**2.)+xo.gt.x(i).and.sqrt(rdx**2.-(y(j+1)-yo)**2.)+xo.lt.x(i+1)) then
                              sign1 = 1.
                           else
                              sign1 =-1.
                           end if
                           a=sign1*(((y(j+1)-yo)*sqrt(rdx**2.-(y(j+1)-yo)**2.)/2.+rdx**2./2.*asin((y(j+1)-yo)/rdx))&
                             -((y(j)-yo)*sqrt(rdx**2.-(y(j)-yo)**2.)/2.+rdx**2./2.*asin((y(j)-yo)/rdx)))&
                             +(y(j+1)-y(j))*(xo-x(i))
                           if (phicon(1).lt.0) then
                              phivf(i,j)=a/(dx*dy)
                           else
                              phivf(i,j)=1.-a/(dx*dy)
                           end if
                        else
                           if (sqrt(rdx**2.-(x(i+1)-xo)**2.)+yo.gt.y(j).and.sqrt(rdx**2.-(x(i+1)-xo)**2.)+yo.lt.y(j+1)) then
                              sign2 = 1.
                           else
                              sign2 =-1.
                           end if
                           a=sign2*(((x(i+1)-xo)*sqrt(rdx**2.-(x(i+1)-xo)**2.)/2.+rdx**2./2.*asin((x(i+1)-xo)/rdx))&
                             -((x(i)-xo)*sqrt(rdx**2.-(x(i)-xo)**2.)/2.+rdx**2./2.*asin((x(i)-xo)/rdx)))&
                             +(x(i+1)-x(i))*(yo-y(j))
                           if (phicon(1).lt.0) then
                              phivf(i,j)=a/(dx*dy)
                           else
                              phivf(i,j)=1.-a/(dx*dy)
                           end if





                        end if
                    end if
                 end if
           end do
        end do
!        open(10,file="data_phivf_incode" // nss // "_" // mum // ".dat")
!           do j=iymn(k1),iymx(k1)
!            do i=ixmn(k1),ixmx(k1)
!              write (10,*) k1,i,j,phivf(i,j,3)
!
!!            do j=iymn+1,iymx-1
!!               do i=ixmn+1,ixmx-1
!!                 write (10,*) ks,i,j,phivf(i,j)
!!               end do
!!             end do



      end subroutine volume_frac








