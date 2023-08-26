cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
        program ENUF-DPD
c            -- simulation of simple solvent with salts using DPD method
c            -- impletation of ENUF algorithm into DPD model 
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
        implicit none
        integer Npmax, Ncharge,maxnab,Ifile,mt,stp,Lx,Ly,Lz,nbin
        integer np_salt, nn_salt

        parameter( Lx       = 30   ) ! length of simulation box
        parameter( Ly       = 30   )
        parameter( Lz       = 30   )
        parameter( Npmax    = 108000 ) ! number density of beads
        parameter( maxnab   = 100*Npmax )
        parameter( Ifile    = 0    )
        parameter( mt       = 3    )
        parameter( stp      = 500  )
        parameter(np_salt   = 20   ) ! AlCl3 salt
        parameter(chp_salt  = 3.0  )
        parameter(nn_salt   = 60   )
        parameter(chn_salt  = -1.0 )
        parameter( Ncharge  = 80   )
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
        double precision charge(Npmax), alpha_chg
        integer idcharge(Ncharge), kspace, n_flag(Npmax)

c       parameters for link-list cell
        integer M, ncell, mapsize ! 
        parameter(M=Lx, ncell=M*M*M, mapsize=13*ncell)
        integer head(ncell), map(mapsize), list(Npmax)

c       params for coordiante, velocity, and force
        double precision
     +  	  rX(Npmax),     rY(Npmax),     rZ(Npmax),
     +  	  rX0(Npmax),    rY0(Npmax),    rZ0(Npmax),
     +  	  vX(Npmax),     vY(Npmax),     vZ(Npmax),
     +  	  voX(Npmax),    voY(Npmax),    voZ(Npmax),
     +  	  FcX(Npmax),    FcY(Npmax),    FcZ(Npmax),
     +  	  FdX(Npmax),    FdY(Npmax),    FdZ(Npmax),
     +  	  FrX(Npmax),    FrY(Npmax),    FrZ(Npmax),
     +  	  FrealX(Npmax), FrealY(Npmax), FrealZ(Npmax),
     +  	  FrecipX(Npmax),FrecipY(Npmax),FrecipZ(Npmax),
     +  	  FbondX(Npmax), FbondY(Npmax), FbondZ(Npmax),
     +  	  vpX(Npmax),    vpY(Npmax),    vpZ(Npmax)

c       params for forces in DPD model
        double precision alpha(3,3), gamma, sigma

        integer Np, Nm, len_chain, n_charge,maxstp,incstp,
     +    	isteps,issstp,msdstp,rstp, n_ci

        double precision 
     +  	  dLx,dLy,dLz,dLxINV,dLyINV,dLzINV,
     +  	  dlen_chain,dlen_chainINV,
     +  	  rho, dNp, dNpINV, rfac,Tfac,
     +  	  dr_gyr_dendrimer,dr_gyr, dr_ee,
     +  	  dt,tmax,vCM,temp,pressX,pressY,pressZ,pres,
     +  	  pm, pmINV,skin, skinSQ,bond, bap

        double precision rmsd(mt),rhoA(Lx,Ly,Lz),rhoB(Lx,Ly,Lz)

        integer iseed
        real R2INIS, dummy

        integer j, ip, in, it,tt,ttt, ix,iy,iz,iii,ijk,nconf
        double precision vol, dtH, dtRTH,dtH2,dtRTH2,aa,bb,cc

        double precision real_eng, recip_eng, nfft_eng,
     +  	  real_refer_eng, recip_refer_eng

cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
        alpha_chg  = 0.35
        kspace     = 7

        rho        = 4.d0        ! total density of all particles

c      len_chain  = 48         ! size of the chain 
c        n_charge = 48
c        n_ci     = 48
c      bond       = 5.30
c      bap        = 7.d0
c      pm         = 1.d0

        alpha(1,1) = 25.d0    ! 1 for water
        alpha(2,2) = 25.d0    ! 2 for cations
        alpha(3,3) = 25.d0    ! 3 for anions
        alpha(1,2) = 25.d0
        alpha(1,3) = 25.d0
        alpha(2,3) = 25.d0

        sigma      = 3.67d0   ! dpd parameter
        iseed      = 612345   ! seed for a random number generator
        dt         = 2.0d-2   ! time step for the integrator

        isteps     = 200000  ! number of steps
        incstp     = 100     ! how often do we calculate quantities
        issstp     = 5000    ! how often do we take snapshot 

cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
        OPEN(3,STATUS='NEW', FILE='data.dat',
     +     	FORM='FORMATTED', ACCESS='SEQUENTIAL')
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
        dLx = dble(Lx)
        dLy = dble(Ly)
        dLz = dble(Lz)
        dLxINV = 1.0/dLx
        dLyINV = 1.0/dLy
        dLzINV = 1.0/dLz

c		     dlen_chain    = dble(len_chain)
c		     dlen_chainINV = 1.d0/dlen_chain

        tmax = dt*dble(maxstp)      ! Time of the production run

        vol    = dLx*dLy*dLz        ! Volume (area) of the system
        Np     = nint( rho*vol )    ! Total number of particles
        dNp    = dble(Np)
        dNpINV = 1.0/dNp            ! Inverse particle number

        pmINV  = 1.0/pm             ! Inverse particle mass
        gamma  = 0.50*sigma*sigma   ! Strength of the dissipative force
                                    ! (using fluctuation-dissipation theorem)
        dtH   = 0.50*dt
        dtRTH = 0.50*dsqrt(dt)

        skin   = 1.250              ! Skin radius for the Verlet neighbor table
        skinSQ = skin*skin

        rfac   = dsqrt(3.0)         ! Scaling of random forces
        Tfac   = 1.d0/(3.0*dNp-3.0) ! Factor for temperature control

        dummy  = R2INIS( iseed )    ! Initialize the random number generator

cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c                                                                      c
c                         Start DPD simulation                         c
c                                                                      c
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

        do it = 1, isteps
            if (it .eq. 1)then
                call MAPS(M, ncell, mapsize, head, map)
                call GENCON(Npmax,Np,len_chain,n_flag,charge,idcharge,
     +  			      Ncharge,rX,rY,rZ,rX0,rY0,rZ0,vX,vY,vZ,dLx,dLy,dLz,
     +  		        dLxINV,dLyINV,dLzINV,dNpINV,pm,pmINV,Tfac,n_charge,
     +              n_ci,np_salt,nn_salt)
                call FElecRecipENUF_INIT(kspace,Lx,Ly,Lz,Ncharge)
            else
                call INTV(Npmax,Np,vX,vY,vZ,FcX,FcY,FcZ,FdX,FdY,FdZ,
     +      		    FrX,FrY,FrZ,pmINV,dtH,dtRTH)
                call INTR(Npmax,Np,rX,rY,rZ,rX0,rY0,rZ0,vX,vY,vZ,dt,
     +      		    dLx,dLy,dLz,dLxINV,dLyINV,dLzINV)
            endif

            call LINKS(Npmax,Np,M,ncell,rX,rY,rZ,dLx,dLy,dLz,head,list)

            call FCDR(Npmax,Np,ncell,mapsize,head,map,list,rX,rY,rZ,
     +          voX,voY,voZ,FcX,FcY,FcZ,FdX,FdY,FdZ,FrX,FrY,FrZ,
     +          vpX,vpY,vpZ,n_flag,alpha,dLx,dLy,dLz,
     +          dLxINV,dLyINV,dLzINV,rfac,gamma,sigma,incstp)
c            call FBondAngle(Npmax,Np,len_chain,rX,rY,rZ,
c     +          FbondX,FbondY,FbondZ,vpX, vpY, vpZ, n_flag,dLx,
c     +          dLy,dLz,dLxINV,dLyINV,dLzINV,bond)
ccc   		  Electrostatic interactions between charges
            call FElecReal(Npmax,Np,Ncharge,rX,rY,rZ,vpX,vpY,vpZ,
     +          charge,idcharge,alpha_chg,dLx,dLy,dLz,dLxINV,dLyINV,
     +          dLzINV,FrealX,FrealY,FrealZ,real_eng)
            call FElecRecipENUF(kspace, Npmax, Ncharge, alpha_chg,
     +  		    dLxINV, dLyINV, dLzINV,idcharge, charge,
     +          rX, rY, rZ, FrecipX, FrecipY, FrecipZ, nfft_eng)
            call ALLFORCE(Npmax,Np,FcX,FcY,FcZ,FbondX, FbondY, FbondZ,
     +          FrealX,FrealY,FrealZ,FrecipX,FrecipY,FrecipZ)

            call INTV(Npmax,Np,vX,vY,vZ,FcX,FcY,FcZ,
     +          FdX,FdY,FdZ,FrX,FrY,FrZ,pmINV,dtH,dtRTH)

            if (mod(it,incstp).eq.0) then
                call VCMTP(Npmax,Np,vX,vY,vZ,dNpINV,pm,pmINV,Tfac,vCM,
     +  	    	    temp,vpX,vpY,vpZ,vol,pressX,pressY,pressZ,pres)
                write(3,9)it,sngl(dble(it)*dt),vCM,sngl(temp),sngl(pres)
                call FLUSH(3)
            endif
 9 		      format(i10, 1x,f10.4, 1x,d12.6, 4(1x,f14.6))

            if (mod(it,issstp).eq.0) then
                call SNAPSHOT(Npmax,Np,it,n_flag,charge,dLx,dLy,dLz,
     +   	    	   dLxINV,dLyINV,dLzINV,rX,rY,rZ,Ifile)
            end if

        end do

        call FElecRecipENUF_FINALIZE()

        close(3)
        stop 'Good luck!'
        end

cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c                                                                      c
c                        End of DPD simulations                        c
c                                                                      c
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
        subroutine GENCON(Npmax,Np,len_chain,n_flag,charge,idcharge,
     +             Ncharge,rX,rY,rZ,rX0,rY0,rZ0,vX,vY,vZ,dLx,dLy,dLz,
     +             dLxINV,dLyINV,dLzINV,dNpINV,pm,pmINV,Tfac,
     +             n_charge,n_ci,np_salt,nn_salt)
        implicit none
        integer Npmax,Np,n_charge,len_chain, n_ci,np_salt,nn_salt,
     +      Ncharge, idcharge(Ncharge), n_flag(Npmax)
        double precision
     +      rX(Npmax), rY(Npmax),rZ(Npmax), rX0(Npmax), rY0(Npmax),
     +      rZ0(Npmax),vX(Npmax), vY(Npmax),vZ(Npmax),charge(Npmax),
     +      dLx,dLy,dLz,dLxINV,dLyINV,dLzINV,dNpINV,pm,pmINV,Tfac
        double precision
     +      drX_tmp, drY_tmp, drZ_tmp, r_dist, d_scale, rrrz,
     +      rXi, rYi, rZi, vCMX, vCMY, vCMZ, vv, tscal, temper,
     +      raaa, rbbb, rccc, rabc, ra1, rb1, rc1, rabc1
        integer iendnumber, igeneration, lipidhead, lipidend, icenter,
     +      iiii, iiiii, iii, in, ip,j,n_every
        real R2S

c ----- Generate random configuration
        do iii = 1, np_salt ! cations
            rx(iii) = (R2S()-0.5)*dLx
            ry(iii) = (R2S()-0.5)*dLy
            rz(iii) = (R2S()-0.5)*dLz
            charge(iii) = 3.0
            n_flag(iii) = 2
        end do

        do iii=np_salt+1, np_salt+nn_salt !anions
            rx(iii) = (R2S()-0.5)*dLx
            ry(iii) = (R2S()-0.5)*dLy
            rz(iii) = (R2S()-0.5)*dLz
            charge(iii) = -1.0
            n_flag(iii) = 3
        end do

        Nsolvent = Npmax - np_salt+nn_salt
        do iii = np_salt+nn_salt+1, Npmax ! solvent
            rx(iii) = (R2S()-0.5)*dLx
            ry(iii) = (R2S()-0.5)*dLy
            rz(iii) = (R2S()-0.5)*dLz
            charge(iii) = 0.0
            n_flag(iii) = 1
        end do

c ----- PBC
        do ip = 1, Np
            rXi = rX(ip)
            rYi = rY(ip)
            rZi = rZ(ip)
            rX(ip) = rXi - dLx*dnint( rXi*dLxINV )
            rY(ip) = rYi - dLy*dnint( rYi*dLyINV )
            rZ(ip) = rZi - dLz*dnint( rZi*dLzINV )
            rX0(ip) = rX(ip)
            rY0(ip) = rY(ip)
            rZ0(ip) = rZ(ip)
        enddo

c ----- Velocities for COM
        vCMX = 0.d0
        vCMY = 0.d0
        vCMZ = 0.d0
        do ip = 1, Np
            vv = dble( R2S() - 0.5 )
            vX(ip) = vv
            vCMX = vCMX + pm*vv
            vv = dble( R2S() - 0.5 )
            vY(ip) = vv
            vCMY = vCMY + pm*vv
            vv = dble( R2S() - 0.5 )
            vZ(ip) = vv
            vCMZ = vCMZ + pm*vv
        end do
  
        vCMX = vCMX*dNpINV*pmINV
        vCMY = vCMY*dNpINV*pmINV
        vCMZ = vCMZ*dNpINV*pmINV
        do ip = 1, Np
            vX(ip) = vX(ip) - vCMX
            vY(ip) = vY(ip) - vCMY
            vZ(ip) = vZ(ip) - vCMZ
        enddo

c ----- Scaling temperatures
        temper = 0.d0             ! Rescale velocities such that 
        do ip = 1, Np             ! Kinetic Temperature equals to one
          temper=temper+pm*(vX(ip)*vX(ip)+vY(ip)*vY(ip)+vZ(ip)*vZ(ip))
        end do
        temper  = Tfac*temper
        tscal = 1.d0/dsqrt( temper )

        do ip = 1, Np
            vX(ip) = vX(ip)*tscal
            vY(ip) = vY(ip)*tscal
            vZ(ip) = vZ(ip)*tscal
        end do

        iii = 0
        do in = 1, Np
            if (charge(in).ne.(0.0)) then
                iii = iii + 1
                idcharge(iii) = in
            end if
        end do

        return
        end
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
        subroutine MAPS(M,ncell,mapsize,head,map)
        implicit none
        integer M,ncell,mapsize,head(ncell),map(mapsize)
        integer ix,iy,iz,imap, icell

        icell(ix,iy,iz) = 1 + mod(ix-1+M,M) + mod(iy-1+M,M)*M
     +                      + mod(iz-1+M,M)*M*M

        do iz = 1, M
            do iy = 1, M
                do ix = 1, M
                    imap = (icell(ix,iy,iz)-1)*13
                    map(imap+1 ) = icell(ix+1, iy  , iz  )
                    map(imap+2 ) = icell(ix+1, iy+1, iz  )
                    map(imap+3 ) = icell(ix  , iy+1, iz  )
                    map(imap+4 ) = icell(ix-1, iy+1, iz  )
                    map(imap+5 ) = icell(ix+1, iy  , iz-1)
                    map(imap+6 ) = icell(ix+1, iy+1, iz-1)
                    map(imap+7 ) = icell(ix  , iy+1, iz-1)
                    map(imap+8 ) = icell(ix-1, iy+1, iz-1)
                    map(imap+9 ) = icell(ix+1, iy  , iz+1)
                    map(imap+10 )= icell(ix+1, iy+1, iz+1)
                    map(imap+11 )= icell(ix  , iy+1, iz+1)
                    map(imap+12 )= icell(ix-1, iy+1, iz+1)
                    map(imap+13 )= icell(ix  , iy  , iz+1)
                end do
            end do
        end do

        return
        end
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
       subroutine LINKS(Npmax,Np,M,ncell,rX,rY,rZ,dLx,dLy,dLz,head,list)
        implicit none
        integer Npmax,Np,M,ncell,head(ncell),list(Npmax),icell,i
        double precision rX(Npmax),rY(Npmax),rZ(Npmax),dLx,dLy,dLz,
     +      celli,cell,cellx,celly,cellz

        do icell = 1,ncell
            head(icell) = 0
        end do

        celli = dble(M)
        cellx = dLx/celli
        celly = dLy/celli
        cellz = dLz/celli
        cell  = dmin1(cellx,celly,cellz)

        if (cell.lt.1.d0) then
            stop 'cell size too small for cutoff'
        endif
      
        do i = 1, Np
            icell = 1+ int((rX(i)/dLx+0.5)*celli)
     +	             + int((rY(i)/dLy+0.5)*celli)*M
     +	             + int((rZ(i)/dLz+0.5)*celli)*M*M
            list(i) = head(icell)
            head(icell) = i
        end do

        return
        end
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
        subroutine INTV(Npmax,Np,vX,vY,vZ,FcX,FcY,FcZ,FdX,FdY,FdZ,
     +	           FrX,FrY,FrZ,pmINV,hcd,hr)
        implicit none
        integer Npmax, Np, ip
        double precision vX(Npmax),vY(Npmax),vZ(Npmax),FcX(Npmax),
     +       FcY(Npmax),FcZ(Npmax),FdX(Npmax),FdY(Npmax),FdZ(Npmax),
     +       FrX(Npmax),FrY(Npmax),FrZ(Npmax),pmINV,hcd,hr

        do ip = 1, Np
            vX(ip) = vX(ip) + pmINV*((FcX(ip)+FdX(ip))*hcd + FrX(ip)*hr)
            vY(ip) = vY(ip) + pmINV*((FcY(ip)+FdY(ip))*hcd + FrY(ip)*hr)
            vZ(ip) = vZ(ip) + pmINV*((FcZ(ip)+FdZ(ip))*hcd + FrZ(ip)*hr)
        end do

        return
        end
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
        subroutine INTR(Npmax,Np,rX,rY,rZ,rX0,rY0,rZ0,vX,vY,vZ,h,
     +             dLx,dLy,dLz,dLxINV,dLyINV,dLzINV)
        implicit none
        integer Npmax, Np, ip
        double precision rX(Npmax), rY(Npmax), rZ(Npmax), rX0(Npmax),
     +         rY0(Npmax),rZ0(Npmax),vX(Npmax), vY(Npmax), vZ(Npmax),h,
     +         dLx,dLy,dLz,dLxINV,dLyINV,dLzINV,rXi,rYi,rZi

        do ip = 1, Np
           rXi = rX(ip) + vX(ip)*h
           rYi = rY(ip) + vY(ip)*h
           rZi = rZ(iP) + vZ(ip)*h
           rX0(ip) = rX0(ip) + vX(ip)*h
           rY0(ip) = rY0(ip) + vY(ip)*h
           rZ0(ip) = rZ0(ip) + vZ(ip)*h
           rX(ip) = rXi - dLx*dnint( rXi*dLxINV )
           rY(ip) = rYi - dLy*dnint( rYi*dLyINV )
           rZ(ip) = rZi - dLz*dnint( rZi*dLzINV )
        end do

        return
        end
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
        subroutine VCMTP(Npmax,Np,vX,vY,vZ,dNpINV,pm,pmINV,Tfac,vCM,
     +             temp,vpX,vpY,vpZ,vol,pressX,pressY,pressZ,pres)
        implicit none
        integer Npmax, Np, ip
        double precision vX(Npmax), vY(Npmax), vZ(Npmax),vpX(Npmax),
     +         vpY(Npmax),vpZ(Npmax),vol,vCM,dNpINV,pm,pmINV,Tfac,temp,
     +         vXi,vYi,vZi,vCMX,vCMY,vCMZ,pressX,pressY,pressZ,pres

        vCMX = 0.d0
        vCMY = 0.d0
        vCMZ = 0.d0
        temp = 0.d0
        pressX = 0.d0
        pressY = 0.d0
        pressZ = 0.d0
      
        do ip = 1, Np
            vXi = vX(ip)
            vYi = vY(ip)
            vZi = vZ(ip)
            vCMX = vCMX + pm*vXi
            vCMY = vCMY + pm*vYi
            vCMZ = vCMZ + pm*vZi
            temp = temp + pm*( vXi*vXi + vYi*vYi + vZi*vZi )
            pressX = pressX + vpX(ip)
            pressY = pressY + vpY(ip)
            pressZ = pressZ + vpZ(ip)
        end do

        vCM = dNpINV*pmINV*dsqrt(vCMX*vCMX + vCMY*vCMY + vCMZ*vCMZ )
        temp = Tfac*temp
        pres = (pressX + pressY +pressZ)/(3.0*vol)+3.0
        
        return
        end
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c     	SUBROUTINE:
c     	The pseudorandom number generator.
c     	This subroutine has been taken as it is.
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
        real function  R2S()
c***********************************************(GB 12/95)
c   Portable long-period (ca. 2.3 * 10^18) random number generator of
c   L'ECUYER [1] with BAYS-DURHAM shuffle [2] and added safeguards as
c   published by PRESS et al. [3] as "ran2" in 1992. In this version
c   (called "R2S" for distinction) no argument is needed, and the
c   initialization can be done by an entry point "R2INIS(iseedS)" with
c   any positive integer "iseedS" as seed. The internal state corres-
c   ponding to  iseedS = 1  is used if no initialization with "R2INIS"
c   is done before the first call of "R2S".
c
c   "R2S" returns a uniform random number in the interval  ]0., 1.[
c   (NOTE: The endpoint values are excluded!)
c
c   *  INITIALIZATION of "R2S":           rdummy = R2INIS(iseedS)
c   *  GENERATION of a random number:     r = R2S()
c
c   NOTE:
c   *  "rdummy" is a dummy variable of type REAL.
c   *  No variable declaractions are necessary in the calling module.
c   *  Parameter RNMX=1.-EPS should approximate the largest floating
c      point value that is less than 1. (EPS for a specific machine
c      can be determined with subroutine "MACHAR" from chapt. 20.1 of
c      ref [3]). (For "IBM RISC System/6000" workstations with "AIX XL
c      Fortran/6000", subroutine MACHAR gives  EPS=1.192092896E-07 )
c
c   REFERENCES:
c   [1]  P. L'ECUYER, Communications of the ACM, vol. 31 (1988) 742-774.
c   [2]  in D.E. KNUTH, "Seminumerical Algorithms" (2nd. ed.), vol. 2 of
c        "The Art of Computer Programming", Addison-Wesley, Reading, MA
c        (1981) paragraphs 3.2-3.3 .
c   [3]  W.H. PRESS, S.A. TEUKOLSKY, W.T. VETTERLING, and B.P. FLANNERY,
c        "Numerical Recipes in FORTRAN" (2nd ed.), Cambridge University
c        Press, Cambridge (1992), chapt. 7.1
c
c   TEST OUTPUT (first 35 numbers for iseed=1, in row-wise sequence):
c   0.285381  0.253358  0.093469  0.608497  0.903420  0.195873  0.462954
c   0.939021  0.127216  0.415931  0.533739  0.107446  0.399671  0.650371
c   0.027072  0.076975  0.068986  0.851946  0.637346  0.573097  0.902278
c   0.887676  0.372177  0.347516  0.207896  0.569131  0.364677  0.392418
c   0.366707  0.432149  0.153942  0.626891  0.749454  0.341041  0.830393
c***********************************************************************
        parameter ( IM1  = 2147483563,
     *              IM2  = 2147483399,
     *              AM   = 1./IM1,
     *              IMM1 = IM1-1,
     *              IA1  = 40014,
     *              IA2  = 40692,
     *              IQ1  = 53668,
     *              IQ2  = 52774,
     *              IR1  = 12211,
     *              IR2  = 3791,
     *              NTAB = 32,
     *              NDIV = 1+IMM1/NTAB,
     *              EPS  = 1.2e-7,
     *              RNMX = 1.-EPS )
        integer ivS(NTAB)
        save ivS, iyS, idum1S, idum2S
        data idum1S/1720212868/, idum2S/1/, iyS/1720212868/
        data ivS/ 1720212868, 1392842846, 1031324961,  718590712,
     *            82237802, 1816996195, 1529538438, 1789446856,
     *            156648835,   52437849, 1441478319,   36906150,
     *            1269685686, 1644535938,  394503142,  310212663,
     *            1596049480,    7553450,  322224693,  445508654,
     *            28884682,  643161691,  407948861,  479214492,
     *            2124954851,  612891482,  112933431, 1814689225,
     *            53445315, 1904850491, 1695805043, 1860990862 /
c***********************************************************************
c ----- Compute "idum1S" by  idum1S = mod( IA1*idum1S, IM1 ) without
c       overflow by SCHRAGE's method (see ref. [3]):
        kS = idum1S / IQ1
        idum1S = IA1*( idum1S - kS*IQ1 ) - kS*IR1
        if  ( idum1S .lt. 0 )  idum1S = idum1S + IM1
c ----- Compute "idum2S" likewise:
        kS = idum2S / IQ2
        idum2S = IA2*( idum2S - kS*IQ2 ) - kS*IR2
        if  ( idum2S .lt. 0 )  idum2S = idum2S + IM2
c ----- "jS" will be in the range [1 (1) NTAB] :
        jS = 1 + iyS/NDIV
c ----- Here "idum1S" is shuffled, and "idum1S" and "idum2S" are combined
c       to produce output:
        iyS = ivS(jS) - idum2S
        ivS(jS) = idum1S
        if ( iyS .lt. 1 )  iyS = iyS + IMM1
c ----- Because users don't expect endpoint values:
        R2S = min( AM*iyS, RNMX )
        return
c >>>>> Initialization: >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
        entry R2INIS (iseedS)
c ----- Be sure to prevent a negative "iseedS" or  iseedS = 0 :
        idum1S = max( iabs(iseedS), 1 )
        idum2S = idum1S
c ----- Load the shuffle table (after 8 warm-ups):
        do jS= NTAB+8, 1, -1
            kS = idum1S / IQ1
            idum1S = IA1*( idum1S - kS*IQ1 ) - kS*IR1
            if  ( idum1S .lt. 0 )  idum1S = idum1S + IM1
            if  ( jS .le. NTAB )  ivS(jS) = idum1S
        end do
        iyS = ivS(1)
        R2INIS = iyS
        return
        end
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
        subroutine FCDR(Npmax,Np,ncell,mapsize,head,map,list,rX,rY,rZ,
     +	           voX,voY,voZ,FcX,FcY,FcZ,FdX,FdY,FdZ,FrX,
     +	           FrY,FrZ,vpX,vpY,vpZ,n_flag,alpha,dLx,dLy,dLz,
     +	           dLxINV,dLyINV,dLzINV,rfac,gamma,sigma,incstp)
        implicit none
        integer Npmax,Np,ncell, mapsize, n_flag(Npmax), incstp
        integer head(ncell),map(mapsize),list(Npmax)
        double precision
     +      rX(Npmax),  rY(Npmax), rZ(Npmax),rX0(Npmax),rY0(Npmax),
     +      rZ0(Npmax),voX(Npmax),voY(Npmax),voZ(Npmax),
     +      FcX(Npmax),FcY(Npmax),FcZ(Npmax),FdX(Npmax),FdY(Npmax),
     +      FdZ(Npmax),FrX(Npmax),FrY(Npmax),FrZ(Npmax),
     +      vpX(Npmax),vpY(Npmax),vpZ(Npmax),alpha(3,3),rfac,
     +      dLx, dLy, dLz, dLxINV,dLyINV, dLzINV,gamma, sigma
        integer ip, jp, im, nfi, nfj,icell,jcell0,nabor,
     +      jcell, iz,izb, iza, isb

        double precision
     +      rXi, rYi, rZi, vXi,vYi,vZi,vpXi,vpYi,vpZi,FcXi,FcYi,FcZi,
     +      FdXi, FdYi, FdZi, FrXi, FrYi, FrZi,rXij, rYij, rZij,
     +      rijSQ,rij,rijINV,omega,eXij,eYij,eZij,vXij,vYij,vZij,
     +      Fcfac, Fdfac, Frfac,FcXij, FcYij, FcZij,
     +      FdXij, FdYij, FdZij,FrXij, FrYij, FrZij,
     +      z1, z2, rZii, rZjj, delq

        real R2S

c ----- Initialization of all forces
        do ip = 1, Np
          FcX(ip) = 0.d0
          FcY(ip) = 0.d0
          FcZ(ip) = 0.d0
          FdX(ip) = 0.d0
          FdY(ip) = 0.d0
          FdZ(ip) = 0.d0
          FrX(ip) = 0.d0
          FrY(ip) = 0.d0
          FrZ(ip) = 0.d0
          vpX(ip) = 0.d0
          vpY(ip) = 0.d0
          vpZ(ip) = 0.d0
        end do

        do icell = 1, ncell
            ip = head(icell)
 1000    	  if (ip.gt.0) then
                rXi =  rX(ip)
                rYi =  rY(ip)
                rZi =  rZ(ip)
                vXi =  voX(ip)
                vYi =  voY(ip)
                vZi =  voZ(ip)
                FcXi = FcX(ip)
                FcYi = FcY(ip)
                FcZi = FcZ(ip)
                FdXi = FdX(ip)
                FdYi = FdY(ip)
                FdZi = FdZ(ip)
                FrXi = FrX(ip)
                FrYi = FrY(ip)
                FrZi = FrZ(ip)
                vpXi = vpX(ip)
                vpYi = vpY(ip)
                vpZi = vpZ(ip)

                nfi = n_flag(ip)
                jp = list(ip)

 2000       	  if (jp.gt.0) then
                    rXij = rXi - rX(jp)
                    rYij = rYi - rY(jp)
                    rZij = rZi - rZ(jp)
                    rXij = rXij - dLx*dnint( rXij*dLxINV )
                    rYij = rYij - dLy*dnint( rYij*dLyINV )
                    rZij = rZij - dLz*dnint( rZij*dLzINV )
                    rijSQ = rXij*rXij + rYij*rYij + rZij*rZij

                    if ( rijSQ .lt. 1.d0 ) then
                        rij    = dsqrt( rijSQ )
                        rijINV = 1.d0/rij
                        omega  = 1.d0 - rij

                        eXij = rXij*rijINV
                        eYij = rYij*rijINV
                        eZij = rZij*rijINV

                        vXij = vXi - voX(jp)
                        vYij = vYi - voY(jp)
                        vZij = vZi - voZ(jp)

                        nfj = n_flag(jp)

                        Fcfac = omega
                        FcXij = Fcfac * eXij * alpha(nfi,nfj)
                        FcYij = Fcfac * eYij * alpha(nfi,nfj)
                        FcZij = Fcfac * eZij * alpha(nfi,nfj)

                       Fdfac=omega*omega*(eXij*vXij+eYij*vYij+eZij*vZij)
                        FdXij = Fdfac*eXij
                        FdYij = Fdfac*eYij
                        FdZij = Fdfac*eZij

                        Frfac=omega * rfac * ( 2.*R2S() - 1.0 )
                        FrXij = Frfac*eXij
                        FrYij = Frfac*eYij
                        FrZij = Frfac*eZij

                        FcXi = FcXi + FcXij
                        FcYi = FcYi + FcYij
                        FcZi = FcZi + FcZij
                        FdXi = FdXi + FdXij
                        FdYi = FdYi + FdYij
                        FdZi = FdZi + FdZij
                        FrXi = FrXi + FrXij
                        FrYi = FrYi + FrYij
                        FrZi = FrZi + FrZij
                        vpXi = vpXi + FcXij*rXij
                        vpYi = vpYi + FcYij*rYij
                        vpZi = vpZi + FcZij*rZij

                        FcX(jp) = FcX(jp) - FcXij
                        FcY(jp) = FcY(jp) - FcYij
                        FcZ(jp) = FcZ(jp) - FcZij
                        FdX(jp) = FdX(jp) - FdXij
                        FdY(jp) = FdY(jp) - FdYij
                        FdZ(jp) = FdZ(jp) - FdZij
                        FrX(jp) = FrX(jp) - FrXij
                        FrY(jp) = FrY(jp) - FrYij
                        FrZ(jp) = FrZ(jp) - FrZij
                    endif

                    jp = list(jp)
                    go to 2000
                endif

                jcell0 = 13*(icell-1)
                do nabor = 1, 13
                    jcell = map(jcell0+nabor)
                    jp = head(jcell)
 3000       	   	  if (jp.gt.0) then
                        rXij = rXi - rX(jp)
                        rYij = rYi - rY(jp)
                        rZij = rZi - rZ(jp)
                        rXij = rXij - dLx*dnint( rXij*dLxINV )
                        rYij = rYij - dLy*dnint( rYij*dLyINV )
                        rZij = rZij - dLz*dnint( rZij*dLzINV )
                        rijSQ = rXij*rXij + rYij*rYij + rZij*rZij

                        if (rijSQ.lt.1.d0) then
                            rij    = dsqrt( rijSQ )
                            rijINV = 1.d0/rij
                            omega  = 1.d0 - rij

                            eXij  = rXij*rijINV
                            eYij  = rYij*rijINV
                            eZij  = rZij*rijINV

                            vXij  = vXi - voX(jp)
                            vYij  = vYi - voY(jp)
                            vZij  = vZi - voZ(jp)

                            nfj   = n_flag(jp)

                            Fcfac = omega
                            FcXij = Fcfac * eXij * alpha(nfi,nfj)
                            FcYij = Fcfac * eYij * alpha(nfi,nfj)
                            FcZij = Fcfac * eZij * alpha(nfi,nfj)

                       Fdfac=omega*omega*(eXij*vXij+eYij*vYij+eZij*vZij)
                            FdXij = Fdfac*eXij
                            FdYij = Fdfac*eYij
                            FdZij = Fdfac*eZij

                            Frfac = omega * rfac*( 2.*R2S() - 1. )
                            FrXij = Frfac*eXij
                            FrYij = Frfac*eYij
                            FrZij = Frfac*eZij
                            FcXi  = FcXi + FcXij
                            FcYi  = FcYi + FcYij
                            FcZi  = FcZi + FcZij
                            FdXi  = FdXi + FdXij
                            FdYi  = FdYi + FdYij
                            FdZi  = FdZi + FdZij
                            FrXi  = FrXi + FrXij
                            FrYi  = FrYi + FrYij
                            FrZi  = FrZi + FrZij
                            vpXi = vpXi + FcXij*rXij
                            vpYi = vpYi + FcYij*rYij
                            vpZi = vpZi + FcZij*rZij

                            FcX(jp) = FcX(jp) - FcXij
                            FcY(jp) = FcY(jp) - FcYij
                            FcZ(jp) = FcZ(jp) - FcZij
                            FdX(jp) = FdX(jp) - FdXij
                            FdY(jp) = FdY(jp) - FdYij
                            FdZ(jp) = FdZ(jp) - FdZij
                            FrX(jp) = FrX(jp) - FrXij
                            FrY(jp) = FrY(jp) - FrYij
                            FrZ(jp) = FrZ(jp) - FrZij

                        endif

                        jp = list(jp)
                        go to 3000

                    endif
                enddo

                FcX(ip) = FcXi
                FcY(ip) = FcYi
                FcZ(ip) = FcZi
                FdX(ip) = FdXi
                FdY(ip) = FdYi
                FdZ(ip) = FdZi
                FrX(ip) = FrXi
                FrY(ip) = FrYi
                FrZ(ip) = FrZi
                vpX(ip) = vpXi
                vpY(ip) = vpYi
                vpZ(ip) = vpZi

                ip = list(ip)
                go to 1000
            endif
        enddo

        do ip = 1, Np
            FdX(ip) = -gamma*FdX(ip)
            FdY(ip) = -gamma*FdY(ip)
            FdZ(ip) = -gamma*FdZ(ip)
            FrX(ip) = sigma*FrX(ip)
            FrY(ip) = sigma*FrY(ip)
            FrZ(ip) = sigma*FrZ(ip)
        end do

        return
        end
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
        subroutine FBondAngle(Npmax,Np,len_chain,rX,rY,rZ,FbondX,FbondY,
     +             FbondZ,vpX, vpY,vpZ,n_flag,dLx,dLy,dLz,dLxINV,
     +             dLyINV,dLzINV,bond) 
        implicit none
        integer Npmax,Np,n_flag(Npmax),iendnumber,icenter,
     +	    ichargechain,inetrochain1,inetrochain2,len_chain,
     +	    lipidhead,lipidend,iiii,minus,n0,n11,n12,n21,n22,iii

        double precision
     +	    rX(Npmax),rY(Npmax),rZ(Npmax),FbondX(Npmax),FbondY(Npmax),
     +	    FbondZ(Npmax),vpX(Npmax),vpY(Npmax),vpZ(Npmax),dLx,dLy,
     +	    dLz,dLxINV,dLyINV,dLzINV,bond,rXij,rYij,rZij,rrr,r,rINV
 
        do iiii = 1, Np
            FbondX(iiii) = 0.d0
            FbondY(iiii) = 0.d0
            FbondZ(iiii) = 0.d0
        end do

        do iii = 2, len_chain
            rXij = rX(iii) - rX(iii-1)
            rYij = rY(iii) - rY(iii-1)
            rZij = rZ(iii) - rZ(iii-1)
            rXij = rXij - dLx*dnint( rXij*dLxINV )
            rYij = rYij - dLy*dnint( rYij*dLyINV )
            rZij = rZij - dLz*dnint( rZij*dLzINV )
            rrr = rXij*rXij + rYij*rYij+ rZij*rZij
            r = dsqrt(rrr)
            rINV = 1.0 / r
            FbondX(iii) = FbondX(iii) - bond * rXij * rINV
            FbondY(iii) = FbondY(iii) - bond * rYij * rINV
            FbondZ(iii) = FbondZ(iii) - bond * rZij * rINV
            FbondX(iii-1) = FbondX(iii-1) + bond * rXij * rINV
            FbondY(iii-1) = FbondY(iii-1) + bond * rYij * rINV
            FbondZ(iii-1) = FbondZ(iii-1) + bond * rZij * rINV
        end do

        return
        end
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
        subroutine SNAPSHOT(Npmax,Np,nconf,n_flag,charge,dLx,dLy,dLz,
     +	           dLxINV,dLyINV,dLzINV,rX,rY,rZ,Ifile)
        implicit none
        integer Npmax, Np,lench,nconf, n_flag(Npmax),Ifile,ip,j,ia
        double precision rX(Npmax), rY(Npmax), rZ(Npmax), charge(Npmax),
     +      drX_dpl(Npmax), drY_dpl(Npmax), drZ_dpl(Npmax),
     +      dLx, dLy, dLz, dLxINV, dLyINV,dLzINV,drX,drY,drZ

        do ip = 1, Np
            drX = rX(ip)
            drY = rY(ip)
            drZ = rZ(ip)
            drX_dpl(ip) = drX - dLx*dnint(drX * dLxINV)
            drY_dpl(ip) = drY - dLy*dnint(drY * dLyINV)
            drZ_dpl(ip) = drZ - dLz*dnint(drZ * dLzINV)
        end do

        do j=1, Np
            write(nconf,925)
     +	        j,n_flag(j),charge(j),drX_dpl(j),drY_dpl(j),drZ_dpl(j)
        end do
 925  	format(1x,1i6,1x,1i3,1x,1f6.3, 3(1x,1f10.6))
        close(nconf)
        return
        end
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
        subroutine FElecReal(Npmax,Np,Ncharge,rX,rY,rZ,vpX,vpY,vpZ,
     +     	       charge,idcharge,alpha_chg,dLx,dLy,dLz,
     +     	       dLxINV,dLyINV,dLzINV,FrealX,FrealY,FrealZ,real_eng)
        implicit none
        integer Npmax,Ncharge,Np,idcharge(Ncharge)
        double precision rX(Npmax), rY(Npmax), rZ(Npmax),charge(Npmax),
     +	    FrealX(Npmax), FrealY(Npmax), FrealZ(Npmax),vpX(Npmax),
     +	    vpY(Npmax), vpZ(Npmax),dLx,dLy,dLz,dLxINV,dLyINV,dLzINV
        integer ip, jp, im, nfi, nfj,ni,nj,icharge,jcharge
        double precision rXi,rYi,rZi,rXj,rYj,rZj,rXij,rYij,rZij,rijSQ,
     +	     rij, rijINV,rrijINV,eXij,eYij,eZij,Freal,FrealXij,FrealYij,
     +	     FrealZij,chgi,chgj,pi,piINV,aa,ab,ac,ad,alpha_chg,
     +	     FrFactor,aaa,abb,g1, g2, g3, g4, real_eng

        pi = 4.0*atan(1.0)
        piINV = 1.0/dsqrt(pi)
        real_eng = 0.0

        do ip = 1, Np
            FrealX(ip) = 0.d0
            FrealY(ip) = 0.d0
            FrealZ(ip) = 0.d0
        end do

        do icharge = 1, Ncharge - 1
            do jcharge = icharge + 1, Ncharge

                ni = idcharge(icharge)
                nj = idcharge(jcharge)

                rXi =  rX(ni)
                rYi =  rY(ni)
                rZi =  rZ(ni)
                rXj =  rX(nj)
                rYj =  rY(nj)
                rZj =  rZ(nj)

                rXij = rXi - rXj
                rYij = rYi - rYj
                rZij = rZi - rZj
                rXij = rXij - dLx*dnint( rXij*dLxINV )
                rYij = rYij - dLy*dnint( rYij*dLyINV )
                rZij = rZij - dLz*dnint( rZij*dLzINV )
                rijSQ = rXij*rXij + rYij*rYij + rZij*rZij

                if  ( rijSQ .lt. 9.d0 )  then ! r_cutoff is tunable
                    chgi = charge(ni)
                    chgj = charge(nj)
                    rij    = dsqrt( rijSQ )
                    rijINV = 1.d0/rij

                    eXij = rXij*rijINV
                    eYij = rYij*rijINV
                    eZij = rZij*rijINV
                    rrijINV = 1.d0/rijSQ

                    g1 = 13.87 / (4*pi)
                    g2 = exp(-2.0*0.929*rij)
                    g3 = 1.0 + 0.929*rij
                    g4 = 2.0*0.929*rij

                    aa = erfc(alpha_chg*rij)*rijINV
                    real_eng = real_eng + chgi*chgj*aa*(g1*(1.0-g2*g3))

                    ab = 2.0*alpha_chg*piINV
                    ac = ab * exp(-alpha_chg*alpha_chg*rijSQ)
                    ad = rrijINV*(g1*(1.0-g2-g2*g3*g4))
                    Freal = chgi*chgj*(ac+aa)*ad

                    FrealXij = Freal*rXij
                    FrealYij = Freal*rYij
                    FrealZij = Freal*rZij
                    FrealX(ni) = FrealX(ni) + FrealXij
                    FrealY(ni) = FrealY(ni) + FrealYij
                    FrealZ(ni) = FrealZ(ni) + FrealZij
                    vpX(ni) = vpX(ni)+ FrealXij*rXij
                    vpY(ni) = vpY(ni)+ FrealYij*rYij
                    vpZ(ni) = vpZ(ni)+ FrealZij*rZij
                    FrealX(nj) = FrealX(nj) - FrealXij
                    FrealY(nj) = FrealY(nj) - FrealYij
                    FrealZ(nj) = FrealZ(nj) - FrealZij

                end if
            end do
        end do
        real_eng = real_eng - alpha_chg*piINV*Ncharge

        return
        end
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
        subroutine ALLFORCE(Npmax,Np,FcX,FcY,FcZ,
     +	           FbondX, FbondY, FbondZ,FrealX,FrealY,
     +	           FrealZ,FrecipX,FrecipY,FrecipZ)
        implicit none
        integer  Npmax,Np,i
        double precision
     +	    FcX(Npmax), FcY(Npmax), FcZ(Npmax),
     +	    FbondX(Npmax), FbondY(Npmax), FbondZ(Npmax),
     +	    FrealX(Npmax), FrealY(Npmax), FrealZ(Npmax),
     +	    FrecipX(Npmax),FrecipY(Npmax),FrecipZ(Npmax)

        do i = 1, Np
            FcX(i)= FcX(i) + FrealX(i) + FrecipX(i) + FbondX(i)
            FcY(i)= FcY(i) + FrealY(i) + FrecipY(i) + FbondY(i)
            FcZ(i)= FcZ(i) + FrealZ(i) + FrecipZ(i) + FbondZ(i)
        end do
        return
        end
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

