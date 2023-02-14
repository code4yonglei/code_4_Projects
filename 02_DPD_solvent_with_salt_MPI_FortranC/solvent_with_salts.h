cccccccccccccccccccccccccccccccccccc main part
        integer
     +      maxstp,              ! Number of time steps after equilibration
     +      incstp,              ! Number of steps between samples
     +      isteps,              ! Number of time steps for equilibration
     +      issstp,              ! Number of steps between snapshots
     +      rho0,
     +      Np, it, n

        double precision
     +      alpha(4,4),          ! Conservative force strengths
     +      gamma,               ! Dissipative force strength
     +      sigma                ! Random force strength

        double precision
     +      vol,
     +      dispXcm, dispYcm,
     +      dt,                  ! Time step for the molecular dynamics
     +      tmax,                ! Time of the simulation after equilibration
     +      vCM,                 ! Center-of-mass velocity
     +      temp,                ! Temperature
     +      d_bond               ! Harmonic force constant for the chain

        double precision dtH,dtRTH
        double precision Ran_Uniform     ! For the random number generator

        integer iseed, nconf, Ifile

        parameter( Ifile = 10 )
  
        common/parameter1/alpha,gamma,sigma,rho0,rfac,Tfac,vol
        common/parameter2/dr_gyr,dr_ee,dispXcm,dispYcm,vCM,temp,d_bond
        common/parameter3/it,dt,dtH,dtRTH

cccccccccccccccccccccccccccccccccccc MPI
   
        integer ierr,myid,numprocs, status, request, rc
        double precision starttime,endtime,elapsetime
        common/mmppii/ierr,myid,numprocs,status,request,rc

cccccccccccccccccccccccccccccccccccc GENCON
 
        integer i,j,k,nxnodes,nynodes,nznodes,num
        parameter ( num =27 )
        integer Npmax,Nc, Nmp,Nm,Nmx,Nmy,Nmz,Ncore,WPVIS,MAXMSG

        parameter ( Npmax = 41472 )
        parameter ( Nmp = 8000 ) ! particle max number of one processor
        parameter ( Nch = 1000 ) ! max number of charge sin one processor
        parameter ( WPVIS = 10 ) ! id,idglob,rX,rY,rZ,vX,vY,vZ,n_flag,chg
        parameter ( MAXMSG = WPVIS*Nmp)

        integer ip,ipp,ipp2,ip2,ip3,ip4
        double precision
     +      dLx,dLy,dLz,dLxINV,dLyINV,dLzINV,
     +      dNm,dNmINV,dNnINV,dNp,dNpINV,
     +      IRVG(MAXMSG),
     +      pm,pmINV,Tfac,rfac,drX_tmp,drY_tmp,Tfacg,
     +      r_dist,d_scale,rXi,rYi,rZi,vv,
     +      vCMX(0:num),vCMY(0:num),vCMZ(0:num),tscal,temper(0:num)

        parameter ( dLx = 24.d0 )
        parameter ( dLy = 24.d0 )
        parameter ( dLz = 24.d0 )

        common/inform/IRVG
        common/num/Ncore,Nm,Nc,Nmx,Nmy,Nmz
        common/para/dLxINV,dLyINV,dLzINV,pm,pmINV,dNp,dNpINV
        common/node/nxnodes,nynodes,nznodes

cccccccccccccccccccccccccccccccccccc EXTVOL

        integer up,down,up1,down1
        common/weizhi/up,down,up1,down1

        integer
     +      nclx,ncly,nclz,ncly1,nclc1,
     +      rX0,rY0,rZ0,l, m,m1,
     +      VPVIS,MAXVSG,CPVIS,MAXCSG,MAXNSG,
     +      msglen1,msglen2,msglen3,msglen4,
     +      msglen11,msglen22,msglen33,msglen44

        parameter ( VPVIS =13) ! id,idglob,rX,rY,rZ,vX,vY,vZ,
            ! n_flag,Cx,Cy or Cz,
            ! vX,vY,vZ

        parameter ( MAXVSG = VPVIS*Nmp)
        parameter ( CPVIS = 4)  ! Cx,Cy,Cz,Cs
        parameter ( MAXCSG = CPVIS*Nmp)

        integer CXYS(MAXCSG)

        double precision
     +      msgbuf1(MAXVSG),msgbuf2(MAXVSG),
     +      msgbuf3(MAXVSG),msgbuf4(MAXVSG),
     +      msgbuf11(MAXVSG),msgbuf22(MAXVSG),
     +      msgbuf33(MAXVSG),msgbuf44(MAXVSG)
       
        common/ncl/nclx,ncly,nclz,ncly1,nclc1
        common/msg1/msgbuf1,msgbuf2,msgbuf3,msgbuf4
        common/msg2/msgbuf11,msgbuf22,msgbuf33,msgbuf44
        common/cell/CXYS

cccccccccccccccccccccccccccccccccccc LINKS

        integer cell,ncell,list(Nmp),
     +        cellix,celliy,celliz,cellx,celly,cellz
        parameter ( ncell = INT(dLz*dLy*dLx) )
        integer head(ncell)
        common/headcell/head,list

cccccccccccccccccccccccccccccccccccc MAPS

        integer ix,iy,iz,imap,icell,mapx
        parameter(mapx=ncell*13)

        integer map(mapx)
        common/mapxs/map

cccccccccccccccccccccccccccccccccccc FCDR + Fenuf(Real+Reciprocal)
 
        integer
     +      FPVIS,MAXFSG,jcell,jcell0,
     +      jp,jp2,jp3,jp4,jp5,jp6,nfi,nfj,nabor,
     +      p,e,f,g,h,t,m2,
     +      ip5,ip6,PPVIS,MAXPSG,
     +      l1,l2,l3,l4,
     +      m3,m4,ZPVIS,MAXZSG

        parameter ( FPVIS = 9)  ! Fcx,Fcy,Fcz,Fdx,Fdy,Fdz,Frx,Fry,Frz
        parameter ( MAXFSG = FPVIS*Nmp)

        parameter ( EPVIS = 8)  ! Real and Reciprocal F and E from ENUF
        parameter ( MAXESG = EPVIS*Nmp)

        parameter (ZPVIS =6)    !Fzx,Fzy,Fzz,Frx,Fry,Frz
        parameter (MAXZSG=ZPVIS*Nmp)
    
        parameter (PPVIS = 3)
        parameter (MAXPSG = PPVIS*Nmp)

        double precision
     +      FCDR(MAXFSG), FenufRF(MAXESG),omega,FZ(MAXZSG),VPXY(MAXPSG),
     +      vX0i, vY0i, vZ0i,vX0j,vY0j,vZ0j,
     +      FcXi, FcYi, FcZi, FdXi, FdYi, FdZi, FrXi, FrYi, FrZi,
     +      rXij, rYij, rZij,rijSQ, rij, rijINV, 
     +      eXij, eYij, eZij,vX0ij, vY0ij, vZ0ij,Fcfac, Fdfac, Frfac,
     +      FcXij, FcYij, FcZij, FdXij, FdYij, FdZij,
     +      FrXij, FrYij, FrZij,
     +      vpXi,vpYi,vpZi, a,b,c,d,w,z

        common/forcez/FZ,FCDR,FenufRF
        common/virial/VPXY

cccccccccccccccccccccccccccccccccccc Velocity COM
        integer VIS,MAXVIS
        parameter( VIS = 3 )
        parameter( MAXVIS=VIS*Nmp)
       
        double precision
     +      vXi,vYi,vZi,tem(0:num),vCMXM,vCMYM,vCMZM,
     +      VXYZ(MAXVIS),press,
     +      pX,pY,pZ,pressX(0:num),pressY(0:num),pressZ(0:num)

        common/velocity/VXYZ
        common/vc/vCMXM,vCMYM,vCMZM

cccccccccccccccccccccccccccccccccccc R_Radius_of_Gyration__End2End_Dist
 
        double precision
     +      dr_gyr, dr_ee,
     +      drX_cm, drY_cm, drX_g, drY_g,
     +      disp_X, disp_Y, drX_tm, drY_tm,
     +      drX_ee, drY_ee
      
        integer lenc,displs(0:num),len(0:num), sm

cccccccccccccccccccccccccccccccccccc INVT & INTR

        double precision hcd0,hr0,hcd,hr

cccccccccccccccccccccccccccccccccccc MIGRATE

        integer myrxyzs(Nmp)

cccccccccccccccccccccccccccccccccccc RESULT

        common/res/nconf

cccccccccccccccccccccccccccccccccccc DENSITY

        integer rhoA(dLx,dLy,dLz),rhoB(dLx,dLy,dLz)

cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
