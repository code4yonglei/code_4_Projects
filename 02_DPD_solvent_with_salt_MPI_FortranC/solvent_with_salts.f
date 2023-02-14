cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

		program  DPD_3D_MPI  ! linear CPU

cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

		implicit none
		include 'mpif.h'
		include 'solvent_with_salts.h'

		rho0       = 3.d0         ! total density of all particles
		pm         = 1.d0         ! particle mass (identical for all particles)

		alpha(1,1) = 25.d0        ! water - water interaction
		alpha(2,2) = 25.d0        ! monomer - monomer
		alpha(3,3) = 25.d0        ! salt-cation
		alpha(4,4) = 25.d0        ! salt-anion

		alpha(1,2) = 50.d0
		alpha(1,3) = 25.d0
		alpha(1,4) = 25.d0

		alpha(2,1) = alpha(1,2)
		alpha(2,3) = 50.00
		alpha(2,4) = 50.00

		alpha(3,4) = 25.d0

		sigma      = 3.d0         ! dpd parameter
		iseed      = 612345       ! seed for a random number generator
		dt         = 5.00d-2      ! time step for the integrator
		isteps     = 500000       ! number of steps for simulation
		incstp     = 100          ! how often do we calculate quantities
		issstp     = 1000         ! how often do we take snapshots

		alpha_chg  = 0.35
		kspace     = 7

cccccccccccccccccccccccccc Fix some constants cccccccccccccccccccccccccc

		dLxINV = 1.d0/dLx         ! Inverse system size
		dLyINV = 1.d0/dLy
		dLzINV = 1.d0/dLz
		vol = dLx*dLy*dLz         ! Volume (area) of the system
		Np = nint( rho0*vol )     ! Total number of particles
		dNpINV=1.d0/dble(Np)

		pmINV = 1.d0/pm           ! Inverse particle mass
		gamma = 0.5d0*sigma*sigma ! strength of Fd (fluctuation-dissipation theorem)
		dtH   = 0.5d0*dt          ! Time steps for the integrator
		dtRTH = 0.5d0*dsqrt(dt)
		rfac  = dsqrt(3.d0)       ! Scaling of random forces
		Tfac  = 1.d0/(3.d0*dble(Np) - 3.0) ! Factor for temperature control

cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c                                                                      c
c                         Start DPD simulation                         c
c                                                                      c
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
		call MPI_INIT(ierr)
		call MPI_COMM_RANK(MPI_COMM_WORLD, myid, ierr)
		call MPI_COMM_SIZE(MPI_COMM_WORLD, numprocs, ierr)
		print*,"Process",myid,"of",numprocs,"is alive"
		if (myid.eq.0) then
			starttime = MPI_WTIME()
			OPEN(3,STATUS='NEW', FILE='data.dat',
	 +     		FORM='FORMATTED', ACCESS='SEQUENTIAL')
			OPEN(5,STATUS='NEW', FILE='time.dat',
	 +     		FORM='FORMATTED', ACCESS='SEQUENTIAL')
		endif


		do it = 1, isteps  ! EQUILIBRATE THE SYSTEM!
			if (it .eq. 1)then
				call Ran_Init(iseed)
				call GENCON  ! including FElecRecipENUF_INIT
				call MAPS
			else
				call INTV
				call INTR
				call MIGRATE   ! tranfer bead info to new process
				call EXTVOL
			endif

			call LINKS
			call F_CDR_EE
c			call FBondAngle
			call INTV

			if ((mod(it,incstp).eq.0).or.(it.eq.1)) call VCMTP
			if ((mod(it,issstp).eq.0).or.(it.eq.1)) call SNAPSHOT

c           output other properties of modelling system
		enddo


		if (myid.eq.0) then
			endtime = MPI_WTIME()
			elapsetime = endtime - starttime
			write(5,50) elapsetime
			call FLUSH(5)
			close(3)
			close(5)
			close(myid+50)
		endif
 50     format(5x,'elapsetime:',g15.6,'s')
		if(myid .eq. 0)
			call FElecRecipENUF_FINALIZE()
		endif
		call MPI_FINALIZE(rc)

		stop 'Good luck!'
		end

cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c                                                                      c
c                        End of DPD simulations                        c
c                                                                      c
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

		subroutine GENCON
		include 'mpif.h'
		include 'solvent_with_salts.h'

		nxnodes = 1
		nynodes = numprocs  ! domain decomposition along Y direction
		nznodes = 1

		if (myid.ne.0) then
			Ncore = Npmax/numprocs
			Nchcore = Nch/numprocs
		else
			Ncore = Npmax - Npmax/numprocs*(numprocs-1)
			Nchcore = Nch - Nch/numprocs*(numprocs-1)
			call FElecRecipENUF_INIT(kspace, dLx, dLy, dLz, Nch)  ! ENUF
		endif

		Nbinary = Ncore / 2 ! binary mixture

		do i = 1, MAXMSG
			IRVG(i)=0
		enddo

		do i = 1,MAXVSG
			msgbuf1(i)=0
			msgbuf2(i)=0
			msgbuf3(i)=0
			msgbuf4(i)=0
		enddo

		vCMX(myid) = 0.d0
		vCMY(myid) = 0.d0
		vCMZ(myid) = 0.d0

cccccccccccccccccccccccccccccccccccccccccccccccccccccccccc 

		do  ip = 1, Nbinary
			ip2 = (ip-1) * WPVIS
			IRVG(ip2+1) = ip
			IRVG(ip2+2) = myid*Nm+ip

			IRVG(ip2+3) =(Ran_Uniform()-0.5)*dLx
			if (myid.ne.numprocs-1) then
				IRVG(ip2+4)=-dLy/2+myid*DINT(dLy/nynodes)
	 +              + Ran_Uniform()*DINT(dLy/nynodes)
			else
				IRVG(ip2+4)=-dLy/2+myid*DINT(dLy/nynodes)
	 +      		+ Ran_Uniform()*(dLy-myid*DINT(dLy/nynodes))
			endif
			IRVG(ip2+5) =(Ran_Uniform()-0.5)*dLz

			vv = dble( Ran_Uniform() - 0.5 )
			IRVG(ip2+6) = vv
			vCMX(myid) = vCMX(myid) + pm*vv

			vv  = dble( Ran_Uniform() - 0.5 )
			IRVG(ip2+7) = vv
			vCMY(myid) = vCMY(myid) + pm*vv

			vv  = dble( Ran_Uniform() - 0.5 )
			IRVG(ip2+8) =vv
			vCMZ(myid) = vCMZ(myid) + pm*vv

			if (Ncore - ip < int(Nchcore/2))then
				IRVG(ip2+9)  = 1   ! bead type
				IRVG(ip2+10) = 0.0 ! neutral beads
			else
				IRVG(ip2+9)  = 3
				IRVG(ip2+10) = 1.0 ! cations
			endif
		enddo


		do  ip = Nbinary+1, Ncore
			ip2 = (ip-1) * WPVIS
			IRVG(ip2+1) = ip
			IRVG(ip2+2) = myid*Nm+ip

			IRVG(ip2+3) =(Ran_Uniform()-0.5)*dLx
			if (myid.ne.numprocs-1) then
				IRVG(ip2+4)=-dLy/2+myid*DINT(dLy/nynodes)
	 +      	+ Ran_Uniform()*DINT(dLy/nynodes)
			else
				IRVG(ip2+4)=-dLy/2+myid*DINT(dLy/nynodes)
	 +      	+ Ran_Uniform()*(dLy-myid*DINT(dLy/nynodes))
			endif
			IRVG(ip2+5) =(Ran_Uniform()-0.5)*dLz

			vv = dble( Ran_Uniform() - 0.5 )
			IRVG(ip2+6) = vv
			vCMX(myid) = vCMX(myid) + pm*vv

			vv  = dble( Ran_Uniform() - 0.5 )
			IRVG(ip2+7) = vv
			vCMY(myid) = vCMY(myid) + pm*vv

			vv  = dble( Ran_Uniform() - 0.5 )
			IRVG(ip2+8) =vv
			vCMZ(myid) = vCMZ(myid) + pm*vv

			if (Ncore - ip < int(Nchcore - Nchcore/2))then
				IRVG(ip2+9)  = 2
				IRVG(ip2+10) = 0.0
			else
				IRVG(ip2+9)  = 4
				IRVG(ip2+10) = -1.0 ! anions
			endif
		enddo

cccccccccccccccccccccccccccccccccccc

		msgbuf1(1) = vCMX(myid)
		msgbuf1(2) = vCMY(myid)
		msgbuf1(3) = vCMZ(myid)

		call MPI_GATHER(msgbuf1(1), 3, MPI_DOUBLE_PRECISION,
	 +       msgbuf2(1), 3, MPI_DOUBLE_PRECISION,
	 +       0, MPI_COMM_WORLD, ierr)

		msgbuf1(1) = 0.d0
		msgbuf1(2) = 0.d0
		msgbuf1(3) = 0.d0

		if (myid.eq.0) then
			vCMXM = 0.d0
			vCMYM = 0.d0
			vCMZM = 0.d0

			l=0
			do i = 0,numprocs-1
				vCMXM = vCMXM + msgbuf2(l+1)
				vCMYM = vCMYM + msgbuf2(l+2)
				vCMZM = vCMZM + msgbuf2(l+3)
				msgbuf2(l+1) = 0.d0
				msgbuf2(l+2) = 0.d0
				msgbuf2(l+3) = 0.d0
				l=l+3
			enddo
			vCMXM = vCMXM * dNpINV * pmINV
			vCMYM = vCMYM * dNpINV * pmINV
			vCMZM = vCMZM * dNpINV * pmINV
		endif

		call MPI_BCAST(vCMXM, 3, MPI_DOUBLE_PRECISION,
	 +                  0, MPI_COMM_WORLD, ierr)

		do ip = 1,Ncore
			ip2 = (ip-1) * WPVIS
			IRVG(ip2+6) = IRVG(ip2+6) - vCMXM
			IRVG(ip2+7) = IRVG(ip2+7) - vCMYM
			IRVG(ip2+8) = IRVG(ip2+8) - vCMZM
		enddo

		temper(myid) = 0.d0    ! rescale velocities

		do ip = 1,Ncore
			ip2 = (ip-1)* WPVIS
			temper(myid) = temper(myid) +
	 +                pm*( IRVG(ip2+6)*IRVG(ip2+6) +
	 +                     IRVG(ip2+7)*IRVG(ip2+7) +
	 +                     IRVG(ip2+8)*IRVG(ip2+8))
		enddo

		msgbuf1(1) = temper(myid)

		call MPI_GATHER(msgbuf1(1), 1, MPI_DOUBLE_PRECISION,
	 +                 msgbuf2(1), 1, MPI_DOUBLE_PRECISION,
	 +                 0, MPI_COMM_WORLD,ierr)
		msgbuf1(1) = 0.d0

		if (myid.eq.0) then
			temp = 0.d0
			l = 0
			do i = 0 , numprocs-1
				temp = temp + msgbuf2(l+1)
				msgbuf2(l+1) = 0.d0
				l = l + 1
			enddo
			temp  = Tfac * temp
			tscal = 1.d0/dsqrt( temp )
		endif

		call MPI_BCAST(tscal, 1, MPI_DOUBLE_PRECISION,
	 +                0, MPI_COMM_WORLD, ierr)

		do ip = 1 , Ncore
			ip2 = (ip-1) * WPVIS
			IRVG(ip2+6) = IRVG(ip2+6) * tscal
			IRVG(ip2+7) = IRVG(ip2+7) * tscal
			IRVG(ip2+8) = IRVG(ip2+8) * tscal
		enddo

		do ip = 1,Ncore
			ip2=(ip-1) * VIS
			ip3=(ip-1) * WPVIS
			vXYZ(ip2+1) = IRVG(ip3+6)
			vXYZ(ip2+2) = IRVG(ip3+7)
			vXYZ(ip2+3) = IRVG(ip3+8)
		enddo

		return
		end
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
c     	SUBROUTINE:
c     	The pseudorandom number generator.
c     	This subroutine has been taken as it is.
c
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
		Subroutine Ran_Init(Ijkl)
		Implicit None
		Integer Ijkl
Cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
C       This Is The Initialization Routine For Ran_Uniform()           C
C       Note: The Seed Variable Should Be In The Range 0<=900 000 000  C
Cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
		Integer I,J,K,L,Ij,Kl
		Integer Ii,Jj,M
		Double Precision S,T
		Double Precision U(97),C,Cd,Cm
		Integer I97,J97
		Logical Initialised
		Common /Raset1/ U,C,Cd,Cm,I97,J97,Initialised

		Initialised = .False.

		If( Ijkl .Lt. 0  .Or.  Ijkl .Gt. 900 000 000) Then
			Write(6,*) 'The Random Number Seed Must Have A Value ',
	 &                 'Between 0 And 900 000 000'
			Stop
		Endif

		Ij = Ijkl / 30082
		Kl = Ijkl - 30082 * Ij
		I = Mod(Ij/177,177) + 2
		J = Mod(Ij    ,177) + 2
		K = Mod(Kl/169,178) + 1
		L = Mod(Kl,   169)

		Do Ii = 1,97
			S = 0.0d0
			T = 0.5d0
			Do Jj = 1,24
				M = Mod(Mod(I*J,179)*K,179)
				I = J
				J = K
				K = M
				L = Mod(53*L+1,169)
				If (Mod(L*M,64) .Ge. 32) S = S + T
				T = 0.5d0 * T
			End Do
			U(Ii) = S
		End Do

		C           = 362436.0d0 / 16777216.0d0
		Cd          = 7654321.0d0 / 16777216.0d0
		Cm          = 16777213.0d0 /16777216.0d0
		I97         = 97
		J97         = 33
		Initialised = .True.

		Return
		End
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
		Function Ran_Uniform()
		Implicit None
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c     Generate A Uniformly Distributed Randomnumber Between 0 And 1    c
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
		double precision U(97),C,Cd,Cm,Uni,Ran_Uniform
		integer I97,J97
		logical Initialised
		common /Raset1/ U,C,Cd,Cm,I97,J97,Initialised

		if (.Not. Initialised) then
			write(6,*)'Ran_Uniform:Initialise Ran_Uniform With Ran_Init'
			stop
		endif

		Uni = U(I97) - U(J97)
		if( Uni .Lt. 0.0d0 ) Uni = Uni + 1.0d0
		U(I97) = Uni
		I97 = I97 - 1
		if(I97 .Eq. 0) I97 = 97
		J97 = J97 - 1
		if(J97 .Eq. 0) J97 = 97
		C = C - Cd
		if( C .Lt. 0.0d0 ) C = C + Cm
		Uni = Uni - C
		if( Uni .Lt. 0.0d0 ) Uni = Uni + 1.0d0
		Ran_Uniform = Uni

		return
		end
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
c                  copy coordinates to upper process                   c
c
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
		subroutine EXTVOL ! 
		implicit none
		include 'mpif.h'
		include 'solvent_with_salts.h'

		nclx = dLx
		nclz = dLz

		if(myid.ne.numprocs-1) then
			ncly = INT(dLy/nynodes)
		else 
			ncly= dLy-(numprocs-1)*INT(dLy/nynodes)
		endif

		ncly1 = ncly+1
		nclc1 = nclx*ncly1*nclz

		do ip = 1, MAXCSG
			CXYS(ip) = 0
		enddo

		do ip = 1,MAXVSG
			msgbuf1(ip) = 0
			msgbuf2(ip) = 0
		enddo

		do ip = 1,Ncore

			ip2=(ip-1)*WPVIS
			ip3=(ip-1)*CPVIS

			rX0 = 0
			if (IRVG(ip2+3).ge.0) then
				rX0 = AINT(IRVG(ip2+3))
			else if (IRVG(ip2+3).lt.0) then
				rX0 = AINT(IRVG(ip2+3))-1
			endif
			CXYS(ip3+1) = rX0 + INT(dLx/2)+1

			rY0 = 0
			if (IRVG(ip2+4).ge.0) then
				rY0 = AINT(IRVG(ip2+4))
			else if (IRVG(ip2+4).lt.0) then
				rY0 = AINT(IRVG(ip2+4))-1
			endif
			if(myid.ne.numprocs-1) then
				CXYS(ip3+2) = rY0 + INT(dLy/2)+1-myid * ncly
			else
				CXYS(ip3+2)= rY0+1-(dLy/2-ncly)
			endif

			rZ0 = 0
			if (IRVG(ip2+5).ge.0) then
				rZ0 = AINT(IRVG(ip2+5))
			else if (IRVG(ip2+5).lt.0) then
				rZ0 = AINT(IRVG(ip2+5))-1
			endif
			CXYS(ip3+3) = rZ0 + INT(dLz/2) + 1

		enddo

		l=0
		do ip = 1 , Ncore
			ip2 = (ip-1) * WPVIS
			ip3 = (ip-1) * CPVIS
			ip4 = (ip-1) * VIS

			if (CXYS(ip3+2).eq.1) then
				msgbuf1(l+1) = IRVG(ip2+2) ! global_id
				msgbuf1(l+2) = IRVG(ip2+3) ! rX
				msgbuf1(l+3) = IRVG(ip2+4) ! rY
				msgbuf1(l+4) = IRVG(ip2+5) ! rZ
				msgbuf1(l+5) = IRVG(ip2+6) ! vX
				msgbuf1(l+6) = IRVG(ip2+7) ! vY
				msgbuf1(l+7) = IRVG(ip2+8) ! vZ
				msgbuf1(l+8) = IRVG(ip2+9) ! n_flag
				msgbuf1(l+9) = IRVG(ip2+10)! charge
				msgbuf1(l+10)= CXYS(ip3+1) ! Cx
				msgbuf1(l+11)= CXYS(ip3+3) ! Cz
				msgbuf1(l+12)= vXYZ(ip4+1) ! vX0
				msgbuf1(l+13)= vXYZ(ip4+2) ! vY0
				msgbuf1(l+14)= vXYZ(ip4+3) ! vZ0
				l=l+14
			endif
		enddo

		msglen1 = l
		if (myid.eq.0) then
			down = numprocs-1
			up = myid+1
		elseif(myid.eq.numprocs-1) then
			down = myid-1
			up = 0
		else
			down = myid-1
			up = myid+1
		endif

		if (down.eq.0) then
			down1=numprocs-1
			up1=up+1
		else if (up.eq.numprocs-1) then
			down1=down-1
			up1=0
		else
			down1=down-1
			up1=up+1
		endif

		call MPI_SENDRECV(msglen1, 1, MPI_INTEGER, down,
	 +                 0, msglen2, 1, MPI_INTEGER, up,
	 +                 0, MPI_COMM_WORLD, status,ierr)
		call MPI_SENDRECV(msgbuf1(1), msglen1, MPI_DOUBLE_PRECISION,
	 +                  down, 1, msgbuf2(1), msglen2,
	 +                  MPI_DOUBLE_PRECISION, up, 1,
	 +                  MPI_COMM_WORLD, status, ierr)

		do i = 1, msglen1
			msgbuf1(i) =0.d0
		enddo
		msglen1 = 0

		l = 0
		do ip = 1, int(msglen2/14)
			ip2 = (Ncore+ip-1) * WPVIS
			ip3 = (Ncore+ip-1) * CPVIS
			ip4 = (Ncore+ip-1) * VIS
			IRVG(ip2+1) = Ncore+ip     ! id_local
			IRVG(ip2+2) = msgbuf2(l+1) ! id_global
			IRVG(ip2+3) = msgbuf2(l+2) ! rX
			IRVG(ip2+4) = msgbuf2(l+3) ! rY
			IRVG(ip2+5) = msgbuf2(l+4) ! rZ
			IRVG(ip2+6) = msgbuf2(l+5) ! vX
			IRVG(ip2+7) = msgbuf2(l+6) ! vY
			IRVG(ip2+8) = msgbuf2(l+7) ! vZ
			IRVG(ip2+9) = msgbuf2(l+8) ! n_flag  
			IRVG(ip2+10)= msgbuf2(l+9) ! charge  
			CXYS(ip3+1) = msgbuf2(l+10)! Cx
			CXYS(ip3+2) = ncly1        ! Cy
			CXYS(ip3+3) = msgbuf2(l+11)! Cz
			vXYZ(ip4+1) = msgbuf2(l+12)! vX0
			vXYZ(ip4+2) = msgbuf2(l+13)! vY0
			vXYZ(ip4+3) = msgbuf2(l+14)! vZ0
			msgbuf2(l+1) = 0.d0
			msgbuf2(l+2) = 0.d0
			msgbuf2(l+3) = 0.d0
			msgbuf2(l+4) = 0.d0
			msgbuf2(l+5) = 0.d0
			msgbuf2(l+6) = 0.d0
			msgbuf2(l+7) = 0.d0
			msgbuf2(l+8) = 0.d0
			msgbuf2(l+9) = 0.d0
			msgbuf2(l+10)= 0.d0
			msgbuf2(l+11)= 0.d0
			msgbuf2(l+12)= 0.d0
			msgbuf2(l+13)= 0.d0
			msgbuf2(l+14)= 0.d0
			l=l+14 
		enddo

		Nmx = Ncore + int(msglen2/14)
		do ip = 1, Nmx
			ip3 = (ip-1) * CPVIS
			ix = CXYS(ip3+1)
			iy = CXYS(ip3+2)
			iz = CXYS(ip3+3)
			CXYS(ip3+4) = 1 
	 +    		+ mod(ix-1+nclx*ncly1*nclz,nclx)
	 +      	+ mod(iy-1+nclx*ncly1*nclz,ncly1)*nclx
	 +      	+ mod(iz-1+nclx*ncly1*nclz,nclz)*nclx*ncly1
		enddo
		return 
		end
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
		subroutine MAPS
		implicit none
		include'mpif.h'
		include'solvent_with_salts.h'

		icell(ix,iy,iz) = 1
	 +       + mod(ix-1+nclx*ncly1*nclz,nclx)
	 +       + mod(iy-1+nclx*ncly1*nclz,ncly1)*nclx
	 +       + mod(iz-1+nclx*ncly1*nclz,nclz)*nclx*ncly1

		do i = 1 ,mapx
			map(i) = 0
		enddo

		do iz = 1,nclz
			do iy = 1,ncly
				do ix = 1,nclx
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
					map(imap+10) = icell(ix+1, iy+1, iz+1)
					map(imap+11) = icell(ix  , iy+1, iz+1)
					map(imap+12) = icell(ix-1, iy+1, iz+1)
					map(imap+13) = icell(ix  , iy  , iz+1)
				enddo
			enddo
		enddo

		return
		end
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
		subroutine LINKS
		implicit none
		include'mpif.h'
		include'solvent_with_salts.h'
		do icell = 1 , nclc1
			head(icell) = 0
		enddo
		do ip = 1, Nmp
			list(ip) = 0
		enddo
		do ip = 1, Nmx
			ip3 = (ip-1) * CPVIS
			icell = CXYS(ip3+4)
			list(ip) = head(icell)
			head(icell) = ip
		enddo

		return
		end
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
		subroutine F_CDR_EE
		implicit none
		include 'mpif.h'
		include 'solvent_with_salts.h'

		do ip = 1, MAXFSG
			FCDR(ip) = 0
		enddo

		do ip = 1, MAXZSG
			FZ(ip) = 0
		enddo

		do ip = 1, MAXPSG
			VPXY(ip) = 0
		enddo

		do icell = 1 , nclc1
			ip = head(icell)
 1000    	if (ip.gt.0) then

				ip2 = (ip-1) * WPVIS
				ip3 = (ip-1) * CPVIS
				ip4 = (ip-1) * FPVIS
				ip5 = (ip-1) * VIS
				ip6 = (ip-1) * PPVIS
				ip7 = (ip-1) * EPVIS
		
				rXi = IRVG(ip2+3)
				rYi = IRVG(ip2+4)
				rZi = IRVG(ip2+5)
		
				nfi = IRVG(ip2+9)
		
				vX0i= vXYZ(ip5+1)
				vY0i= vXYZ(ip5+2)
				vZ0i= vXYZ(ip5+3)
		
				FcXi= FCDR(ip4+1)
				FcYi= FCDR(ip4+2)
				FcZi= FCDR(ip4+3)
		
				FdXi= FCDR(ip4+4)
				FdYi= FCDR(ip4+5)
				FdZi= FCDR(ip4+6)

				FrXi= FCDR(ip4+7)
				FrYi= FCDR(ip4+8)
				FrZi= FCDR(ip4+9)
		
				vpXi= VPXY(ip6+1)
				vpYi= VPXY(ip6+2)
				vpZi= VPXY(ip6+3)

				jp = list(ip)

2000    		if (jp.gt.0) then

					jp2 = (jp-1) * WPVIS
					jp3 = (jp-1) * CPVIS
					jp4 = (jp-1) * FPVIS
					jp5 = (jp-1) * VIS
					jp7 = (jp-1) * EPVIS

					rXij = rXi - IRVG(jp2+3)
					rYij = rYi - IRVG(jp2+4)
					rZij = rZi - IRVG(jp2+5)

					rXij = rXij - dLx*dnint( rXij*dLxINV )
					rYij = rYij - dLy*dnint( rYij*dLyINV )
					rZij = rZij - dLz*dnint( rZij*dLzINV )

					rijSQ = rXij*rXij + rYij*rYij + rZij*rZij

					if (rijSQ.lt.1.d0) then

						rij = dsqrt(rijSQ)
						rijINV = 1.d0/rij
						omega = 1.d0 -rij
			
						eXij = rXij * rijINV
						eYij = rYij * rijINV
						eZij = rZij * rijINV
			
						vX0ij = vX0i - vXYZ(jp5+1)
						vY0ij = vY0i - vXYZ(jp5+2)
						vZ0ij = vZ0i - vXYZ(jp5+3)
			
						nfj = IRVG(jp2+9)

						if ((CXYS(ip3+2).eq.ncly1).and.
	 +                      (CXYS(jp3+2).eq.ncly1)) then
							Fcfac=0
							Fdfac=0
							Frfac=0
						else
							Fcfac = omega ! Concervative force
							Fdfac = omega*omega*( eXij*vX0ij +
	 +              				eYij*vY0ij + eZij*vZ0ij )
							Frfac = omega*rfac*(2.*Ran_Uniform()-1.)
						endif

						FcXij = Fcfac * eXij * alpha(nfi,nfj)
						FcYij = Fcfac * eYij * alpha(nfi,nfj)
						FcZij = Fcfac * eZij * alpha(nfi,nfj)

						vpXi = vpXi + FcXij * rXij
						vpYi = vpYi + FcYij * rYij
						vpZi = vpZi + FcZij * rZij
			
						FdXij = Fdfac * eXij
						FdYij = Fdfac * eYij
						FdZij = Fdfac * eZij
			
						FrXij = Frfac * eXij
						FrYij = Frfac * eYij
						FrZij = Frfac * eZij
			
						FcXi = FcXi + FcXij ! Update forces
						FcYi = FcYi + FcYij
						FcZi = FcZi + FcZij

						FdXi = FdXi + FdXij
						FdYi = FdYi + FdYij
						FdZi = FdZi + FdZij
			
						FrXi = FrXi + FrXij
						FrYi = FrYi + FrYij
						FrZi = FrZi + FrZij

						FCDR(jp4+1) = FCDR(jp4+1) - FcXij
						FCDR(jp4+2) = FCDR(jp4+2) - FcYij
						FCDR(jp4+3) = FCDR(jp4+3) - FcZij
			
						FCDR(jp4+4) = FCDR(jp4+4) - FdXij
						FCDR(jp4+5) = FCDR(jp4+5) - FdYij
						FCDR(jp4+6) = FCDR(jp4+6) - FdZij
			
						FCDR(jp4+7) = FCDR(jp4+7) - FrXij
						FCDR(jp4+8) = FCDR(jp4+8) - FrYij
						FCDR(jp4+9) = FCDR(jp4+9) - FrZij

						! short-range part of Electrostatic force
						if((IRVG(ip2+10).ne.0.0).and.
	 +          			(IRVG(jp2+10).ne.0.0))then
							chgi = charge(ni)
							chgj = charge(nj)
							g1 = 13.87 / (4*pi)
							g2 = exp(-2.0*0.929*rij)
							g3 = 1.0 + 0.929*rij
							g4 = 2.0*0.929*rij
			
							ab = erfc(alpha_chg*rij)*rijINV
								!real_eng = real_eng + chgi*chgj*ab  ! *(g1*(1.0-g2*g3))
							FenufRF(ip7+4)=FenufRF(ip7+4)+chgi*chgj*ab
							FenufRF(jp7+4)=FenufRF(ip7+4)+chgi*chgj*ab

							aa1 = 2.0*alpha_chg*piINV
							aa = aa1 * exp(-alpha_chg*alpha_chg*rijSQ)
							Freal = chgi*chgj*(aa+ab)*rrijINV
	 +                      		*(g1*(1.0-g2-g2*g3*g4))

							FrealXij = Freal*rXij
							FrealYij = Freal*rYij
							FrealZij = Freal*rZij
							FenufRF(ip7+1) = FenufRF(ip7+1) + FrealXij
							FenufRF(ip7+2) = FenufRF(ip7+2) + FrealYij
							FenufRF(ip7+3) = FenufRF(ip7+3) + FrealZij
							FenufRF(jp7+1) = FenufRF(jp7+1) - FrealXij
							FenufRF(jp7+2) = FenufRF(jp7+2) - FrealYij
							FenufRF(jp7+3) = FenufRF(jp7+3) - FrealZij
							VPXY(ip6+1) = VPXY(ip6+1) + FrealXij*rXij
							VPXY(ip6+2) = VPXY(ip6+2) + FrealYij*rYij
							VPXY(ip6+3) = VPXY(ip6+3) + FrealZij*rZij
						endif
					endif

					jp = list(jp)
				go to 2000

			endif

			jcell0 = 13 * (icell-1)

			do nabor = 1 ,13
				jcell = map(jcell0+nabor)
				jp = head(jcell)
 3000     		if (jp.gt.0) then

					jp2 = (jp-1) * WPVIS
					jp3 = (jp-1) * CPVIS
					jp4 = (jp-1) * FPVIS
					jp5 = (jp-1) * VIS
					jp7 = (jp-1) * EPVIS

					rXij = rXi - IRVG(jp2+3)
					rYij = rYi - IRVG(jp2+4)
					rZij = rZi - IRVG(jp2+5)

					rXij = rXij - dLx*dnint( rXij*dLxINV )
					rYij = rYij - dLy*dnint( rYij*dLyINV )
					rZij = rZij - dLz*dnint( rZij*dLzINV )
					rijSQ = rXij*rXij + rYij*rYij + rZij*rZij

					if (rijSQ.lt.1.d0) then
						rij = dsqrt(rijSQ)
						rijINV = 1.d0/rij
						omega = 1.d0 -rij

						eXij = rXij * rijINV
						eYij = rYij * rijINV
						eZij = rZij * rijINV

						vX0ij = vX0i - vXYZ(jp5+1)
						vY0ij = vY0i - vXYZ(jp5+2)
						vZ0ij = vZ0i - vXYZ(jp5+3)

						nfj = IRVG(jp2+9)

						if ((CXYS(ip3+2).eq.ncly1).and.
	 +                      (CXYS(jp3+2).eq.ncly1)) then
							Fcfac=0
							Fdfac=0
							Frfac=0
						else
							Fcfac = omega ! Concervative force
							Fdfac = omega * omega *( eXij*vX0ij +
	 +                      		eYij*vY0ij + eZij*vZ0ij)
							Frfac = omega*rfac*(2.*Ran_Uniform()-1.)
						endif

						FcXij = Fcfac * eXij * alpha(nfi,nfj)
						FcYij = Fcfac * eYij * alpha(nfi,nfj)
						FcZij = Fcfac * eZij * alpha(nfi,nfj)

						vpXi = vpXi + FcXij * rXij
						vpYi = vpYi + FcYij * rYij
						vpZi = vpZi + FcZij * rZij
						
						FdXij = Fdfac * eXij
						FdYij = Fdfac * eYij
						FdZij = Fdfac * eZij

						FrXij = Frfac * eXij
						FrYij = Frfac * eYij
						FrZij = Frfac * eZij

						FcXi = FcXi + FcXij ! Update forces
						FcYi = FcYi + FcYij
						FcZi = FcZi + FcZij

						FdXi = FdXi + FdXij
						FdYi = FdYi + FdYij
						FdZi = FdZi + FdZij

						FrXi = FrXi + FrXij
						FrYi = FrYi + FrYij
						FrZi = FrZi + FrZij
						
						FCDR(jp4+1) = FCDR(jp4+1) - FcXij
						FCDR(jp4+2) = FCDR(jp4+2) - FcYij
						FCDR(jp4+3) = FCDR(jp4+3) - FcZij

						FCDR(jp4+4) = FCDR(jp4+4) - FdXij
						FCDR(jp4+5) = FCDR(jp4+5) - FdYij
						FCDR(jp4+6) = FCDR(jp4+6) - FdZij

						FCDR(jp4+7) = FCDR(jp4+7) - FrXij
						FCDR(jp4+8) = FCDR(jp4+8) - FrYij
						FCDR(jp4+9) = FCDR(jp4+9) - FrZij

c						short-range part of Electrostatic force
						if((IRVG(ip2+10).ne.0.0).and.
	 +          			(IRVG(jp2+10).ne.0.0))then
							chgi = charge(ni)
							chgj = charge(nj)
							g1 = 13.87 / (4*pi)
							g2 = exp(-2.0*0.929*rij)
							g3 = 1.0 + 0.929*rij
							g4 = 2.0*0.929*rij

							aa = erfc(alpha_chg*rij)*rijINV
							ab = aa*(g1*(1.0-g2*g3))
							real_eng = real_eng + chgi*chgj*ab

							ac = 2.0*alpha_chg*piINV
							ad = ac * exp(-alpha_chg*alpha_chg*rijSQ)
							ae = rrijINV*(g1*(1.0-g2-g2*g3*g4))
							Freal = chgi*chgj*(ad+aa)*ae

							FrealXij = Freal*rXij
							FrealYij = Freal*rYij
							FrealZij = Freal*rZij
							FenufRF(ip7+1) = FenufRF(ip7+1) + FrealXij
							FenufRF(ip7+2) = FenufRF(ip7+2) + FrealYij
							FenufRF(ip7+3) = FenufRF(ip7+3) + FrealZij
							FenufRF(jp7+1) = FenufRF(jp7+1) - FrealXij
							FenufRF(jp7+2) = FenufRF(jp7+2) - FrealYij
							FenufRF(jp7+3) = FenufRF(jp7+3) - FrealZij
							VPXY(ip6+1) = VPXY(ip6+1) + FrealXij*rXij
							VPXY(ip6+2) = VPXY(ip6+2) + FrealYij*rYij
							VPXY(ip6+3) = VPXY(ip6+3) + FrealZij*rZij
						endif

					endif
					jp = list(jp)
					goto 3000
				endif
			enddo

			FCDR(ip4+1) = FcXi
			FCDR(ip4+2) = FcYi
			FCDR(ip4+3) = FcZi

			FCDR(ip4+4) = FdXi
			FCDR(ip4+5) = FdYi
			FCDR(ip4+6) = FdZi

			FCDR(ip4+7) = FrXi
			FCDR(ip4+8) = FrYi
			FCDR(ip4+9) = FrZi

			VPXY(ip6+1) = vpXi
			VPXY(ip6+2) = vpYi
			VPXY(ip6+3) = vpZi
			ip = list(ip)
			goto 1000
			endif
		enddo
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
		do ip = 1, Nmx
			ip4 = (ip-1)*FPVIS
			FCDR(ip4+4) = -gamma * FCDR(ip4+4)
			FCDR(ip4+5) = -gamma * FCDR(ip4+5)
			FCDR(ip4+6) = -gamma * FCDR(ip4+6)
			FCDR(ip4+7) =  sigma * FCDR(ip4+7)
			FCDR(ip4+8) =  sigma * FCDR(ip4+8)
			FCDR(ip4+9) =  sigma * FCDR(ip4+9)
		enddo
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

		do i = 1,MAXVSG
			msgbuf1(i)=0
			msgbuf2(i)=0
		enddo

		l1=0
		do ip = Ncore+1 , Nmx

			ip2 = (ip-1) * WPVIS
			ip3 = (ip-1) * CPVIS
			ip4 = (ip-1) * FPVIS
			ip5 = (ip-1) * VIS
			ip6 = (ip-1) * PPVIS
			ip7 = (ip-1) * EPVIS

			if (CXYS(ip3+2).eq.ncly1) then

				msgbuf1(l1+1) = IRVG(ip2+2) ! global_id
				msgbuf1(l1+2) = FCDR(ip4+1)
				msgbuf1(l1+3) = FCDR(ip4+2)
				msgbuf1(l1+4) = FCDR(ip4+3)
				msgbuf1(l1+5) = FCDR(ip4+4)
				msgbuf1(l1+6) = FCDR(ip4+5)
				msgbuf1(l1+7) = FCDR(ip4+6)
				msgbuf1(l1+8) = FCDR(ip4+7)
				msgbuf1(l1+9) = FCDR(ip4+8)
				msgbuf1(l1+10)= FCDR(ip4+9)
				msgbuf1(l1+11)= VPXY(ip6+1)
				msgbuf1(l1+12)= VPXY(ip6+2)
				msgbuf1(l1+13)= VPXY(ip6+3)
				msgbuf1(l1+14)= FenufRF(ip7+1)
				msgbuf1(l1+15)= FenufRF(ip7+2)
				msgbuf1(l1+16)= FenufRF(ip7+3)
				msgbuf1(l1+17)= FenufRF(ip7+4)

				IRVG(ip2+1) = 0
				IRVG(ip2+2) = 0
				IRVG(ip2+3) = 0
				IRVG(ip2+4) = 0
				IRVG(ip2+5) = 0
				IRVG(ip2+6) = 0
				IRVG(ip2+7) = 0
				IRVG(ip2+8) = 0
				IRVG(ip2+9) = 0

				FCDR(ip4+1) = 0
				FCDR(ip4+2) = 0
				FCDR(ip4+3) = 0
				FCDR(ip4+4) = 0
				FCDR(ip4+5) = 0
				FCDR(ip4+6) = 0
				FCDR(ip4+7) = 0
				FCDR(ip4+8) = 0
				FCDR(ip4+9) = 0
		
				CXYS(ip3+1) = 0
				CXYS(ip3+2) = 0
				CXYS(ip3+3) = 0
				CXYS(ip3+4) = 0
		
				vXYZ(ip5+1)  = 0
				vXYZ(ip5+2)  = 0
				vXYZ(ip5+3)  = 0
		
				VPXY(ip6+1) = 0
				VPXY(ip6+2) = 0
				VPXY(ip6+3) = 0
		
				FenufRF(ip7+1) = 0
				FenufRF(ip7+2) = 0
				FenufRF(ip7+3) = 0
				FenufRF(ip7+4) = 0

				l1 = l1 + 17
			endif
		enddo

		Nmx=0
		msglen1 = l1
		call MPI_SENDRECV(msglen1, 1, MPI_INTEGER, up,
	 +                 7, msglen2, 1, MPI_INTEGER, down, 
	 +                 7, MPI_COMM_WORLD, status, ierr)
		call MPI_SENDRECV(msgbuf1, msglen1, MPI_DOUBLE_PRECISION, up,
	 +                 8, msgbuf2, msglen2, MPI_DOUBLE_PRECISION, down,
	 +                 8, MPI_COMM_WORLD, status, ierr)

		do i = 1 , msglen1
			msgbuf1(i) =0
		enddo

		m2 = msglen2/17

		l2 = 0
		do ip = 1 , Ncore
			ip2 = (ip-1) * WPVIS
			ip3 = (ip-1) * CPVIS
			ip4 = (ip-1) * FPVIS
			ip6 = (ip-1) * PPVIS
			ip7 = (ip-1) * EPVIS
			do i = 1 , m2
				if(IRVG(ip2+2).eq.msgbuf2(l2+1)) then
	
					FCDR(ip4+1) = FCDR(ip4+1) + msgbuf2(l2+2)
					FCDR(ip4+2) = FCDR(ip4+2) + msgbuf2(l2+3)
					FCDR(ip4+3) = FCDR(ip4+3) + msgbuf2(l2+4)
					FCDR(ip4+4) = FCDR(ip4+4) + msgbuf2(l2+5)
					FCDR(ip4+5) = FCDR(ip4+5) + msgbuf2(l2+6)
					FCDR(ip4+6) = FCDR(ip4+6) + msgbuf2(l2+7)
					FCDR(ip4+7) = FCDR(ip4+7) + msgbuf2(l2+8)
					FCDR(ip4+8) = FCDR(ip4+8) + msgbuf2(l2+9)
					FCDR(ip4+9) = FCDR(ip4+9) + msgbuf2(l2+10)
		
					VPXY(ip6+1) = VPXY(ip6+1) + msgbuf2(l2+11)
					VPXY(ip6+2) = VPXY(ip6+2) + msgbuf2(l2+12)
					VPXY(ip6+3) = VPXY(ip6+3) + msgbuf2(l2+13)
		
					FenufRF(ip7+1) = FenufRF(ip7+1) + msgbuf2(l2+14)
					FenufRF(ip7+2) = FenufRF(ip7+2) + msgbuf2(l2+15)
					FenufRF(ip7+3) = FenufRF(ip7+3) + msgbuf2(l2+16)
					FenufRF(ip7+4) = FenufRF(ip7+4) + msgbuf2(l2+17)
		
					msgbuf2(l2+1) = 0
					msgbuf2(l2+2) = 0
					msgbuf2(l2+3) = 0
					msgbuf2(l2+4) = 0
					msgbuf2(l2+5) = 0
					msgbuf2(l2+6) = 0
					msgbuf2(l2+7) = 0
					msgbuf2(l2+8) = 0
					msgbuf2(l2+9) = 0
					msgbuf2(l2+10)= 0
					msgbuf2(l2+11)= 0
					msgbuf2(l2+12)= 0
					msgbuf2(l2+13)= 0
					msgbuf2(l2+14)= 0
					msgbuf2(l2+15)= 0
					msgbuf2(l2+16)= 0
					msgbuf2(l2+17)= 0
					l2 = l2 + 17
				endif
			enddo
		enddo
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c                                                                      c
c       Force and Energy from Reciprocal Space using ENUF method       c
c                                                                      c
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

		do ip = 1, MAXVSG
			msgbuf3(ip) = 0
			msgbuf4(ip) = 0
		enddo

		ibuff = 0
		do ip = 1, Ncore
			ip2=(ip-1)*WPVIS
			if (abs(IRVG(ip2+10)).gt.0.0) then
				msgbuf3(ibuff+1) = IRVG(ip2+2)	   !id_global
				msgbuf3(ibuff+2) = IRVG(ip2+3)     !rX
				msgbuf3(ibuff+3) = IRVG(ip2+4)     !rY
				msgbuf3(ibuff+4) = IRVG(ip2+5)     !rZ
				msgbuf3(ibuff+5) = IRVG(ip2+10)    !charge
				ibuff = ibuff+5
			endif
		end do
		msglen3 = ibuff

		call MPI_GATHER(msgbuf3(1), msglen3, MPI_DOUBLE_PRECISION,
	 +                  msgbuf4(1), msglen3, MPI_DOUBLE_PRECISION,
	 +                  0, MPI_COMM_WORLD, ierr)

		ibuff = 0
		do ip = 1, msglen3
			msgbuf3(ip) = 0
			msgbuf3(ip) = 0
			msgbuf3(ip) = 0
			msgbuf3(ip) = 0
			msgbuf3(ip) = 0
		end do

		ibuff = 0
		do ip = 1, Nch
			idchg(ip) = msgbuf4(ibuff+1) !id_global
			rxchg(ip) = msgbuf4(ibuff+2) !rX
			rychg(ip) = msgbuf4(ibuff+3) !rY
			rzchg(ip) = msgbuf4(ibuff+4) !rZ
			eechg(ip) = msgbuf4(ibuff+5) !charge
			msgbuf4(ibuff+1) = 0.0
			msgbuf4(ibuff+2) = 0.0
			msgbuf4(ibuff+3) = 0.0
			msgbuf4(ibuff+4) = 0.0
			msgbuf4(ibuff+5) = 0.0
			ibuff = ibuff + 5
		end do

		if (myid .eq. 0)then
			call FElecRecipENUF(kspace, Nch, alpha_chg, 
	 +  		dLxINV, dLyINV, dLzINV,
	 +  	    idchg, eechg, rxchg, rychg, rzchg,
	 +  	    FrecipX, FrecipY, FrecipZ, nfft_eng)
		endif

		do ip = 1, MAXVSG
			msgbuf3(ip) = 0
			msgbuf4(ip) = 0
		enddo

		do ip = 1, Nch
			msgbuf4((ip-1)*5+1) = idchg(ip)    !id_global
			msgbuf4((ip-1)*5+2) = FrecipX(ip)  !force X
			msgbuf4((ip-1)*5+3) = FrecipY(ip)  !force Y
			msgbuf4((ip-1)*5+4) = FrecipZ(ip)  !force Z
			msgbuf4((ip-1)*5+5) = nfft_eng(ip) !energy
		end do
		
		call MPI_BCAST(msgbuf4(1), Nch*5, MPI_DOUBLE_PRECISION,
	 +                0, MPI_COMM_WORLD, ierr)

		do ip = 1, Ncore
			ip2 = (ip-1) * WPVIS
			ip4 = (ip-1) * EPVIS
			if (abs(IRVG(ip2+10)).gt.0.0) then ! ions
				do ix = 1, Nch
					if (IRVG(ip2+2) .eq. int(msgbuf4((ix-1)*5+1)))then
						FenufRF(ip4+5)= msgbuf4((ix-1)*5+2) ! force
						FenufRF(ip4+6)= msgbuf4((ix-1)*5+3)
						FenufRF(ip4+7)= msgbuf4((ix-1)*5+4)
						FenufRF(ip4+8)= msgbuf4((ix-1)*5+5) ! energy
					endif
				end do
			endif
		enddo
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c                                                                      c
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
		do ip = 1 , Ncore
			ip2 = (ip-1) * WPVIS
			ip4 = (ip-1) * FPVIS
			ip5 = (ip-1) * ZPVIS
			ip7 = (ip-1) * EPVIS

			FZ(ip5+1) = FCDR(ip4+1) + FCDR(ip4+4) + FCDR(ip4+7) 
	 +     		+ FenufRF(ip7+1) + FenufRF(ip7+5)
			FZ(ip5+2) = FCDR(ip4+2) + FCDR(ip4+5) + FCDR(ip4+8) 
	 +     		+ FenufRF(ip7+2) + FenufRF(ip7+6)
			FZ(ip5+3) = FCDR(ip4+3) + FCDR(ip4+6) + FCDR(ip4+9)
	 +     		+ FenufRF(ip7+3) + FenufRF(ip7+7)

			FZ(ip5+4) = FCDR(ip4+7)
			FZ(ip5+5) = FCDR(ip4+8)
			FZ(ip5+6) = FCDR(ip4+9)
		enddo
		return

		end
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
		subroutine VCMTP
		implicit none
		include'mpif.h'
		include'solvent_with_salts.h'
	
		vCMX(myid) = 0.d0
		vCMY(myid) = 0.d0
		vCMZ(myid) = 0.d0
		tem(myid)  = 0.d0
		
		pressX(myid) = 0.d0
		pressY(myid) = 0.d0
		pressZ(myid) = 0.d0
	
		do i = 1, MAXVSG
			msgbuf1(i) =0
			msgbuf2(i) =0 
		enddo

		do ip = 1, Ncore
			ip2 = (ip-1) * WPVIS
			ip3 = (ip-1) * PPVIS
	
			vXi = 0
			vYi = 0
			vZi = 0
	
			vXi = IRVG(ip2+6)
			vYi = IRVG(ip2+7)
			vZi = IRVG(ip2+8)
	
			vCMX(myid) = vCMX(myid) + pm * vXi
			vCMY(myid) = vCMY(myid) + pm * vYi
			vCMZ(myid) = vCMZ(myid) + pm * vZi
	
			tem(myid) = tem(myid) + pm*(vXi*vXi+vYi*vYi+vZi*vZi)
			pressX(myid) = pressX(myid) + VPXY(ip3+1)
			pressY(myid) = pressY(myid) + VPXY(ip3+2)
			pressZ(myid) = pressZ(myid) + VPXY(ip3+3)
	
			VPXY(ip3+1) = 0
			VPXY(ip3+2) = 0
			VPXY(ip3+3) = 0
		enddo

		msgbuf1(1) = vCMX(myid)
		msgbuf1(2) = vCMY(myid)
		msgbuf1(3) = vCMZ(myid)
		msgbuf1(4) = tem(myid)
		msgbuf1(5) = pressX(myid)
		msgbuf1(6) = pressY(myid)
		msgbuf1(7) = pressZ(myid)

		vCMX(myid) = 0
		vCMY(myid) = 0
		vCMZ(myid) = 0
		tem(myid)  = 0
		pressX(myid) = 0
		pressY(myid) = 0
		pressZ(myid) = 0

		call MPI_GATHER(msgbuf1, 7, MPI_DOUBLE_PRECISION,
	 +                	msgbuf2, 7, MPI_DOUBLE_PRECISION,
	 +                  0, MPI_COMM_WORLD, ierr)

		msgbuf1(1) = 0
		msgbuf1(2) = 0
		msgbuf1(3) = 0
		msgbuf1(4) = 0
		msgbuf1(5) = 0
		msgbuf1(6) = 0
		msgbuf1(7) = 0

		if (myid.eq.0) then
			vCMXM = 0.d0
			vCMYM = 0.d0
			vCMZM = 0.d0
			vCM   = 0.d0
			temp  = 0.d0
			pX    = 0.d0
			pY    = 0.d0
			pZ    = 0.d0
			press = 0.d0

			l = 0

			do i = 0 , numprocs-1
				vCMXM = vCMXM + msgbuf2(l+1)
				vCMYM = vCMYM + msgbuf2(l+2)
				vCMZM = vCMZM + msgbuf2(l+3)
				temp  = temp  + msgbuf2(l+4)
				pX    = pX    + msgbuf2(l+5)
				pY    = pY    + msgbuf2(l+6)
				pZ    = pZ    + msgbuf2(l+7)
				msgbuf2(l+1) = 0
				msgbuf2(l+2) = 0
				msgbuf2(l+3) = 0
				msgbuf2(l+4) = 0
				msgbuf2(l+5) = 0
				msgbuf2(l+6) = 0
				msgbuf2(l+7) = 0
				l=l+7
			enddo

			vCMXM = vCMXM * dNpINV * pmINV
			vCMYM = vCMYM * dNpINV * pmINV
			vCMZM = vCMZM * dNpINV * pmINV

			vCM = dNpINV * pmINV
	 +  		* dsqrt(vCMXM*vCMXM+vCMYM*vCMYM+vCMZM*vCMZM)

			temp = Tfac * temp

			press = (pX+pY+pZ)/(3*vol)+3.d0

c      		print*,vCMXM,vCMYM,vCMZM

			write(3,1000) ! Write SNAPSHOTs
	 +           it, sngl(dble(it)*dt),
	 +           vCM, sngl(temp), sngl(press)
			call FLUSH(3)
 1000  		format(i10, 1x,f10.4, 1x,d12.6, 1x,f14.6, 1x,f14.6)
		endif
		return
		end
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
		subroutine INTV
		implicit none
		include 'solvent_with_salts.h'

		do ip = 1, Ncore
			ip2 = (ip-1) * WPVIS
			ip5 = (ip-1) * ZPVIS
			IRVG(ip2+6) = IRVG(ip2+6)
	 + 			+ pmINV*((FZ(ip5+1)-FZ(ip5+4))*dtH + FZ(ip5+4)*dtRTH)
			IRVG(ip2+7) = IRVG(ip2+7)
	 + 			+ pmINV*((FZ(ip5+2)-FZ(ip5+5))*dtH + FZ(ip5+5)*dtRTH)
			IRVG(ip2+8) = IRVG(ip2+8)
	 + 			+ pmINV*((FZ(ip5+3)-FZ(ip5+6))*dtH + FZ(ip5+6)*dtRTH)
		enddo

		return
		end
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
		subroutine INTR
		implicit none
		include 'solvent_with_salts.h'

		do ip = 1, Ncore
			ip2 = (ip-1) * WPVIS
			IRVG(ip2+3) = IRVG(ip2+3) + IRVG(ip2+6) * dt
			IRVG(ip2+4) = IRVG(ip2+4) + IRVG(ip2+7) * dt
			IRVG(ip2+5) = IRVG(ip2+5) + IRVG(ip2+8) * dt
		enddo
	
		do ip = 1 , Ncore
			ip2 = (ip-1) * WPVIS
			rXi = IRVG(ip2+3)
			rYi = IRVG(ip2+4)
			rZi = IRVG(ip2+5)
			IRVG(ip2+3) = rXi - dLx*dnint( rXi*dLxINV )
			IRVG(ip2+4) = rYi - dLy*dnint( rYi*dLyINV )
			IRVG(ip2+5) = rZi - dLz*dnint( rZi*dLzINV )
		enddo

		return
		end
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
		subroutine MIGRATE
		implicit none
		include 'mpif.h'
		include 'solvent_with_salts.h'

		do ip = 1,Nmp
			myrxyzs(ip)=0
		enddo
		do i =1,MAXVSG
			msgbuf1(i)=0
			msgbuf2(i)=0
		enddo

		n=0 
		l1=0
		l2=0

		msglen1=0
		msglen2=0

		do ip = 1, Ncore

			ip2 = (ip-1) * WPVIS
			ip4 = (ip-1) * VIS
	
			rXi = IRVG(ip2+3)
			rYi = IRVG(ip2+4)
			rZi = IRVG(ip2+5)
	
			IRVG(ip2+3) = rXi - dLx*dnint( rXi*dLxINV )
			IRVG(ip2+4) = rYi - dLy*dnint( rYi*dLyINV )
			IRVG(ip2+5) = rZi - dLz*dnint( rZi*dLzINV )

			if(IRVG(ip2+4).lt.(-dLy/2+(numprocs-1)*DINT(dLy/nynodes)))then
				myrxyzs(ip) = DINT((IRVG(ip2+4)+dLy/2)/DINT(dLy/nynodes))
			else
				myrxyzs(ip) =numprocs-1
			endif

			if (myrxyzs(ip).eq.myid) then
				n=n+1

				jp2 = (n-1) * WPVIS
				jp4 = (n-1) * VIS

				IRVG(jp2+1) = n
				IRVG(jp2+2) = IRVG(ip2+2)
				IRVG(jp2+3) = IRVG(ip2+3)
				IRVG(jp2+4) = IRVG(ip2+4)
				IRVG(jp2+5) = IRVG(ip2+5)
				IRVG(jp2+6) = IRVG(ip2+6)
				IRVG(jp2+7) = IRVG(ip2+7)
				IRVG(jp2+8) = IRVG(ip2+8)
				IRVG(jp2+9) = IRVG(ip2+9)
				vXYZ(jp4+1) = vXYZ(ip4+1)
				vXYZ(jp4+2) = vXYZ(ip4+2)
				vXYZ(jp4+3) = vXYZ(ip4+3)

			elseif(myrxyzs(ip).eq.down) then

				msgbuf1(l1+1) = IRVG(ip2+2)
				msgbuf1(l1+2) = IRVG(ip2+3)
				msgbuf1(l1+3) = IRVG(ip2+4)
				msgbuf1(l1+4) = IRVG(ip2+5)
				msgbuf1(l1+5) = IRVG(ip2+6)
				msgbuf1(l1+6) = IRVG(ip2+7)
				msgbuf1(l1+7) = IRVG(ip2+8)
				msgbuf1(l1+8) = IRVG(ip2+9)
				msgbuf1(l1+9) = IRVG(ip2+10)
				msgbuf1(l1+10)= vXYZ(ip4+1)
				msgbuf1(l1+11)= vXYZ(ip4+2)
				msgbuf1(l1+12)= vXYZ(ip4+3)
				l1 = l1 + 12
				IRVG(ip2+1) = 0
				IRVG(ip2+2) = 0
				IRVG(ip2+3) = 0
				IRVG(ip2+4) = 0
				IRVG(ip2+5) = 0
				IRVG(ip2+6) = 0
				IRVG(ip2+7) = 0
				IRVG(ip2+8) = 0
				IRVG(ip2+9) = 0
				IRVG(ip2+10) = 0
				vXYZ(ip4+1) = 0
				vXYZ(ip4+2) = 0
				vXYZ(ip4+3) = 0

			elseif (myrxyzs(ip).eq.up) then

				msgbuf2(l2+1) = IRVG(ip2+2)
				msgbuf2(l2+2) = IRVG(ip2+3)
				msgbuf2(l2+3) = IRVG(ip2+4)
				msgbuf2(l2+4) = IRVG(ip2+5)
				msgbuf2(l2+5) = IRVG(ip2+6)
				msgbuf2(l2+6) = IRVG(ip2+7)
				msgbuf2(l2+7) = IRVG(ip2+8)
				msgbuf2(l2+8) = IRVG(ip2+9)
				msgbuf2(l2+9) = IRVG(ip2+10)
				msgbuf2(l2+10)= vXYZ(ip4+1)
				msgbuf2(l2+11)= vXYZ(ip4+2)
				msgbuf2(l2+12)= vXYZ(ip4+3)

				l2 = l2 + 12

				IRVG(ip2+1) = 0
				IRVG(ip2+2) = 0
				IRVG(ip2+3) = 0
				IRVG(ip2+4) = 0
				IRVG(ip2+5) = 0
				IRVG(ip2+6) = 0
				IRVG(ip2+7) = 0
				IRVG(ip2+8) = 0
				IRVG(ip2+9) = 0
				IRVG(ip2+10) = 0
				vXYZ(ip4+1) = 0
				vXYZ(ip4+2) = 0
				vXYZ(ip4+3) = 0
			endif
		enddo

		Ncore = n
		msglen1=l1
		msglen2=l2
		msglen11=0
		msglen22=0

		call MPI_SENDRECV(msglen1, 1, MPI_INTEGER, down,
	 +      		   1, msglen11,1, MPI_INTEGER, up,
	 +      		   1, MPI_COMM_WORLD, status, ierr)
		call MPI_SENDRECV(msgbuf1(1),msglen1, MPI_DOUBLE_PRECISION,down,
	 +      		  11, msgbuf11(1),msglen11,MPI_DOUBLE_PRECISION,up,
	 +      		  11, MPI_COMM_WORLD, status, ierr)    

		call MPI_SENDRECV(msglen2, 1, MPI_INTEGER, up,
	 +                 2, msglen22,1, MPI_INTEGER, down,
	 +      		   2, MPI_COMM_WORLD, status, ierr)
		call MPI_SENDRECV(msgbuf2(1), msglen2, MPI_DOUBLE_PRECISION, up,
	 +      		 21, msgbuf22(1),msglen22,MPI_DOUBLE_PRECISION,down,
	 +     			 21, MPI_COMM_WORLD, status, ierr)  

		l1=0
		l2=0
		l1 = msglen11/12
		l2 = msglen22/12

		msglen1=0
		msglen2=0
		msglen11=0
		msglen22=0

		if (l1.ne.0) then
			l=0
			do ip = 1, l1
				ip2 = (Ncore+ip-1) * WPVIS
				ip4 = (Ncore+ip-1) * VIS

				IRVG(ip2+1) = Ncore+ip
				IRVG(ip2+2) = msgbuf11(l+1)
				IRVG(ip2+3) = msgbuf11(l+2)
				IRVG(ip2+4) = msgbuf11(l+3)
				IRVG(ip2+5) = msgbuf11(l+4)
				IRVG(ip2+6) = msgbuf11(l+5)
				IRVG(ip2+7) = msgbuf11(l+6)
				IRVG(ip2+8) = msgbuf11(l+7)
				IRVG(ip2+9) = msgbuf11(l+8)
				IRVG(ip2+10)= msgbuf11(l+9)
				vXYZ(ip4+1) = msgbuf11(l+10)
				vXYZ(ip4+2) = msgbuf11(l+11)
				vXYZ(ip4+3) = msgbuf11(l+12)
				l=l+12
			enddo
			Ncore=Ncore+l1
		endif

		if(l2.ne.0) then
			l=0
			do ip = 1,l2
				ip2 = (Ncore+ip-1) * WPVIS
				ip4 = (Ncore+ip-1) * VIS

				IRVG(ip2+1) = Ncore+ip
				IRVG(ip2+2) = msgbuf22(l+1)
				IRVG(ip2+3) = msgbuf22(l+2)
				IRVG(ip2+4) = msgbuf22(l+3)
				IRVG(ip2+5) = msgbuf22(l+4)
				IRVG(ip2+6) = msgbuf22(l+5)
				IRVG(ip2+7) = msgbuf22(l+6)
				IRVG(ip2+8) = msgbuf22(l+7)
				IRVG(ip2+9) = msgbuf22(l+8)
				IRVG(ip2+10)= msgbuf22(l+9)
				vXYZ(ip4+1) = msgbuf22(l+10)
				vXYZ(ip4+2) = msgbuf22(l+11)
				vXYZ(ip4+3) = msgbuf22(l+12)
				l=l+12
			enddo
			Ncore=Ncore+l2
		endif

		l1=0
		l2=0
		l3=0
		l4=0

		do i = 1, MAXVSG
			msgbuf1(i)=0
			msgbuf2(i)=0
			msgbuf11(i)=0
			msgbuf22(i)=0
		enddo

		return
		end
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
		subroutine SNAPSHOT
		implicit none
		include 'mpif.h'
		include 'solvent_with_salts.h'

		open (myid+Ifile+nconf)

		do ip = 1 , Ncore
			ip2 = (ip-1) * WPVIS
			write((myid+Ifile+nconf),100),it,myid,
	+    		INT(IRVG(ip2+2)),IRVG(ip2+3),IRVG(ip2+4),
	+    		IRVG(ip2+5),INT(IRVG(ip2+9))

 100   	format(1x,I10,1x,I5,1x,I10,5x,f25.16,5x,f25.16,5x,f25.16, 5x,I2)
		enddo
		close(myid+Ifile+nconf)

		return
		end
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
