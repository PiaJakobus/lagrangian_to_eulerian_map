MODULE SPH_variables
    !****************************************
  !                                       *
  !  global SPH variables; SKR 27.07.2017 *
  !                                       *
  !****************************************
  IMPLICIT NONE

  INTEGER                                     :: npart     ! particle number
  INTEGER                                     :: n1,n2     ! in principle obsolete,
  INTEGER                                     :: npm       ! but kept to use old
                                                           ! data format
  INTEGER                                     :: cur_dump  ! current data dump number
  DOUBLE PRECISION                            :: t         ! time  
  DOUBLE PRECISION,DIMENSION(:,:),ALLOCATABLE :: pos_u     ! position
  DOUBLE PRECISION,DIMENSION(:,:),ALLOCATABLE :: vel_u     ! velocity
  DOUBLE PRECISION,DIMENSION(:,:),ALLOCATABLE :: S_l       ! canon. momentum, covariant
  DOUBLE PRECISION,DIMENSION(:),  ALLOCATABLE :: u         ! spec. internal energy !!!!NEEDED AS ARRAY???
  DOUBLE PRECISION,DIMENSION(:),  ALLOCATABLE :: h         ! smoothing length
  DOUBLE PRECISION,DIMENSION(:),  ALLOCATABLE :: rhostar   ! comput. frame mass density
  DOUBLE PRECISION,DIMENSION(:),  ALLOCATABLE :: rho       ! mass density in local rest frame
  DOUBLE PRECISION,DIMENSION(:),  ALLOCATABLE :: temp      ! potential future use
  DOUBLE PRECISION,DIMENSION(:),  ALLOCATABLE :: Ye        ! potential future use
  DOUBLE PRECISION,DIMENSION(:),  ALLOCATABLE :: pmass     ! mass per SPH-particle
  DOUBLE PRECISION,DIMENSION(:),  ALLOCATABLE :: ehat      ! canon. energy per baryon
  DOUBLE PRECISION,DIMENSION(:),  ALLOCATABLE :: theta     ! generalized Lorentz factor
  DOUBLE PRECISION,DIMENSION(:),  ALLOCATABLE :: Pr        ! pressure
  DOUBLE PRECISION,DIMENSION(:),  ALLOCATABLE :: cs        ! sound velocity
  DOUBLE PRECISION,DIMENSION(:),  ALLOCATABLE :: av        ! dissipation parameter
  DOUBLE PRECISION,DIMENSION(:),  ALLOCATABLE :: Kent      ! (pseudo-)entropy parameter
  DOUBLE PRECISION,DIMENSION(:),  ALLOCATABLE :: divv      ! velocity divergence
  DOUBLE PRECISION,DIMENSION(:),  ALLOCATABLE :: PU        ! partition of unity
  DOUBLE PRECISION,DIMENSION(:),  ALLOCATABLE :: d_divX    ! rel. deviation of div(X) from 3 
  DOUBLE PRECISION,DIMENSION(:),  ALLOCATABLE :: CD_R      ! shock detector Cullen & Dehnen     
  DOUBLE PRECISION,DIMENSION(:),  ALLOCATABLE :: CD_lim    ! Cullen & Dehnen limiter    
  DOUBLE PRECISION,DIMENSION(:),  ALLOCATABLE :: K_noise   ! Noise detector

  ! in principle obsolete, but kept to use old data format
  DOUBLE PRECISION rstar,mstar,escap,tkin,tgrav,tterm


CONTAINS

  
  SUBROUTINE allocate_SPH_memory
    
    !************************************************
    !                                               *
    ! allocate memory for all SPH variables;        *
    ! also initialize all variables; SKR 22.04.2017 *
    !                                               *
    !************************************************
    
    IMPLICIT NONE
    
    INTEGER allocation_status
    
    PRINT*
    PRINT*,'...allocating memory for SPH variables...'
    PRINT*

    ! positions
    ALLOCATE(pos_u(3,npart),STAT=allocation_status)
    IF(allocation_status > 0)THEN
       PRINT*,'...allocation error for variable pos_u'
       PRINT*
       STOP
    ENDIF
    pos_u(1:3,1:npart)= 0.0D0

    ! velocities
    ALLOCATE(vel_u(3,npart),STAT=allocation_status)
    IF(allocation_status > 0)THEN
       PRINT*,'...allocation error for variable vel_u'
       PRINT*
       STOP
    ENDIF
    vel_u(1:3,1:npart)= 0.0D0

    ! covar. canonical momentum per baryon
    ALLOCATE(S_l(3,npart),STAT=allocation_status)
    IF(allocation_status > 0)THEN
       PRINT*,'...allocation error for variable S_l'
       PRINT*
       STOP
    ENDIF
    S_l(1:3,1:npart)= 0.0D0

    ! specific internal energy
    ALLOCATE(u(npart),STAT=allocation_status)
    IF(allocation_status > 0)THEN
       PRINT*,'...allocation error for variable u'
       PRINT*
       STOP
    ENDIF
    u(1:npart)= 0.0D0

    ! smoothing length
    ALLOCATE(h(npart),STAT=allocation_status)
    IF(allocation_status > 0)THEN
       PRINT*,'...allocation error for variable h'
       PRINT*
       STOP
    ENDIF
    h(1:npart)= 0.0D0

    ! computing frame mass density
    ALLOCATE(rhostar(npart),STAT=allocation_status)
    IF(allocation_status > 0)THEN
       PRINT*,'...allocation error for variable rhostar'
       PRINT*
       STOP
    ENDIF
    rhostar(1:npart)= 0.0D0

    ! local rest frame mass density
    ALLOCATE(rho(npart),STAT=allocation_status)
    IF(allocation_status > 0)THEN
       PRINT*,'...allocation error for variable rho'
       PRINT*
       STOP
    ENDIF
    rho(1:npart)= 0.0D0

    ! temperature
    ALLOCATE(temp(npart),STAT=allocation_status)
    IF(allocation_status > 0)THEN
       PRINT*,'...allocation error for variable temp'
       PRINT*
       STOP
    ENDIF
    temp(1:npart)= 0.0D0

    ! electron fraction
    ALLOCATE(Ye(npart),STAT=allocation_status)
    IF(allocation_status > 0)THEN
       PRINT*,'...allocation error for variable Ye'
       PRINT*
       STOP
    ENDIF
    Ye(1:npart)= 0.0D0

    ! mass per SPH particle
    ALLOCATE(pmass(npart),STAT=allocation_status)
    IF(allocation_status > 0)THEN
       PRINT*,'...allocation error for variable pmass'
       PRINT*
       STOP
    ENDIF
    pmass(1:npart)= 0.0D0

    ! canonical energy per baryon
    ALLOCATE(ehat(npart),STAT=allocation_status)
    IF(allocation_status > 0)THEN
       PRINT*,'...allocation error for variable ehat'
       PRINT*
       STOP
    ENDIF
    ehat(1:npart)= 0.0D0

    ! theta
    ALLOCATE(theta(npart),STAT=allocation_status)
    IF(allocation_status > 0)THEN
       PRINT*,'...allocation error for variable theta'
       PRINT*
       STOP
    ENDIF
    theta(1:npart)= 0.0D0

    ! gas pressure
    ALLOCATE(Pr(npart),STAT=allocation_status)
    IF(allocation_status > 0)THEN
       PRINT*,'...allocation error for variable Pr'
       PRINT*
       STOP
    ENDIF
    Pr(1:npart)= 0.0D0

    ! sound velocity
    ALLOCATE(cs(npart),STAT=allocation_status)
    IF(allocation_status > 0)THEN
       PRINT*,'...allocation error for variable cs'
       PRINT*
       STOP
    ENDIF
    cs(1:npart)= 0.0D0

    ! dissipation parameter
    ALLOCATE(av(npart),STAT=allocation_status)
    IF(allocation_status > 0)THEN
       PRINT*,'...allocation error for variable av'
       PRINT*
       STOP
    ENDIF
    av(1:npart)= 0.0D0

    ! (pseudo-)entropy parameter
    ALLOCATE(Kent(npart),STAT=allocation_status)
    IF(allocation_status > 0)THEN
       PRINT*,'...allocation error for variable Kent'
       PRINT*
       STOP
    ENDIF
    Kent(1:npart)= 0.0D0

    ! divv(v)
    ALLOCATE(divv(npart),STAT=allocation_status)
    IF(allocation_status > 0)THEN
       PRINT*,'...allocation error for variable divv'
       PRINT*
       STOP
    ENDIF
    divv(1:npart)= 0.0D0

    ! partition of unity
    ALLOCATE(PU(npart),STAT=allocation_status)
    IF(allocation_status > 0)THEN
       PRINT*,'...allocation error for variable PU'
       PRINT*
       STOP
    ENDIF
    PU(1:npart)= 1.0D0
    
    !deviation of div(X) from 3
    ALLOCATE(d_divX(npart),STAT=allocation_status)
    IF(allocation_status > 0)THEN
       PRINT*,'...allocation error for variable d_divX'
       PRINT*
       STOP
    ENDIF
    d_divX(1:npart)= 0.0D0

    ! shock detector Cullen & Dehnen
    ALLOCATE(CD_R(npart),STAT=allocation_status)
    IF(allocation_status > 0)THEN
       PRINT*,'...allocation error for variable CD_R'
       PRINT*
       STOP
    ENDIF
    CD_R(1:npart)= 0.0D0

    ! limiter Cullen & Dehnen                                                                                               
    ALLOCATE(CD_lim(npart),STAT=allocation_status)
    IF(allocation_status > 0)THEN
       PRINT*,'...allocation error for variable CD_lim'
       PRINT*
       STOP
    ENDIF
    CD_lim(1:npart)= 1.0D0

    ALLOCATE(K_noise(npart),STAT=allocation_status)
    IF(allocation_status > 0)THEN
       PRINT*,'...allocation error for variable K_noise'
       PRINT*
       STOP
    ENDIF
    K_noise(1:npart)= 0.0D0


  END SUBROUTINE allocate_SPH_memory

  SUBROUTINE deallocate_SPH_memory

    !*******************************************
    !                                          *
    ! deallocate memory for all SPH variables; *
    ! SKR 22.04.2017                           *
    !                                          *
    !*******************************************
    
    IMPLICIT NONE
    
    PRINT*
    PRINT*,'...de-allocating memory for SPH variables...'

    DEALLOCATE(pos_u)
    DEALLOCATE(vel_u)
    DEALLOCATE(S_l)
    DEALLOCATE(u)
    DEALLOCATE(h)
    DEALLOCATE(rhostar)
    DEALLOCATE(rho)
    DEALLOCATE(temp)
    DEALLOCATE(Ye)
    DEALLOCATE(pmass)
    DEALLOCATE(ehat)
    DEALLOCATE(theta)
    DEALLOCATE(Pr)
    DEALLOCATE(cs)
    DEALLOCATE(av)
    DEALLOCATE(Kent)
    DEALLOCATE(divv)
    DEALLOCATE(PU)
    DEALLOCATE(d_divX)
    DEALLOCATE(CD_R)
    DEALLOCATE(CD_lim)
    DEALLOCATE(K_noise)

    PRINT*,'...done...'
    PRINT*
    
  END SUBROUTINE deallocate_SPH_memory


END MODULE SPH_variables
MODULE analyze

  USE constants,     ONLY: huge

  IMPLICIT NONE

  INTEGER          nstep_per_dump,nstep_tot
  DOUBLE PRECISION Wtime_per_dump,Mbar_last

  CONTAINS
  SUBROUTINE COM(npart,pos,pmass,com_x,com_y,com_z,com_d)

    !************************************************************
    !                                                           *
    ! monitor baryonic center of mass;   SKR 28.09.2020         *
    !                                                           *
    !************************************************************

    IMPLICIT NONE

    INTEGER,          INTENT(IN)  :: npart
    DOUBLE PRECISION, INTENT(IN)  :: pos(3,npart),pmass(npart)
    DOUBLE PRECISION, INTENT(OUT) :: com_x,com_y,com_z,com_d
    
    INTEGER           a
    DOUBLE PRECISION  mass

    com_x= 0.0D0
    com_y= 0.0D0
    com_z= 0.0D0
    mass=  0.0D0

    DO a=1,npart
      com_x= com_x + pos(1,a)*pmass(a)
      com_y= com_y + pos(2,a)*pmass(a)
      com_z= com_z + pos(3,a)*pmass(a)
      mass=  mass  + pmass(a)
    ENDDO

    com_x= com_x/mass
    com_y= com_y/mass
    com_z= com_z/mass
    com_d= SQRT(com_x**2 + com_y**2 + com_z**2) 

  END SUBROUTINE COM

END MODULE analyze 



MODULE files

  IMPLICIT NONE
  
  CHARACTER(7)  :: file1,file2,file3,file4
  CHARACTER(10) :: filename
  INTEGER          iprint,inulum,igwaves,icons

END MODULE files


MODULE options

  !**************************************************
  !                                                 *
  ! options to run SPHINCS; SKR 13.01.2017          *
  !                                                 *
  !**************************************************

  IMPLICIT NONE

  CHARACTER(2)     :: coordinate_system 
  CHARACTER(3)     :: unit_system
  CHARACTER(3)     :: metric_type
  CHARACTER(4)     :: eos_type
  CHARACTER(4)     :: eos_str
  CHARACTER(4)     :: data_type
  CHARACTER(5)     :: basename

  INTEGER          :: ibh,ipos,ikernel,isteer,itiming
  INTEGER          :: isort,i_output_nsteps,integ_meth
  INTEGER          :: i_ind,irelax,ider,iweight,ndes
  INTEGER          :: icow
  DOUBLE PRECISION :: damp,av_max,av_min

  ! integration-related, but keep them here
  DOUBLE PRECISION :: dtfac,nout_r,output_Dt,tstop,dt_max

  ! boundary conditions
  INTEGER               iBD,ncore
  DOUBLE PRECISION      xBD1,xBD2
  DOUBLE PRECISION      yBD1,yBD2
  DOUBLE PRECISION      zBD1,zBD2
  DOUBLE PRECISION      rhoBD,rim
  DOUBLE PRECISION      rim_x_1,rim_x_2
  DOUBLE PRECISION      rim_y_1,rim_y_2
  DOUBLE PRECISION      rim_z_1,rim_z_2
  CHARACTER(8)          BDtype
END MODULE options




MODULE pwp_EOS

  !*****************************************************
  !                                                    *
  ! collect here "piecewise polytropic equationd of    *
  ! state (EOS);                                       *
  ! References: Read et al. (2009);                    *
  !                                                    *
  ! a reasonable selection of EOSs: Tab.2.2 in thesis  *
  ! of Tim Dietrich; the exact numbers have been taken *
  ! from Tim's table 2.2                               *
  ! --> already in SPHINCS code units                  *
  !                                                    *
  ! SKR 15.04.2021                                     *
  ! PD & SKR 12/2023: no m0c^2-convention version      *
  !                                                    *
  !*****************************************************

  IMPLICIT NONE
    !------------------!
  !-- for all EOSs --!
  !------------------!
  ! fixed (cold) crust
  DOUBLE PRECISION, PARAMETER :: Gamma0=   1.35692D0
  DOUBLE PRECISION, PARAMETER :: Gamma0_1= Gamma0 - 1.0D0
  DOUBLE PRECISION, PARAMETER :: K0=       8.948185D-2

  ! thermal Gamma
  DOUBLE PRECISION, PARAMETER :: Gamma_th=   1.750D0
  DOUBLE PRECISION, PARAMETER :: Gamma_th_1= Gamma_th - 1.0D0
  
  ! transition densities 
  DOUBLE PRECISION, PARAMETER :: rho_1= 8.12123D-4
  DOUBLE PRECISION, PARAMETER :: rho_2= 1.62040D-3

  ! transition pressures
  DOUBLE PRECISION:: p0, p1, p2

  !-------------------!
  !-- specific EOSs --!
  !-------------------!
  ! SLY
  DOUBLE PRECISION, PARAMETER :: rho_0_SLy=  2.36953D-4
  DOUBLE PRECISION, PARAMETER :: Gamma1_SLy= 3.005D0
  DOUBLE PRECISION, PARAMETER :: Gamma2_SLy= 2.988D0
  DOUBLE PRECISION, PARAMETER :: Gamma3_SLy= 2.851D0

  ! ALF2
  DOUBLE PRECISION, PARAMETER :: rho_0_ALF2=  3.15606D-4
  DOUBLE PRECISION, PARAMETER :: Gamma1_ALF2= 4.070D0
  DOUBLE PRECISION, PARAMETER :: Gamma2_ALF2= 2.411D0
  DOUBLE PRECISION, PARAMETER :: Gamma3_ALF2= 1.890D0

  ! ALF4
  DOUBLE PRECISION, PARAMETER :: rho_0_ALF4=  2.62193D-4
  DOUBLE PRECISION, PARAMETER :: Gamma1_ALF4= 3.009D0
  DOUBLE PRECISION, PARAMETER :: Gamma2_ALF4= 3.438D0
  DOUBLE PRECISION, PARAMETER :: Gamma3_ALF4= 1.803D0

  ! ENG
  DOUBLE PRECISION, PARAMETER :: rho_0_ENG=  2.99450D-4
  DOUBLE PRECISION, PARAMETER :: Gamma1_ENG= 3.514D0
  DOUBLE PRECISION, PARAMETER :: Gamma2_ENG= 3.130D0
  DOUBLE PRECISION, PARAMETER :: Gamma3_ENG= 3.168D0

  ! H4
  DOUBLE PRECISION, PARAMETER :: rho_0_H4=  1.43830D-4
  DOUBLE PRECISION, PARAMETER :: Gamma1_H4= 2.909D0
  DOUBLE PRECISION, PARAMETER :: Gamma2_H4= 2.246D0
  DOUBLE PRECISION, PARAMETER :: Gamma3_H4= 2.144D0

  ! MPA1
  DOUBLE PRECISION, PARAMETER :: rho_0_MPA1=  2.71930D-4
  DOUBLE PRECISION, PARAMETER :: Gamma1_MPA1= 3.446D0
  DOUBLE PRECISION, PARAMETER :: Gamma2_MPA1= 3.572D0
  DOUBLE PRECISION, PARAMETER :: Gamma3_MPA1= 2.887D0

  ! MS1
  DOUBLE PRECISION, PARAMETER :: rho_0_MS1=  1.52594D-4
  DOUBLE PRECISION, PARAMETER :: Gamma1_MS1= 3.224D0
  DOUBLE PRECISION, PARAMETER :: Gamma2_MS1= 3.033D0
  DOUBLE PRECISION, PARAMETER :: Gamma3_MS1= 1.325D0

  ! MS1b
  DOUBLE PRECISION, PARAMETER :: rho_0_MS1b=  1.84169D-4
  DOUBLE PRECISION, PARAMETER :: Gamma1_MS1b= 3.456D0
  DOUBLE PRECISION, PARAMETER :: Gamma2_MS1b= 3.011D0
  DOUBLE PRECISION, PARAMETER :: Gamma3_MS1b= 1.425D0

  ! AP3
  DOUBLE PRECISION, PARAMETER :: rho_0_AP3=  2.61888D-4
  DOUBLE PRECISION, PARAMETER :: Gamma1_AP3= 3.166D0
  DOUBLE PRECISION, PARAMETER :: Gamma2_AP3= 3.573D0
  DOUBLE PRECISION, PARAMETER :: Gamma3_AP3= 3.281D0

  ! WFF1
  DOUBLE PRECISION, PARAMETER :: rho_0_WFF1=  2.85177D-4
  DOUBLE PRECISION, PARAMETER :: Gamma1_WFF1= 2.519D0
  DOUBLE PRECISION, PARAMETER :: Gamma2_WFF1= 3.791D0
  DOUBLE PRECISION, PARAMETER :: Gamma3_WFF1= 3.660D0

  ! WFF2
  DOUBLE PRECISION, PARAMETER :: rho_0_WFF2=  2.70845D-4
  DOUBLE PRECISION, PARAMETER :: Gamma1_WFF2= 2.888D0
  DOUBLE PRECISION, PARAMETER :: Gamma2_WFF2= 3.475D0
  DOUBLE PRECISION, PARAMETER :: Gamma3_WFF2= 3.517D0

  ! GNH3
  DOUBLE PRECISION, PARAMETER :: rho_0_GNH3=  1.08017D-4   
  DOUBLE PRECISION, PARAMETER :: Gamma1_GNH3= 2.664D0
  DOUBLE PRECISION, PARAMETER :: Gamma2_GNH3= 2.194D0
  DOUBLE PRECISION, PARAMETER :: Gamma3_GNH3= 2.304D0

  ! APR4 --> see Hotokezaka et al. (2013b), Tab.1
  DOUBLE PRECISION, PARAMETER :: rho_0_APR4=  2.45191D-4   
  DOUBLE PRECISION, PARAMETER :: Gamma1_APR4= 2.830D0
  DOUBLE PRECISION, PARAMETER :: Gamma2_APR4= 3.445D0
  DOUBLE PRECISION, PARAMETER :: Gamma3_APR4= 3.348D0

  ! invent one: "haso" for "hard --> soft"
  DOUBLE PRECISION, PARAMETER :: rho_0_haso=  1.6119D-4 ! = 10^14 g/ccm   
  DOUBLE PRECISION, PARAMETER :: Gamma1_haso= 2.75D0
  DOUBLE PRECISION, PARAMETER :: Gamma2_haso= 2.75D0
  DOUBLE PRECISION, PARAMETER :: Gamma3_haso= 2.75D0

  ! soft single polytrope
  DOUBLE PRECISION, PARAMETER :: rho_0_soft=  1.6119D-8 ! = 10^10 g/ccm
  DOUBLE PRECISION, PARAMETER :: Gamma1_soft= 2.0D0
  DOUBLE PRECISION, PARAMETER :: Gamma2_soft= 2.0D0
  DOUBLE PRECISION, PARAMETER :: Gamma3_soft= 2.0D0

  ! stiff single polytrope
  DOUBLE PRECISION, PARAMETER :: rho_0_stiff=  1.6119D-8 ! = 10^10 g/ccm
  DOUBLE PRECISION, PARAMETER :: Gamma1_stiff= 2.75D0
  DOUBLE PRECISION, PARAMETER :: Gamma2_stiff= 2.75D0
  DOUBLE PRECISION, PARAMETER :: Gamma3_stiff= 2.75D0
  
  ! for general pwp-EOS
  DOUBLE PRECISION, PRIVATE :: rho_0,Gamma1,Gamma2,Gamma3
  DOUBLE PRECISION, PRIVATE :: K1,K2,K3,Gamma1_1,Gamma2_1
  DOUBLE PRECISION, PRIVATE :: Gamma3_1,a_c(0:3)

  CONTAINS
    SUBROUTINE gen_pwp_eos(rho,u_cold,P_cold,u,P,cs)

    !********************************************
    !                                           *
    ! general piecewise polytropic EOS;         * 
    !                                           *
    ! ATTENTION: a) select_EOS_parameters must  *
    !               have been called first;     *
    !                                           *
    !            b) expects rest-mass density   *
    !               rho= nlrf*m0c2_cu           *
    !               (code units as input)       *
    !                                           *
    ! SKR 12.05.2021                            *
    !                                           *
    ! Made u and P optional. If they are both   *
    ! given, the thermal part is used. If they  *
    ! both not given, the thermal part is set   *
    ! to 0 and the cold u and P are computed.   *
    ! If only one of them is given, the         *
    ! SUBROUTINE stops printing an error        *
    ! message.                                  *
    !                                           *
    ! FT 13.04.2023                             *
    !                                           *
    !********************************************

    IMPLICIT NONE
    
    DOUBLE PRECISION, INTENT(IN)           :: rho
    DOUBLE PRECISION, INTENT(IN),  OPTIONAL:: u
    DOUBLE PRECISION, INTENT(OUT)          :: u_cold, P_cold, cs
    DOUBLE PRECISION, INTENT(OUT), OPTIONAL:: P
    
    DOUBLE PRECISION Gamma,Gamma_1,u_th,P_th
    DOUBLE PRECISION Kap,enth,a_cold

    IF( PRESENT(u).AND.(.NOT.PRESENT(P)) .OR. &
        PRESENT(P).AND.(.NOT.PRESENT(u)) )THEN

      PRINT *, "** ERROR in SUBROUTINE gen_pwp_eos! The optional arguments ", &
               "u and P should be both present (hot system) or both ", &
               "absent (cold system). Please check that this is the case."
      PRINT *
      STOP

    ENDIF

    ! COLD EOS: select correct parameters (Read et al. 2009, Eqs 1-9)
    IF( rho < rho_0 )THEN
       Gamma= Gamma0
       Kap=    K0
       a_cold= a_c(0)
    ELSEIF( rho_0 <= rho .AND. rho < rho_1 )THEN
       Gamma=  Gamma1
       Kap=    K1
       a_cold= a_c(1)
    ELSEIF( rho_1 <= rho .AND. rho < rho_2 )THEN
       Gamma=  Gamma2
       Kap=    K2
       a_cold= a_c(2)
    ELSEIF( rho >= rho_2 )THEN
       Gamma=  Gamma3
       Kap=    K3
       a_cold= a_c(3)
    ELSE
       PRINT*,'...something is wrong in gen_pwp_eos...'
       STOP
    ENDIF
    Gamma_1= Gamma - 1.0D0
    P_cold=  Kap*rho**Gamma
    u_cold=  a_cold + (Kap/Gamma_1)*rho**Gamma_1

    IF( PRESENT(u) )THEN
      u_th= MAX(u - u_cold, 0.0D0)
      P_th= Gamma_th_1*rho*u_th
      P=    P_cold + P_th
    ELSE
      u_th= 0.D0
      P_th= 0.D0
    ENDIF

    enth= 1.0D0 + (u_cold + u_th) + (P_cold + P_th)/rho
    cs=   SQRT((Gamma*P_cold + Gamma_th*P_th)/(rho*enth))
    
  END SUBROUTINE gen_pwp_eos
  


    SUBROUTINE select_EOS_parameters(eos_str)
    
    !***************************************************
    !                                                  *
    ! select parameters for piecewise polytrope;       *
    ! SKR 15.12.2023                                   *
    !                                                  *
    !***************************************************

    IMPLICIT NONE

    CHARACTER(4), INTENT(IN) :: eos_str

    DOUBLE PRECISION eps_0,eps_1,eps_2

    ! parameters according to specified EOS
    SELECT CASE( TRIM(eos_str) )
    CASE( 'SLy' )
       rho_0=  rho_0_SLy
       Gamma1= Gamma1_SLy
       Gamma2= Gamma2_SLy
       Gamma3= Gamma3_SLy
    CASE( 'ALF2' )
       rho_0=  rho_0_ALF2
       Gamma1= Gamma1_ALF2
       Gamma2= Gamma2_ALF2
       Gamma3= Gamma3_ALF2
    CASE( 'ALF4' )
       rho_0=  rho_0_ALF4
       Gamma1= Gamma1_ALF4
       Gamma2= Gamma2_ALF4
       Gamma3= Gamma3_ALF4
    CASE( 'ENG' )
       rho_0=  rho_0_ENG
       Gamma1= Gamma1_ENG
       Gamma2= Gamma2_ENG
       Gamma3= Gamma3_ENG
    CASE( 'H4' )
       rho_0=  rho_0_H4
       Gamma1= Gamma1_H4
       Gamma2= Gamma2_H4
       Gamma3= Gamma3_H4
    CASE( 'MPA1' )
       rho_0=  rho_0_MPA1
       Gamma1= Gamma1_MPA1
       Gamma2= Gamma2_MPA1
       Gamma3= Gamma3_MPA1
    CASE( 'MS1' )
       rho_0=  rho_0_MS1
       Gamma1= Gamma1_MS1
       Gamma2= Gamma2_MS1
       Gamma3= Gamma3_MS1
    CASE( 'MS1b' )
       rho_0=  rho_0_MS1b
       Gamma1= Gamma1_MS1b
       Gamma2= Gamma2_MS1b
       Gamma3= Gamma3_MS1b
    CASE( 'AP3' )
       rho_0=  rho_0_AP3
       Gamma1= Gamma1_AP3
       Gamma2= Gamma2_AP3
       Gamma3= Gamma3_AP3
    CASE( 'WFF1' )
       rho_0=  rho_0_WFF1
       Gamma1= Gamma1_WFF1
       Gamma2= Gamma2_WFF1
       Gamma3= Gamma3_WFF1
    CASE( 'WFF2' )
       rho_0=  rho_0_WFF2
       Gamma1= Gamma1_WFF2
       Gamma2= Gamma2_WFF2
       Gamma3= Gamma3_WFF2
    CASE( 'GNH3' )
       rho_0=  rho_0_GNH3
       Gamma1= Gamma1_GNH3
       Gamma2= Gamma2_GNH3
       Gamma3= Gamma3_GNH3
    CASE( 'APR4' )
       rho_0=  rho_0_APR4
       Gamma1= Gamma1_APR4
       Gamma2= Gamma2_APR4
       Gamma3= Gamma3_APR4
    CASE( 'haso' )
       rho_0=  rho_0_haso
       Gamma1= Gamma1_haso
       Gamma2= Gamma2_haso
       Gamma3= Gamma3_haso
    CASE( 'soft' )
       rho_0=  rho_0_soft
       Gamma1= Gamma1_soft
       Gamma2= Gamma2_soft
       Gamma3= Gamma3_soft
    CASE( 'stif' )
       rho_0=  rho_0_stiff
       Gamma1= Gamma1_stiff
       Gamma2= Gamma2_stiff
       Gamma3= Gamma3_stiff
    CASE DEFAULT
       PRINT*,'...invalid EOS-string: ',eos_str
       PRINT*,'...(keep in mind: 4 letters/numbers) '
       STOP
    END SELECT

    ! Gamma - 1
    Gamma1_1= Gamma1 - 1.0D0
    Gamma2_1= Gamma2 - 1.0D0
    Gamma3_1= Gamma3 - 1.0D0

    ! prefactors from continuity
    K1= K0*rho_0**(Gamma0 - Gamma1)
    K2= K1*rho_1**(Gamma1 - Gamma2)
    K3= K2*rho_2**(Gamma2 - Gamma3)

    ! transition parameters, Eq.(6-7) Read+ 2009
    a_c(0)= 0.0D0
    eps_0=  rho_0 + K0/Gamma0_1*rho_0**Gamma0

    a_c(1)= eps_0/rho_0 - 1.0D0 - K1/Gamma1_1*rho_0**Gamma1_1
    eps_1=  (1.0D0 + a_c(1))*rho_1 + K1/Gamma1_1*rho_1**Gamma1

    a_c(2)= eps_1/rho_1 - 1.0D0 - K2/Gamma2_1*rho_1**Gamma2_1
    eps_2=  (1.0D0 + a_c(2))*rho_2 + K2/Gamma2_1*rho_2**Gamma2

    a_c(3)= eps_2/rho_2 - 1.0D0 - K3/Gamma3_1*rho_2**Gamma3_1

    ! Separating pressures
    p0= K0*rho_0**Gamma0
    p1= K1*rho_1**Gamma1
    p2= K2*rho_2**Gamma2

  END SUBROUTINE select_EOS_parameters




    SUBROUTINE gen_pwp_eos_all(npart,rho,u)

    !********************************************
    !                                           *
    ! general piecewise polytropic EOS;         * 
    !                                           *
    ! ATTENTION: a) select_EOS_parameters must  *
    !               have been called first;     *
    !                                           *
    !            b) expects rest-mass density   *
    !               rho= nlrf*m0c2_cu           *
    !               (code units as input)       *
    !                                           *
    ! version for ALL particles; Pr & cs are    *
    ! passed into module 'SPH_variables'        *
    !                                           *
    ! SKR 14.05.2021                            *
    !                                           *
    !********************************************

    USE SPH_variables, ONLY: Pr,cs

    IMPLICIT NONE

    INTEGER,          INTENT(IN) :: npart
    DOUBLE PRECISION, INTENT(IN) :: rho(npart)
    DOUBLE PRECISION, INTENT(IN) :: u(npart)

    INTEGER          a
    DOUBLE PRECISION Gamma,Gamma_1,P_cold,u_cold,ua,P
    DOUBLE PRECISION u_th,P_th,Kap,enth,a_cold,rhoa

    DO a= 1,npart
       
       rhoa= rho(a)
       ua=   u(a)
       
       ! COLD EOS: select correct parameters (Read et al. 2009, Eqs 1-9)
       IF( rhoa < rho_0 )THEN
          Gamma= Gamma0
          Kap=    K0
          a_cold= a_c(0)
       ELSEIF( rho_0 <= rhoa .AND. rhoa < rho_1 )THEN
          Gamma=  Gamma1
          Kap=    K1
          a_cold= a_c(1)
       ELSEIF( rho_1 <= rhoa .AND. rhoa < rho_2 )THEN
          Gamma=  Gamma2
          Kap=    K2
          a_cold= a_c(2)
       ELSEIF( rhoa >= rho_2 )THEN
          Gamma=  Gamma3
          Kap=    K3
          a_cold= a_c(3)
       ELSE
          PRINT*,'...something is wrong in gen_pwp_eos_all...'
          STOP
       ENDIF
       Gamma_1= Gamma - 1.0D0
       P_cold=  Kap*rhoa**Gamma
       u_cold=  a_cold + (Kap/Gamma_1)*rhoa**Gamma_1

       ! thermal part
       IF( u_cold > ua )THEN
          PRINT*,'...u_cold= ',u_cold
          PRINT*,'...u=      ',ua
          STOP
       ENDIF
       u_th= MAX(ua - u_cold, 0.0D0)
       P_th= Gamma_th_1*rhoa*u_th

       ! deliverables
       Pr(a)= P_cold + P_th
       enth=  1.0D0 + ua + P/rhoa
       cs(a)= SQRT((Gamma*P_cold + Gamma_th*P_th)/(rhoa*enth))
    ENDDO
    
  END SUBROUTINE gen_pwp_eos_all

END MODULE pwp_eos




