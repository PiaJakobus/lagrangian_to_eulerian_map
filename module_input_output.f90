MODULE input_output
  
  !*****************************************************************
  !                                                                *
  !  collect here all subroutines that deal with input and output  *
  !  SKR 08.08.2014                                                *
  !                                                                *
  !*****************************************************************

  USE units,     ONLY: utime,umass,udens,set_units,m0c2_cu
  USE constants, ONLY: amu,Msun
  
  IMPLICIT NONE
  
  INTEGER          :: i,dcount
  INTEGER          :: unit_logfile=  1313 
  INTEGER          :: unit_consfile= 1414
  INTEGER          :: unit_GWfile=   1515
  DOUBLE PRECISION :: t_start,t_end
  
  LOGICAL          :: file_exists
  CHARACTER(17)    :: BH_file='BH_parameters.dat'
  CHARACTER(18)    :: GWfile='GWs_quadrupole.dat'
  CHARACTER(19)    :: logfile='SPHINCS_logfile.txt'
  CHARACTER(20)    :: input_file='SPHINCS_fm_input.dat'
  CHARACTER(24)    :: consfile='SPHINCS_conservation.dat'
 
  
CONTAINS

 
  SUBROUTINE read_SPHINCS_dump(namefile)
    
    !*******************************************
    !                                          *
    ! read SPHINCS data dump; SKR 27.07.2017   *
    ! SKR 27.10.2017                           *
    !                                          *
    !*******************************************
    USE analyze,       ONLY: nstep_tot,COM,Mbar_last
    USE files,         ONLY: filename
    USE sph_variables, ONLY: npart,rstar,mstar,n1,n2,npm,t,&
                             h,escap,tkin,tgrav,tterm,pos_u,&
                             vel_u,u,pmass,rho,temp,Ye,av,divv,&
                             Theta,Pr, allocate_SPH_memory, deallocate_SPH_memory
    USE pwp_EOS,       ONLY: gen_pwp_eos,gen_pwp_eos_all,select_eos_parameters
    USE units
!++++++++++++++++++++++=debugging+++++++++++++++++++++
    USE options,       ONLY: ipos,icow
!+++++++++++++++++++++++++++++++++++++++++++++++++++++

    IMPLICIT NONE

    INTEGER          a
    !DOUBLE PRECISION com_x,com_y,com_z,com_d
    !DOUBLE PRECISION com_vx,com_vy,com_vz,com_v

!++++++++++++++++++debugging++++++++++++++++
    DOUBLE PRECISION adot_0,m1,m2,vr1,vr2!,RCM(3),VCM(3)
!+++++++++++++++++++++++++++++++++++++++++++

    CHARACTER(LEN= *), INTENT(INOUT), OPTIONAL :: namefile
    CHARACTER(LEN= :), ALLOCATABLE             :: finalnamefile

    IF( PRESENT(namefile) )THEN
      finalnamefile= namefile
    ELSE
      finalnamefile= filename
    ENDIF

    ! open file
    OPEN(10,file=finalnamefile,form='UNFORMATTED')
    PRINT*
    WRITE(*,'(a,a)')' ...reading from file ',finalnamefile
    READ(10) npart
    CLOSE(10)
    OPEN(10,file=finalnamefile,form='UNFORMATTED')
    CALL allocate_SPH_memory
    ! read in MAGMA-type format
    READ(10)npart,rstar,mstar,n1,n2,npm,t, & 
            (h(a),a=1,npart),escap,tkin,tgrav,tterm, & 
            (pos_u(1,a),a=1,npart),(pos_u(2,a),a=1,npart),&
            (pos_u(3,a),a=1,npart),(vel_u(1,a),a=1,npart), & 
            (vel_u(2,a),a=1,npart),(vel_u(3,a),a=1,npart),&
            (u(a),a=1,npart),(pmass(a),a=1,npart),           &
            (rho(a),a=1,npart),(temp(a),a=1,npart),      &
            (Ye(a),a=1,npart),(av(a),a=1,npart),          &
            (divv(a),a=1,npart),(Theta(a),a=1,npart),     &
            (Pr(a),a=1,npart)                          !, &
            !
            !-- leave here for potential later use
            ! 
            !(pmasspm(a),a=1,npm),&
            !(pmpos(1,a),a=1,npm),(pmpos(2,a),a=1,npm),&
            !(pmpos(3,a),a=1,npm),(pmvel(1,a),a=1,npm),&
            !(pmvel(2,a),a=1,npm),(pmvel(3,a),a=1,npm),&
            !(pmdvdt(1,a),a=1,npm),(pmdvdt(2,a),a=1,npm),&
            !(pmdvdt(3,a),a=1,npm)
    CLOSE(10)
    PRINT*,'...done...'
    print*, "test", size(h), npart
    !stop 
    ! SKR 15.03.2024:
    ! since the current initial data use our old convention of 
    ! measuring energies in units of m0c^2 and this pressure IS USED
    ! in the first phys2cons, we NEED TO TRANSFORM THE PRESSURE HERE
    Pr= Pr*m0c2_cu

    ! print some numbers for controle
    PRINT*
    PRINT*,'...npart         ',npart
    IF( npart > 0 )THEN
       PRINT*,'...log(rho_max)  ',LOG10(MAXVAL(rho(1:npart)))
       WRITE(*,'(a,E14.3,a)')' ...(rho_max=   ',MAXVAL(rho(1:npart)),')'
       PRINT*,'...log(rho_min)  ',LOG10(MINVAL(rho(1:npart)))
       PRINT*
    ELSE
       icow= -1
    ENDIF
 
    ! monitor total baryonic mass
    Mbar_last= SUM(pmass(1:npart))

    PRINT*,'...Mbar_last = ', Mbar_last

    ! some double-checks
    IF( n1 + n2 .NE. npart )THEN
       PRINT*
       PRINT*,'...n1+n2 do not add up to npart '
       PRINT*,'...n1=    ',n1
       PRINT*,'...n2=    ',n2
       PRINT*,'...n1+n2= ',n1+n2
       PRINT*,'...npart= ',npart
       STOP
    ENDIF

    IF( npart > 0 )then
!+++++++++++++++++++++++++++++debugging++++++++++++++++++++++++
       adot_0= 0.D0 !5.D-2
       IF( adot_0 > 0.D0 )THEN
          PRINT*
          PRINT*,'...WARNING: adding radial velocity adot=',adot_0
          PRINT*
       ENDIF

       start:IF( ipos == 0 )THEN
          m1= SUM(pmass(1:n1))
          m2= SUM(pmass(n1+1:npart))       

          ! add additional radial velocity
          vr1= m2/(m1+m2)*adot_0
          DO a=1,n1
             vel_u(1,a)= vel_u(1,a) + vr1 
          ENDDO
       
          vr2= -m1/(m1+m2)*adot_0
          DO a=n1+1,npart
             vel_u(1,a)= vel_u(1,a) + vr2
          ENDDO
       
          ! reset CM-velocity
!          RCM= 0.D0
!          VCM= 0.D0
!          sum= 0.D0
!          DO a=1,npart
!             RCM= RCM + pmass(a)*pos_u(1:3,a)
!             VCM= VCM + pmass(a)*vel_u(1:3,a)
!             sum= sum + pmass(a)
!          ENDDO
!          RCM= RCM/sum
!          VCM= VCM/sum

!          DO a=1,npart
!             pos_u(1:3,a)= pos_u(1:3,a) - RCM 
!             vel_u(1:3,a)= vel_u(1:3,a) - VCM
!          ENDDO

       ENDIF start
!++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++


       ! baryonic centre of mass
!       CALL COM(npart,pos_u,pmass,com_x,com_y,com_z,com_d)
!       PRINT*
!       PRINT*,'...the baryonic centre of mass is offset' 
!       PRINT*,'   from the origin by',com_d
!       PRINT*
    
       ! baryonic centre of mass velocity 
!       CALL COM(npart,vel_u,pmass,com_vx,com_vy,com_vz,com_v)
!       PRINT*
!       PRINT*,'...the baryonic centre of mass has velocity'
!       PRINT*,'   (code units)      ',com_v
!       PRINT*

    ENDIF

    ! COMMENT SKR 30.10.2017:
    ! metric and Nstar DO NOT NEED to be updated here;
    ! this is done before they are needed in step_RK3
    
    ! initialize total time step counter
    nstep_tot= 0

  END SUBROUTINE read_SPHINCS_dump

END MODULE input_output
 
