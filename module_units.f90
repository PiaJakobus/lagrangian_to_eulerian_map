MODULE units

  USE constants, ONLY: G_grav,MSun,G_Msun,amu,&
                       c2cgs,c3cgs,fm2cm,MeV2erg,&
                       c_light,k_B,m0c2

  IMPLICIT NONE

  DOUBLE PRECISION umass,udist,udens,utime,uerg,uergg
  DOUBLE PRECISION uergcc,upress_nuc_to_code
  DOUBLE PRECISION ueint_nuc_to_code,c2cu,c3cu
  DOUBLE PRECISION MeVnuc_to_code,umass_in_Msun
  DOUBLE PRECISION uvel,c_light_cu,c_light_cu2
  DOUBLE PRECISION c_light_cu3,c_light_cu4
  DOUBLE PRECISION c_light_cu5,Temp_fac
  DOUBLE PRECISION time_cu_to_s,time_cu_to_ms
  DOUBLE PRECISION time_cu_to_min,time_cu_to_hour
  DOUBLE PRECISION m0c2_cu,ud3


CONTAINS


  SUBROUTINE set_units(unit_system)

    !*********************************************
    !                                            *
    !  once the variable 'unit_system' is known  *
    !  calculate unit transformation factors;    *
    !  SKR 18.07.2012                            *
    !                                            *
    !*********************************************
    
    IMPLICIT NONE
    
    CHARACTER(3), INTENT(IN):: unit_system
    
    !-------------------------------------------!
    !-- set mass & distance, rest follows ... --!
    !-------------------------------------------!
    SELECT CASE(unit_system)
    CASE('NSM') ! for neutron star mergers etc.
       umass_in_Msun= 1.0D0
       umass=         umass_in_Msun*MSun
       udist=         147662.50382504018D0 ! so uvel= c 
    CASE('WDS') ! white dwarf stars
       umass_in_Msun= 1.0D0
       umass=         umass_in_Msun*MSun
       udist=         1.0D9     
    CASE('MBH') ! massive BHs, solar type stars etc
       umass_in_Msun= 1.0D0
       umass=         umass_in_Msun*MSun
       udist=         1.0D10    
    END SELECT
    
    !-----------------!
    !-- other units --!
    !-----------------!
    ! density
    ud3=             udist**3
    udens=           umass/ud3
    
    ! time
    utime=           SQRT(ud3/(G_Msun*umass_in_Msun))
    time_cu_to_s=    utime
    time_cu_to_ms=   utime*1.0D3
    time_cu_to_min=  time_cu_to_s/60.D0
    time_cu_to_hour= time_cu_to_s/3600.D0
    
    ! velocity
    uvel= udist/utime

    ! ergs
    uerg= umass*udist**2/utime**2
    
    ! ergs per gram
    uergg= udist**2/utime**2
    
    ! erg per cm**3
    uergcc= umass/(udist*utime**2)
    
    ! pressure conversion factor
    upress_nuc_to_code= (MeV2erg/fm2cm**3)/uergcc
    
    ! specific internal energy density
    ueint_nuc_to_code= (MeV2erg/amu)/uergg
    
    ! c2cu in reactions per code time and target nucleon 
    c2cu= c2cgs*utime
    
    ! c3cu in MeV per code time and target nucleon 
    c3cu= c3cgs*utime/MeV2erg
    
    ! MeV/nucleon in code units
    MeVnuc_to_code= MeV2erg/(amu*uergg)
    
    ! speed of light in code units
    c_light_cu=  c_light/uvel
    c_light_cu2= c_light_cu*c_light_cu
    c_light_cu3= c_light_cu*c_light_cu2
    c_light_cu4= c_light_cu*c_light_cu3
    c_light_cu5= c_light_cu*c_light_cu4

    ! temperature polytropic gas: T_K= (P/rho)*Temp_fac
    Temp_fac=  uergcc*amu/(k_b*udens)

    ! baryon rest mass energy in code units
    m0c2_cu= m0c2/uerg

  END SUBROUTINE set_units
  

END MODULE units

