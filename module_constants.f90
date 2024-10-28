MODULE constants

  !*************************************************
  !                                                *
  ! various useful numbers and conversion factors; *
  ! fundamental constants according to             *
  ! http://pml.nist.gov/                           *
  ! SKR 01.02.2016                                 *
  !                                                *
  !*************************************************

  IMPLICIT NONE

  ! numbers
  DOUBLE PRECISION, PARAMETER :: tiny            = 1.D-16
  DOUBLE PRECISION, PARAMETER :: huge            = 1.D30
  DOUBLE PRECISION, PARAMETER :: ln2             = 0.69314718055995D0
  DOUBLE PRECISION, PARAMETER :: ln2_1           = 1.D0/ln2
  DOUBLE PRECISION, PARAMETER :: sqrt2           = 1.4142136D0
  DOUBLE PRECISION, PARAMETER :: sqrt3           = 1.7320508D0
  DOUBLE PRECISION, PARAMETER :: pi              = 3.14159265358979D0  
  DOUBLE PRECISION, PARAMETER :: two_pi          = 2.D0*pi
  DOUBLE PRECISION, PARAMETER :: four_pi         = 4.D0*pi
  DOUBLE PRECISION, PARAMETER :: third           = 1.0D0/3.0D0
  DOUBLE PRECISION, PARAMETER :: two_thirds      = 2.0D0/3.0D0
  DOUBLE PRECISION, PARAMETER :: three_quarts    = 0.75D0
  DOUBLE PRECISION, PARAMETER :: sixth           = 1.0D0/6.0D0
  DOUBLE PRECISION, PARAMETER :: half            = 0.5D0
  DOUBLE PRECISION, PARAMETER :: quarter         = 0.25D0
  DOUBLE PRECISION, PARAMETER :: fifteenth       = 1.D0/15.D0

  ! constants
  DOUBLE PRECISION, PARAMETER :: RSun            = 6.955D+10            ! in cm
  DOUBLE PRECISION, PARAMETER :: MSun            = 1.9891D+33           ! in g; recommendation IAU
  DOUBLE PRECISION, PARAMETER :: G_grav          = 6.67384D-8           ! in cm^3/g/s^2
  DOUBLE PRECISION, PARAMETER :: G_Msun          = 1.32712440018D26     ! product Msun*G_grav in cgs; --> http://ssd.jpl.nasa.gov/?constants 
  DOUBLE PRECISION, PARAMETER :: c_light         = 2.99792458D10        ! in cm/s
  DOUBLE PRECISION, PARAMETER :: c_light2        = c_light**2
  DOUBLE PRECISION, PARAMETER :: RSS_Sun_cgs     = 2.D0*G_Msun/c_light2 ! Schwarzschild-radius of the Sun (cgs)
  DOUBLE PRECISION, PARAMETER :: AU_cgs          = 1.4959787069D13      ! Astronomical Unit in cgs
  DOUBLE PRECISION, PARAMETER :: amu             = 1.6605402D-24        ! atomic mass unit
  DOUBLE PRECISION, PARAMETER :: m_el            = 9.1093897D-28        ! electron mass
  DOUBLE PRECISION, PARAMETER :: h_Planck        = 6.6260693D-27        ! Planck's constant
  DOUBLE PRECISION, PARAMETER :: k_B             = 1.3806505D-16        ! Boltzmann's constant
  DOUBLE PRECISION, PARAMETER :: sig_SB          = 5.67051D-5           ! Stefan-Boltzmann constant  
  DOUBLE PRECISION, PARAMETER :: e_coul          = 1.60317653D-19       ! elementary charge
  DOUBLE PRECISION, PARAMETER :: c2cgs           = 6.15D-4              ! for neutrino processes  
  DOUBLE PRECISION, PARAMETER :: c3cgs           = 5.04D-10
  DOUBLE PRECISION, PARAMETER :: m0c2            = amu*c_light2         ! atomic mass unit in energy units
  
  ! conversion factors
  DOUBLE PRECISION, PARAMETER :: fm2cm           = 1.D-13
  DOUBLE PRECISION, PARAMETER :: cm2km           = 1.D-5
  DOUBLE PRECISION, PARAMETER :: cm2m            = 1D-2
  DOUBLE PRECISION, PARAMETER :: km2cm           = 1.D+5
  DOUBLE PRECISION, PARAMETER :: MeV2erg         = 1.60217646D-6
  DOUBLE PRECISION, PARAMETER :: MeV2K           = 11598589518.725805D0
  DOUBLE PRECISION, PARAMETER :: K2MeV           = 1.D0/MeV2K
  DOUBLE PRECISION, PARAMETER :: Press_nuc2cgs   = MeV2erg/fm2cm**3
  DOUBLE PRECISION, PARAMETER :: specen_nuc2cgs  = MeV2erg/amu          ! specific energies from (MeV/bar) to (erg/g)
  DOUBLE PRECISION, PARAMETER :: n2rho           = amu/fm2cm**3
  DOUBLE PRECISION, PARAMETER :: press_SI2cgs    = 1.D1                 ! Pascal (SI) to Barye (CGS)
  DOUBLE PRECISION, PARAMETER :: m2cm            = 1.D2
  DOUBLE PRECISION, PARAMETER :: kg2g            = 1.D3
  DOUBLE PRECISION, PARAMETER :: rad2deg         = 180.D0/pi            ! radians to degrees
  DOUBLE PRECISION, PARAMETER :: deg2rad         = pi/180.D0            ! degrees to radians

  ! constants for RCB-tree
  DOUBLE PRECISION, PARAMETER :: safety_factor   = 1.05D0
  DOUBLE PRECISION, PARAMETER :: safety_factor_1 = 1.D0/safety_factor

  
END MODULE constants
