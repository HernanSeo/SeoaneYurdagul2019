module common_par
use nrtype
implicit none

!'-------------------------------------------------------------'
! -------------------model parameters --------------------------
!'-------------------------------------------------------------'
! betay=0.86
real(dp),parameter :: A=1d0!,mu=.205,  rbar= 0.04 ,tempr=1.+rbar , betay= 0.96 ,sig=2.   , theta = 0.33d0, pssig=0d0
!real(dp),parameter :: kappan=0.3235 ,  kappat=kappan, omega= 0.307000693252802, delta=1d0, tau=0d0, thresh=0.969d0 !1.5d0!!0.89 in the min yT for out calibration

! --------------grid  ----------------------------------------
!'-------------------------------------------------------------'

!real(dp) ,parameter :: minb = -0.95d0!-0.95d0!-0.7!-0.85
!real(dp) ,parameter :: largestb=1.75d0!-0.1d0!-0.6
!real(dp) ,parameter :: mind= -1.5d0!-0.5d0!-0.25 !-0.15!-0.1*0.89072920798418365/(1.04-1.0)
!real(dp) ,parameter :: largestd = 0.5d0!this is not needed anymore
!integer  ,parameter :: nb1 = 25
!integer  ,parameter :: nd1 = 25
!real(dp) ,parameter :: gsp = 1! (I have to keep this because of the 2-d interpolation that requires equally spaced grids... original: 1.1

!'-------------------------------------------------------------'
! --------------numerical  ------------------------------
!'-------------------------------------------------------------'

integer ,parameter :: maxitpolicy=300
real(dp),parameter :: tolpol=1.0e-3!
real(dp),parameter :: tolvfi=1.0e-3!1
integer ,parameter :: ntry=2000
integer ,parameter :: maxitbis=200
real(dp),parameter :: xacc=1.e-10
real(dp),parameter :: infinit= HUGE(0.)
integer ,parameter :: maxit=1000
real(dp),parameter :: lambda=.1!1.
real(dp),parameter :: tolzbrent=1.0e-8
real(dp),parameter :: tolbrent=1e-8
!'-------------------------------------------------------------'
! --------------stochastic process ------------------------------
!'-------------------------------------------------------------'

integer, parameter  ::  nip=4 ! no. of T shocks  and nt shocks
integer, parameter  ::  nipp=nip ! no. of Nt shocks       (must equal nvaln defined below)
integer , parameter  :: ny1=nip*nipp

!'-------------------------------------------------------------'
! --------------others ------------------------------
!'-------------------------------------------------------------'

integer, parameter  ::  tcut=1000           ! observation deleted
!integer, parameter  ::  tt=2000      ! number of simulations
real(dp),parameter  :: tollimit =   0.00001
! global
real(dp) :: senn(ny1), se(ny1),p(ny1,ny1)!,ztvect(ny1), znvect(ny1)
!real(dp) :: ktight(nb1,ny1)


end module common_par


module global

 ! supporting module for tauchen hussey

  implicit none

      integer zznvar,zznlag,zzns,zzthet,zzpmat,zziw,zznst,zziw2,  zziw3,zzw  ,nvaln ,nlag,nrquad ,nvar
      parameter(nvaln=4)         ! NUMBER OF DISCRETE POINT FOR EACH VARIABLE
      parameter(nvar=2)
      parameter(nlag=1)
      parameter(nrquad=210)
      parameter(zznvar=2)
      parameter(zznlag=1)
      parameter(zzns=nvaln*nvaln)
      parameter(zzthet=zznvar+(zznvar**2)*(zznlag+2)+1)
      parameter(zzpmat=(zzns**zznlag)*zzns)
      parameter(zziw=2*zzns+2*zznlag)
      parameter(zznst=zzns**zznlag)
      parameter(zzw=4*zzns+2*zznvar+4*zzns*zznvar+2*zznvar*zznlag)
      parameter(zziw2=zznst*zznlag)
      parameter(zziw3=zznst*zzns)

end module global


