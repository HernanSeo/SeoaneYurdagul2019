program solve
    use common_par
implicit none
integer, parameter :: n_b=200,n_z=16,n_zt=4,n_zn=4,n_g=1,n_kappa=1,n_omega=1,n_iter=30000
COMMON /PARAM/ bvect,kappavect,omegavect,r,cmin,bata,kappat,kappan,omega,sigma,eta,&
fixedt,fixedn,coefq,n_iter3,Planner,ztvect,znvect,gvect,pz,pg,pkappa,pomega,allow,gradb
DOUBLE PRECISION bvect(n_b),kappavect(n_kappa),omegavect(n_omega),r,cmin,bata,kappat,kappan,omega&
,sigma,eta,fixedt,fixedn,coefq,pz(n_z,n_z),pg(n_g,n_g),pkappa(n_kappa,n_kappa),pomega(n_omega,n_omega),&
ztvect(n_z),znvect(n_z),gvect(n_g),allow,gradb
INTEGER n_iter3,Planner
SAVE /PARAM/
COMMON /MATS/ b,ct,cn,pn,expmat
DOUBLE PRECISION, dimension(n_b,n_z,n_g,n_kappa,n_omega):: b,ct,cn,pn,expmat
SAVE /MATS/
COMMON /MATS2/ mu,ctlimited,ctunlimited
DOUBLE PRECISION, dimension(n_b,n_z,n_g,n_kappa,n_omega):: mu,ctlimited,ctunlimited
SAVE /MATS2/



double precision, dimension(n_b,n_z,n_g,n_kappa,n_omega)::bnew,pnnew,bguess,blimited,bunlimited,errunlimited,expnext
double precision guesses,choices,shadows,sum0,sum1,sum2,difsol,&
errmat(n_b,n_z,n_g,n_kappa,n_omega),errors,nothing(4),grad,difbrelbest,difprelbest,difsolbest,&
consn,const,uT,difrelbestold,penalty,densdif,difcrit,unlimited(3),ctpr,cnpr,pnpr,mupr,psipr!,eq1,eq2,eq3
!!!! Functions !!!!
double precision intp2,zfun,yfun,Rfun,Rlim,kappatfun,kappanfun,normalcdf!,expfun

!double precision updown,delta_z,z_inf,z_sup,z_right,z_left
double precision bmin,bmax
integer i,lt,ltnext,iter,maxpos(4),inext,load,&
cminhappens(n_b,n_z,n_g,n_kappa,n_omega),ii,iter3soln(n_b,n_z,n_g,n_kappa,n_omega)

double precision phi(n_b,n_z,n_g,n_kappa,n_omega),phinew(n_b,n_z,n_g,n_kappa,n_omega),contr,&
ratio_b,bopt,avgb,avgboy,blim,mupos,avgzt,avgzn,&
avgdifbrel,difbrel(n_b,n_z,n_g,n_kappa,n_omega),difprel(n_b,n_z,n_g,n_kappa,n_omega),limitb


integer index_b,index_kn,i_b,ll,iterdens,h,hnext,x,xnext,tr,trnext,j,n

double precision updown,rhog,mean_g,std_g,g_inf,g_sup,delta_g,g_left,g_right,std_eps,&
pana(n_g,n_g,n_g),pantempa(n_g,n_g,n_g),ptempa1(n_g,n_g,n_g),ptempa2(n_g,n_g,n_g),ptempa3(n_g,n_g,n_g)&
,ptempa4(n_g,n_g,n_g),ppa,qqa,gbar,g0
double precision rhozt,rhozn,mean_z,std_z,z_inf,z_sup,delta_z,z_left,z_right,std_epst,std_epsn,&
ztvecttemp(n_zt),znvecttemp(n_zn),pzt(n_zt,n_zt),pzn(n_zn,n_zn)


     integer t1(8),t2(8)
     real :: f(8) = (/ 0., 0., 0., 0., 3600., 60., 1., .001 /)
     double precision tt,sec,sec0,sec1,sec2

character(len=20) :: myfmt,myfmt_k



!Parameters!
r=0.04d0;bata=0.9d0!! If yearly
!r=0.01d0;bata=0.97d0;delta=0.014d0

kappat=0.32d0;kappan=0.32d0
omega=0.31d0;sigma=2.d0;cmin=0.0001d0;fixedt=0.d0;fixedn=0.d0
difcrit=0.0000001d0;coefq=1.d0


allow=1E-10


!omega=0.27
eta=1/0.83d0-1.d0;
!eta=1.d0
!omega=0.9d0

!!! When solving planner's problem, just set this to 1. !!!!!
Planner=0

      open(unit=300,file="outputs/params_end.txt")
      read (300,*) load
      read (300,*) n_iter3
      read (300,*) r
      read (300,*) bata
      read (300,*) bmin
      read (300,*) bmax
      read (300,*) kappat
      read (300,*) gradb
      read (300,*) Planner
      close(300)
!n_iter3=10
      kappan=kappat
      
!if (n_g==1) then
!gvect=1.d0
!pg=1.d0
!else
     
!pg(1,1)=0.8d0
!pg(2,1)=1-pg(1,1)
!pg(2,2)=0.6d0
!pg(1,2)=1-pg(2,2) 
!
!
!gvect(1)=1.d0
!!gvect(2)=1.15d0
!
!gvect(2)=(2-pg(1,1)-pg(2,2))/(1-pg(1,1))*(1.03-(1-(1-pg(1,1))/(2-pg(1,1)-pg(2,2)))*gvect(1))
!
!
!
!    
!endif    



!!!! RW method !!!!
!      std_g=0.05d0
!      g0=1.10d0
!      ppa=0.9;qqa=0.8
!      rhog=ppa+qqa-1        
!      gbar=(((n_g-1.d0)/(1.d0-rhog**2.d0))**0.5d0)*std_g ! log(0.49d0)!
!     
!gvect(1)=log(g0)-gbar
!do i=2,n_g
!gvect(i)=gvect(i-1)+gbar/((n_g-1.d0)/2.d0)
!enddo
!
!!!! If no types !!!!!!!!!!!!!
!!
!      pana(1,1,1)=1.d0
!
!      do n=2,n_g
!
!      ptempa1=0.d0
!      ptempa2=0.d0
!      ptempa3=0.d0
!      ptempa4=0.d0
!
!      do i=1,n-1                                                        !The transition matrix is generated according to the rule highlighted in KAREN A. KOPECKY notes.
!      do j=1,n-1
!      ptempa1(i,j,n)=pana(i,j,n-1)
!      ptempa2(i,j+1,n)=pana(i,j,n-1)
!      ptempa3(i+1,j,n)=pana(i,j,n-1)
!      ptempa4(i+1,j+1,n)=pana(i,j,n-1)
!      enddo
!      enddo
!      pantempa(:,:,n)=ppa*ptempa1(:,:,n)+(1-ppa)*ptempa2(:,:,n)+(1-qqa)*ptempa3(:,:,n)+qqa*ptempa4(:,:,n)
!
!      pana(1,:,n)=pantempa(1,:,n)
!      pana(n,:,n)=pantempa(n,:,n)
!      do i=2,n-1
!      pana(i,:,n)=pantempa(i,:,n)/2
!      enddo                                                             !Normalizing to make each row sum up to 1
!
!      enddo
!      pg=transpose(pana(:,:,n_g))
!
!gvect=exp(gvect)
!
!!write(*,*) pg(1,:)
!!write(*,*) pg(2,:)
!!write(*,*) gvect
!!pause  
!endif

!!!!!!!!!!!!!!!!!!!!!!
!!! Tauchen method for the g's!!!
updown=3.d0
!! AG !!
!mean_g=log(1.006)
!std_eps=0.03d0*2.d0
!rhog=0.17

!std_eps=0.07d0
!rhog=0.04
!mean_g=log(1.0236d0)-0.5d0*std_eps**2/(1.d0-rhog**2)
!mean_g=log(1.03d0)

!Own calib!
!std_eps=0.04d0
!rhog=0.6d0

!Bayesian est!
!std_eps=0.0282d0
!rhog=0.8512d0

!std_eps=0.0414d0
!rhog=0.2592d0

!std_eps=0.0289d0
!rhog=0.6687d0

!Calibration pc!
mean_g=log(1.01d0)

std_eps=0.0353d0
rhog=0.5499d0

!Alternative persistence!
!std_eps=std_eps*sqrt(1-0.7d0**2.d0)/sqrt(1-rhog**2.d0)
!rhog=0.7d0


!Alternative sgg!

!std_eps=std_eps!*2.d0

!mean_g=0.d0

std_g = std_eps/ SQRT(1 - rhog**2)
g_inf = mean_g - updown*std_g
g_sup = mean_g + updown*std_g
delta_g = updown*std_g / (n_g - 1d+0)


do i=1,n_g
   gvect(i) = exp(g_inf + (g_sup - g_inf) * (i-1) / (n_g - 1))
end do

!write(*,*) log(gvect(1)) + delta_g - rhog*log(gvect(1)),log(gvect(n_g))- delta_g - rhog*log(gvect(n_g))
!write(*,*) log(gvect(1)),log(gvect(n_g)),gvect(1),gvect(n_g)
!pause

do i=1,n_g
   g_left  = (log(gvect(1)) + delta_g - rhog*log(gvect(i))-(1-rhog)*mean_g)/ std_eps
   g_right = (log(gvect(n_g))- delta_g - rhog*log(gvect(i))-(1-rhog)*mean_g)/ std_eps
   pg(1,i) = normalcdf(g_left)
   pg(n_g,i) = 1d+0 - normalcdf(g_right)
   do j=2,n_g-1
     g_right = (log(gvect(j)) + delta_g - rhog*log(gvect(i))-(1-rhog)*mean_g)/ std_eps
     g_left  = (log(gvect(j)) - delta_g - rhog*log(gvect(i))-(1-rhog)*mean_g)/ std_eps
     pg(j,i) = (normalcdf(g_right) - normalcdf(g_left) ) !/prob_man_y
   end do
end do



if (n_g==1) then
gvect=1.d0
pg=1.d0
endif






 
if (n_omega==1) then
omegavect=omega
pomega=1.d0
else
     
omegavect(1)=omega*1.10d0
omegavect(2)=omega*.90d0

pomega(1,1)=0.95d0
pomega(2,1)=1-pomega(1,1)
pomega(2,2)=0.95d0
pomega(1,2)=1-pomega(2,2)     
endif

if (n_kappa==1) then
kappavect=kappat
pkappa=1.d0
else
     
kappavect(1)=kappat*1.1d0
kappavect(2)=kappat*0.9d0

pkappa(1,1)=0.95d0
pkappa(2,1)=1-pkappa(1,1)
pkappa(2,2)=0.95d0
pkappa(1,2)=1-pkappa(2,2)     
endif

      open(unit=260,file="outputs/gvect.txt")
      do i=1,n_g
      write (260,*) gvect(i)
      enddo
      close(260)

      open(unit=240,file="outputs/pg.txt")
      do i=1,n_g
      write (240,*) (pg(i,h),h=1,n_g)
      enddo
      close(240)   
 
 
 
      open(unit=262,file="outputs/omegavect.txt")
      do i=1,n_omega
      write (262,*) omegavect(i)
      enddo
      close(262)

      open(unit=242,file="outputs/pomega.txt")
      do i=1,n_omega
      write (242,*) (pomega(i,h),h=1,n_omega)
      enddo
      close(242)   

      open(unit=262,file="outputs/kappavect.txt")
      do i=1,n_kappa
      write (262,*) kappavect(i)
      enddo
      close(262)

      open(unit=242,file="outputs/pkappa.txt")
      do i=1,n_kappa
      write (242,*) (pkappa(i,h),h=1,n_kappa)
      enddo
      close(242)    



!Indicator of using the previous policy functions as guesses here!
!load

if (n_b>1) then
!bvect(1)=-2.2d0;bvect(n_b)=0.d0
bvect(1)=bmin;bvect(n_b)=bmax

do i=2,n_b-1
bvect(i)=bvect(i-1)+(bvect(n_b)-bvect(1))/(n_b-1)
enddo
else
bvect=0.d0
endif




if (n_z>1) then
call process(n_z)      ! stochastic process
 znvect=senn
 ztvect=se
 pz=p
else 
!!! No shock !!!!
znvect=1.d0
ztvect=1.d0
pz=1.d0
endif
!
!!!!!!!!!!!!!!! zt tauchen !!!!!!!!!!!!!!!!!
!!rhozt=0.d0;std_epst=0.0966d0
!updown=3.d0
!mean_z=0.d0
!rhozt=0.7501d0
!std_epst=0.0532d0
!
!std_z = std_epst/ SQRT(1 - rhozt**2)
!z_inf = mean_z - updown*std_z
!z_sup = mean_z + updown*std_z
!delta_z = updown*std_z / (n_zt - 1d+0)
!
!
!do i=1,n_zt
!   ztvecttemp(i) = exp(z_inf + (z_sup - z_inf) * (i-1) / (n_zt - 1))
!end do
!
!
!
!do i=1,n_zt
!   z_left  = (log(ztvecttemp(1)) + delta_z - rhozt*log(ztvecttemp(i))-(1-rhozt)*mean_z)/ std_epst
!   z_right = (log(ztvecttemp(n_zt))- delta_z - rhozt*log(ztvecttemp(i))-(1-rhozt)*mean_z)/ std_epst
!   pzt(1,i) = normalcdf(z_left)
!   pzt(n_zt,i) = 1d+0 - normalcdf(z_right)
!   do j=2,n_zt-1
!     z_right = (log(ztvecttemp(j)) + delta_z - rhozt*log(ztvecttemp(i))-(1-rhozt)*mean_z)/ std_epst
!     z_left  = (log(ztvecttemp(j)) - delta_z - rhozt*log(ztvecttemp(i))-(1-rhozt)*mean_z)/ std_epst
!     pzt(j,i) = (normalcdf(z_right) - normalcdf(z_left) ) !/prob_man_y
!   end do
!end do
!
!if (n_z==1) then
!ztvecttemp=1.d0
!pzt=1.d0
!endif
!
!!!!!!!!!!!!!!! zn tauchen !!!!!!!!!!!!!!!!!
!
!!rhozn=0.d0;std_epsn=0.05520d0
!
!mean_z=0.d0
!rhozn=0.7963d0
!std_epsn=0.0495d0
!
!
!std_z = std_epsn/ SQRT(1 - rhozn**2)
!z_inf = mean_z - updown*std_z
!z_sup = mean_z + updown*std_z
!delta_z = updown*std_z / (n_zn - 1d+0)
!
!
!do i=1,n_zn
!   znvecttemp(i) = exp(z_inf + (z_sup - z_inf) * (i-1) / (n_zn - 1))
!end do
!
!do i=1,n_zn
!   z_left  = (log(znvecttemp(1)) + delta_z - rhozn*log(znvecttemp(i))-(1-rhozn)*mean_z)/ std_epsn
!   z_right = (log(znvecttemp(n_zn))- delta_z - rhozn*log(znvecttemp(i))-(1-rhozn)*mean_z)/ std_epsn
!   pzn(1,i) = normalcdf(z_left)
!   pzn(n_zn,i) = 1d+0 - normalcdf(z_right)
!   do j=2,n_zn-1
!     z_right = (log(znvecttemp(j)) + delta_z - rhozn*log(znvecttemp(i))-(1-rhozn)*mean_z)/ std_epsn
!     z_left  = (log(znvecttemp(j)) - delta_z - rhozn*log(znvecttemp(i))-(1-rhozn)*mean_z)/ std_epsn
!     pzn(j,i) = (normalcdf(z_right) - normalcdf(z_left) ) !/prob_man_y
!   end do
!end do
!
!if (n_z==1) then
!znvecttemp=1.d0
!pzn=1.d0
!endif
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
!
!do i=1,n_z
!ztvect(i)=ztvecttemp(i-floor((i-1)*1./n_zt)*n_zt)
!znvect(i)=znvecttemp(floor((i-1)*1./n_zn)+1)
!enddo
!
!do i=1,n_z
!do j=1,n_z
!pz(j,i)=pzt(j-floor((j-1)*1./n_zt)*n_zt,i-floor((i-1)*1./n_zt)*n_zt)*pzn(floor((j-1)*1./n_zn)+1,floor((i-1)*1./n_zn)+1)
!enddo
!enddo
!!
!!!write(*,*) ztvecttemp
!!!write(*,*) znvecttemp
!!!write(*,*) ztvect
!!!write(*,*) znvect
!!!pause
!!
!!!write(*,*) sum(pzt),sum(pzn),sum(pz(:,2)),sum(pz(:,4))
!!
!!


      open(unit=262,file="outputs/bvect.txt")
      do i=1,n_b
      write (262,*) bvect(i)
      enddo
      close(262)



      open(unit=262,file="outputs/ztvect.txt")
      do i=1,n_z
      write (262,*) ztvect(i)
      enddo
      close(262)
      
      open(unit=262,file="outputs/znvect.txt")
      do i=1,n_z
      write (262,*) znvect(i)
      enddo
      close(262)

      open(unit=242,file="outputs/pz.txt")
      do i=1,n_z
      write (242,*) (pz(i,j),j=1,n_z)
      enddo
      close(242)

   !   open(unit=262,file="outputs/ztvect.txt")
   !   do i=1,n_zt
   !   write (262,*) ztvect(i)
   !   enddo
   !   close(262)

   !   open(unit=242,file="outputs/pzt.txt")
   !   do i=1,n_zt
   !   write (242,*) (pzt(i,j),j=1,n_zt)
   !   enddo
   !   close(242)

   !   open(unit=262,file="outputs/znvect.txt")
   !   do i=1,n_zn
   !   write (262,*) znvect(i)
   !   enddo
   !   close(262)

   !   open(unit=242,file="outputs/pzn.txt")
   !   do i=1,n_zn
   !   write (242,*) (pzn(i,j),j=1,n_zn)
   !   enddo
   !   close(242)

!!!!! Simple guesses for the policy functions!!!!

b=0.d0
!expectb=0.d0

 !$OMP PARALLEL private(i,lt,tr,x,h)
 !$OMP DO SCHEDULE(RUNTIME)

do i=1,n_b
do lt=1,n_z
do tr=1,n_g
do x=1,n_kappa
do h=1,n_omega

ct(i,lt,tr,x,h)=max(cmin,fixedt+yfun(lt,tr,1)+bvect(i)*Rfun(i,lt,tr,x,h)-b(i,lt,tr,x,h)*gvect(tr))
cn(i,lt,tr,x,h)=max(cmin,fixedn+yfun(lt,tr,2))
!lamda(i,lt,tr,x,h)=(omega*ct(i,lt,tr,x,h)**(-eta)+&
!                   (1-omega)*cn(i,lt,tr,x,h)**(-eta))**(-1/eta-1)*omega*ct(i,lt,tr,x,h)**(-eta-1)&
!                 *((omega*ct(i,lt,tr,x,h)**(-eta)+(1-omega)*cn(i,lt,tr,x,h)**(-eta))**(-1/eta))**(-sigma)
enddo
enddo
enddo
enddo
enddo

 !$OMP END DO
 !$OMP END PARALLEL
bguess=b

 !$OMP PARALLEL private(i,lt,tr,x,h)
 !$OMP DO SCHEDULE(RUNTIME)

do i=1,n_b
do lt=1,n_z
do tr=1,n_g
do x=1,n_kappa
do h=1,n_omega
pnnew(i,lt,tr,x,h)=(1-omega)/omegavect(h)*(ct(i,lt,tr,x,h)/cn(i,lt,tr,x,h))**(eta+1)
if (pnnew(i,lt,tr,x,h)<0.001d0.and.bvect(i)>-0.3d0) then
write(*,*) 'pnlow',ct(i,lt,tr,x,h),bvect(i)
pause
endif
enddo
enddo
enddo
enddo
enddo

 !$OMP END DO
 !$OMP END PARALLEL

pn=min(1.d0,pnnew)


mu=0.d0

 !!!! Simple guess for density !!!!!

phi=1.d0/n_b/n_z/n_g/n_kappa/n_omega

 !!! Guesses over!!!

 !!! If load==1, use the previously-exported text file as guesses!!!
write(*,*) load

if (load==1) then
  open (unit=24,file ="outputs/b.txt")
  print*,'model baseline'

do i=1,n_b
do lt=1,n_z
do tr=1,n_g
do x=1,n_kappa
do h=1,n_omega
!read(204,*) nothing(1),nothing(2),nothing(3),b(i,lt,tr,x,h)&
!,mu(i,lt,tr,x,h),ct(i,lt,tr,x,h),cn(i,lt,tr,x,h)&
!,pn(i,lt,tr,x,h),errmat(i,lt,tr,x,h),ii


!read(24,*) nothing(1),nothing(2),nothing(3),nothing(4),b(i,lt,tr,x,h),nothing(1)&
!,mu(i,lt,tr,x,h),nothing(1),nothing(2),ct(i,lt,tr,x,h),cn(i,lt,tr,x,h)&
!,pn(i,lt,tr,x,h),nothing(1),nothing(1),errmat(i,lt,tr,x,h),nothing(1),ii

read(24,*) nothing(1),nothing(2),nothing(3),nothing(4),b(i,lt,tr,x,h)&
,mu(i,lt,tr,x,h),ct(i,lt,tr,x,h),cn(i,lt,tr,x,h)&
,pn(i,lt,tr,x,h),errmat(i,lt,tr,x,h),ii,nothing(1),nothing(2),nothing(1),nothing(2)&
,ctlimited(i,lt,tr,x,h),ctunlimited(i,lt,tr,x,h),expmat(i,lt,tr,x,h)



enddo
enddo
enddo
enddo
enddo
close(24)


open (unit=15,file ="outputs/guesses.txt")
do i=1,n_b
do lt=1,n_z
do tr=1,n_g
do x=1,n_kappa
do h=1,n_omega
read(15,*) bguess(i,lt,tr,x,h)
enddo
enddo
enddo
enddo
enddo
close(15)

endif


!write(*,*) 'after reading',sum(bguess),sum(b),sum(pn)




!!! Will start the main iteration here. difsol: simple max abs val distance in policies (b),difbrelbest is relative measures of distance. !!!
iter=1;difsol=1000.d0;difbrelbest=1000.d0;difrelbestold=1100.d0;difbrel=1000.d0
do while (iter<n_iter.and.maxval(difbrel)>difcrit)
!write(*,*) 'iter', iter
call date_and_time(VALUES=t1)
sec0=f(5)*t1(5)+f(6)*t1(6)+f(7)*t1(7)+f(8)*t1(8)

cminhappens=0
 !$OMP PARALLEL private(i,lt,tr,x,h,guesses,choices,shadows,errors,const,consn,uT,unlimited)
 !$OMP DO ORDERED SCHEDULE(RUNTIME)
do i=1,n_b
do lt=1,n_z
do tr=1,n_g
do x=1,n_kappa
do h=1,n_omega

!!! Use previous policies for b' as guesses in each solution. !!!

!$OMP ORDERED


blimited(i,lt,tr,x,h)=-(kappatfun(lt,x)*coefq*yfun(lt,tr,1)+kappanfun(lt,x)*coefq*pn(i,lt,tr,x,h)*&
yfun(lt,tr,2))/(Rlim(i,lt,tr,x,h)*gvect(tr))


if (yfun(lt,tr,1)+bvect(i)*Rfun(i,lt,tr,x,h)-blimited(i,lt,tr,x,h)*gvect(tr)<cmin) then

bnew(i,lt,tr,x,h)=blimited(i,lt,tr,x,h) !-(kappatfun(lt,x)*coefq*yfun(lt,tr,1)+kappanfun(lt,x)*coefq*pn(i,lt,tr,x,h)*&
 !yfun(lt,tr,2))/(Rlim(i,lt,tr,x,h)*gvect(tr))


errmat(i,lt,tr,x,h)=0.d0
cminhappens(i,lt,tr,x,h)=1
mu(i,lt,tr,x,h)=0.d0

bunlimited(i,lt,tr,x,h)=bnew(i,lt,tr,x,h)
errunlimited(i,lt,tr,x,h)=0.d0
ctunlimited(i,lt,tr,x,h)=cmin



else
!guesses=min(max(bguess(i,lt,tr,x,h),blimited(i,lt,tr,x,h)),bvect(n_b))
guesses=bguess(i,lt,tr,x,h)


!!! Subroutine solve3 generates policies, shadow values, and the errors in equations as output!!!
call solve3(i,lt,tr,x,h,guesses,choices,shadows,errors,iter3soln(i,lt,tr,x,h),unlimited)

bnew(i,lt,tr,x,h)=choices
mu(i,lt,tr,x,h)=shadows
errmat(i,lt,tr,x,h)=errors

bunlimited(i,lt,tr,x,h)=unlimited(1)
errunlimited(i,lt,tr,x,h)=unlimited(2)
ctunlimited(i,lt,tr,x,h)=unlimited(3)



endif


!$OMP END ORDERED
enddo
enddo
enddo
enddo
enddo
 !$OMP END DO
 !$OMP END PARALLEL

!!! Difsol is a simple sum of max abs val distance measure. !!!
difsol=maxval(abs(bnew-b))

call date_and_time(VALUES=t1)
sec=f(5)*t1(5)+f(6)*t1(6)+f(7)*t1(7)+f(8)*t1(8)

!write(*,*) 'sol3 takes',sec-sec0

!!! These are more sophisticated distance measures, adjusting by the level of the policy in each point.!!!

difbrel=0.d0
difprel=0.d0
do i=1,n_b
do lt=1,n_z
do tr=1,n_g
do x=1,n_kappa
do h=1,n_omega

difbrel(i,lt,tr,x,h)=abs((b(i,lt,tr,x,h)-bnew(i,lt,tr,x,h))/&
(abs(b(i,lt,tr,x,h)+abs(bnew(i,lt,tr,x,h))*0.5d0+0.01d0)))

enddo
enddo
enddo
enddo
enddo


!bguess=bnew
bguess=bunlimited


!!! Relaxation parameter in updating. !!!
grad=0.d0
b=grad*b+(1-grad)*bnew;



!write(*,*) 'b updated'

!!! Find the c's and lambda corresponding to our policies !!!
 !$OMP PARALLEL private(i,lt,tr,x,h)
 !$OMP DO SCHEDULE(RUNTIME)

do i=1,n_b
do lt=1,n_z
do tr=1,n_g
do x=1,n_kappa
do h=1,n_omega

ct(i,lt,tr,x,h)=max(cmin,fixedt+yfun(lt,tr,1)+bvect(i)*Rfun(i,lt,tr,x,h)-b(i,lt,tr,x,h)*gvect(tr))
cn(i,lt,tr,x,h)=max(cmin,fixedn+yfun(lt,tr,2))

ctlimited(i,lt,tr,x,h)=max(cmin,fixedt+yfun(lt,tr,1)+bvect(i)*Rfun(i,lt,tr,x,h)-blimited(i,lt,tr,x,h)*gvect(tr))


enddo
enddo
enddo
enddo
enddo

 !$OMP END DO
 !$OMP END PARALLEL
 


 
 
 

 !!! Find the prices pn !!!

 !$OMP PARALLEL private(i,lt,tr,x,h)
 !$OMP DO SCHEDULE(RUNTIME)

do i=1,n_b
do lt=1,n_z
do tr=1,n_g
do x=1,n_kappa
do h=1,n_omega
pnnew(i,lt,tr,x,h)=(1-omega)/omegavect(h)*(ct(i,lt,tr,x,h)/cn(i,lt,tr,x,h))**(eta+1)
enddo
enddo
enddo
enddo
enddo

 !$OMP END DO
 !$OMP END PARALLEL
 
difprel=0.d0
do i=1,n_b
do lt=1,n_z
do tr=1,n_g
do x=1,n_kappa
do h=1,n_omega

difprel(i,lt,tr,x,h)=abs((pn(i,lt,tr,x,h)-pnnew(i,lt,tr,x,h))/&
(abs(pn(i,lt,tr,x,h)+abs(pnnew(i,lt,tr,x,h))*0.5d0+0.01d0)))

enddo
enddo
enddo
enddo
enddo

if (maxval(difbrel)<difbrelbest.and.iter.ne.1) then
difbrelbest=maxval(difbrel)
difsolbest=difsol
difprelbest=maxval(difprel)
endif

!i=12;lt=2;tr=9;x=1;h=1
!
!write(*,*) 'pns start'
!!i=maxloc(difprel,1)!;lt=maxloc(difprel,2);tr=maxloc(difprel,3);x=1;h=1
!
!write(*,*) maxloc(difprel)!,'==',i,lt,tr!,maxval(difprel),maxloc(difbrel)
!write(*,*) pn(i,lt,tr,x,h),pnnew(i,lt,tr,x,h)
!!write(*,*) ct(i,lt,tr,x,h),ctlimited(i,lt,tr,x,h),b(i,lt,tr,x,h),blimited(i,lt,tr,x,h)
!
!write(*,*) 'pns end'
!write(*,*) 'in between',maxval(difbrel),maxval(difprel)


!!! Update the prices !!!
grad=0.!995d0
pn=min(100.d0,grad*pn+(1.d0-grad)*pnnew)
!grad=0.d0
!lamda=grad*lamda+(1.d0-grad)*lamdanew

!write(*,*) 'after solving',sum(bguess),sum(b),sum(pn)

!$OMP  PARALLEL DO SCHEDULE(RUNTIME)  PRIVATE(i,ltnext,trnext,xnext,hnext,ctpr,cnpr,penalty,psipr,mupr,pnpr) 


expnext=0.d0
do i=1,n_b
do ltnext=1,n_z
do trnext=1,n_g
do xnext=1,n_kappa
do hnext=1,n_omega


ctpr=max(cmin,min(ctlimited(i,ltnext,trnext,xnext,hnext),ctunlimited(i,ltnext,trnext,xnext,hnext)))
cnpr=max(cmin,fixedn+yfun(ltnext,trnext,2))

penalty=0.d0

psipr=0.d0;mupr=0.d0;pnpr=0.d0
if (Planner==1) then
mupr=intp2(bvect(i),mu(:,ltnext,trnext,xnext,hnext))
pnpr=min(max(0.d0,(1-omega)/omegavect(hnext)*(ctpr/cnpr)**(eta+1)),100.d0) !Planner New!
psipr=kappanfun(ltnext,xnext)*pnpr*cnpr/ctpr*(1+eta) !Planner New!
endif

expnext(i,ltnext,trnext,xnext,hnext)=&
(omegavect(hnext)*ctpr**(-eta)+(1-omega)*cnpr**(-eta))**(-1/eta-1)*omegavect(hnext)*ctpr**(-eta-1)&
                *((omegavect(hnext)*ctpr**(-eta)+(1-omega)*cnpr**(-eta))**(-1/eta))**(-sigma)+mupr*psipr !Planner New!

enddo
enddo
enddo
enddo
enddo

!if (Planner==1) then
!expfun=expfun-mu(i,lt,tr,x,h)*(kappanfun(lt,x)*pn(i,lt,tr,x,h)*cn(i,lt,tr,x,h)/ct(i,lt,tr,x,h)*(1+eta)-1.d0)
!endif

 !$OMP END PARALLEL DO
 
 
 
!$OMP  PARALLEL DO SCHEDULE(RUNTIME) private(i,lt,tr,x,h,ltnext,trnext,xnext,hnext,sum0)


do i=1,n_b
do lt=1,n_z
do tr=1,n_g
do x=1,n_kappa
do h=1,n_omega
sum0=0.d0
  !$OMP  PARALLEL DO &
  !$OMP  PRIVATE(ltnext,trnext,xnext,hnext) SHARED(i,lt,tr,x,h)      &
  !$OMP  REDUCTION(+ : sum0) 

do ltnext=1,n_z
do trnext=1,n_g
do xnext=1,n_kappa
do hnext=1,n_omega

sum0=sum0+pg(trnext,tr)*pz(ltnext,lt)*pkappa(xnext,x)*pomega(hnext,h)*Rfun(1,1,ltnext)*expnext(i,ltnext,trnext,xnext,hnext)

enddo
enddo
enddo
enddo
  !$OMP END PARALLEL DO

expmat(i,lt,tr,x,h)=bata/gvect(tr)**sigma*sum0

!if (Planner==1) then
!expfun=expfun-mu(i,lt,tr,x,h)*(kappanfun(lt,x)*pn(i,lt,tr,x,h)*cn(i,lt,tr,x,h)/ct(i,lt,tr,x,h)*(1+eta)-1.d0)
!endif

enddo
enddo
enddo
enddo
enddo

 !$OMP END PARALLEL DO







call date_and_time(VALUES=t1)
sec=f(5)*t1(5)+f(6)*t1(6)+f(7)*t1(7)+f(8)*t1(8)

!write(*,*) 'updated prices at',sec-sec0


if (1==2.and.(iter-25*int(iter/25)==0.or.iter==n_iter.or.maxval(difbrel)<=difcrit)) then

call date_and_time(VALUES=t1)
sec=f(5)*t1(5)+f(6)*t1(6)+f(7)*t1(7)+f(8)*t1(8)


!!! Every 25 iteration, compute the ergodic distribution and obtain debt/output levels and the frequency of binding constraint!!!
!!! This is not crucial, just for us to see where we are in the calibration before running the simulations (in separate file)!!!

iterdens=0;densdif=100.d0
do while (iterdens<100.and.densdif>1E-14)

phinew=0.d0

do i=1,n_b
do lt=1,n_z
do tr=1,n_g
do x=1,n_kappa
do h=1,n_omega

bopt=b(i,lt,tr,x,h)
index_b=min(max(floor((bopt-bvect(1))/(bvect(2)-bvect(1))),1),n_b-1)

do while (bvect(index_b)<bopt.and.index_b<n_b)
index_b=index_b+1
enddo
index_b=max(index_b-1,1)
!write(*,*) 'in int0',index_b


ratio_b = (bopt - bvect(index_b)) / (bvect(index_b+1) - bvect(index_b))



   do i_b=0,1
        contr = ((1-i_b)*(1 - ratio_b) + i_b*ratio_b)*phi(i,lt,tr,x,h)

   if (contr.ne.contr) then
   write(*,*) 'contr',(1-i_b)*(1 - ratio_b) + i_b*ratio_b,phi(i,lt,tr,x,h)
   pause
   endif

   do ll=1,n_z
   do trnext=1,n_g
   do xnext=1,n_kappa
   do hnext=1,n_omega
   phinew(index_b + i_b,ll,trnext,xnext,hnext)=phinew(index_b + i_b,ll,trnext,xnext,hnext)&
+pg(trnext,tr)*pz(ll,lt)*pkappa(x,xnext)*pomega(hnext,h)*contr
   enddo
   enddo
   enddo
   enddo

   end do

enddo
enddo
enddo
enddo
enddo

if (abs(sum(phinew)-1.d0)>1E-5) then
write*, 'phi new not density', sum(phinew)
pause
endif

densdif=sum(abs(phinew-phi))/n_b/n_z/n_g/n_kappa/n_omega
!write(*,*) densdif

phi=phinew

iterdens=iterdens+1
end do


avgb=0.d0;avgboy=0.d0;avgzt=0.d0;avgzn=0.d0
blim=0.d0;mupos=0.d0
avgdifbrel=0.d0

do i=1,n_b
do lt=1,n_z
do tr=1,n_g
do x=1,n_kappa
do h=1,n_omega
avgb=avgb+phi(i,lt,tr,x,h)*bvect(i)
avgboy=avgboy+phi(i,lt,tr,x,h)*bvect(i)/(coefq*pn(i,lt,tr,x,h)*yfun(lt,tr,2)+coefq*yfun(lt,tr,1))
avgzt=avgzt+phi(i,lt,tr,x,h)*ztvect(lt)
avgzn=avgzn+phi(i,lt,tr,x,h)*znvect(lt)


if (mu(i,lt,tr,x,h)>0.d0) mupos=mupos+phi(i,lt,tr,x,h)


limitb=-(kappanfun(lt,x)*coefq*pn(i,lt,tr,x,h)*yfun(lt,tr,2)+kappatfun(lt,x)*coefq*yfun(lt,tr,1))/&
(Rlim(i,lt,tr,x,h)*gvect(tr))
if (b(i,lt,tr,x,h)<=limitb+allow) blim=blim+phi(i,lt,tr,x,h)

!if (mu(i,lt,tr,x,h)>0.d0)  then 
!write(*,*) mu(i,lt,tr,x,h),limitb,b(i,lt,tr,x,h)
!pause
!endif

avgdifbrel=avgdifbrel+phi(i,lt,tr,x,h)*difbrel(i,lt,tr,x,h)

enddo
enddo
enddo
enddo
enddo

call date_and_time(VALUES=t1)
sec0=f(5)*t1(5)+f(6)*t1(6)+f(7)*t1(7)+f(8)*t1(8)

write(*,*) 'density took',sec0-sec,'secs'
write(*,*) 'avg(b)=',avgb,'avg(boy)=',avgboy,'mu pos',mupos

open (unit=15,file ="outputs/moments.txt")
write(15,*) mupos
write(15,*) avgboy
write(15,*) avgb
write(15,*) maxval(difbrel)
write(15,*) sum(difbrel)/(n_b*n_z*n_g*n_kappa*n_omega)
write(15,*) avgdifbrel
write(15,*) sum(abs(errmat))/(n_b*n_z*n_g*n_kappa*n_omega)
close(15)


endif

!if (iter<50.or.iter-int(iter/100)*100==0) then
write(*,*) '--------------------------------'
write(*,*) iter,difsol,maxval(difbrel),maxval(difprel)
!write(*,*) maxloc(difbrel)
write(*,*) sum(abs(errmat))/(n_b*n_z*n_g*n_kappa*n_omega),maxval(abs(errmat))
write(*,*) count(b<bvect(1)+1E-5),count(b>bvect(n_b)-1E-5),sum(iter3soln)*1.d0/(n_b*n_z*n_g*n_kappa*n_omega),maxval(iter3soln)
write(*,*) '--------------------------------'
!endif
iter=iter+1

open (unit=204,file ="outputs/dif.txt")
write(204,*) maxval(difbrel)
write(204,*) maxval(difprel)
write(204,*) difsol
write(204,*) iter
close(204)

nothing=0.d0

!!!! These are the policy functions.
if ((iter-int(iter/10)*10==0.or.iter==n_iter.or.n_iter3>=100).and.difbrelbest<difrelbestold.and.iter.ne.2) then
!if (1==1) then

write(*,*) 'writing b and guess'
!if (difbrelbest<difbrelbestold) then
!if (relaxkn==1d0 .and. relaxkt==1d0 .and. kappat .lt.90d0) then
  open (unit=204,file ="outputs/b.txt")
!elseif (relaxkn.lt.1d0 .and. relaxkt==1d0 .and. kappat .lt.90d0) then
!  open (unit=204,file ="outputs/b_nokn.txt")
!elseif (relaxkn==1d0 .and. relaxkt.lt.1d0 .and. kappat .lt.90d0) then
!  open (unit=204,file ="outputs/b_nokt.txt")
!elseif (relaxkn==1d0 .and. relaxkt==1d0 .and. kappat .gt.90d0) then
!  open (unit=204,file ="outputs/b_noblim.txt")
!elseif (relaxkn.lt.1d0 .and. relaxkt.lt.1d0 .and. kappat .lt.90d0) then
!  open (unit=204,file ="outputs/b_noirrev.txt")
!elseif (relaxkn.lt.1d0 .and. relaxkt.lt.1d0 .and. kappat .gt.90d0) then
!  open (unit=204,file ="outputs/b_noconstraints.txt")
!endif
!open (unit=204,file ="outputs/b.txt")
do i=1,n_b
do lt=1,n_z
do tr=1,n_g
do x=1,n_kappa
do h=1,n_omega
!write(204,'(F14.9,A,F14.9,A,F14.9,A,F14.9&
!&,A,F18.5,A,F14.9,A,F14.9&
!&,A,F18.5,A,F19.11,A,I2)') &
!bvect(i),',',ztvect(lt),',',znvect(lt),',',b(i,lt,tr,x,h),',',&
!mu(i,lt,tr,x,h),',',ct(i,lt,tr,x,h),',',cn(i,lt,tr,x,h),',',&
!pn(i,lt,tr,x,h),',',min(1000.d0,errmat(i,lt,tr,x,h)),',',cminhappens(i,lt,tr,x,h)

write(204,'(F14.9,A,F14.9,A,F14.9,A,F14.9,A,F14.9&
&,A,F18.5,A,F18.5,A,F18.5&
&,A,F18.5,A,F20.13,A,I2,A,F12.8,A,F12.8,A,F14.9,A,F14.9,A,F14.9,A,F14.9,A,F24.17)') &
bvect(i),',',ztvect(lt),',',znvect(lt),',',gvect(tr),',',b(i,lt,tr,x,h),',',&
mu(i,lt,tr,x,h),',',ct(i,lt,tr,x,h),',',cn(i,lt,tr,x,h),',',&
pn(i,lt,tr,x,h),',',min(1000.d0,errmat(i,lt,tr,x,h)),',',cminhappens(i,lt,tr,x,h)&
,',',difbrel(i,lt,tr,x,h),',',difprel(i,lt,tr,x,h),',',blimited(i,lt,tr,x,h),',',bunlimited(i,lt,tr,x,h)&
,',',ctlimited(i,lt,tr,x,h),',',ctunlimited(i,lt,tr,x,h),',',min(100000.d0,expmat(i,lt,tr,x,h))


enddo
enddo
enddo
enddo
enddo
close(204)


open (unit=15,file ="outputs/guesses.txt")
do i=1,n_b
do lt=1,n_z
do tr=1,n_g
do x=1,n_kappa
do h=1,n_omega
write(15,'(F14.9)') bguess(i,lt,tr,x,h)
enddo
enddo
enddo
enddo
enddo
close(15)

!write(*,*) 'just written',sum(bguess),sum(b),sum(pn)





open (unit=204,file ="outputs/difbest.txt")
write(204,*) difbrelbest
write(204,*) difprelbest
write(204,*) difsolbest
write(204,*) iter
close(204)

difrelbestold=difbrelbest
endif

enddo
end program

!!!! The rhs of the euler equation, the part that uses expectations!!!
!double precision function expfun(bpr,lt,i,tr,x,h,ind)
!implicit none
!integer, parameter :: n_b=200,n_z=16,n_zt=4,n_zn=4,n_g=1,n_kappa=1,n_omega=1,n_iter=30000
!COMMON /PARAM/ bvect,kappavect,omegavect,r,cmin,bata,kappat,kappan,omega,sigma,eta,&
!fixedt,fixedn,coefq,n_iter3,Planner,ztvect,znvect,gvect,pz,pg,pkappa,pomega,allow,gradb
!DOUBLE PRECISION bvect(n_b),kappavect(n_kappa),omegavect(n_omega),r,cmin,bata,kappat,kappan,omega&
!,sigma,eta,fixedt,fixedn,coefq,pz(n_z,n_z),pg(n_g,n_g),pkappa(n_kappa,n_kappa),pomega(n_omega,n_omega),&
!ztvect(n_z),znvect(n_z),gvect(n_g),allow,gradb
!INTEGER n_iter3,Planner
!SAVE /PARAM/
!COMMON /MATS/ b,ct,cn,pn,expmat
!DOUBLE PRECISION, dimension(n_b,n_z,n_g,n_kappa,n_omega):: b,ct,cn,pn,expmat
!SAVE /MATS/
!COMMON /MATS2/ mu,ctlimited,ctunlimited
!DOUBLE PRECISION, dimension(n_b,n_z,n_g,n_kappa,n_omega):: mu,ctlimited,ctunlimited
!SAVE /MATS2/
!
!
!integer ind,lt,ltnext,i,h,hnext,x,xnext,tr,trnext
!double precision bpr,yfun,yderfun,ctpr,cnpr,intp2,lamdapr,sum0,pnpr,b2pr,flow,ctprAGG,cnprAGG,&
!ctraw,cnraw,penalty,Rfun,mupr,psipr,kappanfun
!
!
!
!if (ind==1) then
!
!sum0=0.d0
!
!do ltnext=1,n_z
!do trnext=1,n_g
!do xnext=1,n_kappa
!do hnext=1,n_omega
!
!
!
!!ctpr=max(cmin,intp2(bpr,ktpr,ct(:,:,ltnext)))
!!cnpr=max(cmin,fixedn+yfun(ltnext,ktot-ktpr,2))
!
!ctraw=min(intp2(bpr,ctlimited(:,ltnext,trnext,xnext,hnext)),intp2(bpr,ctunlimited(:,ltnext,trnext,xnext,hnext)))
!cnraw=fixedn+yfun(ltnext,trnext,2)
!
!ctpr=max(cmin,ctraw)
!cnpr=max(cmin,cnraw)
!
!penalty=0.d0
!
!psipr=0.d0;mupr=0.d0;pnpr=0.d0
!if (Planner==1) then
!mupr=intp2(bpr,mu(:,ltnext,trnext,xnext,hnext))
!pnpr=min(max(0.d0,(1-omega)/omegavect(hnext)*(ctpr/cnpr)**(eta+1)),100.d0) !Planner New!
!psipr=kappanfun(ltnext,xnext)*pnpr*cnpr/ctpr*(1+eta) !Planner New!
!endif
!
!
!
!lamdapr=(omegavect(hnext)*ctpr**(-eta)+(1-omega)*cnpr**(-eta))**(-1/eta-1)*omegavect(hnext)*ctpr**(-eta-1)&
!                *((omegavect(hnext)*ctpr**(-eta)+(1-omega)*cnpr**(-eta))**(-1/eta))**(-sigma)+mupr*psipr !Planner New!
!
!sum0=sum0+pg(trnext,tr)*pz(ltnext,lt)*pkappa(xnext,x)*pomega(hnext,h)*Rfun(1,1,ltnext)*(lamdapr+penalty)
!enddo
!enddo
!enddo
!enddo
!expfun=bata/gvect(tr)**sigma*sum0
!
!if (Planner==1) then
!expfun=expfun-mu(i,lt,tr,x,h)*(kappanfun(lt,x)*pn(i,lt,tr,x,h)*cn(i,lt,tr,x,h)/ct(i,lt,tr,x,h)*(1+eta)-1.d0)
!endif
!
!
!endif
!end function



!!! Solve for the three variables (b',kt',kn') and spit out the lagrange mults (mu,gammat,gamman) and errors. !!!
subroutine solve3(i,lt,tr,x,h,guesses,choices,shadows,errors,iterb,unlimited)
implicit none
integer, parameter :: n_b=200,n_z=16,n_zt=4,n_zn=4,n_g=1,n_kappa=1,n_omega=1,n_iter=30000
COMMON /PARAM/ bvect,kappavect,omegavect,r,cmin,bata,kappat,kappan,omega,sigma,eta,&
fixedt,fixedn,coefq,n_iter3,Planner,ztvect,znvect,gvect,pz,pg,pkappa,pomega,allow,gradb
DOUBLE PRECISION bvect(n_b),kappavect(n_kappa),omegavect(n_omega),r,cmin,bata,kappat,kappan,omega&
,sigma,eta,fixedt,fixedn,coefq,pz(n_z,n_z),pg(n_g,n_g),pkappa(n_kappa,n_kappa),pomega(n_omega,n_omega),&
ztvect(n_z),znvect(n_z),gvect(n_g),allow,gradb
INTEGER n_iter3,Planner
SAVE /PARAM/
COMMON /MATS/ b,ct,cn,pn,expmat
DOUBLE PRECISION, dimension(n_b,n_z,n_g,n_kappa,n_omega):: b,ct,cn,pn,expmat
SAVE /MATS/



double precision guesses,choices,shadows,errors,bpr,dif3,consn,const,uT,limitb,eq1,&
adjb,critdifb,bprmin,bprmax,bprbest,muprbest,eq1best,&
grad,zfun,bprold,difbinbest,penalty,yfun,Rlim,Rfun,kappatfun,kappanfun,errb,unlimited(3),intp2
integer i,lt,iterb,n_iterb,iterbbest,timesimproved,triedb,h,x,tr,ltnext,trnext,hnext,xnext

adjb=10.d0
critdifb=1E-15
grad=0.d0
!iter3best=0
!relaxkt=1d0
!relaxkn=1d0
!adjkn=1.d0;adjkt=1.d0
bpr=guesses
bprbest=bpr;

limitb=max(bvect(1)-0.d0,&
-(kappanfun(lt,x)*coefq*pn(i,lt,tr,x,h)*yfun(lt,tr,2)+kappatfun(lt,x)*coefq*yfun(lt,tr,1))/&
(Rlim(i,lt,tr,x,h)*gvect(tr)))

!dif3=10000.d0!;difbinbest=10000.d0
!timesimproved=0

!do while (iter3<n_iter3.and.dif3>critdif3)

!!!!!!! Update bpr !!!!!!!!!!!!
iterb=1;eq1=100.d0;n_iterb=1000000
bprmin=bvect(1)-0.d0!limitb
!;bprmax=bvect(n_b)!+15.d0


bprmax=max(bprmin,min((fixedt+yfun(lt,tr,1)+bvect(i)*Rfun(i,lt,tr,x,h))/gvect(tr),bvect(n_b)))

if (bprmin>bprmax-1E-15) bpr=bprmin

!mupr=0.d0!;bpr=limitb
bprold=bpr;triedb=1;errb=10000.d0
!do while (iterb<=n_iterb.and.abs(eq1)>critdifb.and.bprmin<bprmax-1E-12)
do while (iterb<=n_iterb.and.errb>critdifb.and.bprmin<bprmax-1E-15)

const=fixedt+yfun(lt,tr,1)+bvect(i)*Rfun(i,lt,tr,x,h)-bpr*gvect(tr)
consn=fixedn+yfun(lt,tr,2)
penalty=0.d0
if (const<cmin) then
penalty=penalty+1000.d0*(cmin-const)**2
const=cmin
endif

if (consn<cmin) then
penalty=penalty+1000.d0*(cmin-consn)**2
consn=cmin
endif
uT=(omegavect(h)*const**(-eta)+(1-omega)*consn**(-eta))**(-1/eta-1)*omegavect(h)*const**(-eta-1)&
    *((omegavect(h)*const**(-eta)+(1-omega)*consn**(-eta))**(-1/eta))**(-sigma)+penalty   
eq1=uT-intp2(bpr,expmat(:,lt,tr,x,h))

errb=abs(eq1)
!if (iterb==1) errb=errb/10.d0

if (errb>critdifb.and.iterb<n_iterb) then !a!


if (eq1<0.d0) then !b!
bprmin=bpr
if (bpr-eq1*adjb>=bprmax) then !c!
bpr=0.5*bpr+0.5*bprmax
else !c!
bpr=bpr-eq1*adjb
endif !c!
else if (eq1>0.d0) then !b!
bprmax=bpr
if (bpr-eq1*adjb<=bprmin) then !c!
bpr=0.5*bpr+0.5*bprmin
else !c!
bpr=bpr-eq1*adjb
endif !c!
endif !b!

endif !a!
iterb=iterb+1
enddo

unlimited(1)=bpr
unlimited(2)=eq1
unlimited(3)=const

!if (errb>0.00001d0) then
!write(*,*) 'in solve',errb,eq1,iterb,bpr,bprmin,bprmax
!write(*,*) i,lt,tr
!write(*,*) 'guess and limit',guesses,limitb
!pause
!endif

bprbest=bpr
muprbest=0.d0
if (bprbest<=limitb+allow) then

!gradb=0.95d0
!write(*,*) bpr,bprold,iterb

!bpr=gradb*bprold+(1.d0-gradb)*bpr
!!!!!!!!!!! Updated bpr !!!!!!!!!!!!!!!!!!!

bpr=limitb
bprbest=limitb

const=fixedt+yfun(lt,tr,1)+bvect(i)*Rfun(i,lt,tr,x,h)-bpr*gvect(tr)
consn=fixedn+yfun(lt,tr,2)

penalty=0.d0
if (const<cmin) then
penalty=penalty+1000.d0*(cmin-const)**2
const=cmin
endif

if (consn<cmin) then
penalty=penalty+1000.d0*(cmin-consn)**2
consn=cmin
endif


uT=(omegavect(h)*const**(-eta)+(1-omega)*consn**(-eta))**(-1/eta-1)*omegavect(h)*const**(-eta-1)&
    *((omegavect(h)*const**(-eta)+(1-omega)*consn**(-eta))**(-1/eta))**(-sigma)+penalty

eq1=uT-intp2(bpr,expmat(:,lt,tr,x,h))       !expfun(bpr,lt,i,tr,x,h,1) 
muprbest=eq1


!dif3=0.d0
!if (bpr>limitb+allow.or.eq1<0.d0) dif3=abs(eq1)

endif

eq1best=eq1-muprbest

!iterbbest=iterb
!timesimproved=timesimproved+1
!iter3best=iter3
!endif

!iter3=iter3+1
!
!enddo

!if (bpr<=limitb+allow.and.eq1<=0.d0) then
!write(*,*) bpr,bprbest,limitb,eq1,iterb,muprbest
!pause
!endif


choices=bprbest
shadows=muprbest
errors=eq1best

end subroutine

!!!! Linear interpolation in 1 dimension !!!!
double precision function intp2(bpr,mat)
implicit none
integer, parameter :: n_b=200,n_z=16,n_zt=4,n_zn=4,n_g=1,n_kappa=1,n_omega=1,n_iter=30000
COMMON /PARAM/ bvect,kappavect,omegavect,r,cmin,bata,kappat,kappan,omega,sigma,eta,&
fixedt,fixedn,coefq,n_iter3,Planner,ztvect,znvect,gvect,pz,pg,pkappa,pomega,allow,gradb
DOUBLE PRECISION bvect(n_b),kappavect(n_kappa),omegavect(n_omega),r,cmin,bata,kappat,kappan,omega&
,sigma,eta,fixedt,fixedn,coefq,pz(n_z,n_z),pg(n_g,n_g),pkappa(n_kappa,n_kappa),pomega(n_omega,n_omega),&
ztvect(n_z),znvect(n_z),gvect(n_g),allow,gradb
INTEGER n_iter3,Planner
SAVE /PARAM/

double precision bpr,mat(n_b),ratio_b,weight,acum
integer index_b,i_b


!relaxkt=1d0
!relaxkn=1d0

!index_b=1
index_b=min(max(floor((bpr-bvect(1))/(bvect(2)-bvect(1))),1),n_b-1)
if (bvect(index_b)>bpr.and.bpr>bvect(1)) then
write(*,*) 'problem in intp3 b',index_b,bpr
pause
endif

do while (bvect(index_b)<bpr.and.index_b<n_b)
index_b=index_b+1
enddo
index_b=max(index_b-1,1)
!write(*,*) 'in int0',index_b

ratio_b = (bpr - bvect(index_b)) / (bvect(index_b+1) - bvect(index_b))


acum=0
   do i_b=0,1
        weight = ((1-i_b)*(1 - ratio_b) + i_b*ratio_b)
        acum = acum + mat(index_b + i_b)*weight
   end do


intp2=acum



end function


!!! Given the shock, we define what will be the productivity in each state. This way, we allow for one sector to have constant z.
double precision function zfun(ll,sect)
implicit none
integer, parameter :: n_b=200,n_z=16,n_zt=4,n_zn=4,n_g=1,n_kappa=1,n_omega=1,n_iter=30000
COMMON /PARAM/ bvect,kappavect,omegavect,r,cmin,bata,kappat,kappan,omega,sigma,eta,&
fixedt,fixedn,coefq,n_iter3,Planner,ztvect,znvect,gvect,pz,pg,pkappa,pomega,allow,gradb
DOUBLE PRECISION bvect(n_b),kappavect(n_kappa),omegavect(n_omega),r,cmin,bata,kappat,kappan,omega&
,sigma,eta,fixedt,fixedn,coefq,pz(n_z,n_z),pg(n_g,n_g),pkappa(n_kappa,n_kappa),pomega(n_omega,n_omega),&
ztvect(n_z),znvect(n_z),gvect(n_g),allow,gradb
INTEGER n_iter3,Planner
SAVE /PARAM/

integer ll,sect


!relaxkt=1d0
!relaxkn=1d0

if (sect==1) then
zfun=ztvect(ll)
else
zfun=znvect(ll)
endif

end function

!!! Output function for either of the sectors !!!
double precision function yfun(ll,tr,sect)
implicit none
integer, parameter :: n_b=200,n_z=16,n_zt=4,n_zn=4,n_g=1,n_kappa=1,n_omega=1,n_iter=30000
COMMON /PARAM/ bvect,kappavect,omegavect,r,cmin,bata,kappat,kappan,omega,sigma,eta,&
fixedt,fixedn,coefq,n_iter3,Planner,ztvect,znvect,gvect,pz,pg,pkappa,pomega,allow,gradb
DOUBLE PRECISION bvect(n_b),kappavect(n_kappa),omegavect(n_omega),r,cmin,bata,kappat,kappan,omega&
,sigma,eta,fixedt,fixedn,coefq,pz(n_z,n_z),pg(n_g,n_g),pkappa(n_kappa,n_kappa),pomega(n_omega,n_omega),&
ztvect(n_z),znvect(n_z),gvect(n_g),allow,gradb
INTEGER n_iter3,Planner
SAVE /PARAM/

integer ll,sect,tr


if (sect==1) then
yfun=gvect(tr)*ztvect(ll)
else
yfun=gvect(tr)*znvect(ll)
endif

end function


!!!! In case we want to make the interest rate a function of the states, we can use this function!!!!
double precision function Rfun(i,lt,tr,x,h)
implicit none
integer, parameter :: n_b=200,n_z=16,n_zt=4,n_zn=4,n_g=1,n_kappa=1,n_omega=1,n_iter=30000
COMMON /PARAM/ bvect,kappavect,omegavect,r,cmin,bata,kappat,kappan,omega,sigma,eta,&
fixedt,fixedn,coefq,n_iter3,Planner,ztvect,znvect,gvect,pz,pg,pkappa,pomega,allow,gradb
DOUBLE PRECISION bvect(n_b),kappavect(n_kappa),omegavect(n_omega),r,cmin,bata,kappat,kappan,omega&
,sigma,eta,fixedt,fixedn,coefq,pz(n_z,n_z),pg(n_g,n_g),pkappa(n_kappa,n_kappa),pomega(n_omega,n_omega),&
ztvect(n_z),znvect(n_z),gvect(n_g),allow,gradb
INTEGER n_iter3,Planner
SAVE /PARAM/

integer i,lt,tr,x,h

!rfun=r*exp(maxval(ztvect)-ztvect(lt))
Rfun=(1+r)!/ztvect(lt)**2.d0

!if (ztvect(lt)<=minval(ztvect)+1E-4) Rfun=(1+r)*1.25d0

!if (ztvect(lt)<1.d0) Rfun=(1+r)/ztvect(lt)**2.d0

end function

!!!! In case we want to use an interest rate in the borrowing limit a function, we can just change this !!!!
double precision function Rlim(i,lt,tr,x,h)
implicit none
integer, parameter :: n_b=200,n_z=16,n_zt=4,n_zn=4,n_g=1,n_kappa=1,n_omega=1,n_iter=30000
COMMON /PARAM/ bvect,kappavect,omegavect,r,cmin,bata,kappat,kappan,omega,sigma,eta,&
fixedt,fixedn,coefq,n_iter3,Planner,ztvect,znvect,gvect,pz,pg,pkappa,pomega,allow,gradb
DOUBLE PRECISION bvect(n_b),kappavect(n_kappa),omegavect(n_omega),r,cmin,bata,kappat,kappan,omega&
,sigma,eta,fixedt,fixedn,coefq,pz(n_z,n_z),pg(n_g,n_g),pkappa(n_kappa,n_kappa),pomega(n_omega,n_omega),&
ztvect(n_z),znvect(n_z),gvect(n_g),allow,gradb
INTEGER n_iter3,Planner
SAVE /PARAM/

integer i,lt,tr,x,h

!rlim=r*exp(maxval(ztvect)-ztvect(lt))
Rlim=1.d0 !(1+r)/ztvect(lt)**2.d0
!Rlim=1.d0/gvect(tr)

end function

!!!! We can make the tightness parameter of the constraint a function with this !!!!
double precision function kappatfun(lt,x)
implicit none
integer, parameter :: n_b=200,n_z=16,n_zt=4,n_zn=4,n_g=1,n_kappa=1,n_omega=1,n_iter=30000
COMMON /PARAM/ bvect,kappavect,omegavect,r,cmin,bata,kappat,kappan,omega,sigma,eta,&
fixedt,fixedn,coefq,n_iter3,Planner,ztvect,znvect,gvect,pz,pg,pkappa,pomega,allow,gradb
DOUBLE PRECISION bvect(n_b),kappavect(n_kappa),omegavect(n_omega),r,cmin,bata,kappat,kappan,omega&
,sigma,eta,fixedt,fixedn,coefq,pz(n_z,n_z),pg(n_g,n_g),pkappa(n_kappa,n_kappa),pomega(n_omega,n_omega),&
ztvect(n_z),znvect(n_z),gvect(n_g),allow,gradb
INTEGER n_iter3,Planner
SAVE /PARAM/

integer i,lt,x

kappatfun=kappavect(x)!*ztvect(lt)**2.d0

end function

!!!! We can make the tightness parameter of the constraint a function with this !!!!
double precision function kappanfun(lt,x)
implicit none
integer, parameter :: n_b=200,n_z=16,n_zt=4,n_zn=4,n_g=1,n_kappa=1,n_omega=1,n_iter=30000
COMMON /PARAM/ bvect,kappavect,omegavect,r,cmin,bata,kappat,kappan,omega,sigma,eta,&
fixedt,fixedn,coefq,n_iter3,Planner,ztvect,znvect,gvect,pz,pg,pkappa,pomega,allow,gradb
DOUBLE PRECISION bvect(n_b),kappavect(n_kappa),omegavect(n_omega),r,cmin,bata,kappat,kappan,omega&
,sigma,eta,fixedt,fixedn,coefq,pz(n_z,n_z),pg(n_g,n_g),pkappa(n_kappa,n_kappa),pomega(n_omega,n_omega),&
ztvect(n_z),znvect(n_z),gvect(n_g),allow,gradb
INTEGER n_iter3,Planner
SAVE /PARAM/

integer i,lt,x

kappanfun=kappavect(x)!*ztvect(lt)**2.d0

end function

!!!! Only used for the Tauchen.!!!
double precision function normalcdf(x)
implicit none

double precision x,pi,inside,erf,a1,a2,a3,a4,a5,rat,t,p


p = 0.3275911d0
a1 = 0.254829592d0
a2 = -0.284496736d0
a3 = 1.421413741d0
a4 = -1.453152027d0
a5 = 1.061405429d0


inside=(abs(x)-0.d0)/(sqrt(2.d0))
rat=1.d0/(1.d0+p*inside)


erf=1-(a1*rat+a2*rat**2+a3*rat**3+a4*rat**4+a5*rat**5)*exp(-inside**2)
if (x<0) erf=-erf
!if (x>=0) then
normalcdf=0.5d0*(1+erf)
!else
!normalcdf=1-0.5d0*(1+erf)
!endif


end function
