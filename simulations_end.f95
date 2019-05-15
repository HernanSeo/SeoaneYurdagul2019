program sim
implicit none

integer, parameter :: n_b=400,n_z=121,n_g=25,n_kappa=1,n_omega=1,n_id=2000,n_time=500
COMMON /PARAM/ bvect,kappavect,omegavect,r,cmin,bata,kappat,kappan,omega,sigma,eta,&
fixedt,fixedn,coefq,n_iter3,ztvect,znvect,gvect,pz,pg,pkappa,pomega,allow,gradb
DOUBLE PRECISION bvect(n_b),kappavect(n_kappa),omegavect(n_omega),r,cmin,bata,kappat,kappan,omega&
,sigma,eta,fixedt,fixedn,coefq,pz(n_z,n_z),pg(n_g,n_g),pkappa(n_kappa,n_kappa),pomega(n_omega,n_omega),&
ztvect(n_z),znvect(n_z),gvect(n_g),allow,gradb
INTEGER n_iter3
SAVE /PARAM/
COMMON /MATS/ b,ct,cn,mu,pn,blimited,bunlimited
DOUBLE PRECISION, dimension(n_b,n_z,n_g,n_kappa,n_omega):: b,ct,cn,mu,pn,blimited,bunlimited
SAVE /MATS/

      double precision, dimension(n_id,n_time):: b_sim,ct_sim,cn_sim,yt_sim,yn_sim,c_sim,randomzt,randomzn,&
mu_sim, pn_sim,randomg,tax_sim,blim_sim!,randomkappa,randomomega
      integer, dimension(n_id,n_time):: z_sim,kappa_sim,omega_sim,g_sim

double precision nothing(4),expnow,ctnow,cnnow,lamdanow,psinow,&
tax(n_b,n_z,n_g,n_kappa,n_omega),sum0,ctraw,cnraw,ctpr,cnpr,psipr,mupr,bpr,pnpr
!!!! Functions !!!!
double precision intp2,yfun,Rfun,Rlim,kappatfun,kappanfun
integer id,t,i,j,lt,ii,h,x,tr,ltnext,hnext,xnext,trnext,Planner


Planner=0
!Parameters!
r=0.01d0;bata=0.97d0;kappat=0.32d0;kappan=0.32d0;


omega=0.31d0;sigma=2.d0;eta=1/0.83d0-1.d0;cmin=0.0001d0;fixedt=0.d0;fixedn=0.d0

      open(unit=300,file="outputs/params_end.txt")
      read (300,*) j
      read (300,*) i
      read (300,*) r
      read (300,*) bata
      read (300,*) nothing(1)
      read (300,*) nothing(1)
      read (300,*) kappat
      read (300,*) gradb
      read (300,*) Planner
      close(300)

      kappan=kappat


      open(unit=262,file="outputs/bvect.txt")
      do i=1,n_b
      read (262,*) bvect(i)
      enddo
      close(262)


      open(unit=262,file="outputs/kappavect.txt")
      do i=1,n_kappa
      read (262,*) kappavect(i)
      enddo
      close(262)

      open(unit=242,file="outputs/pkappa.txt")
      do i=1,n_kappa
      read (242,*) (pkappa(i,j),j=1,n_kappa)
      enddo
      close(242)

      open(unit=262,file="outputs/omegavect.txt")
      do i=1,n_omega
      read (262,*) omegavect(i)
      enddo
      close(262)

      open(unit=242,file="outputs/pomega.txt")
      do i=1,n_omega
      read (242,*) (pomega(i,j),j=1,n_omega)
      enddo
      close(242) 

      open(unit=260,file="outputs/gvect.txt")
      do i=1,n_g
      read (260,*) gvect(i)
      enddo
      close(260)

      open(unit=240,file="outputs/pg.txt")
      do i=1,n_g
      read (240,*) (pg(i,j),j=1,n_g)
      enddo
      close(240)   

!call process(n_z)      ! stochastic process
! znvect=senn
! ztvect=se
! pz=p

      open(unit=262,file="outputs/ztvect.txt")
      do i=1,n_z
      read (262,*) ztvect(i)
      enddo
      close(262)
      
      open(unit=262,file="outputs/znvect.txt")
      do i=1,n_z
      read (262,*) znvect(i)
      enddo
      close(262)

      open(unit=242,file="outputs/pz.txt")
      do i=1,n_z
      read (242,*) (pz(i,j),j=1,n_z)
      enddo
      close(242)

!do lt=1,n_z
!write(*,*) lt,ztvect(lt),znvect(lt)
!enddo
!pause

open (unit=204,file ="outputs/b.txt")

do i=1,n_b
do lt=1,n_z
do tr=1,n_g
do x=1,n_kappa
do h=1,n_omega
!read(204,*) nothing(1),nothing(2),nothing(3),b(i,lt,tr,x,h)&
!,mu(i,lt,tr,x,h),ct(i,lt,tr,x,h),cn(i,lt,tr,x,h)&
!,pn(i,lt,tr,x,h),nothing(1),ii

!read(204,*) bvect(i),nothing(1),ztvect(lt),znvect(lt),b(i,lt,tr,x,h),nothing(1),&
!mu(i,lt,tr,x,h),nothing(1),nothing(1),ct(i,lt,tr,x,h),cn(i,lt,tr,x,h),&
!pn(i,lt,tr,x,h),nothing(1),nothing(1),nothing(1),&
!nothing(1),ii

read(204,*) bvect(i),ztvect(lt),znvect(lt),gvect(tr),b(i,lt,tr,x,h),&
mu(i,lt,tr,x,h),ct(i,lt,tr,x,h),cn(i,lt,tr,x,h),&
pn(i,lt,tr,x,h),nothing(1),ii,&
nothing(1),nothing(1),blimited(i,lt,tr,x,h),bunlimited(i,lt,tr,x,h),nothing(2),nothing(2),nothing(2)

!bvect(i),',',ztvect(lt),',',znvect(lt),',',gvect(tr),',',b(i,lt,tr,x,h),',',&
!mu(i,lt,tr,x,h),',',ct(i,lt,tr,x,h),',',cn(i,lt,tr,x,h),',',&
!pn(i,lt,tr,x,h),',',min(1000.d0,errmat(i,lt,tr,x,h)),',',cminhappens(i,lt,tr,x,h)&
!,',',difbrel(i,lt,tr,x,h),',',difprel(i,lt,tr,x,h),',',blimited(i,lt,tr,x,h),',',bunlimited(i,lt,tr,x,h)&
!,',',ctlimited(i,lt,tr,x,h),',',ctunlimited(i,lt,tr,x,h),',',min(100000.d0,expmat(i,lt,tr,x,h))



enddo
enddo
enddo
enddo
enddo
close(204)


if (Planner==1) then
!!!! Compute the optimal tax !!!!!
do i=1,n_b
do lt=1,n_z
do tr=1,n_g
do x=1,n_kappa
do h=1,n_omega
sum0=0.d0

do ltnext=1,n_z
do trnext=1,n_g
do xnext=1,n_kappa
do hnext=1,n_omega




bpr=b(i,lt,tr,x,h)

ctraw=intp2(bpr,ct(:,ltnext,trnext,xnext,hnext))
cnraw=fixedn+yfun(ltnext,trnext,2)

ctpr=max(cmin,ctraw)
cnpr=max(cmin,cnraw)

psipr=0.d0;mupr=0.d0;pnpr=0.d0
mupr=max(0.d0,intp2(bpr,mu(:,ltnext,trnext,xnext,hnext)))
pnpr=min(max(0.d0,(1-omega)/omegavect(hnext)*(ctpr/cnpr)**(eta+1)),100.d0) !Planner New!
psipr=kappanfun(ltnext,xnext)*pnpr*cnpr/ctpr*(1+eta) !Planner New!



sum0=sum0+pg(trnext,tr)*pz(ltnext,lt)*pkappa(xnext,x)*pomega(hnext,h)*Rfun(1,1,ltnext)*mupr*psipr 
enddo
enddo
enddo
enddo
expnow=bata/gvect(tr)**sigma*sum0


ctnow=ct(i,lt,tr,x,h)
cnnow=cn(i,lt,tr,x,h)
lamdanow=(omegavect(h)*ctnow**(-eta)+(1-omega)*cnnow**(-eta))**(-1/eta-1)*omegavect(h)*ctnow**(-eta-1)&
                *((omegavect(h)*ctnow**(-eta)+(1-omega)*cnnow**(-eta))**(-1/eta))**(-sigma)

psinow=kappanfun(lt,x)*pn(i,lt,tr,x,h)*cnnow/ctnow*(1+eta)

tax(i,lt,tr,x,h)=(expnow-mu(i,lt,tr,x,h)*psinow)/lamdanow

enddo
enddo
enddo
enddo
enddo

open (unit=24,file ="outputs/tax.txt")
do i=1,n_b
do lt=1,n_z
do tr=1,n_g
do x=1,n_kappa
do h=1,n_omega
write(24,'(F12.8,A,I2,A,F12.8,A,I2,A,I2,A,F12.8)') &
bvect(i),',',lt,',',gvect(tr),',',x,',',h,',',tax(i,lt,tr,x,h)
enddo
enddo
enddo
enddo
enddo
close(24)

endif
!!!!!!!!!!!!!!!!


!     ii=max(floor((bvect(n_b)+0.2-bvect(1))/(bvect(2)-bvect(1))),1)
!     write(*,*) ii,bvect(ii:min(n_b,ii+2))
!     pause
!
! i=10;jt=n_kt;jn=10;l=2
!write(*,*) intp3(bvect(i),ktvect(jt)+1.,knvect(jn),b(:,:,:,l)),b(i,jt-1,jn,l),b(i,jn,l)
!pause
!              open (unit=201,file = 'outputs/randomzt.txt')
!              open (unit=202,file = 'outputs/randomkappa.txt')
!              open (unit=203,file = 'outputs/randomomega.txt')
              open (unit=201,file = 'outputs/randomzt5m.txt')
!              open (unit=202,file = 'outputs/randomkappa_1m.txt')
!              open (unit=203,file = 'outputs/randomomega_1m.txt')
              open (unit=204,file = 'outputs/randomg5m.txt')
              do i=1,n_id
              do t=1,n_time-1
              read(201,'(F11.8)') randomzt(i,t)
!              read(202,'(F11.8)') randomkappa(i,t)
!              read(203,'(F11.8)') randomomega(i,t)
              read(204,'(F11.8)') randomg(i,t)
              enddo
              enddo
              close(201)
!              close(202)
!              close(203)
              close(204)
!print*,randomzt
        write(*,*) 'starting sim'


         !!! Variables in simulations !!!

        !z_sim: Income shock
        !b_sim: Debt
        !cj_sim: Consumption

        do id=1,n_id
        z_sim(id,1)=int(n_z/2)+1!;zn_sim(id,1)=int(n_zn/2)+1
        b_sim(id,1)=bvect(int(n_b/2))
        kappa_sim(id,1)=1;omega_sim(id,1)=1;g_sim(id,1)=1
        enddo


        do id=1,n_id
        do t=1,n_time-1
               
        yt_sim(id,t)=yfun(z_sim(id,t),g_sim(id,t),1)
        yn_sim(id,t)=yfun(z_sim(id,t),g_sim(id,t),2)    

        pn_sim(id,t)=intp2(b_sim(id,t),pn(:,z_sim(id,t),g_sim(id,t),kappa_sim(id,t),omega_sim(id,t)))
        
        if (pn_sim(id,t)<0.001d0.and.b_sim(id,t)>-0.3d0) then
        write(*,*) b_sim(id,t),z_sim(id,t),g_sim(id,t),kappa_sim(id,t),omega_sim(id,t)
        pause
        endif
        
        
!        b_sim(id,t+1)=max(-(kappanfun(z_sim(id,t),kappa_sim(id,t))*pn_sim(id,t)*yn_sim(id,t)+&
!        kappatfun(z_sim(id,t),kappa_sim(id,t))*yt_sim(id,t))/(Rlim(1,1,z_sim(id,t),g_sim(id,t),1,1)*gvect(g_sim(id,t)))&
!,intp2(b_sim(id,t),b(:,z_sim(id,t),g_sim(id,t),kappa_sim(id,t),omega_sim(id,t))))
        
        b_sim(id,t+1)=max(intp2(b_sim(id,t),blimited(:,z_sim(id,t),g_sim(id,t),kappa_sim(id,t),omega_sim(id,t))),&
                          intp2(b_sim(id,t),bunlimited(:,z_sim(id,t),g_sim(id,t),kappa_sim(id,t),omega_sim(id,t))))
                          
        blim_sim(id,t+1)=intp2(b_sim(id,t),blimited(:,z_sim(id,t),g_sim(id,t),kappa_sim(id,t),omega_sim(id,t)))                 

        mu_sim(id,t)=intp2(b_sim(id,t),mu(:,z_sim(id,t),g_sim(id,t),kappa_sim(id,t),omega_sim(id,t)))
        tax_sim(id,t)=intp2(b_sim(id,t),tax(:,z_sim(id,t),g_sim(id,t),kappa_sim(id,t),omega_sim(id,t)))

        ct_sim(id,t)=max(cmin,fixedt+yt_sim(id,t)+b_sim(id,t)*Rfun(1,1,z_sim(id,t))-b_sim(id,t+1)*gvect(g_sim(id,t)))
        cn_sim(id,t)=max(cmin,fixedn+yn_sim(id,t))

        c_sim(id,t)=(omegavect(omega_sim(id,t))*ct_sim(id,t)**(-eta)+(1-omega)*cn_sim(id,t)**(-eta))**(-1/eta)


       
    

        z_sim(id,t+1)=1
        do while (sum(pz(1:z_sim(id,t+1),z_sim(id,t)))<randomzt(id,t).and.z_sim(id,t+1)<n_z)
        z_sim(id,t+1)=z_sim(id,t+1)+1
        enddo
!         z_sim(id,t+1)=8

        
     
        kappa_sim(id,t+1)=1
!        do while (sum(pkappa(1:kappa_sim(id,t+1),kappa_sim(id,t)))<randomkappa(id,t).and.kappa_sim(id,t+1)<n_kappa)
!        kappa_sim(id,t+1)=kappa_sim(id,t+1)+1
!        enddo
        
        omega_sim(id,t+1)=1
!        do while (sum(pomega(1:omega_sim(id,t+1),omega_sim(id,t)))<randomomega(id,t).and.omega_sim(id,t+1)<n_omega)
!        omega_sim(id,t+1)=omega_sim(id,t+1)+1
!        enddo
        
        g_sim(id,t+1)=1
        do while (sum(pg(1:g_sim(id,t+1),g_sim(id,t)))<randomg(id,t).and.g_sim(id,t+1)<n_g)
        g_sim(id,t+1)=g_sim(id,t+1)+1
        enddo
!        
!        if (t<100) then 
!!        kappa_sim(id,t+1)=1;omega_sim(id,t+1)=2!;z_sim(id,t+1)=7
!        z_sim(id,t+1)=7
!        else if (t==100) then
!!        kappa_sim(id,t+1)=2;omega_sim(id,t+1)=2!;z_sim(id,t+1)=7
!        z_sim(id,t+1)=6
!!        else if (t<103) then
!!!        kappa_sim(id,t+1)=2;omega_sim(id,t+1)=2!;z_sim(id,t+1)=7
!!        z_sim(id,t+1)=6
!!        else 
!!        omega_sim(id,t+1)=2!kappa_sim(id,t+1)=2
!!        z_sim(id,t+1)=6
!        endif



!        write(*,*) z_sim(id,t),z_sim(id,t+1),pz(z_sim(id,t),:),randomz(id,t)

        enddo
        enddo

print*,'simulations done'
!       write(*,*) b_sim(10,:)
!       pause


          open (UNIT=21, FILE="outputs\data_sim.txt", status = 'replace')


              do id=1,n_id
              do t=1,n_time
WRITE(21,'(I5,A,I5,A,F12.8,A,F12.8,A,F12.8,A,F12.8,A,F12.8,A,F12.8,A,F12.8,A,F12.8,A,F12.8,A,&
&F18.5,A,F18.5,A,F14.5,A,F12.8)') &
id,',',t,',',b_sim(id,t),',',ztvect(z_sim(id,t)),',',znvect(z_sim(id,t))&
,',',c_sim(id,t),',',ct_sim(id,t),',',cn_sim(id,t),',',yt_sim(id,t),',',yn_sim(id,t),',',pn_sim(id,t)&
,',',mu_sim(id,t),',',gvect(g_sim(id,t)),',',tax_sim(id,t),',',blim_sim(id,t)
              enddo
              enddo
         close(21)


        end program

!!!! Linear interpolation in three dimensions !!!!
double precision function intp2(bpr,mat)
implicit none
integer, parameter :: n_b=400,n_z=121,n_g=25,n_kappa=1,n_omega=1,n_id=2000,n_time=500
COMMON /PARAM/ bvect,kappavect,omegavect,r,cmin,bata,kappat,kappan,omega,sigma,eta,&
fixedt,fixedn,coefq,n_iter3,ztvect,znvect,gvect,pz,pg,pkappa,pomega,allow,gradb
DOUBLE PRECISION bvect(n_b),kappavect(n_kappa),omegavect(n_omega),r,cmin,bata,kappat,kappan,omega&
,sigma,eta,fixedt,fixedn,coefq,pz(n_z,n_z),pg(n_g,n_g),pkappa(n_kappa,n_kappa),pomega(n_omega,n_omega),&
ztvect(n_z),znvect(n_z),gvect(n_g),allow,gradb
INTEGER n_iter3
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

!write(*,*) 'in int',intp3
!pause

end function

!!! Output function for either of the sectors !!!
double precision function yfun(ll,tr,sect)
implicit none
integer, parameter :: n_b=400,n_z=121,n_g=25,n_kappa=1,n_omega=1,n_id=2000,n_time=500
COMMON /PARAM/ bvect,kappavect,omegavect,r,cmin,bata,kappat,kappan,omega,sigma,eta,&
fixedt,fixedn,coefq,n_iter3,ztvect,znvect,gvect,pz,pg,pkappa,pomega,allow,gradb
DOUBLE PRECISION bvect(n_b),kappavect(n_kappa),omegavect(n_omega),r,cmin,bata,kappat,kappan,omega&
,sigma,eta,fixedt,fixedn,coefq,pz(n_z,n_z),pg(n_g,n_g),pkappa(n_kappa,n_kappa),pomega(n_omega,n_omega),&
ztvect(n_z),znvect(n_z),gvect(n_g),allow,gradb
INTEGER n_iter3
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
integer, parameter :: n_b=400,n_z=121,n_g=25,n_kappa=1,n_omega=1,n_id=2000,n_time=500
COMMON /PARAM/ bvect,kappavect,omegavect,r,cmin,bata,kappat,kappan,omega,sigma,eta,&
fixedt,fixedn,coefq,n_iter3,ztvect,znvect,gvect,pz,pg,pkappa,pomega,allow,gradb
DOUBLE PRECISION bvect(n_b),kappavect(n_kappa),omegavect(n_omega),r,cmin,bata,kappat,kappan,omega&
,sigma,eta,fixedt,fixedn,coefq,pz(n_z,n_z),pg(n_g,n_g),pkappa(n_kappa,n_kappa),pomega(n_omega,n_omega),&
ztvect(n_z),znvect(n_z),gvect(n_g),allow,gradb
INTEGER n_iter3
SAVE /PARAM/

integer i,lt,tr,x,h

!rlim=r*exp(maxval(ztvect)-ztvect(lt))
Rfun=(1+r)!/ztvect(lt)**2.d0

!if (ztvect(lt)<=minval(ztvect)+1E-4) Rfun=(1+r)*1.25d0

!if (ztvect(lt)<1.d0) Rfun=(1+r)/ztvect(lt)!**2.d0

end function

!!!! In case we want to use an interest rate in the borrowing limit a function, we can just change this !!!!
double precision function Rlim(i,lt,tr,x,h)
implicit none
integer, parameter :: n_b=400,n_z=121,n_g=25,n_kappa=1,n_omega=1,n_id=2000,n_time=500
COMMON /PARAM/ bvect,kappavect,omegavect,r,cmin,bata,kappat,kappan,omega,sigma,eta,&
fixedt,fixedn,coefq,n_iter3,ztvect,znvect,gvect,pz,pg,pkappa,pomega,allow,gradb
DOUBLE PRECISION bvect(n_b),kappavect(n_kappa),omegavect(n_omega),r,cmin,bata,kappat,kappan,omega&
,sigma,eta,fixedt,fixedn,coefq,pz(n_z,n_z),pg(n_g,n_g),pkappa(n_kappa,n_kappa),pomega(n_omega,n_omega),&
ztvect(n_z),znvect(n_z),gvect(n_g),allow,gradb
INTEGER n_iter3
SAVE /PARAM/

integer i,lt,tr,x,h

!rlim=r*exp(maxval(ztvect)-ztvect(lt))
Rlim=1.d0!(1+r)/ztvect(lt)**2.d0
!Rlim=1.d0/gvect(tr)


end function

!!!! We can make the tightness parameter of the constraint a function with this !!!!
double precision function kappatfun(lt,x)
implicit none
integer, parameter :: n_b=400,n_z=121,n_g=25,n_kappa=1,n_omega=1,n_id=2000,n_time=500
COMMON /PARAM/ bvect,kappavect,omegavect,r,cmin,bata,kappat,kappan,omega,sigma,eta,&
fixedt,fixedn,coefq,n_iter3,ztvect,znvect,gvect,pz,pg,pkappa,pomega,allow,gradb
DOUBLE PRECISION bvect(n_b),kappavect(n_kappa),omegavect(n_omega),r,cmin,bata,kappat,kappan,omega&
,sigma,eta,fixedt,fixedn,coefq,pz(n_z,n_z),pg(n_g,n_g),pkappa(n_kappa,n_kappa),pomega(n_omega,n_omega),&
ztvect(n_z),znvect(n_z),gvect(n_g),allow,gradb
INTEGER n_iter3
SAVE /PARAM/
integer i,lt,tr,x,h

kappatfun=kappavect(x)!*ztvect(lt)**2.d0

end function


!!!! We can make the tightness parameter of the constraint a function with this !!!!
double precision function kappanfun(lt,x)
implicit none
integer, parameter :: n_b=400,n_z=121,n_g=25,n_kappa=1,n_omega=1,n_id=2000,n_time=500
COMMON /PARAM/ bvect,kappavect,omegavect,r,cmin,bata,kappat,kappan,omega,sigma,eta,&
fixedt,fixedn,coefq,n_iter3,ztvect,znvect,gvect,pz,pg,pkappa,pomega,allow,gradb
DOUBLE PRECISION bvect(n_b),kappavect(n_kappa),omegavect(n_omega),r,cmin,bata,kappat,kappan,omega&
,sigma,eta,fixedt,fixedn,coefq,pz(n_z,n_z),pg(n_g,n_g),pkappa(n_kappa,n_kappa),pomega(n_omega,n_omega),&
ztvect(n_z),znvect(n_z),gvect(n_g),allow,gradb
INTEGER n_iter3
SAVE /PARAM/

integer i,lt,tr,x,h

kappanfun=kappavect(x)!*ztvect(lt)**2.d0

end function
