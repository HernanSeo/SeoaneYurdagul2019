*cd "C:\Users\eyurdagu\Dropbox\project - overborr_missallocation\codes\alternative4\outputs"

*cd "\outputs"


insheet id t b zt zn c ct cn yt yn pn mu g tax blim using "data_sim.txt", comma clear 


destring tax, replace force


gen kappa=0.42
gen omega=0.31

tsset id t


gen eta=1/0.83-1

gen xr=(omega^(1/(1+eta))+(1-omega)^(1/(1+eta))*pn^(eta/(1+eta)))^((1+eta)/eta)

gen gdpint=yt+pn*yn
gen gdp=(yt+pn*yn)/xr
gen ctilde=(ct+pn*cn)/xr

gen boy=b/xr/gdp
gen boyt=b/xr/yt
gen shareyt=yt/xr/gdp
gen tboy=(yt-ct)/xr/gdp
gen ytoyn=yt/yn


sum t 
gen tup=r(max)
drop if t>=tup


gen debt=-b
gen debtoy=-boy

foreach x in debt debtoy zt zn c ct cn yt yn gdp xr pn g {
gen log`x'=log(`x')
}

foreach x in yt yn ct cn gdp c {
gen growth`x'=`x'/L1.`x'*L1.g-1
}

*gen ca=F1.b-b
gen ca=(g*F1.b-b)/xr
gen caoy=ca/gdp


gen blimhat=-kappa*(yt+pn*yn)/g

gen bbind=(F1.b<=blim+1E-10)
gen mupos=(mu>=1E-10)

gen bpr=F1.b

gen goodtimes=1 if t>100
sum growthgdp if t>100, d
replace goodtimes=0 if growthgdp<r(p50) & t>100


gen goodtimesdetr=1 if t>100
sum gdp if t>100, d
replace goodtimesdetr=0 if gdp<r(p50) & t>100

gen goodtimesintdetr=1 if t>100
sum gdpint if t>100, d
replace goodtimesintdetr=0 if gdpint<r(p50) & t>100

gen goodtimesg=1 if t>100
sum g if t>100, d
replace goodtimesg=0 if g<r(p50) & t>100


/*
sum mupos if goodtimes==1, d
gen muposgood=r(mean)*100
sum mupos if goodtimes==0, d
gen muposbad=r(mean)*100

sum mupos if goodtimesdetr==1, d
gen muposgooddetr=r(mean)*100
sum mupos if goodtimesdetr==0, d
gen muposbaddetr=r(mean)*100

sum mupos if goodtimesintdetr==1, d
gen muposgoodintdetr=r(mean)*100
sum mupos if goodtimesintdetr==0, d
gen muposbadintdetr=r(mean)*100

sum mupos if goodtimesg==1, d
gen muposgoodg=r(mean)*100
sum mupos if goodtimesg==0, d
gen muposbadg=r(mean)*100


noisily sum muposgood muposbad muposgoodintdetr muposbadintdetr muposgooddetr muposbaddetr muposgoodg muposbadg, separator(0)*/
****Comment *****

/**/

foreach x of varlist caoy tboy logpn logxr {
corr `x' L1.`x' if t>100
gen rho`x'=r(rho)
}

corr logpn ytoyn if t>100
gen rhopnytoyn=r(rho)
corr logxr ytoyn if t>100
gen rhoxrytoyn=r(rho)

corr growthyt growthyn if t>100
gen rhogytgyn=r(rho)
corr growthyt L1.growthyt if t>100
gen rhogyt=r(rho)
corr growthyn L1.growthyn if t>100
gen rhogyn=r(rho)
corr growthyn L1.growthyt if t>100
gen rhogynlgyt=r(rho)
*corr growthgdp L1.growthgdp if t>100
*gen rhoggdp=r(rho)

corr growthyt tboy if t>100
gen rhogyttboy=r(rho)
corr growthyt caoy if t>100
gen rhogytcaoy=r(rho)
corr growthyt pn if t>100
gen rhogytpn=r(rho)

foreach x of varlist caoy tboy logxr growthc {
corr `x' growthgdp if t>100
gen rhogy`x'=r(rho)
}

gen logyint=log(gdp*xr)
foreach x of varlist caoy tboy logxr logc {
corr `x' logyint if t>100
gen rhoyint`x'=r(rho)
}

gen blimdif=F1.b-blim


gen nx=yt-ct
gen nxoy=nx/xr/gdp

*** Define the Sudden Stop ***

*gen ssnew=         (nxoy-L1.nxoy>0.02 & gdp/L1.gdp*L1.g-1<-0.02 & t>100)
gen ssnew=         (nxoy-L1.nxoy>0.02 & gdp/L1.gdp*L1.g-1<-0.02)

gen ssnew2=ssnew
replace ssnew2=0 if ssnew==1 & ( L1.ssnew==1 | L2.ssnew==1 | L3.ssnew==1 | L4.ssnew==1 | L5.ssnew==1)

rename ssnew ss
rename ssnew2 ssnew


*** This part creates the ultimate statistics which include those in Table 3 and the model statistics mentioned in Section 5. ***

preserve

collapse (mean) ss ssnew rho* boy caoy tboy shareyt mupos* growthyt growthyn growthgdp ///
(sd) sdlogc=logc sdgyt=growthyt sdgyn=growthyn sdgct=growthct sdgcn=growthcn sdggdp=growthgdp sdgc=growthc sdcaoy=caoy sdtboy=tboy sdpn=logpn sdxr=logxr if t>100

noisily sum sdgyt rhogyt sdgyn rhogyn growthyt growthyn rhogytgyn rhogynlgyt sdgct sdgcn sdpn rhologpn rhologxr boy mupos ss ssnew caoy  rhocaoy tboy rhotboy shareyt ///
rhopnytoyn rhoxrytoyn rhogytcaoy rhogyttboy rhogytpn sdgc sdxr sdcaoy sdtboy rhogygrowthc rhogylogxr rhogycaoy rhogytboy sdlogc rhoyintlogc rhoyintlogxr rhoyintcaoy rhoyinttboy ///
, separator(0)

restore

**** IF not Solving Social Planner's problem, comment out the line below, which creates the tax numbers in Table 4 ****
run "tax.do"

