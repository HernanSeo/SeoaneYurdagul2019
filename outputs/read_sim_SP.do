/*
*cd "C:\Users\yurdagul\Dropbox\Overborrowing\endowment\7_calib_ar_mg300_sg070_rho040_ng11_nz1_kappa100_r6"
cd "C:\Users\eyurdagu\Dropbox\Overborrowing\endowment\8_calib_baseline"
insheet id t b kt kn zt zn c ct cn yt yn pn mu kappa omega g using "data_sim.txt", comma clear 
save "data_sim.dta", replace

cd "C:\Users\eyurdagu\Dropbox\Overborrowing\endowment\8_calib_baseline\SP"
insheet id t b kt kn zt zn c ct cn yt yn pn mu kappa omega g using "data_sim.txt", comma clear 
save "data_sim.dta", replace

gen model="BchiSP"

append using "C:\Users\eyurdagu\Dropbox\Overborrowing\endowment\8_calib_baseline\data_sim.dta"
replace model="Bchi" if model==""
*/
*cd "C:\Users\eyurdagu\Dropbox\Overborrowing\endowment\20_calib_ar_40012125_mg100_sg0353_rho5499_rhot7501rhon7963sgt0532sgn0495_kappa36_r4_beta91_zgupdown3"
*cd "C:\Users\eyurdagu\Dropbox\Overborrowing\endowment\20_calib_ar_40012125_mg100_sg0353_rho5499_rhot7501rhon7963sgt0532sgn0495_kappa36_r4_beta91_zgupdown3\nz64"
*cd "C:\Users\eyurdagu\Dropbox\Overborrowing\endowment\20_calib_ar_40012125_mg100_sg0353_rho5499_rhot7501rhon7963sgt0532sgn0495_kappa36_r4_beta91_zgupdown3\ng13"

////

//// First read the Decentralized Equilibrium Simulations //////
insheet id t b zt zn c ct cn yt yn pn mu g tax blim using "data_sim.txt", comma clear 
gen model="BM"
save "data_sim.dta", replace

//// Then read the Social Planner Simulations //////

*cd "C:\Users\yurdagul\Dropbox\Overborrowing\endowment\10_calib_ar_nb100_mg300_sg0289_rho6687_ng11_kappa38_r9_bup05\SP"
insheet id t b zt zn c ct cn yt yn pn mu g tax blim using "SP\data_sim.txt", comma clear 
save "SP\data_sim.dta", replace

gen model="BMSP"
append using "data_sim.dta", force


egen idgp=group(id model)
tsset idgp t

gen gdp=yt+pn*yn
gen boy=b/gdp
gen boyt=b/yt
gen shareyt=yt/gdp
gen tboy=(yt-ct)/gdp
gen ytoyn=yt/yn

gen kappa=0.42
gen omega=0.31
gen eta=1/0.83-1

gen sigma=2
gen utility=(c^(1-sigma)-1)/(1-sigma)


sum t 
gen tup=r(max)
drop if t>=tup

gen xr=(omega^(1/(1+eta))+(1-omega)^(1/(1+eta))*pn^(eta/(1+eta)))^((1+eta)/eta)

gen debt=-b
gen debtoy=-boy

foreach x in debt debtoy zt zn c ct cn yt yn gdp {
gen log`x'=log(`x')
}

foreach x in yt yn ct cn gdp c {
gen growth`x'=`x'/L1.`x'*L1.g-1
}

gen ca=F1.b-b
gen caoy=ca/gdp


gen bbind=(F1.b<=blim+1E-10)
gen mupos=(mu>=1E-10)

gen bpr=F1.b

gen bg=b*L1.g
gen bog=b/L1.g

**** This creates the histogram in Figure 6 ****

twoway (hist bg if model=="BM", bin(20) lcolor(red) fcolor(red) fintensity(20)) ///
(hist bg if model=="BMSP", bin(20) lcolor(blue) lpattern(dash) fcolor(none) fintensity(40)) if t>100 & bg<-0.6 & bg>-1.2, /// & bg>-1.1, ///
xtitle("Assets") ytitle("Density")  ylabel(,nogrid) legend(order(1 "Decentralized" 2 "Social Planner") rows(1)  region(lcolor(white)) size(medium)) graphregion(fcolor(white) lcolor(white)) plotregion(margin(large))
graph export "SP\histb.pdf", replace

/*
sum bg if model=="BM", d
*sum bg if model=="Bchi", d

gen bgminBM=r(min)
gen bgp10BM=r(p10)
gen bgp25BM=r(p25)
gen bgp50BM=r(p50)
gen bgmeanBM=r(mean)

sum utility if model=="BM", d
gen up50BM=r(p50)
gen umeanBM=r(mean)

sum bg if model=="BMSP", d
*sum bg if model=="BchiSP", d

gen bgminSP=r(min) 
gen bgp10SP=r(p10)
gen bgp25SP=r(p25)
gen bgp50SP=r(p50)
gen bgmeanSP=r(mean)

sum utility if model=="BMSP", d
gen up50SP=r(p50)
gen umeanSP=r(mean)

gen bgminOB=(bgminBM-bgminSP)/bgminSP*100
gen bgp10OB=(bgp10BM-bgp10SP)/bgp10SP*100
gen bgp25OB=(bgp25BM-bgp25SP)/bgp25SP*100
gen bgp50OB=(bgp50BM-bgp50SP)/bgp50SP*100
gen bgmeanOB=(bgmeanBM-bgmeanSP)/bgmeanSP*100

gen up50OB=(up50BM-up50SP)/up50SP*100
gen umeanOB=(umeanBM-umeanSP)/umeanSP*100

*/




/*
twoway (hist bg if model=="BM", bin(50) lcolor(red) fcolor(red) fintensity(20)) ///
(hist bg if model=="BMSP", bin(50) lcolor(blue) lpattern(dash) fcolor(none) fintensity(40)) if t>100 & bg<-0.4, /// & bg>-1.1, ///
xtitle("Assets") ytitle("Density")  ylabel(,nogrid) legend(order(1 "Decentralized" 2 "Social Planner") rows(1)  region(lcolor(white)) size(medium)) graphregion(fcolor(white) lcolor(white)) plotregion(margin(large))
graph export "SP\histb_50bin.pdf", replace
/*
