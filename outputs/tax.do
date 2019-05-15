replace tax=0 if tax<0
corr growthgdp tax if t>100
gen rhogytax=r(rho)
corr debtoy tax if t>100
gen rhodebtoytax=r(rho)
corr loggdp tax if t>100
gen rhogdptax=r(rho)

sum tax if goodtimes==1, d
gen taxgood=r(mean)*100
sum tax if goodtimes==0, d
gen taxbad=r(mean)*100

sum tax if t>100, d
gen taxmean=r(mean)*100
gen taxmedian=r(p50)*100


sum tax if goodtimesdetr==1, d
gen taxgooddetr=r(mean)*100
sum tax if goodtimesdetr==0, d
gen taxbaddetr=r(mean)*100

sum tax if goodtimesintdetr==1, d
gen taxgoodintdetr=r(mean)*100
sum tax if goodtimesintdetr==0, d
gen taxbadintdetr=r(mean)*100

sum tax if goodtimesg==1, d
gen taxgoodg=r(mean)*100
sum tax if goodtimesg==0, d
gen taxbadg=r(mean)*100

//////

gen reversal=nxoy-L1.nxoy
gen outputg=gdp/L1.gdp*L1.g-1

sum tax if t>100 & ssnew, d
gen taxmeanssl0=r(mean)*100

sum tax if t>100 & F1.ssnew, d
gen taxmeanssl1=r(mean)*100

sum tax if t>100 & F2.ssnew, d
gen taxmeanssl2=r(mean)*100

sum tax if t>100 & F3.ssnew, d
gen taxmeanssl3=r(mean)*100

sum tax if t>100 & mupos, d
gen taxmeanmuposl0=r(mean)*100

sum tax if t>100 & F1.mupos, d
gen taxmeanmuposl1=r(mean)*100

sum tax if t>100 & F2.mupos, d
gen taxmeanmuposl2=r(mean)*100

sum tax if t>100 & F3.mupos, d
gen taxmeanmuposl3=r(mean)*100

gen logtax=log(tax)
gen intloggdpbaddetr=loggdp*(1-goodtimesdetr)
reg logtax loggdp intloggdpbaddetr if t>100
gen coefintloggdpbaddetr=_b[intloggdpbaddetr]

gen intloggdpgooddetr=loggdp*goodtimesdetr
reg logtax loggdp intloggdpgooddetr if t>100
gen coefintloggdpgooddetr=_b[intloggdpgooddetr]


reg logtax loggdp if t>100
gen coefloggdp=_b[loggdp]

reg logtax loggdp if t>100 & goodtimesdetr==1
gen coefloggdpgooddetr=_b[loggdp]

reg logtax loggdp if t>100 & goodtimesdetr==0
gen coefloggdpbaddetr=_b[loggdp]

noisily sum taxmean taxgooddetr taxbaddetr taxgood taxbad  taxgoodg taxbadg rhogytax  ///
taxmedian taxgoodintdetr taxbadintdetr rhodebtoytax rhogdptax coefintloggdpgooddetr coefintloggdpbaddetr coefloggdp coefloggdpgooddetr coefloggdpbaddetr, separator(0)
