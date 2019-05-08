clear all;  clc;
warning('off')
j=1;j1=1;j2=1;j3=1;j4=1;

country    ={'Argentina','Bolivia' ,'Brazil'   ,'Chile'   ,'Colombia','CostaRica','Ecuador' ,'Egypt'    ,'Guatemala','Honduras','India'   ,'Indonesia','Korea'   ,'Malaysia','Mexico'  ,'Panama'  ,'Peru'    ,'Philipines', 'SouthAfrica', 'Uruguay'  };
country_num=[ 1         ; 2        ; 3         ; 4        ; 5        ; 6         ; 7        ; 8         ; 9         ; 10       ; 11       ; 12        ; 13       ; 14       ; 15       ; 16       ; 17       ; 18         ; 19            ; 20        ];

def_usa=xlsread('t_nt_crosscountry.xlsx','USA','e18:bi18'); % 1960-2016
tt=1960:2016;

for country_num=1:20
    data = xlsread('t_nt_crosscountry.xlsx',country{country_num},'e2:bi28');
    expr=['agric=data(1,:);'];  eval(expr);
    expr=['expo=data(3,:);'];  eval(expr);
    expr=['gdp=data(7,:);'];  eval(expr);
    expr=['impo=data(10,:);']; eval(expr);
    expr=['manuf=data(14,:);']; eval(expr);
    expr=['gdp_dol=data(22,:);']; eval(expr);
    expr=['def=data(23,:);']; eval(expr);
    expr=['ner=data(24,:);']; eval(expr);
    if country_num==4|| country_num==13|| country_num==20
        debt=gdp*NaN;
    else
        expr=['debt=data(25,:);']; eval(expr);
    end
    if isempty(ner)~=1
        rer=ner.*def_usa./def;
    end
    nxy=(expo-impo)./gdp;
    
    
    T=agric+manuf;
    NT=gdp-T;
    
    % growth rates, we lose the first observation
    gY=diff(log(gdp)); gY=[gY(1) gY];
    gT=diff(log(T)); gT=[gT(1) gT];
    gNT=diff(log(NT)); gNT=[gNT(1) gNT];
    gNX=diff(nxy); gNX=[gNX(1) gNX];
    
    % identify SS
    ssepi=find(gNX>.02 & gY<-.02);
    
    % Pìck the SS dates
    ss_t=tt(ssepi);
    
    hini=5;
    hend=1;
    wdw=5;
    for m=1:length(ssepi)
        i=ssepi(m);
        if m==1 || ( m>1 && ssepi(m)-ssepi(m-1)>5)
            if i-wdw>=1 && i+5<=length(gT)
                YGTaux=gT(i-wdw:i+5);
                YGNaux=gNT(i-wdw:i+5);
                YGaux=gY(i-wdw:i+5);
                if isempty(rer)~=1
                    RERaux=rer(i-wdw:i+5)./rer(i);
                end
                NXYaux=nxy(i-wdw:i+5);
                if i+5<=length(debt)
                    DEBTYaux=debt(i-wdw:i+5)./gdp_dol(i-wdw:i+5);
                end
 
                YG (:,j)=YGaux;
                NXY(:,j)=NXYaux;
                j=j+1;
                
                if any(isnan(YGTaux))~=1 && any(isnan(YGNaux))~=1 
                    YGT(:,j1)=YGTaux;
                    YGN(:,j1)=YGNaux;
                    j1=j1+1;
                end
                
                if any(isnan(RERaux))~=1
                    RER(:,j2)=RERaux;
                    j2=j2+1;
                end           
                
                if any(isnan(DEBTYaux))~=1
                    DEBTY(:,j3)=DEBTYaux;
                    j3=j3+1;
                end


            end
            
            
        end
    end
end



range=-wdw:1:5;
figure(1)
subplot(3,2,1);plot(-wdw:5,mean(YGT,2)*100,'LineWidth',2);title('Tradable output growth');axis('tight');grid('on');;xticks([range])
subplot(3,2,2);plot(-wdw:5,mean(YGN,2)*100,'LineWidth',2);title('Non-tradable output growth');axis('tight');grid('on');;xticks([range])
subplot(3,2,3);plot(-wdw:5,mean(YG,2)*100,'LineWidth',2);title('Output growth');axis('tight');grid('on');;xticks([range])
subplot(3,2,4);plot(-wdw:5,mean(1./RER,2),'LineWidth',2);title('Real exchange rate');axis('tight');grid('on');;xticks([range])
subplot(3,2,5);plot(-wdw:5,mean(NXY,2)*100,'LineWidth',2);title('Net exports to output ratio');axis('tight');grid('on');;xticks([range])
subplot(3,2,6);plot(-wdw:5,mean(DEBTY,2),'LineWidth',2);title('Debt to output ratio');axis('tight');grid('on');;xticks([range])
