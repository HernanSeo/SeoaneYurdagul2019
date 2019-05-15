clear all
nid=2000;
nt=500;
eta=1/0.83-1;

ss_criteria=0; % baseline in the paper

for m=1:2
    clear t    
    if m==1
        q={'id' 't' 'b' 'zt' 'zn' 'c' 'ct' 'cn' 'yt' 'yn' 'pn' 'mu' 'g' 'tax' 'blim'}; 
        kappa=0.36;
        omega=0.31;
        polf=csvread('data_sim.txt'); % ours
    elseif m==2
        clear kappa omega
        kappa=0.32;
        omega=0.31;
        nid=2000;
        nt=500;
        q={'id' 't' 'b' 'zt' 'zn' 'c' 'ct' 'cn' 'yt' 'yn' 'pn' 'mu' 'g' 'tax' 'blim'}; 
        polf=csvread('data_sim_section3.txt');
    end
    
    r = 0.04;
    
    ygobs=NaN(500,nid);
    y=NaN(500,nid);
    yobs=NaN(500,nid);
    
    for j=1:nid
        for i=1:size(polf,2)
            expr=[char(q(i)) '(:,j)=polf(nt*(j-1)+1:nt*j,i);'];
            eval(expr);
        end
        y(:,j)=yt(:,j)+pn(:,j).*yn(:,j);
        rer(:,j) = (omega.^(1/(1+eta))+(1-omega).^(1/(1+eta)).*pn(:,j).^(eta/(1+eta))).^((1+eta)/eta); % increase appreciation
        
        yobs(:,j)=y(:,j)./rer(:,j); %rewrite output in terms of consumption package
        ygobsaux=diff(log(yobs(:,j)));
        ygobs(:,j)=[ygobsaux(1)+log(g(1,j)); ygobsaux+log(g(1:end-1,j))];
        
        ygtaux=diff(log(yt(:,j)));
        ygt(:,j)=[ygtaux(1)+log(g(1,j)); ygtaux+log(g(1:end-1,j))];
        
        ygnaux=diff(log(yn(:,j)));
        ygn(:,j)=[ygnaux(1)+log(g(1,j)); ygnaux+log(g(1:end-1,j))];   
        
        cgaux=diff(log(c(:,j)));
        cg(:,j)=[cgaux(1)+log(g(1,j)); cgaux+log(g(1:end-1,j))];
        
        ztgaux=diff(log(zt(:,j)));
        ztg(:,j)=[ztgaux(1); ztgaux];
        
        zngaux=diff(log(zn(:,j)));
        zng(:,j)=[zngaux(1); zngaux];
        
        rergnaux=diff(log(rer(:,j)));
        rerg(:,j)=[rergnaux(1); rergnaux];
    end
    
    
    clear expr i j
    nxy = (yt-ct)./y;
    
    borrlim= kappa.*yt+kappa.*pn.*yn;
    
    
    by=b(2:end,:).*g(1:end-1,:)./y(1:end-1,:); by=[by; by(end,:)];
    rby=r*b./y;
    cay=nxy+rby;
    reversal = diff(nxy);
    reversal = [reversal(1,:); reversal];
    
    tini=-5;
    tend=5;
    BY    = NaN(110000,-tini+tend+1);
    RBY   = NaN(110000,-tini+tend+1);
    NXY   = NaN(110000,-tini+tend+1);
    CAY   = NaN(110000,-tini+tend+1);
    PN    = NaN(110000,-tini+tend+1);
    RER   = NaN(110000,-tini+tend+1);
    YGOBS = NaN(110000,-tini+tend+1);
    YGT   = NaN(110000,-tini+tend+1);
    YGN   = NaN(110000,-tini+tend+1);
    YOBS  = NaN(110000,-tini+tend+1);
    CT    = NaN(110000,-tini+tend+1);
    YT    = NaN(110000,-tini+tend+1);
    B     = NaN(110000,-tini+tend+1);
    MU    = NaN(110000,-tini+tend+1);
    CG    = NaN(110000,-tini+tend+1);
    ZTG   = NaN(110000,-tini+tend+1);
    ZNG   = NaN(110000,-tini+tend+1);
    GG   = NaN(110000,-tini+tend+1);
    
    itv=5;
    jj=0;
    event =NaN(110000,1);
    AA=0;
    for i=1:nid
        if m==1
            aa=find(reversal(:,i)>0.02 & ygobs(:,i)<-0.02);
        else
            aa=find(reversal(:,i)>0.02 & ygobs(:,i)<-0.02); 
        end
        AA=AA+size(aa,1);
        
        if isempty(aa)~=1
            for j=1:size(aa,1)
                if j==1
                    if aa(j)>100 && aa(j)<500-itv && aa(j)>4
                        jj = jj+1;
                        tss = aa(j);                        
                        event(jj)  = 1;
                        BY(jj,:)   = by(tss+tini:tss+tend,i);
                        RBY(jj,:)  = rby(tss+tini:tss+tend,i);
                        PN(jj,:)   = pn(tss+tini:tss+tend,i);
                        RER(jj,:)  = rer(tss+tini:tss+tend,i);
                        NXY(jj,:)  = nxy(tss+tini:tss+tend,i);
                        CAY(jj,:)  = cay(tss+tini:tss+tend,i);
                        YGOBS(jj,:)= ygobs(tss+tini:tss+tend,i);
                        YGT(jj,:)  = ygt(tss+tini:tss+tend,i);
                        YGN(jj,:)  = ygn(tss+tini:tss+tend,i);
                        YOBS(jj,:) = yobs(tss+tini:tss+tend,i);
                        CT(jj,:)   = ct(tss+tini:tss+tend,i);
                        YT(jj,:)   = yt(tss+tini:tss+tend,i);
                        B(jj,:)    = b(tss+tini:tss+tend,i);
                        MU(jj,:)   = mu(tss+tini:tss+tend,i);
                        CG(jj,:)   = cg(tss+tini:tss+tend,i);
                        ZTG(jj,:)  = ztg(tss+tini:tss+tend,i);
                        ZNG(jj,:)  = zng(tss+tini:tss+tend,i);
                        GG(jj,:)   = log(g(tss+tini:tss+tend,i));
                    end
                else
                    if aa(j)>100 && aa(j)<500-itv && aa(j)-aa(j-1)>5
                        jj = jj+1;
                        tss = aa(j);
                        event(jj)  = 1;
                        BY(jj,:)   = by(tss+tini:tss+tend,i);
                        RBY(jj,:)  = rby(tss+tini:tss+tend,i);
                        PN(jj,:)   = pn(tss+tini:tss+tend,i);
                        RER(jj,:)  = rer(tss+tini:tss+tend,i);
                        NXY(jj,:)  = nxy(tss+tini:tss+tend,i);
                        CAY(jj,:)  = cay(tss+tini:tss+tend,i);
                        YGOBS(jj,:)= ygobs(tss+tini:tss+tend,i);
                        YGT(jj,:)  = ygt(tss+tini:tss+tend,i);
                        YGN(jj,:)  = ygn(tss+tini:tss+tend,i);
                        YOBS(jj,:) = yobs(tss+tini:tss+tend,i);
                        CT(jj,:)   = ct(tss+tini:tss+tend,i);
                        YT(jj,:)   = yt(tss+tini:tss+tend,i);
                        B(jj,:)    = b(tss+tini:tss+tend,i);
                        MU(jj,:)   = mu(tss+tini:tss+tend,i);
                        CG(jj,:)   = cg(tss+tini:tss+tend,i);
                        ZTG(jj,:)  = ztg(tss+tini:tss+tend,i);
                        ZNG(jj,:)  = zng(tss+tini:tss+tend,i);
                        GG(jj,:)   = log(g(tss+tini:tss+tend,i));
                    end
                end
            end
        end
    end
    disp([m AA])
    if jj<110000
        event=event(1:jj,:);
        BY=BY(1:jj,:);
        MU=MU(1:jj,:);
        CAY=CAY(1:jj,:);
        RBY=RBY(1:jj,:);
        NXY=NXY(1:jj,:);
        YGOBS=YGOBS(1:jj,:);
        YGT=YGT(1:jj,:);
        YGN=YGN(1:jj,:);
        YOBS=YOBS(1:jj,:);
        PN=PN(1:jj,:);
        RER=RER(1:jj,:);
        CT=CT(1:jj,:);
        CG=CG(1:jj,:);
        YT=YT(1:jj,:);
        B=B(1:jj,:);
        ZTG=ZTG(1:jj,:);
        ZNG=ZNG(1:jj,:);
        GG=GG(1:jj,:);        
    end
    disp('number of episodes')
    disp(jj)
    disp('frequency in %')
    disp(jj/(2000*(500-(itv+1)-101+1))*100) % We discard any ss in the first 100 and in the last 5 periods, because we need to construct -5:5 intervals
    events{m} = event;
    rbyss{m}  = mean(RBY,1);
    cayss{m}  = mean(CAY,1);
    byss{m}   = mean(BY,1);
    nxyss{m}  = mean(NXY,1)*100;
    yobsss{m} = mean(YOBS,1);
    ygobsss{m}= mean(YGOBS,1)*100;
    ygtss{m}  = mean(YGT,1)*100;
    ygnss{m}  = mean(YGN,1)*100;
    pnss{m}   = mean(PN,1);
    rerss{m}  = mean(RER,1);
    ctss{m}   = mean(CT,1);
    cgss{m}   = mean(CG,1); %final consumption growth
    ytss{m}   = mean(YT,1);
    bss{m}    = mean(B,1);
    muss{m}    = mean(MU,1);
    
    zngss{m}=mean(ZNG,1)*100;
    ztgss{m}=mean(ZTG,1)'*100;
    ggss{m}=mean(GG,1)'*100;
    
    clear BY MU CAY RBY NXY YGOBS YGT YGN YOBD PN RER CT CG YT B ZTG ZTN GG
end

ttt=tini:tend;
load('ss_data_arg.mat')
    % SS episodes used in the file
    %     1.9140    1.9140   -0.0444    0.0946   -0.0114       NaN    9.0000   13.0000
    %     1.9310    1.9310    0.0216    0.0097   -0.0059       NaN   26.0000   30.0000
    %     1.9590    1.9590    0.0425   -0.0440    0.0141    0.0321   54.0000   58.0000
    %     1.9820    1.9820    0.0488   -0.0359   -0.0169   -0.0113   77.0000   81.0000
    %     1.9950    1.9950   -0.0221    0.0116    0.0199    0.0582   90.0000   94.0000


%  1  2  3  4  5 6 7 8 9 10 11
% -5 -4 -3 -2 -1 0 1 2 3 4 5
range=6+tini:11;
nxyss{3}=nanmean(NXY(range,:),2)';
ygobsss{3}=nanmean(YG(range,:),2)'*100;
ygtss{3}=nanmean(YGT(range,:),2)'*100;
ygnss{3}=nanmean(YGN(range,:),2)'*100;
rerss{3}=nanmean(1./RER(range,:),2)';
cayss{3}=nanmean(CAY(range,:),2)';
cgss{3}=nanmean(CG(range,3:end),2)';% first 2 SS we dont have data for consumption


zngss{3}=nanmean(ZNG(range,:),2)'*100;
ztgss{3}=nanmean(ZTG(range,:),2)'*100;
ggss{3}=nanmean(G(range,:),2)'*100;


% keyboard
rerss0data=rerss{3};
rerss0data=rerss0data(6); % normalize to 1 in the SS
range=tini:1:tend;
for m=1:2
    figure(m)
    rerss0=rerss{m};
    rerss0=rerss0(6);% normalize to 1 in the SS
    subplot(3,2,4);plot(ttt,nxyss{m},'LineWidth',2);hold on; plot(ttt,nxyss{3},'-.','LineWidth',2);hold off; title('Net exports to output ratio');grid on ;axis('tight');xticks([range])
    subplot(3,2,1);plot(ttt,ygobsss{m},'LineWidth',2);hold on; plot(ttt,ygobsss{3},'-.','LineWidth',2);hold off;title('Output growth');grid on ;axis('tight');xticks([range])
    subplot(3,2,3);plot(ttt,rerss{m}/rerss0,'LineWidth',2);hold on; plot(ttt,rerss{3}/rerss0data,'-.','LineWidth',2);hold off;title('Real exchange rate');grid on ;axis('tight');xticks([range])
    subplot(3,2,2);plot(ttt,cgss{m}*100,'LineWidth',2);hold on; plot(ttt,cgss{3}*100,'-.','LineWidth',2);hold off;title('Consumption growth');grid on ;axis('tight');xticks([range])
    subplot(3,2,5);plot(ttt,ygtss{m},'LineWidth',2);hold on; plot(ttt,ygtss{3},'-.','LineWidth',2);hold off;title('Tradable output growth');grid on ;xticks([range]);axis('tight')
    subplot(3,2,6);plot(ttt,ygnss{m},'LineWidth',2);hold on; plot(ttt,ygnss{3},'-.','LineWidth',2);hold off;title('Non-tradable output growth');grid on ;xticks([range]);axis('tight')
    if m==1
        legend('Trend shocks model','Data');legend boxoff
    elseif m==2
        legend('Transitory shocks model','Data');legend boxoff
    end
    figure(10*m)
    subplot(3,1,1);plot(ttt,ztgss{m},'LineWidth',2);hold on; plot(ttt,ztgss{3},'-.','LineWidth',2);hold off; title('z_T growth');grid on ;axis('tight');xticks([range])
    subplot(3,1,2);plot(ttt,zngss{m},'LineWidth',2);hold on; plot(ttt,zngss{3},'-.','LineWidth',2);hold off;title('z_N growth');grid on ;axis('tight');xticks([range])
    subplot(3,1,3);plot(ttt,ggss{m},'LineWidth',2);hold on; plot(ttt,1+ggss{3},'-.','LineWidth',2);hold off;title('Trend growth');grid on ;axis('tight');xticks([range])
    if m==1
        legend('Trend shocks model','Smoothed estimates of unobserved processes');legend boxoff
    elseif m==2
        legend('Transitory shocks model','Smoothed estimates of unobserved processes');legend boxoff
    end
    
end
