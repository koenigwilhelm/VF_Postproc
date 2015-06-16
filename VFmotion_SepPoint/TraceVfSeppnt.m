close all
clear all

%% User option flags: (1: yes; 0: no)
flagSaveVFanim=0;
flagNormalizeTime=1;
flagNormalizeForceX=0;
flagNormalizeForceY=0;
flagNormalizeEnergy=0;
flagNormalizeBernoulli=0;
% start and end times for VF animation
tsttAnim=531;
tendAnim=1014;
% start and end times for VF SVD
tsttSVD=6000;
tendSVD=9000;
% start and end times for interested VF vibration cycle
tsttVFcycle=.1773;
tendVFcycle=.1816;

%% get unique datasets (restart of simulation often leads to duplicate time steps)
seppnt=importdata('poa_seppnt.txt',' ');
sbctrc=importdata('poa_sbctrace.txt',' ');

tt_orig=sbctrc(:,1);
xsbc_orig=sbctrc(:,2:2:end-1);
ysbc_orig=sbctrc(:,3:2:end);
[tt,iforw,iback]=unique(tt_orig);
xsbc=xsbc_orig(iforw,:);
ysbc=ysbc_orig(iforw,:);
xsep_orig=seppnt(:,2);
ysep_orig=seppnt(:,3);
Qsep_orig=seppnt(:,4);
xsep=xsep_orig(iforw);
ysep=ysep_orig(iforw);
Qsep=Qsep_orig(iforw);

%% re-order solid boundary nodes so that it goes clockwise
xsbc=xsbc(:,(ysbc_orig(1,:)~=0)|(xsbc_orig(1,:)==min(xsbc_orig(1,:)))|(xsbc_orig(1,:)==max(xsbc_orig(1,:))));
ysbc=ysbc(:,(ysbc_orig(1,:)~=0)|(xsbc_orig(1,:)==min(xsbc_orig(1,:)))|(xsbc_orig(1,:)==max(xsbc_orig(1,:))));
nnvfsurf=length(xsbc(1,:));
nt=length(tt);
Arsort=[xsbc' ysbc'];
Brsort=sortrows(Arsort,[1 -(nt+1)]);
xsbc=Brsort(:,1:nt)';
ysbc=Brsort(:,nt+1:end)';
clear Arsort Brsort seppnt sbctrc xsbc_orig ysbc_orig xsep_orig ysep_orig
clear tt_orig

%% parameters for figures
nstep=7;
xshift=min(xsbc(1,:));
xlim_min=min(xsbc(1,:))-.1;
xlim_max=max(xsbc(1,:))+.3;

%% save animation
if (flagSaveVFanim==1)
    figure()
    filename='AnimVFseppntMov.gif';
    flagAppend=0;
    % for iit=1:nt
    for iit=tsttAnim:tendAnim
        clf
        hold on
        plot(xsbc(iit,:)-xshift,ysbc(iit,:))
        for ix=1:nstep:nnvfsurf
            plot(xsbc(tsttAnim:iit,ix)-xshift,ysbc(tsttAnim:iit,ix),'m')
        end
        plot(xsep(iit-10:iit)-xshift,ysep(iit-10:iit),'k')
        plot(xsep(iit)-xshift,ysep(iit),'-o')
        xlabel('x (cm)')
        ylabel('y (cm)')
        axis([xlim_min-xshift xlim_max-xshift 0 1.4])
        hold off
        %pause(.05)
        drawnow
        frame=getframe(1);
        im=frame2im(frame);
        [A,map]=rgb2ind(im,256);
        if flagAppend==0
            imwrite(A,map,filename,'gif','LoopCount',Inf,'DelayTime',.02);
            flagAppend=1;
        else
            imwrite(A,map,filename,'gif','WriteMode','Append','DelayTime',.02);
        end
    end
end

%% plot results
figure()
hold on
for ix=1:nstep:nnvfsurf
    plot(xsbc(:,ix)-xshift,ysbc(:,ix),'m')
end
plot(xsep-xshift,ysep,'k')
plot(mean(xsbc,1)-xshift,mean(ysbc,1))
set(gca,'FontSize',14)
xlabel('x (cm)')
ylabel('y (cm)')
axis([xlim_min-xshift xlim_max-xshift 0 1.4])
hold off

%% SVD of VF motion
Xsbc=[xsbc(tsttSVD:tendSVD,:) ysbc(tsttSVD:tendSVD,:)];
Xsbcmean=mean(Xsbc);
Xsbc=Xsbc-repmat(Xsbcmean,tendSVD-tsttSVD+1,1);
[Usbc,Ssbc,Vsbc]=svd(Xsbc,'econ');
figure()
hold on
for ix=1:nstep:nnvfsurf
    plot(xsbc(tsttSVD:tendSVD,ix)-xshift,ysbc(tsttSVD:tendSVD,ix),'m')
end
plot(mean(xsbc,1)-xshift,mean(ysbc,1))
plot(mean(xsbc,1)+Vsbc(1:nnvfsurf,1)'*Ssbc(1,1)*min(Usbc(:,1))-xshift,mean(ysbc,1)+Vsbc(nnvfsurf+1:end,1)'*Ssbc(1,1)*min(Usbc(:,1)),'r')
plot(mean(xsbc,1)+Vsbc(1:nnvfsurf,1)'*Ssbc(1,1)*max(Usbc(:,1))-xshift,mean(ysbc,1)+Vsbc(nnvfsurf+1:end,1)'*Ssbc(1,1)*max(Usbc(:,1)),'g')
set(gca,'FontSize',14)
axis([xlim_min-xshift xlim_max-xshift 0 1.4])
hold off
clear Xsbc Usbc Ssbc Vsbc

%% VF vibration frequency
ymaxVSt=max(ysbc,[],2);
hglottis=1.397-ymaxVSt;
hglottisVar=hglottis-mean(hglottis);
hglottisReal=(hglottis-0.0114).*(hglottis-0.0114>0);
Nfft=2^(nextpow2(nt)+6);
Yhg=fft(hglottisVar,Nfft)/nt;
freq=linspace(0,1,Nfft/2+1)/(2*(tt(2)-tt(1)));
figure()
subplot(2,1,1)
plot(freq,abs(Yhg(1:Nfft/2+1))/max(abs(Yhg(1:Nfft/2+1))),'LineWidth',2);set(gca,'FontSize',14)
axis([0 500 0 1.1])
xlabel('frequency (Hz)')
ylabel('relative amplitude')
subplot(2,1,2)
plot(tt,hglottis,'LineWidth',2);set(gca,'FontSize',14)
xlabel('t (s)')
ylabel('h (cm)')

%% Bernoulli's and Control Volume Analyses
poamomentumX=importdata('poa_momentumX.txt',' ');
poamomentumY=importdata('poa_momentumY.txt',' ');
poaenergy=importdata('poa_energy.txt',' ');
poaberups=importdata('poa_bernoulliUps.txt',' ');
poaberdwn=importdata('poa_bernoulliDwn.txt',' ');
poamassrt=importdata('poa_massCons.txt',' ');

rateMX=poamomentumX(iforw,2)+poamomentumX(iforw,3);
rateMflowX=poamomentumX(iforw,8);
rateMdensX=poamomentumX(iforw,9);
fFluidX=poamomentumX(iforw,10);
drivPX=poamomentumX(iforw,4);
dragFX=poamomentumX(iforw,5);
resdMX=poamomentumX(iforw,6);
dragF2X=poamomentumX(iforw,7);

rateMY=poamomentumY(iforw,2)+poamomentumY(iforw,3);
rateMflowY=poamomentumY(iforw,8);
rateMdensY=poamomentumY(iforw,9);
fFluidY=poamomentumY(iforw,10);
drivPY=poamomentumY(iforw,4);
dragFY=poamomentumY(iforw,5);
resdMY=poamomentumY(iforw,6);
dragF2Y=poamomentumY(iforw,7);

rateKE=poaenergy(iforw,2)+poaenergy(iforw,3);
powrPp=poaenergy(iforw,4);
powrFr=poaenergy(iforw,5);
powrDs=poaenergy(iforw,6);
resdEn=poaenergy(iforw,7);
volint = poaenergy(iforw,8);
rateKEflow=poaenergy(iforw,9);
rateKEdens=poaenergy(iforw,10);
powrFluid=poaenergy(iforw,11);
powrCompr=poaenergy(iforw,12);

massInt=poamassrt(iforw,2);
massOut=poamassrt(iforw,3);
massRat=poamassrt(iforw,4);

if (flagNormalizeTime==1)
    tx=(tt-tsttVFcycle)/(tendVFcycle-tsttVFcycle);
else
    tx=tt;
end
if (flagNormalizeForceX==1)
    DenomForceX=max(abs(drivPX.*(tt>=tsttVFcycle).*(tt<=tendVFcycle)));
else
    DenomForceX=1;
end
if (flagNormalizeForceY==1)
    DenomForceY=max(abs(drivPY.*(tt>=tsttVFcycle).*(tt<=tendVFcycle)));
else
    DenomForceY=1;
end
if (flagNormalizeEnergy==1)
    DenomEnergy=max(abs((powrPp+powrFluid+powrCompr).*(tt>=tsttVFcycle).*(tt<=tendVFcycle)));
else
    DenomEnergy=1;
end
if (flagNormalizeBernoulli==1)
    DenomBernoulliDwn=max(abs(poaberdwn(iforw,2).*(tt>=tsttVFcycle).*(tt<=tendVFcycle)));
    DenomBernoulliUps=max(abs(poaberups(iforw,2).*(tt>=tsttVFcycle).*(tt<=tendVFcycle)));
else
    DenomBernoulliDwn=1;
    DenomBernoulliUps=1;
end

figure()
plot(tx,-drivPX/DenomForceX,'g',...
    tx,dragF2X/DenomForceX,'c',...
    tx,-rateMX/DenomForceX,'b',...
    tx,(rateMX+drivPX-dragF2X)/DenomForceX,'-.k',...
    'LineWidth',2);set(gca,'FontSize',14)
if (flagNormalizeTime==1)
    xlabel('t/T')
    xlim([0 1])
else
    xlabel('t (s)')
end
if (flagNormalizeForceX==1)
    ylabel('Normalized Force Terms')
else
    ylabel('Force Terms (g\cdotcm/s)')
end
hlg=legend('$F_{\Delta p}+F_{VF}$','$F_f$','$-\dot{P}_{CV}$','residual');
set(hlg,'interpreter','latex','location','eastoutside')

figure()
plot(tx,fFluidX/DenomForceX,'g',...
    tx,(-drivPX-fFluidX)/DenomForceX,'r',...
    tx,dragF2X/DenomForceX,'c',...
    tx,-(rateMX+rateMflowX+rateMdensX)/DenomForceX,'b',...
    tx,rateMflowX/DenomForceX,'--b',...
    tx,rateMdensX/DenomForceX,'-.b',...
    tx,(rateMX+drivPX-dragF2X)/DenomForceX,'-.k',...
    'LineWidth',2);set(gca,'FontSize',14)
if (flagNormalizeTime==1)
    xlabel('t/T')
    xlim([0 1])
else
    xlabel('t (s)')
end
if (flagNormalizeForceX==1)
    ylabel('Normalized Force Terms')
else
    ylabel('Force Terms (g\cdotcm/s)')
end
hlg=legend('$F_{\Delta p}$','$F_{VF}$','$F_f$',...
    '$-\dot{P}_{CV,acc}$','$\dot{P}_{CV,inf}$','$\dot{P}_{CV,den}$',...
    'residual');
set(hlg,'interpreter','latex','location','eastoutside')

figure()
plot(tx,-drivPX/DenomForceX,'g',...
    tx,dragF2X/DenomForceX,'c',...
    tx,-(rateMX+rateMflowX+rateMdensX)/DenomForceX,'b',...
    tx,rateMflowX/DenomForceX,'--b',...
    tx,rateMdensX/DenomForceX,'-.b',...
    tx,(rateMX+drivPX-dragF2X)/DenomForceX,'-.k',...
    'LineWidth',2);set(gca,'FontSize',14)
if (flagNormalizeTime==1)
    xlabel('t/T')
    xlim([0 1])
else
    xlabel('t (s)')
end
if (flagNormalizeForceX==1)
    ylabel('Normalized Force Terms')
else
    ylabel('Force Terms (g\cdotcm/s)')
end
hlg=legend('$F_{\Delta p}+F_{VF}$','$F_f$',...
    '$-\dot{P}_{CV,acc}$','$\dot{P}_{CV,inf}$','$\dot{P}_{CV,den}$',...
    'residual');
set(hlg,'interpreter','latex','location','eastoutside')

figure()
plot(tx,fFluidY/DenomForceY,'g',...
    tx,(-drivPY-fFluidY)/DenomForceY,'r',...
    tx,dragF2Y/DenomForceY,'c',...
    tx,-(rateMY+rateMflowY+rateMdensY)/DenomForceY,'b',...
    tx,rateMflowY/DenomForceY,'--b',...
    tx,rateMdensY/DenomForceY,'-.b',...
    tx,(rateMY+drivPY-dragF2Y)/DenomForceY,'-.k',...
    'LineWidth',2);set(gca,'FontSize',14)
title('Y-momentum')
if (flagNormalizeTime==1)
    xlabel('t/T')
    xlim([0 1])
else
    xlabel('t (s)')
    xlim([tsttVFcycle tendVFcycle])
end
if (flagNormalizeForceY==1)
    ylabel('Normalized Force Terms')
else
    ylabel('Force Terms (g\cdotcm/s)')
end
hlg=legend('$F_{\Delta p}$','$F_{VF}$','$F_f$',...
    '$-\dot{P}_{CV,acc}$','$\dot{P}_{CV,inf}$','$\dot{P}_{CV,den}$',...
    'residual');
set(hlg,'interpreter','latex','location','eastoutside')

figure()
plot(tx,-powrPp/DenomEnergy,'g',...
    tx,(volint+powrDs)/DenomEnergy,'c',...
    tx,-rateKE/DenomEnergy,'b',...
    tx,-powrDs/DenomEnergy,'m',...
    tx,(rateKE+powrPp-volint)/DenomEnergy,'-.k',...
    'LineWidth',2);set(gca,'FontSize',14)
if (flagNormalizeTime==1)
    xlabel('t/T')
    xlim([0 1])
else
    xlabel('t (s)')
end
if (flagNormalizeForceX==1)
    ylabel('Normalized Energy Terms')
else
    ylabel('Energy Terms (g\cdotcm^2/s^2)')
end
hlg=legend('$\dot{W}_{\Delta p}+\dot{W}_{VF}-\dot{W}_{compression}$',...
    '$\dot{W}_{f}$','$-\dot{KE}_{CV}$','$-\dot{W}_{dissipation}$',...
    'residual');
set(hlg,'interpreter','latex','location','eastoutside')

figure()
plot(tx,powrFluid/DenomEnergy,'g',...
    tx,(-powrPp-powrFluid-powrCompr)/DenomEnergy,'r',...
    tx,(volint+powrDs)/DenomEnergy,'c',...
    tx,-(rateKE+rateKEflow+rateKEdens)/DenomEnergy,'b',...
    tx,rateKEflow/DenomEnergy,'--b',...
    tx,rateKEdens/DenomEnergy,'-.b',...
    tx,powrCompr/DenomEnergy,'y',...
    tx,-powrDs/DenomEnergy,'m',...
    tx,(rateKE+powrPp-volint)/DenomEnergy,'-.k',...
    'LineWidth',2);set(gca,'FontSize',14)
if (flagNormalizeTime==1)
    xlabel('t/T')
    xlim([0 1])
else
    xlabel('t (s)')
end
if (flagNormalizeForceX==1)
    ylabel('Normalized Energy Terms')
else
    ylabel('Energy Terms (g\cdotcm^2/s^2)')
end
hlg=legend('$\dot{W}_{\Delta p}$','$\dot{W}_{VF}$','$\dot{W}_{f}$',...
    '$-\dot{KE}_{CV,acc}$','$\dot{KE}_{CV,inf}$','$\dot{KE}_{CV,den}$',...
    '$-\dot{W}_{compression}$',...
    '$-\dot{W}_{dissipation}$',...
    'residual');
set(hlg,'interpreter','latex','location','eastoutside')

figure()
plot(tx,poaberups(iforw,2)/DenomBernoulliUps,'b',...
    tx,poaberups(iforw,3)/DenomBernoulliUps,'r',...
    tx,poaberups(iforw,4)/DenomBernoulliUps,'g',...
    tx,poaberups(iforw,5)/DenomBernoulliUps,'m',...
    tx,-poaberups(iforw,6)/DenomBernoulliUps,'c',...
    tx,poaberups(iforw,7)/DenomBernoulliUps,'-.k',...
    'LineWidth',2);set(gca,'FontSize',14)
if (flagNormalizeTime==1)
    xlabel('t/T')
    xlim([0 1])
else
    xlabel('t (s)')
end
if (flagNormalizeForceX==1)
    ylabel('Normalized Bernoulli Terms')
else
    ylabel('Bernoulli Terms (cm^2/s^2)')
end
title('Contraction Region')
hlg=legend('$E_{con}$',...
    '$E_{pre}$',...
    '$A_{uns}$',...
    '$E_{den}$',...
    '$E_{vis}$',...
    'residual');
set(hlg,'interpreter','latex','location','eastoutside')

figure()
plot(tx,poaberdwn(iforw,2)/DenomBernoulliDwn,'b',...
    tx,poaberdwn(iforw,3)/DenomBernoulliDwn,'r',...
    tx,poaberdwn(iforw,4)/DenomBernoulliDwn,'g',...
    tx,poaberdwn(iforw,5)/DenomBernoulliDwn,'m',...
    tx,-poaberdwn(iforw,6)/DenomBernoulliDwn,'c',...
    tx,poaberdwn(iforw,7)/DenomBernoulliDwn,'-.k',...
    'LineWidth',2);set(gca,'FontSize',14)
if (flagNormalizeTime==1)
    xlabel('t/T')
    xlim([0 1])
else
    xlabel('t (s)')
end
if (flagNormalizeForceX==1)
    ylabel('Normalized Bernoulli Terms')
else
    ylabel('Bernoulli Terms (cm^2/s^2)')
end
title('Jet Region')
hlg=legend('$E_{con}$',...
    '$E_{pre}$',...
    '$A_{uns}$',...
    '$E_{den}$',...
    '$E_{vis}$',...
    'residual');
set(hlg,'interpreter','latex','location','eastoutside')

figure()
plot(tx,massInt,...
    tx,massOut,...
    tx,massRat,...
    tx,poamassrt(iforw,5),...
    'LineWidth',2);set(gca,'FontSize',14)
if (flagNormalizeTime==1)
    xlabel('t/T')
    xlim([0 1])
else
    xlabel('t (s)')
end
ylabel('Mass flowrate (g/s)')
hlg=legend('mass intake','mass outflow','rate of mass change','residual');
set(hlg,'interpreter','latex','location','eastoutside')

figure()
subplot(3,1,1)
plot(tt,Qsep,'LineWidth',2);set(gca,'FontSize',14)
% axis([0 .01 0 70])
xlabel('t (s)')
ylabel('Q_{sep} (cm^2/s)')
subplot(3,1,2)
plot(tt,1.397-ysep,'LineWidth',2);set(gca,'FontSize',14)
% axis([0 .01 0 .11])
xlabel('t (s)')
ylabel('S_{sep} (cm)')
subplot(3,1,3)
plot(tt,xsep,'LineWidth',2);set(gca,'FontSize',14)
% axis([0 .01 3.2 3.35])
xlabel('t (s)')
ylabel('x_{sep} (cm)')

figure()
plot(tt,hglottisReal,'LineWidth',2);set(gca,'FontSize',14)
% xlim([0 1])
xlabel('t (s)')
ylabel('h (cm)')