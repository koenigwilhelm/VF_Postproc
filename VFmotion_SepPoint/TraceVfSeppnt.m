close all
clear all

seppnt=importdata('poa_seppnt.txt',' ');
sbctrc=importdata('poa_sbctrace.txt',' ');

%% get unique datasets (restart of simulation often leads to duplicate time steps)
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
xlim_min=min(xsbc(1,:))-.1;
xlim_max=max(xsbc(1,:))+.3;

%% save animation
% figure()
% filename='AnimVFseppntMov.gif';
% flagAppend=0;
% tstt=531;
% tend=1014;
% % for iit=1:nt
% for iit=tstt:tend
%     clf
%     hold on
%     plot(xsbc(iit,:),ysbc(iit,:))
%     for ix=1:nstep:nnvfsurf
%         plot(xsbc(tstt:iit,ix),ysbc(tstt:iit,ix),'m')
%     end
%     plot(xsep(iit-10:iit),ysep(iit-10:iit),'k')
%     plot(xsep(iit),ysep(iit),'-o')
%     xlabel('x (cm)')
%     ylabel('y (cm)')
%     axis([xlim_min xlim_max 0 1.4])
%     hold off
%     %pause(.05)
%     drawnow
%     frame=getframe(1);
%     im=frame2im(frame);
%     [A,map]=rgb2ind(im,256);
%     if flagAppend==0
%         imwrite(A,map,filename,'gif','LoopCount',Inf,'DelayTime',.02);
%         flagAppend=1;
%     else
%         imwrite(A,map,filename,'gif','WriteMode','Append','DelayTime',.02);
%     end
% end

%% plot results
figure()
hold on
for ix=1:nstep:nnvfsurf
    plot(xsbc(:,ix),ysbc(:,ix),'m')
end
plot(xsep,ysep,'k')
plot(mean(xsbc,1),mean(ysbc,1))
set(gca,'FontSize',14)
xlabel('x (cm)')
ylabel('y (cm)')
axis([xlim_min xlim_max 0 1.4])
hold off

%% SVD of VF motion
tstt=1;
tend=500;
Xsbc=[xsbc(tstt:tend,:) ysbc(tstt:tend,:)];
Xsbcmean=mean(Xsbc);
Xsbc=Xsbc-repmat(Xsbcmean,tend-tstt+1,1);
[Usbc,Ssbc,Vsbc]=svd(Xsbc,'econ');
figure()
hold on
for ix=1:nstep:nnvfsurf
    plot(xsbc(tstt:tend,ix),ysbc(tstt:tend,ix),'m')
end
plot(mean(xsbc,1),mean(ysbc,1))
plot(mean(xsbc,1)+Vsbc(1:nnvfsurf,1)'*Ssbc(1,1)*min(Usbc(:,1)),mean(ysbc,1)+Vsbc(nnvfsurf+1:end,1)'*Ssbc(1,1)*min(Usbc(:,1)),'r')
plot(mean(xsbc,1)+Vsbc(1:nnvfsurf,1)'*Ssbc(1,1)*max(Usbc(:,1)),mean(ysbc,1)+Vsbc(nnvfsurf+1:end,1)'*Ssbc(1,1)*max(Usbc(:,1)),'g')
set(gca,'FontSize',14)
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

%% Bernoulli's and Control Volume Analysis
poamomentumX=importdata('poa_momentumX.txt',' ');
poaenergy=importdata('poa_energy.txt',' ');
poaberups=importdata('poa_bernoulliUps.txt',' ');
poaberdwn=importdata('poa_bernoulliDwn.txt',' ');
poamassrt=importdata('poa_massCons.txt',' ');

t=poamomentumX(:,1);
rateM=poamomentumX(:,2)+poamomentumX(:,3);
rateMflow=poamomentumX(:,8);
rateMdens=poamomentumX(:,9);
fFluid=poamomentumX(:,10);
drivP=poamomentumX(:,4);
dragF=poamomentumX(:,5);
resdM=poamomentumX(:,6);

dragF2=poamomentumX(:,7);

rateKE=poaenergy(:,2)+poaenergy(:,3);
powrPp=poaenergy(:,4);
powrFr=poaenergy(:,5);
powrDs=poaenergy(:,6);
resdEn=poaenergy(:,7);
volint = poaenergy(:,8);
rateKEflow=poaenergy(:,9);
rateKEdens=poaenergy(:,10);
powrFluid=poaenergy(:,11);
powrCompr=poaenergy(:,12);

massInt=poamassrt(:,2);
massOut=poamassrt(:,3);
massRat=poamassrt(:,4);

tx=(t-0)/(.01-0);
elwidth=.005;

figure()
plot(tx,-drivP,'g',...
    tx,dragF2,'c',...
    tx,-rateM,'b',...
    tx,(rateM+drivP-dragF2),'-.k',...
    'LineWidth',2);set(gca,'FontSize',14)
% xlim([0 1])
xlabel('t/T')
ylabel('force terms')
hlg=legend('$F_{\Delta p}+F_{VF}$','$F_f$','$-\dot{P}_{CV}$','residual');
set(hlg,'interpreter','latex','location','eastoutside')

figure()
plot(tx,fFluid,'g',...
    tx,-drivP-fFluid,'r',...
    tx,dragF2,'c',...
    tx,-(rateM+rateMflow+rateMdens),'b',...
    tx,rateMflow,'--b',...
    tx,rateMdens,'-.b',...
    tx,(rateM+drivP-dragF2),'-.k',...
    'LineWidth',2);set(gca,'FontSize',14)
% xlim([0 1])
xlabel('t/T')
ylabel('force terms')
hlg=legend('$F_{\Delta p}$','$F_{VF}$','$F_f$',...
    '$-\dot{P}_{CV,acc}$','$\dot{P}_{CV,inf}$','$\dot{P}_{CV,den}$',...
    'residual');
set(hlg,'interpreter','latex','location','eastoutside')

figure()
plot(tx,-drivP,'g',...
    tx,dragF2,'c',...
    tx,-(rateM+rateMflow+rateMdens),'b',...
    tx,rateMflow,'--b',...
    tx,rateMdens,'-.b',...
    tx,(rateM+drivP-dragF2),'-.k',...
    'LineWidth',2);set(gca,'FontSize',14)
% xlim([0 1])
xlabel('t/T')
ylabel('force terms')
hlg=legend('$F_{\Delta p}+F_{VF}$','$F_f$',...
    '$-\dot{P}_{CV,acc}$','$\dot{P}_{CV,inf}$','$\dot{P}_{CV,den}$',...
    'residual');
set(hlg,'interpreter','latex','location','eastoutside')

figure()
plot(tx,-powrPp,'g',...
    tx,volint+powrDs,'c',...
    tx,-rateKE,'b',...
    tx,-powrDs,'m',...
    tx,(rateKE+powrPp-volint),'-.k',...
    'LineWidth',2);set(gca,'FontSize',14)
% xlim([0 1])
xlabel('t/T')
ylabel('energy terms')
hlg=legend('$\dot{W}_{\Delta p}+\dot{W}_{VF}-\dot{W}_{compression}$',...
    '$\dot{W}_{f}$','$-\dot{KE}_{CV}$','$-\dot{W}_{dissipation}$',...
    'residual');
set(hlg,'interpreter','latex','location','eastoutside')

figure()
plot(tx,powrFluid,'g',...
    tx,-powrPp-powrFluid-powrCompr,'r',...
    tx,volint+powrDs,'c',...
    tx,-(rateKE+rateKEflow+rateKEdens),'b',...
    tx,rateKEflow,'--b',...
    tx,rateKEdens,'-.b',...
    tx,powrCompr,'y',...
    tx,-powrDs,'m',...
    tx,(rateKE+powrPp-volint),'-.k',...
    'LineWidth',2);set(gca,'FontSize',14)
% xlim([0 1])
xlabel('t/T')
ylabel('energy terms')
hlg=legend('$\dot{W}_{\Delta p}$','$\dot{W}_{VF}$','$\dot{W}_{f}$',...
    '$-\dot{KE}_{CV,acc}$','$\dot{KE}_{CV,inf}$','$\dot{KE}_{CV,den}$',...
    '$-\dot{W}_{compression}$',...
    '$-\dot{W}_{dissipation}$',...
    'residual');
set(hlg,'interpreter','latex','location','eastoutside')

figure()
plot(tx,poaberups(:,2)/elwidth,'b',...
    tx,poaberups(:,3)/elwidth,'g',...
    tx,poaberups(:,4)/elwidth,'c',...
    tx,poaberups(:,5)/elwidth,'m',...
    tx,-poaberups(:,6)/elwidth,'r',...
    tx,poaberups(:,7)/elwidth,'-.k',...
    'LineWidth',2);set(gca,'FontSize',14)
% xlim([0 1])
xlabel('t/T')
ylabel('Bernoulli terms (cm^2/s^2)')
title('Contraction Region')
hlg=legend('$E_{con}$',...
    '$E_{pre}$',...
    '$A_{uns}$',...
    '$E_{den}$',...
    '$E_{vis}$',...
    'residual');
set(hlg,'interpreter','latex','location','eastoutside')

figure()
plot(tx,poaberdwn(:,2)/elwidth,'b',...
    tx,poaberdwn(:,3)/elwidth,'g',...
    tx,poaberdwn(:,4)/elwidth,'c',...
    tx,poaberdwn(:,5)/elwidth,'m',...
    tx,-poaberdwn(:,6)/elwidth,'r',...
    tx,poaberdwn(:,7)/elwidth,'-.k',...
    'LineWidth',2);set(gca,'FontSize',14)
% xlim([0 1])
xlabel('t/T')
ylabel('Bernoulli terms (cm^2/s^2)')
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
    tx,poamassrt(:,5),...
    'LineWidth',2)
legend('mass intake','mass outflow','rate of mass change','residual')

figure()
subplot(3,1,1)
plot(tt,Qsep,'LineWidth',2);set(gca,'FontSize',14)
% axis([0 .01 0 70])
xlabel('t/T')
ylabel('Q_{sep} (cm^2/s)')
subplot(3,1,2)
plot(tt,1.397-ysep,'LineWidth',2);set(gca,'FontSize',14)
% axis([0 .01 0 .11])
xlabel('t/T')
ylabel('S_{sep} (cm)')
subplot(3,1,3)
plot(tt,xsep,'LineWidth',2);set(gca,'FontSize',14)
% axis([0 .01 3.2 3.35])
xlabel('t/T')
ylabel('x_{sep} (cm)')

figure()
plot(tt,hglottisReal,'LineWidth',2);set(gca,'FontSize',14)
% xlim([0 1])
xlabel('t/T')
ylabel('h (cm)')