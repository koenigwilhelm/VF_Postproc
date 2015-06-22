close all
clear all

% Halfspace(0)/Fullspace(1)
flagHalfFullSpace=0;

%% Options for looks of figures
flagFiguresVersion=1;
flagShowResidual=0;

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
if (flagHalfFullSpace==0)
    tsttSVD=7000;
    tendSVD=9000;
else
    tsttSVD=3000;
    tendSVD=4000;
end
% start and end times for interested VF vibration cycle
if (flagHalfFullSpace==0)
    tsttVFcycle=.1773;
    tendVFcycle=.1816;
else
    tsttVFcycle=.09240;
    tendVFcycle=.09666;
end

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
xsbc=xsbc(:,(ysbc_orig(1,:)<2.974)&(ysbc_orig(1,:)>0)|(xsbc_orig(1,:)==min(xsbc_orig(1,:)))|(xsbc_orig(1,:)==max(xsbc_orig(1,:))));
ysbc=ysbc(:,(ysbc_orig(1,:)<2.974)&(ysbc_orig(1,:)>0)|(xsbc_orig(1,:)==min(xsbc_orig(1,:)))|(xsbc_orig(1,:)==max(xsbc_orig(1,:))));
nnvfsurf=length(xsbc(1,:));
nt=length(tt);
Arsort=[xsbc' ysbc'];
Brsort=sortrows(Arsort,[1 -(nt+1)]);
xsbc=Brsort(:,1:nt)';
ysbc=Brsort(:,nt+1:end)';
clear Arsort Brsort seppnt sbctrc xsbc_orig ysbc_orig xsep_orig ysep_orig
clear tt_orig

%% parameters for figures
nstep=10;
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
        axis([xlim_min-xshift xlim_max-xshift 0 1.45])
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
if (flagFiguresVersion==0)
    figure()
    hold on
    for ix=1:nstep:nnvfsurf
        plot(xsbc(:,ix)-xshift,ysbc(:,ix),'m')
    end
    plot(xsep-xshift,ysep,'k')
    plot(mean(xsbc,1)-xshift,mean(ysbc,1))
    set(gca,'FontSize',14)
    xh=xlabel('x (cm)');
    yh=ylabel('y (cm)');
    set([xh,yh],'FontWeight','bold')
    axis([xlim_min-xshift xlim_max-xshift 0 1.45])
    hold off
end

if (flagHalfFullSpace==0)
    %% SVD of VF motion
    Xsbc=[xsbc(tsttSVD:tendSVD,:) ysbc(tsttSVD:tendSVD,:)];
    Xsbcmean=mean(Xsbc);
    Xsbc=Xsbc-repmat(Xsbcmean,tendSVD-tsttSVD+1,1);
    [Usbc,Ssbc,Vsbc]=svd(Xsbc,'econ');
    figure()
    hold on
    if (flagFiguresVersion==1)
        for ix=1:nstep:nnvfsurf
            plot(xsbc(tsttSVD:tendSVD,ix)-xshift,ysbc(tsttSVD:tendSVD,ix),'m',...
                'LineWidth',2)
        end
    %     plot(mean(xsbc,1)-xshift,mean(ysbc,1),'g',...
    %         'LineWidth',2)
        plot(mean(xsbc,1)+Vsbc(1:nnvfsurf,1)'*Ssbc(1,1)*min(Usbc(:,1))-xshift,...
            mean(ysbc,1)+Vsbc(nnvfsurf+1:end,1)'*Ssbc(1,1)*min(Usbc(:,1)),'--r',...
            'LineWidth',2)
        plot(mean(xsbc,1)+Vsbc(1:nnvfsurf,1)'*Ssbc(1,1)*max(Usbc(:,1))-xshift,...
            mean(ysbc,1)+Vsbc(nnvfsurf+1:end,1)'*Ssbc(1,1)*max(Usbc(:,1)),'-.b',...
            'LineWidth',2)
    else
        for ix=1:nstep:nnvfsurf
            plot(xsbc(tsttSVD:tendSVD,ix)-xshift,ysbc(tsttSVD:tendSVD,ix),'m')
        end
        plot(mean(xsbc,1)-xshift,mean(ysbc,1))
        plot(mean(xsbc,1)+Vsbc(1:nnvfsurf,1)'*Ssbc(1,1)*min(Usbc(:,1))-xshift,...
            mean(ysbc,1)+Vsbc(nnvfsurf+1:end,1)'*Ssbc(1,1)*min(Usbc(:,1)),'r')
        plot(mean(xsbc,1)+Vsbc(1:nnvfsurf,1)'*Ssbc(1,1)*max(Usbc(:,1))-xshift,...
            mean(ysbc,1)+Vsbc(nnvfsurf+1:end,1)'*Ssbc(1,1)*max(Usbc(:,1)),'g')
    end
    set(gca,'FontSize',14)
    axis([xlim_min-xshift xlim_max-xshift 0 1.45])
    xlabel('x (cm)')
    ylabel('y (cm)')
    hold off
    clear Xsbc Usbc Ssbc Vsbc

    %% VF vibration frequency
    ymaxVSt=max(ysbc,[],2);
    hglottis=1.397-ymaxVSt;
    hglottisVar=hglottis-mean(hglottis);
    % hglottisReal=(hglottis-0.0114).*(hglottis-0.0114>0);
    hglottisReal=(hglottis-0.0073).*(hglottis-0.0073>0);
    Nfft=2^(nextpow2(nt)+6);
    Yhg=fft(hglottisVar,Nfft)/nt;
    freq=linspace(0,1,Nfft/2+1)/(2*(tt(2)-tt(1)));
    figure()
    subplot(2,1,1)
    if (flagFiguresVersion==1)
        plot(freq,abs(Yhg(1:Nfft/2+1))/max(abs(Yhg(1:Nfft/2+1))),'Color','b','LineWidth',2)
    else
        plot(freq,abs(Yhg(1:Nfft/2+1))/max(abs(Yhg(1:Nfft/2+1))),'LineWidth',2)
    end
    set(gca,'FontSize',14)
    axis([0 500 0 1.1])
    xlabel('frequency (Hz)')
    ylabel('relative amplitude')
    subplot(2,1,2)
    if (flagFiguresVersion==1)
        plot(tt,hglottisReal,'Color','b','LineWidth',2)
    else
        plot(tt,hglottis,'LineWidth',2)
    end
    ylim([-.01 .06])
    set(gca,'FontSize',14)
    xlabel('t (s)')
    ylabel('h (cm)')
    if (flagFiguresVersion==1)
        figure()
        plot(tt(tsttSVD:tendSVD),hglottis(tsttSVD:tendSVD),'Color','b','LineWidth',2)
        axis([tt(tsttSVD) tt(tendSVD) -.01 .07])
        set(gca,'FontSize',14)
        xlabel('t (s)')
        ylabel('h (cm)')
    end
end

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

if (flagFiguresVersion==1)
    % all terms separate, except decompostion of momentum flux
    figure()
    hold on
    plot(tx,fFluidX/DenomForceX,'g',...
        tx,(-drivPX-fFluidX)/DenomForceX,'r',...
        tx,dragF2X/DenomForceX,'c',...
        tx,-(rateMX)/DenomForceX,'b',...
        'LineWidth',2);set(gca,'FontSize',14)
    if (flagShowResidual==1)
        plot(tx,(rateMX+drivPX-dragF2X)/DenomForceX,'-.k',...
            'LineWidth',2)
    end
    hold off
    if (flagNormalizeTime==1)
        xlabel('t/T')
        if (flagHalfFullSpace==0)
            xlim([0 1])
        else
            xlim([0 3])
        end
    else
        xlabel('t (s)')
    end
    if (flagNormalizeForceX==1)
        ylabel('Normalized Force Terms')
    else
        ylabel('Force Terms (g\cdotcm/s)')
    end
    hlg=legend('$F_{\Delta p}$','$F_{VF}$','$F_f$',...
        '$-\dot{P}_{CV}$',...
        'residual');
    set(hlg,'interpreter','latex','location','eastoutside')
    % pressure drive and drag combined as a net pressure force
    figure()
    hold on
    plot(tx,-drivPX/DenomForceX,'g',...
        tx,dragF2X/DenomForceX,'c',...
        tx,-rateMX/DenomForceX,'b',...
        'LineWidth',2);set(gca,'FontSize',14)
    if (flagShowResidual==1)
        plot(tx,(rateMX+drivPX-dragF2X)/DenomForceX,'-.k',...
            'LineWidth',2)
    end
    hold off
    if (flagNormalizeTime==1)
        xlabel('t/T')
        if (flagHalfFullSpace==0)
            xlim([0 1])
        else
            xlim([0 3])
        end
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
    % decomposition of the momentum flux term
    figure()
    plot(tx,(rateMX)/DenomForceX,'b',...
        tx,(rateMX+rateMflowX+rateMdensX)/DenomForceX,'-.r',...
        tx,-rateMflowX/DenomForceX,'--m',...
        tx,-rateMdensX/DenomForceX,':g',...
        'LineWidth',2);set(gca,'FontSize',14)
    if (flagNormalizeTime==1)
        xlabel('t/T')
        if (flagHalfFullSpace==0)
            xlim([0 1])
        else
            xlim([0 3])
        end
    else
        xlabel('t (s)')
    end
    if (flagNormalizeForceX==1)
        ylabel('Normalized Force Terms')
    else
        ylabel('Force Terms (g\cdotcm/s)')
    end
    hlg=legend('$\dot{P}_{CV}$',...
        '$\dot{P}_{CV,acc}$','$-\dot{P}_{CV,inf}$','$-\dot{P}_{CV,den}$');
    set(hlg,'interpreter','latex','location','eastoutside')
else
    % pressure drive and drag combined as a net pressure force
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
    % all terms plotted separately, even 3 momemtum flux terms
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
    % pressure drive and drag combined as a net pressure force;
    % 3 momemtum flux terms plotted separately
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
end

% Y-momentum
figure()
if (flagFiguresVersion==1)
    hold on
    plot(tx,(-drivPY-fFluidY)/DenomForceY,'r',...
            tx,-(rateMY)/DenomForceY,'b',...
        'LineWidth',2);set(gca,'FontSize',14)
    if (flagShowResidual==1)
        plot(tx,fFluidY/DenomForceY,'g',...
            tx,dragF2Y/DenomForceY,'c',...
            tx,(rateMY+drivPY-dragF2Y)/DenomForceY,'-.k',...
            'LineWidth',2)
    end
    hold off
    hlg=legend('$F_{VF}$','$-\dot{P}_{CV}$','$F_{\Delta p}$','$F_f$',...
        'residual');
else
    plot(tx,fFluidY/DenomForceY,'g',...
        tx,(-drivPY-fFluidY)/DenomForceY,'r',...
        tx,dragF2Y/DenomForceY,'c',...
        tx,-(rateMY+rateMflowY+rateMdensY)/DenomForceY,'b',...
        tx,rateMflowY/DenomForceY,'--b',...
        tx,rateMdensY/DenomForceY,'-.b',...
        tx,(rateMY+drivPY-dragF2Y)/DenomForceY,'-.k',...
        'LineWidth',2);set(gca,'FontSize',14)
    title('Y-momentum')
    hlg=legend('$F_{\Delta p}$','$F_{VF}$','$F_f$',...
    '$-\dot{P}_{CV,acc}$','$\dot{P}_{CV,inf}$','$\dot{P}_{CV,den}$',...
    'residual');
end
if (flagNormalizeTime==1)
    xlabel('t/T')
    if (flagHalfFullSpace==0)
        xlim([0 1])
    else
        xlim([0 3])
    end
else
    xlabel('t (s)')
    xlim([tsttVFcycle tendVFcycle])
end
if (flagNormalizeForceY==1)
    ylabel('Normalized Force Terms')
else
    ylabel('Force Terms (g\cdotcm/s)')
end
set(hlg,'interpreter','latex','location','eastoutside')

% pressure work on boundary and compression work combined
figure()
hold on
plot(tx,-powrPp/DenomEnergy,'g',...
    tx,(volint+powrDs)/DenomEnergy,'c',...
    tx,-rateKE/DenomEnergy,'b',...
    tx,-powrDs/DenomEnergy,'m',...
    'LineWidth',2);set(gca,'FontSize',14)
if (flagShowResidual==1)
    plot(tx,(rateKE+powrPp-volint)/DenomEnergy,'-.k',...
        'LineWidth',2)
end
hold off
if (flagNormalizeTime==1)
    xlabel('t/T')
    if (flagHalfFullSpace==0)
        xlim([0 1])
    else
        xlim([0 3])
    end
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
if (flagFiguresVersion==1)
    % all terms plotted separately
    figure()
    hold on
    plot(tx,powrFluid/DenomEnergy,'g',...
        tx,(-powrPp-powrFluid-powrCompr)/DenomEnergy,'r',...
        tx,(volint+powrDs)/DenomEnergy,'c',...
        tx,-(rateKE)/DenomEnergy,'b',...
        tx,powrCompr/DenomEnergy,'y',...
        tx,-powrDs/DenomEnergy,'m',...
        'LineWidth',2);set(gca,'FontSize',14)
    if (flagShowResidual==1)
        plot(tx,(rateKE+powrPp-volint)/DenomEnergy,'-.k',...
            'LineWidth',2)
    end
    hold off
    if (flagNormalizeTime==1)
        xlabel('t/T')
        if (flagHalfFullSpace==0)
            xlim([0 1])
        else
            xlim([0 3])
        end
    else
        xlabel('t (s)')
    end
    if (flagNormalizeForceX==1)
        ylabel('Normalized Energy Terms')
    else
        ylabel('Energy Terms (g\cdotcm^2/s^2)')
    end
    hlg=legend('$\dot{W}_{\Delta p}$','$\dot{W}_{VF}$','$\dot{W}_{f}$',...
        '$-\dot{KE}_{CV}$',...
        '$-\dot{W}_{compression}$',...
        '$-\dot{W}_{dissipation}$',...
        'residual');
    set(hlg,'interpreter','latex','location','eastoutside')
    % KE rate decomposed
    figure()
    plot(tx,(rateKE)/DenomEnergy,'b',...
        tx,(rateKE+rateKEflow+rateKEdens)/DenomEnergy,'-.r',...
        tx,-rateKEflow/DenomEnergy,'--m',...
        tx,-rateKEdens/DenomEnergy,':g',...
        'LineWidth',2);set(gca,'FontSize',14)
    if (flagNormalizeTime==1)
        xlabel('t/T')
        if (flagHalfFullSpace==0)
            xlim([0 1])
        else
            xlim([0 3])
        end
    else
        xlabel('t (s)')
    end
    if (flagNormalizeForceX==1)
        ylabel('Normalized Energy Terms')
    else
        ylabel('Energy Terms (g\cdotcm^2/s^2)')
    end
    hlg=legend('$\dot{KE}_{CV}$',...
        '$\dot{KE}_{CV,acc}$','$-\dot{KE}_{CV,inf}$','$-\dot{KE}_{CV,den}$');
    set(hlg,'interpreter','latex','location','eastoutside')
else
    % all terms plotted separately, even KE rate decomposed
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
end

if (flagHalfFullSpace==0)
    figure()
    hold on
    plot(tx,poaberups(iforw,2)/DenomBernoulliUps,'b',...
        tx,poaberups(iforw,3)/DenomBernoulliUps,'r',...
        tx,poaberups(iforw,4)/DenomBernoulliUps,'g',...
        tx,poaberups(iforw,5)/DenomBernoulliUps,'m',...
        tx,-poaberups(iforw,6)/DenomBernoulliUps,'c',...
        'LineWidth',2);set(gca,'FontSize',14)
    if (flagShowResidual==1)
        plot(tx,poaberups(iforw,7)/DenomBernoulliUps,'-.k',...
            'LineWidth',2)
    end
    hold off
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
    if (flagFiguresVersion==0)
        title('Contraction Region')
    end
    hlg=legend('$E_{con}$',...
        '$E_{pre}$',...
        '$A_{uns}$',...
        '$E_{den}$',...
        '$E_{vis}$',...
        'residual');
    set(hlg,'interpreter','latex','location','eastoutside')

    figure()
    hold on
    plot(tx,poaberdwn(iforw,2)/DenomBernoulliDwn,'b',...
        tx,poaberdwn(iforw,3)/DenomBernoulliDwn,'r',...
        tx,poaberdwn(iforw,4)/DenomBernoulliDwn,'g',...
        tx,poaberdwn(iforw,5)/DenomBernoulliDwn,'m',...
        tx,-poaberdwn(iforw,6)/DenomBernoulliDwn,'c',...
        'LineWidth',2);set(gca,'FontSize',14)
    if (flagShowResidual==1)
        plot(tx,poaberdwn(iforw,7)/DenomBernoulliDwn,'-.k',...
            'LineWidth',2)
    end
    hold off
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
    if (flagFiguresVersion==0)
        title('Jet Region')
    end
    hlg=legend('$E_{con}$',...
        '$E_{pre}$',...
        '$A_{uns}$',...
        '$E_{den}$',...
        '$E_{vis}$',...
        'residual');
    set(hlg,'interpreter','latex','location','eastoutside')

    figure()
    hold on
    plot(tx,massInt,'g',...
        tx,massOut,'r',...
        tx,massRat,'b',...
        'LineWidth',2);set(gca,'FontSize',14)
    if (flagShowResidual==1)
        plot(tx,poamassrt(iforw,5),'-.k',...
            tx,tx*0,'--m',...
            'LineWidth',2)
    end
    hold off
    if (flagNormalizeTime==1)
        xlabel('t/T')
        xlim([0 1])
    else
        xlabel('t (s)')
    end
    ylabel('Mass flowrate (g/s)')
    hlg=legend('mass intake','mass outflow','rate of mass change','residual');
    set(hlg,'interpreter','latex','location','eastoutside')

    if (flagFiguresVersion==0)
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
    end

    figure()
    plot(tx,hglottisReal,'LineWidth',2);set(gca,'FontSize',14)
    axis([0 1 -.01 .06])
    xlabel('t (s)')
    ylabel('h (cm)')
end