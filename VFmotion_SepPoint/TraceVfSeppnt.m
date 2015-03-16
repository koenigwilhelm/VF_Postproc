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
xsep_orig=seppnt(:,3);
ysep_orig=seppnt(:,4);
xsep=xsep_orig(iforw);
ysep=ysep_orig(iforw);

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
figure()
filename='AnimVFseppntMov.gif';
flagAppend=0;
tstt=6790;
tend=7290;
% for iit=1:nt
for iit=tstt:tend
    clf
    hold on
    plot(xsbc(iit,:),ysbc(iit,:))
    for ix=1:nstep:nnvfsurf
        plot(xsbc(tstt:iit,ix),ysbc(tstt:iit,ix),'m')
    end
    plot(xsep(iit-10:iit),ysep(iit-10:iit),'k')
    plot(xsep(iit),ysep(iit),'-o')
    xlabel('x (cm)')
    ylabel('y (cm)')
    axis([xlim_min xlim_max 0 1.4])
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

%% plot results
figure()
hold on
for ix=1:nstep:nnvfsurf
    plot(xsbc(:,ix),ysbc(:,ix),'m')
end
plot(xsep,ysep,'k')
plot(mean(xsbc,1),mean(ysbc,1))
xlabel('x (cm)')
ylabel('y (cm)')
axis([xlim_min xlim_max 0 1.4])
hold off

%% SVD of VF motion
tstt=4410;
tend=7290;
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
hold off
clear Xsbc Usbc Ssbc Vsbc

%% VF vibration frequency
ymaxVSt=max(ysbc,[],2);
hglottis=1.397-ymaxVSt;
hglottis=hglottis-mean(hglottis);
Nfft=2^(nextpow2(nt)+6);
Yhg=fft(hglottis,Nfft)/nt;
freq=linspace(0,1,Nfft/2+1)/(2*(tt(2)-tt(1)));
figure()
plot(freq,abs(Yhg(1:Nfft/2+1))/max(abs(Yhg(1:Nfft/2+1))))
axis([0 500 0 1.1])
xlabel('frequency (Hz)')
ylabel('relative amplitude')