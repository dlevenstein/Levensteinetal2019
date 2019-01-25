%%WCadapt_DwellStats
dropboxroot = '/Users/dlevenstein/Dropbox/Research/';
%dropboxroot = '/mnt/data1/Dropbox/research/';
figfolder = fullfile(dropboxroot,'Current Projects/SlowOscillation/Modeling/Figures');
simdatafolder =  fullfile(dropboxroot,'Current Projects/SlowOscillation/AnalysisScripts/simdata/');

%% Color maps
DOWNcolor = makeColorMap([0.5 0.5 0.5],[0 0 0.8]);
UPcolor = makeColorMap([0.5 0.5 0.5],[0.8 0 0]);
UPDOWNcolor = [flipud(DOWNcolor);UPcolor];
%UP/DOWN stats   
CVcolor = [1 1 1; ...
    makeColorMap(0.8.*[1 1 1],0.8.*[1 1 1]);...
    makeColorMap(0.8.*[1 1 1],[0 0.5 0]);...
    makeColorMap([0 0.5 0],[0.7 0.5 0])];



%% I-W space

simdatafilename = 'dwellmatch6.mat';
simdatafullfile = fullfile(simdatafolder,simdatafilename);

resim=false;
if ~exist(simdatafullfile,'file') || resim
  % UP/DOWN Duration Map
    numI = 50;
    numW = 50;
    testI = linspace(-3.3,-1,numI);
    testbw = linspace(3.5,7.5,numW);
    testI = repmat(testI,numW,1);testI = testI(:);
    testbw = repmat(testbw,1,numI);testbw = testbw(:);
    numsims = length(testI);

    for nn = 1:numsims
        modelparms(nn).N_neurons = 1;
        modelparms(nn).I_in = testI(nn);
        modelparms(nn).W = testbw(nn);
        modelparms(nn).beta = 1;
        modelparms(nn).tau_r = 1;
        modelparms(nn).tau_a = 25;
        modelparms(nn).A0 = 0.5;
        modelparms(nn).Ak = 15;
        modelparms(nn).noiseamp =0.25;
        modelparms(nn).noisefreq = 0.05;
    end

    simparms.simtime = 60000;
    simparms.dt = 1;

    [ dwelltimes_sim,ratehist ] = BatchSimulate_WCadapt( modelparms,simparms );
    save(simdatafullfile,'dwelltimes_sim','ratehist','modelparms','simparms');
else 
    load(simdatafullfile);
end

%Extract Parameters
modelparmsvec = CollapseStruct(modelparms);
ratehistvec = CollapseStruct(ratehist);

I_in = unique(modelparmsvec.I_in);
W = unique(modelparmsvec.W);
meanrate = reshape(ratehistvec.mean,length(I_in),length(W));

%Load Bif'n Diagram
bifnlocation = fullfile(dropboxroot,'Current Projects/SlowOscillation/Figures/IWbifn.dat');
[ IWbifn ] = BifnFromXPP( bifnlocation );
hopfline = [IWbifn(:,1),IWbifn(:,6)];
snline = [IWbifn(:,1),IWbifn(:,13)];

jumpthresh = 0.05;
jumpsize = diff(snline);
snline(jumpsize>jumpthresh)=nan;
jumpsize = diff(hopfline);
hopfline(jumpsize>jumpthresh)=nan;

%% IW: Dwell time stats
clear dwelltimemeans; clear dwelltimestd
numsims = length(dwelltimes_sim);
for nn = 1:numsims
    dwelltimemeans.UP(nn) = nanmean(dwelltimes_sim(nn).UP);
    dwelltimemeans.DOWN(nn) = nanmean(dwelltimes_sim(nn).DOWN);
    
    dwelltimestd.UP(nn) = nanstd(dwelltimes_sim(nn).UP);
    dwelltimestd.DOWN(nn) = nanstd(dwelltimes_sim(nn).DOWN);
    
    dwelltimemeans.numUPs(nn) = length(dwelltimes_sim(nn).UP);
end

%Calculate mean ratio, CVs
dwelltimemeans.ratio = dwelltimemeans.UP./dwelltimemeans.DOWN;
dwelltimeCV.UP = dwelltimestd.UP./ dwelltimemeans.UP;
dwelltimeCV.DOWN = dwelltimestd.DOWN./ dwelltimemeans.DOWN;

%Reshape for image
dwelltimemeans.UP = reshape(dwelltimemeans.UP,length(I_in),length(W));
dwelltimemeans.DOWN = reshape(dwelltimemeans.DOWN,length(I_in),length(W));
dwelltimemeans.numUPs = reshape(dwelltimemeans.numUPs,length(I_in),length(W));
dwelltimeCV.UP = reshape(dwelltimeCV.UP,length(I_in),length(W));
dwelltimeCV.DOWN = reshape(dwelltimeCV.DOWN,length(I_in),length(W));

dwelltimemeans.ratio = dwelltimemeans.UP./dwelltimemeans.DOWN;


%% Figure : IW
ratecolor = makeColorMap([1 1 1],[0.5 0.5 0.8],[0.8 0.5 0.5]);
figure
colormap(ratecolor)
subplot(3,3,9)
    imagesc(I_in,W,meanrate)
    hold on
    plot(snline(:,1),snline(:,2),'k','Linewidth',2)
    plot(hopfline(:,1),hopfline(:,2),'--','color','k','Linewidth',1)
    colorbar
    caxis([0 1])
    %plot(snline(:,1),snline(:,2),'k','Linewidth',2)
    %plot(hopfline(:,1),hopfline(:,2),'--','color',hopfcolor,'Linewidth',1)
    axis xy
    %xlim([-3.5 0]);ylim([0 7.5])
   
    
   

subplot(3,3,3)
colormap(gca,UPcolor)
    h = imagesc(I_in,W,log10(dwelltimemeans.UP));
    set(h, 'AlphaData', dwelltimemeans.numUPs>8);
    hold on
    plot(snline(:,1),snline(:,2),'k','Linewidth',2)
    plot(hopfline(:,1),hopfline(:,2),'--','color','k','Linewidth',1)
    colorbar
    axis xy
    caxis([1 3])
    
subplot(3,3,6)
colormap(gca,CVcolor)
    h = imagesc(I_in,W,(dwelltimeCV.UP));
    set(h, 'AlphaData', dwelltimemeans.numUPs>8);
    hold on
    plot(snline(:,1),snline(:,2),'k','Linewidth',2)
    plot(hopfline(:,1),hopfline(:,2),'k--','Linewidth',1)
    colorbar
    axis xy
    caxis([0.1 1])

    
 
subplot(3,3,7)
colormap(gca,DOWNcolor)
    h = imagesc(I_in,W,log10(dwelltimemeans.DOWN));
    set(h, 'AlphaData', dwelltimemeans.numUPs>8);
    hold on
    plot(snline(:,1),snline(:,2),'k','Linewidth',2)
    plot(hopfline(:,1),hopfline(:,2),'k--','Linewidth',1)
    colorbar
    axis xy
    caxis([1 3])
    
subplot(3,3,8)
colormap(gca,CVcolor)
    h = imagesc(I_in,W,(dwelltimeCV.DOWN));
    set(h, 'AlphaData', dwelltimemeans.numUPs>8);
    hold on
    plot(snline(:,1),snline(:,2),'k','Linewidth',2)
    plot(hopfline(:,1),hopfline(:,2),'k--','Linewidth',1)
    colorbar
    axis xy
    caxis([0.1 1])
    
subplot(2,2,1)
colormap(gca,UPDOWNcolor)
    h = imagesc(I_in,W,log10(dwelltimemeans.ratio));
    set(h, 'AlphaData', dwelltimemeans.numUPs>8);
    hold on
    plot(snline(:,1),snline(:,2),'k','Linewidth',2)
    plot(hopfline(:,1),hopfline(:,2),'k--','Linewidth',1)
    colorbar
    caxis([-2 2])
    axis xy


NiceSave('IWdwellstats',figfolder,'WCadapt')  



%% I-B space

simdatafilename = 'dwellmatch_IB2.mat';
simdatafullfile = fullfile(simdatafolder,simdatafilename);

resim=false;
if ~exist(simdatafullfile,'file') || resim
  % UP/DOWN Duration Map
    numI = 50;
    numB = 50;
    testI = linspace(-5,0,numI);
    testb = linspace(0,3.5,numB);
    testI = repmat(testI,numB,1);testI = testI(:);
    testb = repmat(testb,1,numI);testb = testb(:);
    numsims = length(testI);

    for nn = 1:numsims
        modelparms(nn).N_neurons = 1;
        modelparms(nn).I_in = testI(nn);
        modelparms(nn).W = 6;
        modelparms(nn).beta = testb(nn);
        modelparms(nn).tau_r = 1;
        modelparms(nn).tau_a = 25;
        modelparms(nn).A0 = 0.5;
        modelparms(nn).Ak = 15;
        modelparms(nn).noiseamp =0.25;
        modelparms(nn).noisefreq = 0.05;
    end

    simparms.simtime = 60000;
    simparms.dt = 1;

    [ dwelltimes_sim,ratehist ] = BatchSimulate_WCadapt( modelparms,simparms );
    save(simdatafullfile,'dwelltimes_sim','ratehist','modelparms','simparms');
else 
    load(simdatafullfile);
end

%Extract Parameters
modelparmsvec = CollapseStruct(modelparms);
ratehistvec = CollapseStruct(ratehist);

I_in = unique(modelparmsvec.I_in);
B = unique(modelparmsvec.beta);
meanrate = reshape(ratehistvec.mean,length(I_in),length(B));

%Load Bif'n Diagram
bifnlocation = fullfile(dropboxroot,'Current Projects/SlowOscillation/Figures/IBbifn.dat');
[ IBbifn ] = BifnFromXPP( bifnlocation );
hopfline = [IBbifn(:,1),IBbifn(:,6)];
hopfline2 = [IBbifn(:,1),IBbifn(:,7)];
snline = [IBbifn(:,1),IBbifn(:,14)];

jumpthresh = 0.05;
jumpsize = diff(snline);
snline(jumpsize>jumpthresh)=nan;
jumpsize = diff(hopfline);
hopfline(jumpsize>jumpthresh)=nan;

%% IB: Dwell time stats
clear dwelltimemeans; clear dwelltimestd
numsims = length(dwelltimes_sim);
for nn = 1:numsims
    dwelltimemeans.UP(nn) = nanmean(dwelltimes_sim(nn).UP);
    dwelltimemeans.DOWN(nn) = nanmean(dwelltimes_sim(nn).DOWN);
    
    dwelltimestd.UP(nn) = nanstd(dwelltimes_sim(nn).UP);
    dwelltimestd.DOWN(nn) = nanstd(dwelltimes_sim(nn).DOWN);
    
    dwelltimemeans.numUPs(nn) = length(dwelltimes_sim(nn).UP);
end

%Calculate mean ratio, CVs
dwelltimemeans.ratio = dwelltimemeans.UP./dwelltimemeans.DOWN;
dwelltimeCV.UP = dwelltimestd.UP./ dwelltimemeans.UP;
dwelltimeCV.DOWN = dwelltimestd.DOWN./ dwelltimemeans.DOWN;

%Reshape for image
dwelltimemeans.UP = reshape(dwelltimemeans.UP,length(I_in),length(B));
dwelltimemeans.DOWN = reshape(dwelltimemeans.DOWN,length(I_in),length(B));
dwelltimemeans.numUPs = reshape(dwelltimemeans.numUPs,length(I_in),length(B));
dwelltimeCV.UP = reshape(dwelltimeCV.UP,length(I_in),length(B));
dwelltimeCV.DOWN = reshape(dwelltimeCV.DOWN,length(I_in),length(B));

dwelltimemeans.ratio = dwelltimemeans.UP./dwelltimemeans.DOWN;

%% Figure : IB
ratecolor = makeColorMap([1 1 1],[0.5 0.5 0.8],[0.8 0.5 0.5]);
figure
colormap(ratecolor)
subplot(3,3,9)
    imagesc(I_in,B,meanrate)
    hold on
    plot(snline(:,1),snline(:,2),'k','Linewidth',2)
    plot(hopfline(:,1),hopfline(:,2),'--','color','k','Linewidth',1)
    colorbar
    caxis([0 1])
    %plot(snline(:,1),snline(:,2),'k','Linewidth',2)
    %plot(hopfline(:,1),hopfline(:,2),'--','color',hopfcolor,'Linewidth',1)
    axis xy
    %xlim([-3.5 0]);ylim([0 7.5])
   
    
   

subplot(3,3,3)
colormap(gca,UPcolor)
    h = imagesc(I_in,B,log10(dwelltimemeans.UP));
    set(h, 'AlphaData', dwelltimemeans.numUPs>8);
    hold on
    plot(snline(:,1),snline(:,2),'k','Linewidth',2)
    plot(hopfline(:,1),hopfline(:,2),'k:','Linewidth',1)
    plot(hopfline2(:,1),hopfline2(:,2),'k:','Linewidth',1)
    colorbar
    xlim([-4 0]);ylim([0 3.5])
    axis xy
    caxis([1 3])
    
subplot(3,3,6)
colormap(gca,CVcolor)
    h = imagesc(I_in,B,(dwelltimeCV.UP));
    set(h, 'AlphaData', dwelltimemeans.numUPs>8);
    hold on
    plot(snline(:,1),snline(:,2),'k','Linewidth',2)
    plot(hopfline(:,1),hopfline(:,2),'k:','Linewidth',1)
    plot(hopfline2(:,1),hopfline2(:,2),'k:','Linewidth',1)
    colorbar
    xlim([-4 0]);ylim([0 3.5])
    axis xy
    caxis([0.1 1])

    
 
subplot(3,3,7)
colormap(gca,DOWNcolor)
    h = imagesc(I_in,B,log10(dwelltimemeans.DOWN));
    set(h, 'AlphaData', dwelltimemeans.numUPs>8);
    hold on
    plot(snline(:,1),snline(:,2),'k','Linewidth',2)
    plot(hopfline(:,1),hopfline(:,2),'k:','Linewidth',1)
    plot(hopfline2(:,1),hopfline2(:,2),'k:','Linewidth',1)
    colorbar
    xlim([-4 0]);ylim([0 3.5])
    axis xy
    caxis([1 3])
    
subplot(3,3,8)
colormap(gca,CVcolor)
    h = imagesc(I_in,B,(dwelltimeCV.DOWN));
    set(h, 'AlphaData', dwelltimemeans.numUPs>8);
    hold on
    plot(snline(:,1),snline(:,2),'k','Linewidth',2)
    plot(hopfline(:,1),hopfline(:,2),'k:','Linewidth',1)
    plot(hopfline2(:,1),hopfline2(:,2),'k:','Linewidth',1)
    colorbar
    xlim([-4 0]);ylim([0 3.5])
    colorbar
    axis xy
    caxis([0.1 1])
    
subplot(2,2,1)
colormap(gca,UPDOWNcolor)
    h = imagesc(I_in,B,log10(dwelltimemeans.ratio));
    set(h, 'AlphaData', dwelltimemeans.numUPs>8);
    hold on
    plot(snline(:,1),snline(:,2),'k','Linewidth',2)
    plot(hopfline(:,1),hopfline(:,2),'k:','Linewidth',1)
    plot(hopfline2(:,1),hopfline2(:,2),'k:','Linewidth',1)
    colorbar
    caxis([-2 2])
    xlim([-4 0]);ylim([0 3.5])
    axis xy


NiceSave('IBdwellstats',figfolder,'WCadapt')  




   



%% I-W space: b=0

simdatafilename = 'dwellmatch_IWb0.mat';
simdatafullfile = fullfile(simdatafolder,simdatafilename);

resim=false;
if ~exist(simdatafullfile,'file') || resim
  % UP/DOWN Duration Map
    numI = 10;
    numW = 10;
    testI = linspace(-4,-1,numI);
    testw = linspace(3.5,7.5,numW);
    testI = repmat(testI,numW,1);testI = testI(:);
    testw = repmat(testw,1,numI);testw = testw(:);
    numsims = length(testI);

    for nn = 1:numsims
        modelparms(nn).N_neurons = 1;
        modelparms(nn).I_in = testI(nn);
        modelparms(nn).W = testw(nn);
        modelparms(nn).beta = 0;
        modelparms(nn).tau_r = 1;
        modelparms(nn).tau_a = 20;
        modelparms(nn).A0 = 0.5;
        modelparms(nn).Ak = 15;
        modelparms(nn).noiseamp =0.8;
        modelparms(nn).noisefreq = 0.05;
    end

    simparms.simtime = 30000;
    simparms.dt = 1;

    [ dwelltimes_sim,ratehist ] = BatchSimulate_WCadapt( modelparms,simparms );
    save(simdatafullfile,'dwelltimes_sim','ratehist','modelparms','simparms');
else 
    load(simdatafullfile);
end

%Extract Parameters
modelparmsvec = CollapseStruct(modelparms);
ratehistvec = CollapseStruct(ratehist);

I_in = unique(modelparmsvec.I_in);
W = unique(modelparmsvec.W);
meanrate = reshape(ratehistvec.mean,length(I_in),length(W));

%Load Bif'n Diagram
%bifnlocation = fullfile(simfolder,'Current Projects/SlowOscillation/Modeling/simresults/FICurveNullclines_B0/WIbifnB0.dat' );
% [ WIbifnlines ] = BifnFromXPP( bifnlocation);
%snline = [WIbifnlines(:,1),WIbifnlines(:,5)];
%% IW-B0: Dwell time stats
clear dwelltimemeans; clear dwelltimestd
numsims = length(dwelltimes_sim);
for nn = 1:numsims
    dwelltimemeans.UP(nn) = nanmean(dwelltimes_sim(nn).UP);
    dwelltimemeans.DOWN(nn) = nanmean(dwelltimes_sim(nn).DOWN);
    
    dwelltimestd.UP(nn) = nanstd(dwelltimes_sim(nn).UP);
    dwelltimestd.DOWN(nn) = nanstd(dwelltimes_sim(nn).DOWN);
    
    dwelltimemeans.numUPs(nn) = length(dwelltimes_sim(nn).UP);
end

%Calculate mean ratio, CVs
dwelltimemeans.ratio = dwelltimemeans.UP./dwelltimemeans.DOWN;
dwelltimeCV.UP = dwelltimestd.UP./ dwelltimemeans.UP;
dwelltimeCV.DOWN = dwelltimestd.DOWN./ dwelltimemeans.DOWN;

%Reshape for image
dwelltimemeans.UP = reshape(dwelltimemeans.UP,length(I_in),length(W));
dwelltimemeans.DOWN = reshape(dwelltimemeans.DOWN,length(I_in),length(W));
dwelltimemeans.numUPs = reshape(dwelltimemeans.numUPs,length(I_in),length(W));
dwelltimeCV.UP = reshape(dwelltimeCV.UP,length(I_in),length(W));
dwelltimeCV.DOWN = reshape(dwelltimeCV.DOWN,length(I_in),length(W));

dwelltimemeans.ratio = dwelltimemeans.UP./dwelltimemeans.DOWN;

%% Figure : IW-B0
ratecolor = makeColorMap([1 1 1],[0.5 0.5 0.8],[0.8 0.5 0.5]);
figure
colormap(ratecolor)
subplot(3,3,9)
    imagesc(I_in,W,meanrate)
    hold on
    plot(snline(:,1),snline(:,2),'k','Linewidth',2)
    plot(hopfline(:,1),hopfline(:,2),'--','color','k','Linewidth',1)
    colorbar
    caxis([0 1])
    %plot(snline(:,1),snline(:,2),'k','Linewidth',2)
    %plot(hopfline(:,1),hopfline(:,2),'--','color',hopfcolor,'Linewidth',1)
    axis xy
    %xlim([-3.5 0]);ylim([0 7.5])
   
    
   

subplot(3,3,3)
colormap(gca,UPcolor)
    h = imagesc(I_in,W,log10(dwelltimemeans.UP));
    set(h, 'AlphaData', dwelltimemeans.numUPs>8);
    hold on
    plot(snline(:,1),snline(:,2),'k','Linewidth',2)
    plot(hopfline(:,1),hopfline(:,2),'--','color','k','Linewidth',1)
    colorbar
    axis xy
    caxis([1 3])
    
subplot(3,3,6)
colormap(gca,CVcolor)
    h = imagesc(I_in,W,(dwelltimeCV.UP));
    set(h, 'AlphaData', dwelltimemeans.numUPs>8);
    hold on
    plot(snline(:,1),snline(:,2),'k','Linewidth',2)
    plot(hopfline(:,1),hopfline(:,2),'k--','Linewidth',1)
    colorbar
    axis xy
    caxis([0.1 1])

    
 
subplot(3,3,7)
colormap(gca,DOWNcolor)
    h = imagesc(I_in,W,log10(dwelltimemeans.DOWN));
    set(h, 'AlphaData', dwelltimemeans.numUPs>8);
    hold on
    plot(snline(:,1),snline(:,2),'k','Linewidth',2)
    plot(hopfline(:,1),hopfline(:,2),'k--','Linewidth',1)
    colorbar
    axis xy
    caxis([1 3])
    
subplot(3,3,8)
colormap(gca,CVcolor)
    h = imagesc(I_in,W,(dwelltimeCV.DOWN));
    set(h, 'AlphaData', dwelltimemeans.numUPs>8);
    hold on
    plot(snline(:,1),snline(:,2),'k','Linewidth',2)
    plot(hopfline(:,1),hopfline(:,2),'k--','Linewidth',1)
    colorbar
    axis xy
    caxis([0.1 1])
    
subplot(2,2,1)
colormap(gca,UPDOWNcolor)
    h = imagesc(I_in,W,log10(dwelltimemeans.ratio));
    set(h, 'AlphaData', dwelltimemeans.numUPs>8);
    hold on
    plot(snline(:,1),snline(:,2),'k','Linewidth',2)
    plot(hopfline(:,1),hopfline(:,2),'k--','Linewidth',1)
    colorbar
    caxis([-2.5 2.5])
    axis xy


NiceSave('IWB0dwellstats',figfolder,'WCadapt')  









%% Rescaled
simdatafilename = 'dwellmatch_rescale2.mat';
simdatafullfile = fullfile(simdatafolder,simdatafilename);

resim=false;
if ~exist(simdatafullfile,'file') || resim
  % UP/DOWN Duration Map
    numI = 50;
    numbw = 50;
    testI_rs = linspace(-1,1,numI);
    testbw = linspace(0,6,numbw);
    testI_rs = repmat(testI_rs,numbw,1);testI_rs = testI_rs(:);
    testbw = repmat(testbw,1,numI);testbw = testbw(:);
    numsims = length(testI_rs);
    
    %Resale parameters to "raw" WC model
    w_rs = 2;
    k = 15;
    tau = 20;
    I0 = 0;
    testb_rs = testbw.*w_rs;
    [ testI,w,testb ] = WCadapt_RescaleParms( testI_rs,w_rs,testb_rs,k,I0,tau );

    for nn = 1:numsims
        modelparms(nn).N_neurons = 1;
        modelparms(nn).I_in = testI(nn);
        modelparms(nn).I_rs = testI_rs(nn);
        modelparms(nn).W = w;
        modelparms(nn).W_rs = w_rs;
        modelparms(nn).beta = testb(nn);
        modelparms(nn).b_rs = testb_rs(nn);
        modelparms(nn).tau_r = 1;
        modelparms(nn).tau_a = 20;
        modelparms(nn).A0 = 0.5;
        modelparms(nn).Ak = 15;
        modelparms(nn).noiseamp =0.8;
        modelparms(nn).noisefreq = 0.05;
    end

    simparms.simtime = 50000;
    simparms.dt = 1;

    [ dwelltimes_sim,ratehist ] = BatchSimulate_WCadapt( modelparms,simparms );
    save(simdatafullfile,'dwelltimes_sim','ratehist','modelparms','simparms');
else 
    load(simdatafullfile);
end

%Extract Sim parameters
modelparmsvec = CollapseStruct(modelparms);
ratehistvec = CollapseStruct(ratehist);

I_in = unique(modelparmsvec.I_rs);
BW = unique(modelparmsvec.b_rs./modelparmsvec.W_rs);
meanrate = reshape(ratehistvec.mean,length(I_in),length(BW));

%Load Bifurcation Lines
bifnlocation = fullfile(dropboxroot,'Current Projects/SlowOscillation/Modeling/simresults/RescaledNullclinesBifns/bifn_Ibeta.dat');
[ Ibeta_bifnlines ] = BifnFromXPP(bifnlocation);

%% Rescaled: Dwell time stats

clear dwelltimemeans; clear dwelltimestd; clear dwelltimetotals
numsims = length(dwelltimes_sim);
for nn = 1:numsims
    dwelltimemeans.UP(nn) = nanmean(dwelltimes_sim(nn).UP);
    dwelltimemeans.DOWN(nn) = nanmean(dwelltimes_sim(nn).DOWN);
    
    dwelltimetotals.UP(nn) = nansum(dwelltimes_sim(nn).UP);
    dwelltimetotals.DOWN(nn) = nansum(dwelltimes_sim(nn).DOWN);
    
    dwelltimestd.UP(nn) = nanstd(dwelltimes_sim(nn).UP);
    dwelltimestd.DOWN(nn) = nanstd(dwelltimes_sim(nn).DOWN);
    
    dwelltimemeans.numUPs(nn) = length(dwelltimes_sim(nn).UP);
end

%Calculate mean ratio, CVs
dwelltimeCV.UP = dwelltimestd.UP./ dwelltimemeans.UP;
dwelltimeCV.DOWN = dwelltimestd.DOWN./ dwelltimemeans.DOWN;

%Calculate Reltive durations
dwelltimetotals.pUP = dwelltimetotals.UP./(dwelltimetotals.UP+dwelltimetotals.DOWN);
dwelltimetotals.pDOWN = dwelltimetotals.DOWN./(dwelltimetotals.UP+dwelltimetotals.DOWN);
dwelltimetotals.Prat = dwelltimetotals.pUP./dwelltimetotals.pDOWN;

%Reshape for image
dwelltimemeans.UP = reshape(dwelltimemeans.UP,length(I_in),length(BW));
dwelltimemeans.DOWN = reshape(dwelltimemeans.DOWN,length(I_in),length(BW));
dwelltimemeans.numUPs = reshape(dwelltimemeans.numUPs,length(I_in),length(BW));
dwelltimeCV.UP = reshape(dwelltimeCV.UP,length(I_in),length(BW));
dwelltimeCV.DOWN = reshape(dwelltimeCV.DOWN,length(I_in),length(BW));
dwelltimetotals.Prat = reshape(dwelltimetotals.Prat,length(I_in),length(BW));

dwelltimemeans.ratio = dwelltimemeans.UP./dwelltimemeans.DOWN;

dwelltimemeans.UPDOWNonly = isinf(dwelltimemeans.ratio)-isnan(dwelltimemeans.ratio);

%% Figure: rescaled

figure
subplot(2,2,1)
colormap(gca,UPDOWNcolor)
    h = imagesc(I_in,BW,log10(dwelltimetotals.Prat));
    set(h, 'AlphaData', dwelltimemeans.numUPs>8);
    %h2 = imagesc(I_in,BW,dwelltimemeans.UPDOWNonly);
    %set(h2, 'AlphaData', dwelltimemeans.numUPs<8);
    hold on
        plot(Ibeta_bifnlines(:,1),Ibeta_bifnlines(:,6),'k.')
        hold on
        plot(Ibeta_bifnlines(:,1),Ibeta_bifnlines(:,11),'k--')
        plot(Ibeta_bifnlines(:,1),Ibeta_bifnlines(:,12),'k--')
        xlim([-1 1]);ylim([0 6])
        %LogScale('y',2)
        xlabel('I^*');ylabel('b^*/w^*')
    colorbar
    caxis([-2 2])
    axis xy

subplot(3,3,3)
colormap(gca,UPcolor)
    h = imagesc(I_in,BW,log10(dwelltimemeans.UP));
    set(h, 'AlphaData', dwelltimemeans.numUPs>8);
    hold on
        plot(Ibeta_bifnlines(:,1),Ibeta_bifnlines(:,6),'k.')
        hold on
        plot(Ibeta_bifnlines(:,1),Ibeta_bifnlines(:,11),'k--')
        plot(Ibeta_bifnlines(:,1),Ibeta_bifnlines(:,12),'k--')
        xlim([-1 1]);ylim([0 6])
    colorbar
    axis xy
    caxis([1 3])
    
subplot(3,3,6)
colormap(gca,CVcolor)
    h = imagesc(I_in,BW,(dwelltimeCV.UP));
    set(h, 'AlphaData', dwelltimemeans.numUPs>8);
    hold on
        plot(Ibeta_bifnlines(:,1),Ibeta_bifnlines(:,6),'k.')
        hold on
        plot(Ibeta_bifnlines(:,1),Ibeta_bifnlines(:,11),'k--')
        plot(Ibeta_bifnlines(:,1),Ibeta_bifnlines(:,12),'k--')
        xlim([-1 1]);ylim([0 6])
    colorbar
    axis xy
    caxis([0 1])

    
 
subplot(3,3,7)
colormap(gca,DOWNcolor)
    h = imagesc(I_in,BW,log10(dwelltimemeans.DOWN));
    set(h, 'AlphaData', dwelltimemeans.numUPs>8);
    hold on
        plot(Ibeta_bifnlines(:,1),Ibeta_bifnlines(:,6),'k.')
        hold on
        plot(Ibeta_bifnlines(:,1),Ibeta_bifnlines(:,11),'k--')
        plot(Ibeta_bifnlines(:,1),Ibeta_bifnlines(:,12),'k--')
        xlim([-1 1]);ylim([0 6])
    colorbar
    axis xy
    caxis([1 3])
    
subplot(3,3,8)
colormap(gca,CVcolor)
    h = imagesc(I_in,BW,(dwelltimeCV.DOWN));
    set(h, 'AlphaData', dwelltimemeans.numUPs>8);
    hold on
        plot(Ibeta_bifnlines(:,1),Ibeta_bifnlines(:,6),'k.')
        hold on
        plot(Ibeta_bifnlines(:,1),Ibeta_bifnlines(:,11),'k--')
        plot(Ibeta_bifnlines(:,1),Ibeta_bifnlines(:,12),'k--')
        xlim([-1 1]);ylim([0 6])
    colorbar
    axis xy
    caxis([0 1])

%NiceSave('Resacleddwellstats',figfolder,'WCadapt')  
   


