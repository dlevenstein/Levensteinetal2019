analysisfolder = '/Users/dlevenstein/Dropbox/Research/Current Projects/SlowOscillation/AnalysisScripts/AnalysisFigs/SPWRDurationAnalysis';
DwellMatchAll = GetMatResults(analysisfolder,'_SPWRDurationAnalysis');
figfolder = '/Users/dlevenstein/Dropbox/Research/Current Projects/SlowOscillation/AnalysisScripts/AnalysisFigs/SPWRDurationAnalysis';
dwellhistbins = DwellMatchAll(1).dwellhistbins;
numrecs = length(DwellMatchAll);
DwellMatchAll = CollapseStruct(DwellMatchAll);
%%
analysisfolder = '/Users/dlevenstein/Dropbox/Research/Current Projects/SlowOscillation/AnalysisScripts/AnalysisFigs/SPWRDurationStatsAnalysis';
DwellStatsAll = GetMatResults(analysisfolder,'_SPWRDurationStatsAnalysis');
DwellStatsAll = CollapseStruct(DwellStatsAll,1,'justcat',true);

%%
figure
subplot(2,4,1)
    [DwellStatsAll.CVmean,DwellStatsAll.CVstd] = BoxAndScatterPlot(...
        {DwellStatsAll.dwelltimestats.CV.SWR,DwellStatsAll.dwelltimestats.CV.iSWR},...
        'colors',[[0.8 0 0];[0 0 0.8]],'labels',{'SWR','iSWR'});
    box off
    ylabel('CV')
    ylim([0 1.5])
    
subplot(2,4,2)
    [DwellStatsAll.ratmean,DwellStatsAll.ratstd] = BoxAndScatterPlot(...
        {log10(1./DwellStatsAll.dwelltimestats.meanrat),log10(1./DwellStatsAll.dwelltimestats.Prat)},...
        'labels',{'mean','P(state)'});
    hold on
    plot([0.5 2.5],[0 0],'k--')
    ylim([-2 2])
    LogScale('y',10)
    box off
    ylabel('iSWR:SWR')
    
    
subplot(2,4,3)
    BoxAndScatterPlot(...
        {log10(DwellStatsAll.dwelltimestats.mean.SWR),log10(DwellStatsAll.dwelltimestats.mean.iSWR)},...
        'colors',[[0.8 0 0];[0 0 0.8]],'labels',{'SWR','iSWR'});
    ylim([-1.5 0.5])
    LogScale('y',10)
    box off
    ylabel('Mean (s)')
    
subplot(2,4,4)
    [DwellStatsAll.meanmean,DwellStatsAll.meanstd] = BoxAndScatterPlot(...
        {(DwellStatsAll.dwelltimestats.mean.SWR),(DwellStatsAll.dwelltimestats.mean.iSWR)},...
        'colors',[[0.8 0 0];[0 0 0.8]],'labels',{'SWR','iSWR'});
    %ylim([-1 0.7])
    %LogScale('y',10)
    box off
    ylabel('Mean (s)')
    
 NiceSave('HPCDurStats',figfolder,[]) 

%% Get Simulation data
dropboxroot = '/Users/dlevenstein/Dropbox/Research/';
%dropboxroot = '/mnt/data1/Dropbox/research/';
simdatafolder =  fullfile(dropboxroot,'Current Projects/SlowOscillation/AnalysisScripts/simdata/');
simdatafilename = 'dwellmatch6.mat';
simdatafullfile = fullfile(simdatafolder,simdatafilename);
load(simdatafullfile);

%%
modelparmsvec = CollapseStruct(modelparms);
ratehistvec = CollapseStruct(ratehist);

I_in = unique(modelparmsvec.I_in);
W = unique(modelparmsvec.W);
meanrate = reshape(ratehistvec.mean,length(I_in),length(W));

% Dwell time stats
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


%% Conditions
CVrange = [DwellStatsAll.CVmean' - 2.*DwellStatsAll.CVstd' ...
    DwellStatsAll.CVmean' + 2.*DwellStatsAll.CVstd'];

%log10(DOWN:UP)
ratiorange = [DwellStatsAll.ratmean(1)' - 2.*DwellStatsAll.ratstd(1)' ...
    DwellStatsAll.ratmean(1)' + 2.*DwellStatsAll.ratstd(1)'];


fitconditions = dwelltimemeans.numUPs>8 &...
    dwelltimeCV.UP > CVrange(1,1) &  dwelltimeCV.UP < CVrange(1,2) &...
    dwelltimeCV.DOWN > CVrange(2,1) &  dwelltimeCV.DOWN < CVrange(2,2) &...
    log10(1./dwelltimemeans.ratio) > ratiorange(1) & log10(1./dwelltimemeans.ratio) < ratiorange(2) ;

B = bwboundaries(fitconditions);

figure
imagesc(fitconditions)
hold on
for k = 1:length(B)
   boundary = B{k};
   plot(boundary(:,2), boundary(:,1), 'k', 'LineWidth', 1)
end
axis xy
 
 
%%

modelparms = DwellMatchAll.modelparmsvec(1);
DwellFits = CollapseStruct(DwellMatchAll.dwellfitvec_all,1);
bestfitparms = DwellMatchAll.bestparmfit;
DwellTimes = CollapseStruct(DwellMatchAll.dwelltimes,1);

%% Calculate good fit over all recordings
%nanzeroselectionparm = DwellFits.selectionparm;
%nanzeroselectionparm(isnan(nanzeroselectionparm))=0;
meansimilarity = nanmean(DwellFits.selectionparm,1);
meantimescale = nanmean(DwellFits.bestsf,1);
bestfits_I = modelparms.I_in(bestfitparms);
bestfits_W = modelparms.W(bestfitparms);

taur_idx = sub2ind(size(DwellFits.selectionparm), [1:numrecs], bestfitparms);
bestfits_taur = DwellFits.bestsf(taur_idx);
bestfits_sim = DwellFits.selectionparm(taur_idx);
[~,sortsel] = sort(bestfits_sim);

[~,allbestfit] = max(meansimilarity);
allbest_I = modelparms.I_in(allbestfit);
allbest_W = modelparms.W(allbestfit);
allbest_sf = meantimescale(allbestfit);

%% Dwell Times Statistics

DwellTimes.UPhist_norm = DwellTimes.UPhist./diff(dwellhistbins(1:2));
DwellTimes.DOWNhist_norm = DwellTimes.DOWNhist./diff(dwellhistbins(1:2));

UPhistmean = mean(DwellTimes.UPhist_norm,1);
DOWNhistmean = mean(DwellTimes.DOWNhist_norm,1);
UPhiststd = std(DwellTimes.UPhist_norm,[],1);
DOWNhiststd = std(DwellTimes.DOWNhist_norm,[],1);

%% For Plotting Birfurcation Diagram
IWbifn = importdata('/Users/dlevenstein/Dropbox/Research/Current Projects/SlowOscillation/Figures/IWbifn.dat');
jumpthresh = 0.05;
snline = [IWbifn(:,1) IWbifn(:,3)];
snline(IWbifn(:,5)~=3,:) =nan;
jumpsize = diff(snline);
snline(jumpsize>jumpthresh)=nan;

hopfline = [IWbifn(:,1) IWbifn(:,3)];
hopfline(IWbifn(:,5)~=2,:) =nan;
jumpsize = diff(hopfline);
hopfline(jumpsize>jumpthresh)=nan;

%% Simulate best fit

simtime = 100000;
dt = 1;
modelparms_sim = structfun(@mean,modelparms,'UniformOutput',false);
modelparms_sim.I_in = allbest_I;
modelparms_sim.W = allbest_W;

[ T, Y_sol ] = WCadapt_run(simtime,dt,modelparms_sim);
T = T.*allbest_sf;

%%
simdwellhistbins = linspace(dwellhistbins(1) ,dwellhistbins(end),30);
[thresh,cross,~] = BimodalThresh(Y_sol(:,1),'Schmidt');
simUPdwell = diff(cross.upints,1,2).*allbest_sf;
simDOWNdwell = diff(cross.downints,1,2).*allbest_sf;
simUPdwellhist = hist(log10(simUPdwell),simdwellhistbins);
simUPdwellhist = simUPdwellhist./sum(simUPdwellhist);
simUPdwellhist = simUPdwellhist./diff(simdwellhistbins(1:2));
simDOWNdwellhist = hist(log10(simDOWNdwell),simdwellhistbins);
simDOWNdwellhist = simDOWNdwellhist./sum(simDOWNdwellhist);
simDOWNdwellhist = simDOWNdwellhist./diff(simdwellhistbins(1:2));


%%
%DwellStatsAll.dwelltimestats.meanrat)
simratio = log10((mean(simDOWNdwell)./mean(simUPdwell)))
simUPCV = std(simUPdwell)./mean(simUPdwell)
simDOWNCV = std(simDOWNdwell)./mean(simDOWNdwell)


%% Figure
msize = 60;
hopfcolor = [0 0 0];
fitcolormap = [1 1 1; makeColorMap(1*[1 1 1],[0 0.6 0],[0.8 0.5 0])];
DOWNcolormap = makeColorMap([1 1 1],[0 0 0.8]);
UPcolormap = makeColorMap([1 1 1],[0.8 0 0]);
taucolormap = [1 1 1; makeColorMap([0 0 0.8],[0.8 0 0])];

figure

    subplot(6,2,2)
        colormap(gca,UPcolormap)
        imagesc(dwellhistbins,[1 numrecs],DwellTimes.UPhist(sortsel,:))
    subplot(6,2,4)
        colormap(gca,DOWNcolormap)
        imagesc(dwellhistbins,[1 numrecs],DwellTimes.DOWNhist(sortsel,:))
     subplot(6,2,6)
       bar(simdwellhistbins,simUPdwellhist,'facecolor','none','edgecolor','r','linewidth',1)
        hold on
        bar(simdwellhistbins,simDOWNdwellhist,'facecolor','none','edgecolor','b','linewidth',1)
        errorshade(dwellhistbins,DOWNhistmean,DOWNhiststd,DOWNhiststd,'b','scalar')
        errorshade(dwellhistbins,UPhistmean,UPhiststd,UPhiststd,'r','scalar')
        plot(dwellhistbins,DOWNhistmean,'b','linewidth',3)
        plot(dwellhistbins,UPhistmean,'r','linewidth',3)
        axis tight
        box off
       % ylim([0 0.12])
    subplot(6,2,7)
        histogram(bestfits_taur,3,'facecolor','r')
        xlim([0 0.025])
    subplot(3,10,21)
       % histogram(bestfits_sim,7,'facecolor','r')

           [meansim,stdsim] = BoxAndScatterPlot(...
        {bestfits_sim},'labels',{'Sim.'});
     ylim([0.5 1])
     box off
        xlim([0.8 1.2])
        
        
    subplot(2,2,1)
    I_vals = unique(modelparms.I_in);w_vals=unique(modelparms.W);
    I_vals = I_vals-0.5.*diff(I_vals(1:2));
    w_vals = w_vals-0.5.*diff(w_vals(1:2));
     colormap(gca,fitcolormap)
         hold on
        %scatter(modelparms.I_in,modelparms.W,msize,meansimilarity,'s','filled','MarkerEdgeColor',0.7.*[1 1 1],'linewidth',0.001)
        h =pcolor(I_vals,w_vals,reshape(meansimilarity,sqrt(length(meansimilarity)),sqrt(length(meansimilarity))))
        set(h,'EdgeColor',0.5.*[1 1 1])
        %imagesc(unique(modelparms.I_in),unique(modelparms.W),reshape(meansimilarity,sqrt(length(meansimilarity)),sqrt(length(meansimilarity))))
        plot(bestfits_I,bestfits_W,'r*','markersize',1)
        plot(allbest_I,allbest_W,'ro','markersize',5)
        colorbar
        box on
        caxis([0 0.8])
        axis xy
        plot(snline(:,1),snline(:,2),'k','Linewidth',1)
        plot(hopfline(:,1),hopfline(:,2),'--','color',hopfcolor,'Linewidth',1)
        %legend('Saddle-Node','Hopf','location','southwest')
        xlim([-3.3 -1]);ylim([3.5 7.5])
        set(gca,'ytick',get(gca,'ylim'))
        set(gca,'xtick',get(gca,'xlim'))
        ylabel('Recurrent Weight, W');xlabel('Tonic Input, I')
        
        
%     subplot(2,2,4)
%         plot(bestfits_I,bestfits_taur,'.')
     subplot(2,2,4)
     colormap(gca,taucolormap)
     plotsf =meantimescale;plotsf(meansimilarity<0.4)=0;
         hold on
         %imagesc(unique(modelparms.I_in),unique(modelparms.W),reshape(plotsf,sqrt(length(plotsf)),sqrt(length(plotsf))))
         h = pcolor(I_vals,w_vals,reshape(plotsf,sqrt(length(plotsf)),sqrt(length(plotsf))))
                 set(h,'EdgeColor',0.5.*[1 1 1])
        plot(allbest_I,allbest_W,'ro','markersize',5)
        colorbar
        caxis([0 0.02])
        axis xy
        plot(snline(:,1),snline(:,2),'k','Linewidth',1)
        plot(hopfline(:,1),hopfline(:,2),'--','color',hopfcolor,'Linewidth',1)
        %legend('Saddle-Node','Hopf','location','southwest')
        xlim([-3.3 -1]);ylim([3.5 7.5])
        set(gca,'ytick',get(gca,'ylim'))
        set(gca,'xtick',get(gca,'xlim'))
        ylabel('Recurrent Weight, W');xlabel('Tonic Input, I')



NiceSave('SPWDwellTimeMatch',figfolder,[])

%%
figure
    % subplot(6,2,6)
       bar(simdwellhistbins,simUPdwellhist,'facecolor','none','edgecolor','r','linewidth',1)
        hold on
        bar(simdwellhistbins,simDOWNdwellhist,'facecolor','none','edgecolor','b','linewidth',1)
        errorshade(dwellhistbins,DOWNhistmean,DOWNhiststd,DOWNhiststd,'b','scalar')
        errorshade(dwellhistbins,UPhistmean,UPhiststd,UPhiststd,'r','scalar')
        plot(dwellhistbins,DOWNhistmean,'b','linewidth',3)
        plot(dwellhistbins,UPhistmean,'r','linewidth',3)
        axis tight
        box off
        
NiceSave('SPWDistribution',figfolder,[])

%%
figure
subplot(3,1,1)
plot(T,Y_sol(:,1),'k','linewidth',2)
hold on
plot(T,Y_sol(:,2)-1.2,'k','linewidth',1)
xlim([5 25])
box off
        y = get(gca,'Ylim');
        plot(cross.upints'.*allbest_sf,ones(size(cross.upints))','r','linewidth',5)
%         patch([cross.downints(:,1) cross.downints(:,2) cross.downints(:,2) cross.downints(:,1)]'.*allbest_sf,...
%             repmat([y(1) y(1) y(2) y(2)],length(cross.downints(:,1)),1)',...
%             'g','FaceAlpha',0.2,'EdgeColor','none');

NiceSave('SPW_BestFitExample',figfolder,[])


%%
simfolder = '/Users/dlevenstein/Dropbox/Research/Current Projects/SlowOscillation/Modeling/simresults/';
%'bestfitnull.dat'

[ null1,null2 ] = NullclinesFromXPP(fullfile(simfolder,'HPCbestnullclines.dat'));
%%
DOWNcolor = [0 0 0.8];
UPcolor = [0.8 0 0];
figure
subplot(3,1,1)
plot(T,Y_sol(:,1),'k','linewidth',2)
hold on
plot(T,Y_sol(:,2)-1.2,'k','linewidth',1)
xlim([5 20])
box off
        y = get(gca,'Ylim');
        patch([cross.upints(:,1) cross.upints(:,2) cross.upints(:,2) cross.upints(:,1)]'.*allbest_sf,...
            repmat([y(1) y(1) y(2) y(2)],length(cross.upints(:,1)),1)',...
            UPcolor,'FaceAlpha',0.2,'EdgeColor','none');

    subplot(4,4,9)
        plot(null2(:,1),null2(:,2),'.','color',DOWNcolor,'markersize',7)
        hold on
        plot(null1(:,1),null1(:,2),'k.','markersize',7)
        plot([0 0],get(gca,'ylim'),'k');
        plot(get(gca,'xlim'),[0 0],'k');
        xlabel('Pop. Rate');ylabel('Adapt')
        box on
        xlim([0 1]);ylim([0 1])
        
        
        
NiceSave('HPCBestFitExample',figfolder,[])

%% VIDEO
savevid = fullfile(figfolder,['HPCbestfitvid']);
bestfitvid = VideoWriter(savevid);
bestfitvid.FrameRate = 125;
open(bestfitvid);

xwin = [1800 6265];
tail =50;
tt = xwin(1);



vidfig = figure;
set(vidfig,'color','w');
    subplot(3,3,5)
        plot(null2(:,1),null2(:,2),'.','color',DOWNcolor,'markersize',7)
        hold on
        plot(null1(:,1),null1(:,2),'k.','markersize',7)
        point = plot(Y_sol(tt,1),Y_sol(tt,2),'k.','markersize',20);
        pointtail = plot(Y_sol(max(tt-tail,xwin(1)):tt,1),Y_sol(max(tt-tail,xwin(1)):tt,2),'color',[0.7 0.7 0.7],'markersize',20);
        plot([0 0],get(gca,'ylim'),'k');
        set(gca,'xtick',[0 1]);set(gca,'ytick',[0 1])
        plot(get(gca,'xlim'),[0 0],'k');
        xlabel('Pop. Rate');ylabel('Adaptation')
        box on

for tt = xwin(1)+1:1:xwin(2)
    
subplot(3,1,1)
plot(T(xwin(1):tt),Y_sol(xwin(1):tt,1),'k','linewidth',2)
hold on
plot(T(xwin(1):tt),Y_sol(xwin(1):tt,2)-1.2,'k','linewidth',1)
xlim(T(xwin));ylim([-1.2 1])
box off
        set(gca,'XTick',[])
        set(gca,'YTick',[-0.5 0.5])
        set(gca,'YTickLabels',{'Adaptation','Pop. Rate'})
        
        hold off
        box off
        drawnow
    
subplot(3,3,5)
    delete(point);delete(pointtail);
    point = plot(Y_sol(tt,1),Y_sol(tt,2),'k.','markersize',20);
    pointtail = plot(Y_sol(max(tt-tail,xwin(1)):tt,1),Y_sol(max(tt-tail,xwin(1)):tt,2),'color',[0.7 0.7 0.7],'markersize',20);
    drawnow
    

    
       imgFrame = getframe(vidfig);
       writeVideo(bestfitvid,imgFrame.cdata);
    
end
close(bestfitvid);


