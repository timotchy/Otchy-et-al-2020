%Script to create figure for paper

figure(7054557)
%paths = paths_0514to0510;
%paths = paths_0531to0526;
%Show baseline template to orient the figure
subplot(4,1,1)
imagesc(-1*template)
axis xy, axis tight
set(gca,'XTick',[],'YTick',[])
set(gca,'TickDir','out')

%Show a gradient plot of the stretch over the whole figure
subplot(4,1,2)
%imagesc(diff(paths.pathMean),[-1,1]) %the is the mean path that results from warping the exp day to the baseline day
spread = (basepoints(2)*44):(basepoints(end-1)*44);
HVClinPath = (LLDTWpath(spread,2));
%imagesc(diff(HVClinPath-(1:length(HVClinPath))),[-1,1])
imagesc((diff((HVClinPath')-(1:length(HVClinPath)))),[-1,1])
%colormap('bone')
set(gca,'XTick',[],'YTick',[])
set(gca,'TickDir','out')

%Show the HVC power traces, warped to baseline day, overlay
subplot(4,1,3)
hold on
% trim1 = mNeuro_0510CH3(44:(end-44));
% trim2 = mNeuro_0514CH3(44:(end-44));
% trim1 = mNeuro_0526CH3(44:(end-44));
% trim2 = mNeuro_0531CH3(44:(end-44));
%  trim1 = smooth(baseline_mNeuro(44:(end-44)),220);
%  trim2 = smooth(future_mNeuro_LLDTWWNB(44:(end-44)),220);

% trim1 = smooth(baseline_mNeuro(44:(end-44)),88)';
% trim2 = smooth(future_mNeuro_LLDTWWNB(44:(end-44)),88)';

trim1 = smooth(baseline_mNeuro(spread),88)'; %2ms smoothing
trim2 = smooth(future_mNeuro_LLDTW(spread),88)';

timeVect = (1:length(trim1))/44;
plot(timeVect,trim1)
plot(timeVect,trim2,'r')
axis tight
set(gca,'XTick',[])
hold off
ylabel('Power (V^2)')
set(gca,'TickDir','out')

%Show the running pearson correlation for the two traces, smoothed
subplot(4,1,4)
rho = LocalPearson(trim1,trim2);
rho_sm = interp1(1:length(rho),rho,linspace(1,length(rho),(basepoints(end-1)-basepoints(2)+1)));
plot(rho_sm)
axis tight
ylim([0,1])
xlabel('Template Time (ms)')
ylabel('Correlation')
set(gca,'TickDir','out')

%Get Pearson values for the segments
target = (basepoints(6):basepoints(8))-basepoints(2);
preTarget = 1:(target(1)-100);
postTarget = target(end)+100:length(rho_sm);


% target = templatesyllBreaks(3,1):templatesyllBreaks(4,1);
% preTarget = (templatesyllBreaks(3,1):templatesyllBreaks(4,1))-length(target);
% postTarget = (templatesyllBreaks(3,1):templatesyllBreaks(4,1))+length(target);
% if postTarget(end)>length(rho_sm)
%     postTarget = templatesyllBreaks(4,2):length(rho_sm);
% end

Rtarget = mean(rho_sm(target))
Rpre = mean(rho_sm(preTarget))
Rpost = mean(rho_sm(postTarget))


