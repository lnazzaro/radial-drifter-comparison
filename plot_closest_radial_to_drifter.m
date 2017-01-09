close all

printDir='....';

% load bathymetry data
if(site_loc(2)>38&site_loc(2)<40.7)
    load njbathymetry_3sec
    mapoffset=.2;
else
    load bathymetry_USeastcoast
    mapoffset=.5;
end

% define spatial quality limit you want to count as "bad"
espclim=20;

if(~exist([printDir '/' drifterid '/' site '/'],'dir'))
    mkdir([printDir '/' drifterid '/' site '/']);
end

%% image w/ drifter track, histogram of spatial quality, and scatterplot of radial vs. drifter velocity (with correlation stats)

figure

% drifter track
subplot(2,2,1)
% plot coastline and 20m, 50m isobaths
contour(loni,lati,depthi,[0 0],'k');
hold on
contour(loni,lati,depthi,[-20 -50],'color',[.5 .5 .5]);
% set image boundaries
xlim([min([site_loc(1) Lon])-mapoffset max([site_loc(1) Lon])+mapoffset])
ylim([min([site_loc(2) Lat])-mapoffset max([site_loc(2) Lat])+mapoffset])
project_mercator(gca)
% plot codar site location
scatter(site_loc(1),site_loc(2),100,'m','filled','marker','^');
% plot entire drifter track in black
plot(Lon,Lat,'color','k')
% plot drifter track only where radial measurement and drifter track are
% within maxdist of each other, colored by date
ind_dist=find(DIST<=maxdist);
scatter(Lon(ind_dist),Lat(ind_dist),3,time(ind_dist),'filled')
cb=colorbar;
caxis([min(time) max(time)])
dint=ceil(range(time)/8);
set(cb,'ytick',ceil(min(time)):dint:ceil(min(time))+dint*8,...
    'yticklabel',datestr(ceil(min(time)):dint:ceil(min(time))+dint*8,'mm/dd'));
title({['Site: ' site];['Drifter: ' drifterid]})

% histogram of ESPC (spatial quality)
subplot(2,2,3)
hist(ESPC(ESPC<999))
xlabel('ESPC')
ylabel('Frequency')
title({'Spatial Quality Distribution';['No SD Count = ' int2str(sum(ESPC==999))]});

% scatterplot of drifter and radial velocities; legend with correlation
% statistics
subplot(2,2,[2 4])

% scatter all points in blue
ind=find(~isnan(VELO));
scatter(VELO(ind),VeloDrift(ind),'b.')
set(gca,'xtick',-90:15:90,'ytick',-90:15:90)
% define axis limits as +/- maximum velocity, with a 10% buffer both ends
axlim=max(abs([VELO(ind) VeloDrift(ind)]))*1.1;
xlim(axlim*[-1 1])
ylim(axlim*[-1 1])
axis square
hold on
% get best-fit line for all data
mb=polyfit(VELO(ind),VeloDrift(ind),1);
% calculate r (correlation coefficient) and p (p-value of correlation) for
% all data
[r,p]=corrcoef(VELO(ind),VeloDrift(ind));
r=sprintf('%.2f',r(2));
r(1)='r';
p=sprintf('%.2f',p(2));
p(1)='p';
% calculate RMSE for all data
rmse=sqrt(sum((VELO(ind)-VeloDrift(ind)).^2)/length(ind));
rmse=sprintf('%.2f',rmse);

% scatter points where no spatial quality is calculated in red (will cover
% up blue points)
ind=find(ESPC==999);
scatter(VELO(ind),VeloDrift(ind),'r.')
if(isempty(ind))
    scatter(axlim+2,axlim+2,'r.');
end
ind=find(ESPC~=999&~isnan(VELO));
% best-fit line for all data WITH a valid spatial quality measurement
mb999=polyfit(VELO(ind),VeloDrift(ind),1);
% r and p for all data with a valid spatial quality measurement
[r999,p999]=corrcoef(VELO(ind),VeloDrift(ind));
r999=sprintf('%.2f',r999(2));
r999(1)='r';
p999=sprintf('%.2f',p999(2));
p999(1)='p';
% RMSE for all data with a valid spatial quality measurement
rmse999=sqrt(sum((VELO(ind)-VeloDrift(ind)).^2)/length(ind));
rmse999=sprintf('%.2f',rmse999);

% scatter points where spatial quality exists and is more than 'espclim' in
% cyan (will cover up blue points)
ind=find(ESPC~=999&ESPC>=espclim);
scatter(VELO(ind),VeloDrift(ind),'c.')
if(isempty(ind))
    scatter(axlim+2,axlim+2,'c.');
end
% best-fit line for all data with a valid spatial quality measurement LESS
% than 'espclim'
ind=find(ESPC<espclim);
mbespc=polyfit(VELO(ind),VeloDrift(ind),1);
% r and p for all data with a valid spatial quality measurement <espclim
[respc,pespc]=corrcoef(VELO(ind),VeloDrift(ind));
respc=sprintf('%.2f',respc(2));
respc(1)='r';
pespc=sprintf('%.2f',pespc(2));
pespc(1)='p';
% RMSE for all data with valid spatial quality measurement <espclim
rmseespc=sqrt(sum((VELO(ind)-VeloDrift(ind)).^2)/length(ind));
rmseespc=sprintf('%.2f',rmseespc);

% plot best-fit line for all data
y=mb(1)*axlim*[-1 1]+mb(2);
plot(axlim*[-1 1],y,'color',[.5 0 .5]);
% plot best-fit line for all data with valid ESPC
y=mb999(1)*axlim*[-1 1]+mb999(2);
plot(axlim*[-1 1],y,'color',[.5 0 .5],'linestyle','--');
% plot best-fit line for all data with valid ESPC <espclim
y=mbespc(1)*axlim*[-1 1]+mbespc(2);
plot(axlim*[-1 1],y,'color',[.5 0 .5],'linestyle',':','linewidth',1.5);
% plot 1:1 line
plot(axlim*[-1 1],axlim*[-1 1],'k','linewidth',1.25);

xlabel('CODAR Velocity')
ylabel('Drifter Velocity')

% add legend with correlation statistics
legend(['ESPC<' num2str(espclim) ', N=' int2str(sum(ESPC<espclim))],...
    ['ESPC=999 (DNE), N=' int2str(sum(ESPC==999))],...
    ['ESPC>' num2str(espclim) ', N=' int2str(sum(ESPC>=espclim&ESPC~=999))],...
    ['BF All: ' r ',' p ',rmse' rmse],...
    ['BF ESPC~=999: ' r999 ',' p999 ',rmse' rmse999],...
    ['BF ESPC<' num2str(espclim) ': ' respc ',' pespc ',rmse' rmseespc],...
    '1:1',...
    'location','northoutside');

print(figure(1),[printDir '/' drifterid '/' site '/drifter_comp_' name_str '_' site '.png'],'-dpng','-r300');


%% image with time-series of distance between drifter and radial measurement, drifter velocity, and radial velocity, and time-series of spatial and temporal quality

figure

% change 999s for spatial and temporal quality to nans
ETMP(ETMP==999)=nan;
ESPC(ESPC==999)=nan;

% determine spacing of x tick-marks
xtickint=ceil(range(time)/8);

% time series of spatial separation and velocities
subplot(2,1,1)

% plot time-series of separation distance between drifter location and
% radial measurement in shaded red
ind=find(~isnan(DIST));
patch([time(ind) time(ind(end:-1:1)) time(ind(1))],[DIST(ind) -DIST(ind(end:-1:1)) DIST(ind(1))],'r','facealpha',.25)
hold on
% plot radial velocity time-series in blue
plot(time,VELO,'color','b','linewidth',1.5,'marker','.')
% plot drifter velocity (rotated to radial direction) in green
plot(time,VeloDrift,'color','g','linewidth',1.5,'marker','.')
axlim=max([axlim DIST]);
s1=[];
% place dotted vertical lines at 00:00 each day
for t=floor(min(time)):ceil(max(time))
    s1=[s1,plot([t t],axlim*[-1 1],'color',[.5 .5 .5],'linestyle',':')];
end
set(gca,'xtick',ceil(min(time)):xtickint:ceil(max(time)),...
    'xticklabel',datestr(ceil(min(time)):xtickint:ceil(max(time)),'mm/dd'))
xlim([min(time) max(time)])
ylim(axlim*[-1 1])
xlabel('Time')
ylabel({'Radial Vector Component (cm/s)';'Distance (km)'})
legend('Distance to Closest Radial Measurement',[site ' Velocity'],['Drifter ' drifterid ' Velocity'],...
    'location','northoutside','orientation','horizontal');

% time-series of spatial and temporal quality
subplot(2,1,2)

% plot spatial quality (when it exists) in blue
plot(time,ESPC,'color','b','marker','.')
hold on
% plot dots at times where spatial quality measurement doesn't exist at
% y=-1
ind=find(isnan(ESPC));
scatter(time(ind),-ones(size(ind)),10,'c','filled')
% plot temporal quality (when it exists) in red
plot(time,ETMP,'color','r','marker','.')
% plot dots at times where temporal quality measurement doesn't exist at
% y=-1.5
ind=find(isnan(ETMP));
scatter(time(ind),-1.5*ones(size(ind)),10,'m','filled')
% place dotted vertical lines at 00:00 each day
s2=[];
for t=floor(min(time)):ceil(max(time))
    s2=[s2,plot([t t],axlim*[-1 1],'color',[.5 .5 .5],'linestyle',':')];
end
set(gca,'xtick',ceil(min(time)):xtickint:ceil(max(time)),...
    'xticklabel',datestr(ceil(min(time)):xtickint:ceil(max(time)),'mm/dd'))
xlim([min(time) max(time)])
ylim([-2 max([ESPC ETMP])+2])
xlabel('Time')
ylabel({'Radial Quality'; '(Standard Devation, cm/s)'})
legend('Spatial Quality (ESPC)','No Valid ESPC','Temporal Quality (ETMP)','No Valid ETMP',...
    'location','northoutside','orientation','horizontal');

print(figure(2),[printDir '/' drifterid '/' site '/drifter_comp_' name_str '_' site '_time.png'],'-dpng','-r300');


% assign tick mark locations (every 6 hours) and labels (daily) for time
% subsets
tm=floor(min(time)):1/4:ceil(max(time));
tl=cellstr(datestr(tm,'mm/dd'));
tl(mod(tm,1)>0)={''};

delete(s1,s2)
for s=1:2
    subplot(2,1,s)
    set(gca,'xtick',tm,...
        'xticklabel',tl,...
        'xgrid','on');
end

% loop through subsets of time to look at finer detail - add print command 
% at end of while loop to save each subset image
t1=floor(min(time));
while(t1<max(time))
    t2=t1+7;
    if(t1<min(time))
        t1=min(time);
    end
    if(t2+3>max(time))
        t2=max(time);
    end
    for s=1:2
        subplot(2,1,s)
        xlim([t1 t2])
    end
    t1=t2;
    
    print(figure(2),[printDir '/' drifterid '/' site '/drifter_comp_' name_str '_' site '_time_wk.png'],'-dpng','-r300');

end


%% image with scatterplots of difference between drifter and radial velocities vs. ESPC, ETMP, distance between measurements, and drifter distance to codar site

figure

% load hsv colormap (pretty circular) to plot bearing - repeat 1st color at
% end of colormap
hsvmap=colormap('hsv');
hsvmap(end+1,:)=hsvmap(1,:);

% change nans for spatial and temporal quality to -1s so they can be
% plotted with data
ESPC(isnan(ESPC))=-1;
ETMP(isnan(ETMP))=-1;

% calculate difference in velocities between drifter and radial measurement
veldiff=abs(VELO-VeloDrift);

% calculate difference in bearing to site between drifter location and
% radial measurement
xd=(Lon-site_loc(1))*111.12*cosd(site_loc(2));
yd=(Lat-site_loc(2))*111.12;
bear=atan2(yd,xd);
bear=bear*180/pi;
bear=-(bear-90);
bear(bear<0)=bear(bear<0)+360;
beardiff=abs(BEAR-bear);
beardiff(beardiff>180)=beardiff(beardiff>180)-360;
beardiff=abs(beardiff);
% assign ranges for different colors to plot for bearing difference
beardiffcats=[0 10 .8
    10 20 .7
    20 30 .6
    30 50 .5
    50 75 .4
    75 100 .3
    100 125 .2
    125 150 .1
    150 181 0];
bearlegend=[];
for c=1:size(beardiffcats,1)
    bearlegend=[bearlegend,{sprintf('%d-%d',beardiffcats(c,1),beardiffcats(c,2))}];
end

% plot difference in velocity vs. spatial quality (ESPC=nan will be at -1)
subplot(2,2,1)
scatter(ESPC,veldiff,'b.')
ind=find(ESPC>-1);
% calculate and plot best-fit line
mb=polyfit(ESPC(ind),veldiff(ind),1);
hold on
x=[0 max(ESPC)*1.1];
y=mb(1)*x+mb(2);
plot(x,y,'k');
xlabel('Spatial Quality')
ylabel('Velocity Difference')
% calculate r and p and add to title
[r,p]=corrcoef(ESPC(ind),veldiff(ind));
title({[site ' ' type ' vs. ' drifterid];['r=' num2str(r(2),'%0.2f') ', p=' num2str(p(2),'%0.2f')]})
xint=ceil(max(ESPC)*1.1/7);
set(gca,'xlim',[-2 max(ESPC)*1.1],...
    'ylim',[0 max(veldiff)*1.1],...
    'xtick',[-1 xint:xint:xint*6],...
    'xticklabel',[{'NA'},cellstr(int2str([xint:xint:xint*6]'))']);
box on
grid on

% plot difference in velocity vs. temporal quality (ETMP=nan will be at -1)
subplot(2,2,2)
scatter(ETMP,veldiff,'b.')
ind=find(ETMP>-1);
% calculate and plot best-fit line
mb=polyfit(ETMP(ind),veldiff(ind),1);
hold on
x=[0 max(ETMP)*1.1];
y=mb(1)*x+mb(2);
plot(x,y,'k');
xlabel('Temporal Quality')
ylabel('Velocity Difference')
% calculate r and p and add to title
[r,p]=corrcoef(ETMP(ind),veldiff(ind));
title({[site ' ' type ' vs. ' drifterid];['r=' num2str(r(2),'%0.2f') ', p=' num2str(p(2),'%0.2f')]})
xint=ceil(max(ETMP)*1.1/7);
set(gca,'xlim',[-2 max(ETMP)*1.1],...
    'ylim',[0 max(veldiff)*1.1],...
    'xtick',[-1 xint:xint:xint*6],...
    'xticklabel',[{'NA'},cellstr(int2str([xint:xint:xint*6]'))']);
box on
grid on

% plot difference in velocity vs. distance between drifter and radial
% measurement (difference in bearing to site plotted in shades of gray)
subplot(2,2,3)
% loop through each "range" of bearing differences and plot velocity vs.
% distance in appropriate shade of gray
for c=1:size(beardiffcats,1)
    ind=find(beardiff>=beardiffcats(c,1)&beardiff<beardiffcats(c,2));
    scatter(DIST(ind),veldiff(ind),'markeredgecolor',beardiffcats(c,3)*[1 1 1],'markerfacecolor',beardiffcats(c,3)*[1 1 1],'marker','.')
    if(isempty(ind))
        scatter(-5,-5,'markeredgecolor',beardiffcats(c,3)*[1 1 1],'markerfacecolor',beardiffcats(c,3)*[1 1 1],'marker','.')
    end
    hold on
end
% calculate and plot best-fit line
ind=find(~isnan(veldiff));
mb=polyfit(DIST(ind),veldiff(ind),1);
hold on
xmax=min([max(DIST)*1.1 maxdist]);
x=[0 xmax];
y=mb(1)*x+mb(2);
plot(x,y,'k');
xlabel({'Distance to Radial Measurement';'Color: Bearing Difference'})
ylabel('Velocity Difference')
% calculate r and p and add to title
[r,p]=corrcoef(DIST(ind),veldiff(ind));
title({[site ' ' type ' vs. ' drifterid];['r=' num2str(r(2),'%0.2f') ', p=' num2str(p(2),'%0.2f')]})
xint=ceil(xmax/8);
if(xmax<8)
    xint=ceil(xmax*10/8)/10;
end
set(gca,'xlim',[0 xmax],...
    'ylim',[0 max(veldiff)*1.3],...
    'xtick',[0:xint:xint*7]);
ind=find(beardiffcats(:,1)<=max(beardiff));
legend(bearlegend(ind),'location','north','orientation','horizontal')
box on
grid on

% plot difference in velocity vs. drifter distance to site, colored by
% bearing at closest radial measurement
subplot(2,2,4)
DIST_site=dist_ref(site_loc(1),site_loc(2),Lon,Lat);
scatter(DIST_site,veldiff,10,BEAR,'filled')
% calculate and plot best-fit line
ind=find(~isnan(veldiff));
mb=polyfit(DIST_site(ind),veldiff(ind),1);
hold on
x=[0 max(DIST_site)*1.1];
y=mb(1)*x+mb(2);
plot(x,y,'k');
xlabel('Distance to Site')
ylabel('Velocity Difference')
% calculate r and p and add to title
[r,p]=corrcoef(DIST_site(ind),veldiff(ind));
title({[site ' ' type ' vs. ' drifterid];['r=' num2str(r(2),'%0.2f') ', p=' num2str(p(2),'%0.2f')]})
xint=ceil(max(DIST_site)*1.1/8);
set(gca,'xlim',[0 max(DIST_site)*1.1],...
    'ylim',[0 max(veldiff)*1.1],...
    'xtick',[0:xint:xint*7]);
colorbar
colormap(hsvmap)
ylabel(colorbar,'Bearing')
caxis([0 360])
box on
grid on

print(figure(3),[printDir '/' drifterid '/' site '/drifter_comp_' name_str '_' site '_diffscatter.png'],'-dpng','-r300');

