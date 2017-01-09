clear all
close all

drifterid='43411';

site='BRNT';
type='RDLm'; % RDLm for measured, RDLi for ideal
subdir=true; % are all radial files in one directory, or subdirectories .../yyyy_mm/...
cut_dist=true; % do you want to only include radial measurements within some distance of drifter location? (6km for 5MHz, 2.5km for 13&25MHz)

rdlDir='/home/codaradm/data/radials/'; % directory containing radial files
rext='.ruv'; % file extension

% make sure HFR Progs toolbox is in matlab path
addpath(genpath('/home/codaradm/HFR_Progs-2_1_3beta/matlab/'))

% which fields to match from radial file (if don't exist, will be NaN)
savetypes={'LOND','LATD','VELU','VELV','VFLG','ESPC','ETMP',...
    'RNGE','BEAR','VELO','HEAD','SPRC',...
    'GRNG','LRNG','SPKE','TRND','STCK','GRAD'};

% which fields to switch to NaN if closest radial measurement to drifter is
% too far away (if cut_dist is set to true)
adjustdist={'VELU','VELV','ESPC','ETMP',...
    'RNGE','BEAR','VELO','HEAD','SPRC',...
    'GRNG','LRNG','SPKE','TRND','STCK','GRAD',...
    'U','V','VeloDrift'};

% list long-range sites
LR_sites={'AMAG','ASSA','BLCK','BRIG','CEDR','CORE','CStM','DUCK','FARO',...
    'GMNB','GrnI','HATY','HEMP','HOOK','LISL','LOVE','MABO','MRCH','MVCO',...
    'NANT','NAUS','PYFC','WILD'};

% list medium-range sites
MR_sites={'BRAD','BRMR','BRNT','CDDO','FURA','RATH','SEAB','SPRK','WOOD'};

% list short-range sites
SR_sites={'BISL','CMPT','CPHN','GCAP','HLPN','MISQ','MNTK','OLDB','PORT',...
    'SILD','STLI','SUNS','VIEW'};

% load drifter data - add path if not saved in working directory
load(['USCG_' drifterid '_2016_05_10.mat'])

% drifter 43372 flat-lines after May 31st - remove that data
if(strcmp(drifterid,'43372'))
    time(time>datenum(2016,5,31))=nan;
end

% initialize variables (set to same size as time variable from drifter)
% being matched from radial data
for n=1:length(savetypes)
    eval([savetypes{n} '=nan(size(time));']);
end

DIST=nan(size(time)); % initialize distance variable (how far is closest radial to drifter location
VeloDrift=nan(size(time)); % initialize drifter-velocity variable, rotated to radial direction


for n=1:length(time) % loop through drifter timestamps
    t=time(n);
    if(isnan(t)) % skip if timestamp is nan
        continue;
    end
    
    % load radial file
    if(subdir)
        rdl=loadRDLFile([rdlDir,site,'/',datestr(t,'yyyy_mm/'),...
            type '_' site '_' datestr(t,'yyyy_mm_dd_HH00') rext],true);
    else
        rdl=loadRDLFile([rdlDir,site,'/',...
            type '_' site '_' datestr(t,'yyyy_mm_dd_HH00') rext],true);
    end
    
    if(~isempty(rdl.LonLat))
        site_loc=rdl.SiteOrigin; % location of CODAR site
        
        % find index of table-start row in header
        tablestart_ind=find(strncmp('%TableStart',rdl.OtherMetadata.Header,length('%TableStart')));
        tablestart_ind=tablestart_ind(1);
        
        % divide into pre-data and post-data subsets of header
        subset1=rdl.OtherMetadata.Header(1:tablestart_ind+2);
        subset3=rdl.OtherMetadata.Header(tablestart_ind+3:end);
        
        % get number of columns (in header description)
        colnum_ind=find(strncmp('%TableColumns',subset1,length('%TableColumns')));
        colnum=str2num(subset1{colnum_ind}(length('%TableColumns:  '):end));
        
        % get 4-char variable IDs
        coltype_ind=strncmp('%TableColumnTypes',subset1,length('%TableColumnTypes'));
        coltype=subset1{coltype_ind};
        ind=find(coltype==':');
        coltype=coltype(ind+1:end);
        coltype(coltype==' ')='';
        if(length(coltype)~=colnum*4)
            warning('column type labels do not match up with column count')
        end
        
        % get long-name header to each variable
        header=cell(1,colnum);
        for k=1:colnum
            header{k}=coltype(k*4-3:k*4);
        end
        
        % remove any data flagged in VFLG (mostly/all AngSeg
        % flags)
        indvflg=find(strcmp(header,'VFLG'));
        indvflg_good=find(rdl.OtherMetadata.RawData(:,indvflg)==0);
        rdl.OtherMetadata.RawData=rdl.OtherMetadata.RawData(indvflg_good,:);
        
        % find closest lon/lat radial measurement to current
        % drifter location
        indlon=find(strcmp(header,'LOND'));
        indlat=find(strcmp(header,'LATD'));
        d=dist_ref(Lon(n),Lat(n),rdl.OtherMetadata.RawData(:,indlon),rdl.OtherMetadata.RawData(:,indlat));
        indd=find(d==min(d));
        DIST(n)=d(indd);
        
        % loop through all variabless in 'savetypes' and assign
        % radial data closest to drifter location
        for k=1:length(savetypes)
            ind=find(strcmp(header,savetypes{k})); % find which column has data for that variable
            if(~isempty(ind))
                eval([savetypes{k} '(n)=rdl.OtherMetadata.RawData(indd,ind);']);
            end
        end
        
        % rotate drifter velocity according to bearing at
        % closest radial measurement
        [~,VeloDrift(n)]=rot(U(n),V(n),BEAR(n));
        
        clear rdl header col* ind* subset* table* d k t
    end
end

% reverse rotated drifter velocity to match direction to/from
% CODAR site
VeloDrift=-VeloDrift;

clear n

% if cut_dist is true, set max separation between drifter location and
% radial measurement to 6km for LR sites or 2.5km for MR and SR sites; if
% cut_dist is false, set max separation to 999km
if(cut_dist)
    if(ismember(site,LR_sites))
        maxdist=6;
    elseif(ismember(site,MR_sites))
        maxdist=2.5;
    elseif(ismember(site,SR_sites))
        maxdist=2.5;
    end
    cut_dist_lab=['_maxseparation' num2str(maxdist) 'km'];
else
    cut_dist_lab='';
    maxdist=999;
end

% find any radial matches that are separated from the drifter track by more
% than maxdist, and replace matching data with NaNs.
ind_dist=find(DIST>maxdist);
for nd=1:length(adjustdist)
    eval([adjustdist{nd} '(ind_dist)=nan;']);
end
