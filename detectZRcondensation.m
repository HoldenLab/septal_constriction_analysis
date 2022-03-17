function [condensationState switchPt muInfo]= detectZRcondensation(widthSeries, varargin)

minSpatialDist = 50;
maxNumChanges=1;
stat='std';
for ii=1:numel(varargin)
    if strcmp(varargin{ii},'MinSpatialDist')
        minSpatialDist=varargin{ii+1};
        ii=ii+2;
    elseif strcmp(varargin{ii},'Statistic')
        stat=varargin{ii+1};
        ii=ii+2;
    else
        ii=ii+1;
    end
end

switchPt =findchangepts(widthSeries,'MaxNumChanges',1,'Statistic',stat);
%order the states (lowest is most condensed, implicitly assumes condensed state always detected )

condensationState = ones(size(widthSeries));
didCondense=false;
if ~isempty(switchPt)
    mu1=mean(widthSeries(1:switchPt));
    mu2=mean(widthSeries(switchPt:end));
    if abs(mu1-mu2)>minSpatialDist
        muInfo=[mu1, mu2, abs(mu1-mu2)];
        if mu1>=mu2
            didCondense=true;
            condensationState(1:switchPt)=2;
        else
            %enforce temporal ordering decondensed-> condensed
            %note need to treat any decondensation phenotypes as special case
            didCondense=false;
        end
    else
        %two states separated by <minSpatialDist
        didCondense=false;
    end
else
    %only one state detected
    didCondense=false;
end

if ~didCondense
    switchPt=-1;
    mu=mean(widthSeries);
    muInfo=[mu, mu, 0];
end

