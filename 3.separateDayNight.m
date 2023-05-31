function [index_sunrad_day,index_sunrad_night,index_day,index_night,index_rad_day,index_rad_night,index_sun_day,index_sun_night] = separateDayNight(Timestamp,~,RGlobal,Latitude,Longitude,time_offset)


% function [dayData, nightData] = separateDayNight(Timestamp,InputData,RGlobal,Latitude,Longitude,time_offset)
%
% The function separateDayNight will separate the input data into daytime data and nighttime data
% according to incoming shortwave radiation and calculated times for sunrise and sunset.
%
% Input Data:
% Timestamp	= timestamp vector for each point in the input data in local time
% InputData	= input data vector for which seperation will be done
% RGlobal	= total incoming shortwave radiation vector (global radiation) for each input time period
% Latitude	= scalar, latitude of the location of measurements in degrees, minutes should be
% 		  represented as a decimal value, N=positive, S=negative.
% Longitude	= scalar, longitude of the location of measurements in degrees, minutes should be
% 		  represented as a decimal value, E=positive, W=negative.
% time_offset	= scalar, difference between GMT and local time, E=positive, W=negative.
%
% Output Data:
% index_sunrad_day	= index for InputData where daytime data is located, and both radiation as
% 			  well as sunset/sunrise indicate daytime
% index_sunrad_night	= index for InputData where nighttime data is located, and both radiation as
% 			  well as sunset/sunrise indicate nighttime
% index_day		= index for InputData where daytime data is located, indicated by either
% 			  radiation or sunset/sunrise showing daytime
% index_night		= index for InputData where nighttime data is located, indicated by either
% 			  radiation or sunset/sunrise showing nighttime
%


% set the limit for incoming shortwave radiation in W/m^2 below which data will be treated as
% nighttime data
RadiationLimit = 20;

% find the locations where incoming shortwave is below the threshold
index_rad_night = find(RGlobal<RadiationLimit);

% find the locations where incoming shortwave is above the threshold
index_rad_day = find(RGlobal>=RadiationLimit);

% get the input data for the times where global radiation indicates nighttime
% radiation_night = InputData(index_rad_night);

% get the input data for the times where global radiation indicates daytime
% radiation_day = InputData(index_rad_day);

% crosscheck radiation nighttime with sunset and sunrise data

% get a list of the days for which there are data.
days = unique(floor(Timestamp-1/1440)); % the 1/1440 adjusts the timestamps beacuse 
                             		% midnight in matlab is 00:00:00.

% caulculate the sun rise and sun set times                               
[rise_time, set_time] = sunrise_set_times(days, Latitude, Longitude,...
    time_offset);

% convert the outputs of sunrise_set_times.m function call to matlab
% datenumbers.
rise_time = datenum([zeros(size(days,1),3) rise_time(:,1) rise_time(:,2) zeros(size(days,1),1)]) + days;
set_time = datenum([zeros(size(days,1),3) set_time(:,1) set_time(:,2) zeros(size(days,1),1)]) + days;

% create empty index variables
index_sun_night = [];
index_sun_day = [];

% loop through each day.
for i = 1:length(days)

    % get the data for the ith day.
    t_i = Timestamp(floor(Timestamp-1/1440) == days(i));
    % data_i = InputData(floor(Timestamp-1/1440) == days(i),:);

    % get the indices of the 'daytime' data
    idxday = find(t_i > rise_time(i) & t_i <= set_time(i));

    % next do night time. get indices of 'nighttime' data.
    idxnite = find(t_i <= rise_time(i) | t_i > set_time(i));

    % find the indices for nighttime and daytime data in the yearly timevector
    [C,index_sun_night_i,ib] = intersect(Timestamp,t_i(idxnite));
    [D,index_sun_day_i,id] = intersect(Timestamp,t_i(idxday));
    index_sun_night = [index_sun_night; index_sun_night_i];
    index_sun_day = [index_sun_day; index_sun_day_i];

end

% combine the indices

% the variable index_sunrad_night/day indicates the indices for which both criteria are fullfilled,
% radiation is below/above the threshold AND sunrise/sunset data indicates night/day
index_sunrad_night = intersect(index_rad_night,index_sun_night);
index_sunrad_day = intersect(index_rad_day, index_sun_day);

% the variable index_night/day indicates the indices for which either of the criteria is fullfilled,
% radiation is below/above the threshold OR sunrise/sunset data indicates night/day

index_night = unique([index_rad_night; index_sun_night]);
index_day = unique([index_rad_day; index_sun_day]);





%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [rise_time, set_time] = sunrise_set_times(timestamps, lat, lon,...
    time_offset)

%uses NOAA algorithms to calculate the sunrise and set times for the dates
%specified in the timestamps input arument. Details about the NOAA calculations 
%are available at
%
%           http://www.esrl.noaa.gov/gmd/grad/solcalc/calcdetails.html [verified
%           2013-Mar].
%
%   inputs: timestamps - nx1 vector of timestamps as matlab datenumbers.
%           lat - scalar, the latitude of the location in degrees. (minutes 
%               should be represented as a decimal value. [N=positive; S=negative].
%           lon - scalar, the longitude of the location in degrees. minutes
%               should be represented as a decimal value. [E=positive; W=negative].
%           time_offset - scalar, the difference between GMT and local time
%               [E=positive; W=negative].
%           
%    outputs: rise_time - nx2 matrix of sun rise times. column 1 is hours
%                   and column 2 is minutes.
%             set_times - nx2 matrix of sun set times. column 1 is hours
%                   and column 2 is minutes.
%
%   NOTE - I checked agreement between the matlab function and the NOAA
%   spreadsheet at absolute latitudes < 60 degrees, agreement was within
%   +/- 4 minutes (due to rounding errors). If higher precision is required 
%   a closer look at the algorithm is needed.
%
% written by jeff wood

%convert the matlab datenumbers to julian date and calculate the julian
%century.
jd = timestamps - 693961 + 2415018.5 - time_offset/24;
julcent = (jd - 2451545)./36525;

%solar calculations

%Geom mean long sun angle(degrees)
GMLS = mod(280.46646 + julcent.*(36000.76983 + julcent.*0.0003032), 360);

%geom mean anom sun (degrees)
GMAS = 357.52911 + julcent.*(35999.05029 - 0.0001537.*julcent);

%eccent earth orbit
EEO = 0.016708634 - julcent.*(0.000042037 + 0.0000001267.*julcent);

%sun equation of center
SEC = sin(pi/180.*GMAS).*(1.914602 - julcent.*(0.004817 + 0.000014.*julcent))...
    + sin(pi/180.*(2.*GMAS)).*(0.019993 - 0.000101.*julcent) + sin(pi/180.*(3.*GMAS)).*0.000289;

%sun true long (degrees)
STL = GMLS + SEC;

%sun true anom (degrees)
STA = GMAS + SEC;

%sun rad veector (AUs)
SRV = (1.000001018.*(1 - EEO.^2))./(1 + EEO.*cos(pi/180.*STA));

%sun app long (degrees)
SAL = STL - 0.00569 - 0.00478.*sin(pi/180.*(125.04 - 1934.136.*julcent));

%mean oblique ecliptic (degrees)
MOE = 23 + (26 + ((21.448 - julcent.*(46.815 + julcent.*(0.00059 - julcent.*0.001813))))/60)/60;

%Oblique correction
OC = MOE + 0.00256.*cos(pi/180.*(125.04 - 1934.136.*julcent));

%sun right accension (deg)
SRtAc = (pi/2 - atan2(cosd(SAL), cosd(OC).*sind(SAL))).*180./pi; %needed to adjust to conform to excel

%sun declination (degrees)
SD = asind(sind(OC).*sind(SAL));

%var y
vary = tand((OC/2)).*tand((OC/2));

%equation of time (min)
ET = 4.*(vary.*sind(2.*GMLS) - 2.*EEO.*sind(GMAS) + 4.*EEO.*vary.*...
    sind(GMAS).*cosd(2.*GMLS) - 0.5.*vary.^2.*sind(4.*GMLS) - ...
    1.25.*EEO.^2.*sind(2.*GMAS)).*180./pi;

%HA sunrise (degrees)
HAS = acosd(cosd(90.833)./(cosd(lat).*cosd(SD)) - tand(lat).*tand(SD));

%solar noon (LST)
SN   = (720 - 4.*lon - ET + time_offset.*60)./1440;

%sunrise time (LST)
rise_time_hrs = (SN - HAS.*4./1440).*24;

%sunset time (LST)
set_time_hrs = (SN + HAS.*4./1440).*24;

%put sunrise and set time data into output vectors
rise_time = [floor(rise_time_hrs), floor((rise_time_hrs - floor(rise_time_hrs)).*60)];
set_time = [floor(set_time_hrs), floor((set_time_hrs - floor(set_time_hrs)).*60)];
