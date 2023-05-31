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
