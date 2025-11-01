clear
clf
%
%
%           1) compute power production (in kW) for solar panel
%           2) compare to measured values
%
%

%===============================================
%        set variables and file names
%===============================================
file_in = "SolarArrayProduction.xlsx"; % power production:6/1/2018 to 6/30/2018, 
% header: Timestamp, AH3(power), LSParking, 
% rows: 2880 (30*24*4) rows (15 min. intervals)
% time range for day time hours
tmin = 5.5; %hr
tmax = 20; %hr
dt   = 0.25; %hr
day  = 26; %between 1-30, day in June to compare model and obs., use 26

%Earth tilte angle
% see e.g.: https://www.pveducation.org/pvcdrom/properties-of-sunlight/declination-angle
dec = 23.45;       % Summer solstice
lat = 42 + 17/60;  % 42 degree, 17 minutes

%===============================================
%        compute theoretical power production
%===============================================
dec = deg2rad(dec);
lat = deg2rad(lat);

% create time vector for day time hours
t = tmin:dt:tmax; 

% Calculate local solar time (LST)
LST = t - 1 + 14.6/60; 

% amount of solar irradiance on the solar panels
sunangle = sin(dec)*sin(lat) + cos(dec)*cos(lat)*cosd(15*(LST - 12));
% sun incidence (peak enery for 90 degree incidence)
S_inc = 1.4883*0.7.^(sunangle.^(-0.678)); 

% compute power production during daytime
production_theory = 270*S_inc.*sunangle;
%Calculate final  production with inverter limit (207 kW).
% note! min returns min of array, but here second argument caps vector to
% 207
production_theory = min(production_theory, 207);

%===============================================
%         load spreadsheet with measurements
%===============================================
% load data from spreadsheet with 'readtable'
production     = readtable(file_in);  
% Reorganize data into a matrix of times by days (4*24 = 96 times, by 30 days)
dailyPower     = reshape(production.AH3, 96, 30);
% select specfic day
singleDayPower = dailyPower(:, day);

%===============================================
%          plots
%===============================================
fig = figure(1);
hold on;
tfullday = 0:0.25:23.75;
plot(tfullday,singleDayPower,'.-')
%plot(production.Timestamp,production.AH3)
plot(t,production_theory)
legend( strcat('Obs.: 6/', int2str(day)), 'Model')
xlabel('Time of day')
ylabel('Power Production (kW)')
title('Theoretical max power production for AH3 solar array')
saveas(fig, strcat('SolarPower_6_', int2str(day),'.png'))