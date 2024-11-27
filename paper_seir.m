%% data import
clearvars;close all;clc;
[tableConfirmed,tableDeaths,tableRecovered,~] = getDataCOVID();


data_v = readtable("owid-covid-data.csv");
est_day = 7; %estimation length
%% Country information

Location = 'Korea, South';
Iso_code = "KOR";
Npop= 51710000; % population
%fit = [7.0329 413.9947]; %korea
fit = [9.1832 75.0307]; %omicron in korea
tIndex = ["2020-11-01", "2020-11-23", "2020-12-22",...
     "2021-01-02", "2021-03-01", "2021-04-19", "2021-05-06", "2021-07-01", "2021-07-15",...
     "2021-08-06", "2021-09-06",...
     "2021-09-20", "2021-11-01", "2021-11-22", "2021-12-14", "2022-01-04", "2022-01-28"];
pIndex = [0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0];
%% preprocessing

day_start = tIndex(1);
day_end = tIndex(end);
guess = [1 0.25];
tStart1 = datetime(day_start); % Beginning of the first wave
tEnd = datetime(day_end); % End of simulation
try
 indR = find(contains(tableRecovered.CountryRegion,Location)==1);
 indV = find(contains(data_v.iso_code,Iso_code)==1);
 
catch exception
 searchLoc = strfind(tableRecovered.CountryRegion,Location);
 indR = find(~cellfun(@isempty,searchLoc)) ;
end
time = tStart1:tEnd;

Recovered = table2array(tableRecovered(indR,5:end));
%% SEIRD model

 Deaths_time = data_v(indV, [4,8]);
 Daily_time = data_v(indV, [4,6]);
 
 Vaccination1_time = data_v(indV, [4,36]);
 Vaccination2_time = data_v(indV, [4,37]);
 Vaccination3_time = data_v(indV, [4,38]);
 
 Deaths_time((Deaths_time.date<tStart1 | Deaths_time.date>tEnd), :) = [];
 Daily_time((Daily_time.date<tStart1 | Daily_time.date>tEnd), :) = [];
 
 Vaccination1_time((Vaccination1_time.date<tStart1 | Vaccination1_time.date>tEnd), :) = [];
 Vaccination2_time((Vaccination2_time.date<tStart1 | Vaccination2_time.date>tEnd), :) = [];
 Vaccination3_time((Vaccination3_time.date<tStart1 | Vaccination3_time.date>tEnd), :) = [];
 
 
 
 
 Recovered(time<datetime(tStart1) | time>datetime(tEnd)) = [];
 time(time<datetime(tStart1) | time>datetime(tEnd)) = [];
 Deaths = table2array(Deaths_time(:, 2));
 Daily = table2array(Daily_time(:, 2));
 
 Daily(Daily<=0) = NaN;
 nanx = isnan(Daily);
 t = 1:numel(Daily);
 Daily(nanx) = interp1(t(~nanx), Daily(~nanx), t(nanx));
 
 Vaccination1 = table2array(Vaccination1_time(:, 2));
 Vaccination2 = table2array(Vaccination2_time(:, 2));
 Vaccination3 = table2array(Vaccination3_time(:, 2));
 
 Vaccination11 = fillmissing(Vaccination1, "linear", "EndValues", "nearest");
 Vaccination22 = fillmissing(Vaccination2, "linear", "EndValues", "nearest");
 Vaccination33 = fillmissing(Vaccination3, "linear", "EndValues", "nearest");
 
 Vaccination111 = diff(Vaccination11);
 Vaccination222 = diff(Vaccination22);
 Vaccination333 = diff(Vaccination33);
 
 vac_len1 = length(Daily)-length(Vaccination111);
 vac_len2 = length(Daily)-length(Vaccination222);
 vac_len3 = length(Daily)-length(Vaccination333);
 
 Vaccination1111 = [repelem(0,vac_len1), transpose(Vaccination111)];
 Vaccination2222 = [repelem(0,vac_len1), transpose(Vaccination222)];
 Vaccination3333 = [repelem(0,vac_len1), transpose(Vaccination333)];
 
 Confirmed = cumsum(Daily);
 
 [I,Q,R,D,newT] = getMultipleWaves_7(guess,Npop,time,Confirmed,Recovered,Deaths,...
     repelem(0, length(Daily)),tStart1,tIndex,pIndex,fit,est_day);
 
 Daily_predicted = I+Q+R+D;
 
  for kk = 2:length(tIndex)
     
     dup_index0 = find(newT==datetime(tIndex(kk)));
     dup_index = dup_index0(1);
     newT(dup_index) = [];
     Daily_predicted(dup_index) = [];
     D(dup_index) = [];
  end 
 
 % 7 days prediction
 
 
 
 Daily_predicted_0 = Daily_predicted(:,1:(end-(est_day*24+1)*2));
 Daily_predicted_1 = Daily_predicted(:,[1:(end-(est_day*24+1)*3) (end-(est_day*24+1)*2+1):(end-(est_day*24+1))]);
 Daily_predicted_2 = Daily_predicted(:,[1:(end-(est_day*24+1)*3) (end-(est_day*24+1)+1):(end)]);
 
 D_0 = D(:,1:(end-(est_day*24+1)*2));
 dD = gradient(D_0)./gradient(datenum(newT(1:end-(est_day*24+1)*2)));
 
 dy = gradient(Daily_predicted_0)./gradient(datenum(newT(1:end-(est_day*24+1)*2)));
 dy_tot = dy';
 dy_smoothed = smoothdata(dy, 'movmean', est_day*24);
 dy1 = gradient(Daily_predicted_1)./gradient(datenum(newT(1:end-(est_day*24+1)*2)));
 dy1_smoothed = smoothdata(dy1, 'movmean', est_day*24);
 dy2 = gradient(Daily_predicted_2)./gradient(datenum(newT(1:end-(est_day*24+1)*2)));
 dy2_smoothed = smoothdata(dy2, 'movmean', est_day*24);
 
 % SEIR Figure
%  f = figure;
%  f.Position = [500 300 1200 600];
%  semilogy(time,Daily,'b.');
%  hold on
%  semilogy(newT(1:end-(est_day*24+1)*3), dy(1:end-(est_day*24+1)), 'k');
%  hold on
%  semilogy(newT(end-(est_day*24+1)*3+1:end-(est_day*24+1)*2), dy(end-(est_day*24):end), 'r','LineWidth',2);
%  hold on
%  semilogy(newT(end-(est_day*24+1)*2+1:end-(est_day*24+1)), dy1(end-(est_day*24):end), 'g','LineWidth',2);
%  hold on
%  semilogy(newT(end-(est_day*24+1)+1:end), dy2(end-(est_day*24):end), 'y','LineWidth',2);
%  hold on
%  xline(datetime(tIndex(datetime(tIndex)>=datetime(tIndex(end))-364)), ":")
%  ylabel('Number of cases')
%  xlabel('time (days)')
%  leg = {'Reported', 'Fitted', 'Predicted (Maintained)',...
%      'Predicted (reinforced)','Predicted (relaxed)'};
%  legend(leg{:},'location','eastoutside')
%  set(gcf,'color','w')
%  set(gca,'yscale','lin')
%  title(['Location: ', "Korea"])
%  hold off
%  shg

Daily_time_test = data_v(indV, [4,6]);
Daily_time_test((Daily_time_test.date<=tEnd | Daily_time_test.date>tEnd+est_day), :) = [];
Daily_test = table2array(Daily_time_test(:, 2));
Deaths_time_test = data_v(indV, [4,9]);
Deaths_time_test((Deaths_time_test.date<tStart1 | Deaths_time_test.date>tEnd+est_day), :) = [];
Deaths_test = table2array(Deaths_time_test(:, 2));

from = tStart1;
to = tEnd+est_day;

% data_cases = [[Daily; Daily_test], round(dy(end-hours(to-from):24:end)'),...
%     Deaths_test, round(dD(end-hours(to-from):24:end)')];
% data_date = time(end-days(to-from)+est_day):time(end)+est_day;

% data_cases = [[Daily; Daily_test], round(dy(end-hours(to-from):24:end)'),...
%     Deaths_test, round(dD(end-hours(to-from):24:end)')];
% data_date = time(end-days(to-from)+est_day):time(end)+est_day;
% data_date = datetime(data_date, "Format", "uuuu-MM-dd");

data_cases = [[Daily; Daily_test], round(dy(end-hours(to-from):24:end)'),...
    Deaths_test, round(dD(end-hours(to-from):24:end)')];
data_date = time(end-days(to-from)+est_day):time(end)+est_day;
data_date = datetime(data_date, "Format", "uuuu-MM-dd");
%% SEIQRD model

 
 
 Deaths_time = data_v(indV, [4,8]); % total deaths
 Daily_time = data_v(indV, [4,6]); % new cases
 
 Vaccination1_time = data_v(indV, [4,36]); % people vaccinated
 Vaccination2_time = data_v(indV, [4,37]); % people fully vaccinated
 Vaccination3_time = data_v(indV, [4,38]); % total boosters
 
 Deaths_time((Deaths_time.date<tStart1 | Deaths_time.date>tEnd), :) = [];
 Daily_time((Daily_time.date<tStart1 | Daily_time.date>tEnd), :) = [];
 
 Vaccination1_time((Vaccination1_time.date<tStart1 | Vaccination1_time.date>tEnd), :) = [];
 Vaccination2_time((Vaccination2_time.date<tStart1 | Vaccination2_time.date>tEnd), :) = [];
 Vaccination3_time((Vaccination3_time.date<tStart1 | Vaccination3_time.date>tEnd), :) = [];
 
 Recovered(time<datetime(tStart1) | time>datetime(tEnd)) = [];
 time(time<datetime(tStart1) | time>datetime(tEnd)) = [];
 Deaths = table2array(Deaths_time(:, 2));
 Daily = table2array(Daily_time(:, 2));
 
 Daily(Daily<=0) = NaN;
 nanx = isnan(Daily);
 t = 1:numel(Daily);
 Daily(nanx) = interp1(t(~nanx), Daily(~nanx), t(nanx));
 
 Vaccination1 = table2array(Vaccination1_time(:, 2));
 Vaccination2 = table2array(Vaccination2_time(:, 2));
 Vaccination3 = table2array(Vaccination3_time(:, 2));
 
 Vaccination11 = fillmissing(Vaccination1, "linear", "EndValues", "nearest");
 Vaccination22 = fillmissing(Vaccination2, "linear", "EndValues", "nearest");
 Vaccination33 = fillmissing(Vaccination3, "linear", "EndValues", "nearest");
 
 Vaccination111 = diff(Vaccination11);
 Vaccination222 = diff(Vaccination22);
 Vaccination333 = diff(Vaccination33);
 
 vac_len1 = length(Daily)-length(Vaccination111);
 vac_len2 = length(Daily)-length(Vaccination222);
 vac_len3 = length(Daily)-length(Vaccination333);
 
 Vaccination1111 = [repelem(0,vac_len1), transpose(Vaccination111)];
 Vaccination2222 = [repelem(0,vac_len1), transpose(Vaccination222)];
 Vaccination3333 = [repelem(0,vac_len1), transpose(Vaccination333)];
 
 Vaccination1111(Vaccination1111 ~= 0) = 0;
 Vaccination2222(Vaccination2222 ~= 0) = 0;
 Vaccination3333(Vaccination3333 ~= 0) = 0;
 
 Confirmed = cumsum(Daily);
 
 [Q,R,D,newT] = getMultipleWaves_5(guess,Npop,time,Confirmed,Recovered,Deaths,...
     Vaccination1111,Vaccination2222,Vaccination3333,tStart1,tIndex,pIndex,fit,est_day);
 
 Daily_predicted = Q+R+D;
 
  for kk = 2:length(tIndex)
     
     dup_index0 = find(newT==datetime(tIndex(kk)));
     dup_index = dup_index0(1);
     newT(dup_index) = [];
     Daily_predicted(dup_index) = [];
     D(dup_index) = [];
  end 
 
 % 7 days prediction
 
 
 
 Daily_predicted_0 = Daily_predicted(:,1:(end-(est_day*24+1)*2));
 Daily_predicted_1 = Daily_predicted(:,[1:(end-(est_day*24+1)*3) (end-(est_day*24+1)*2+1):(end-(est_day*24+1))]);
 Daily_predicted_2 = Daily_predicted(:,[1:(end-(est_day*24+1)*3) (end-(est_day*24+1)+1):(end)]);
 
 D_0 = D(:,1:(end-(est_day*24+1)*2));
 dD = gradient(D_0)./gradient(datenum(newT(1:end-(est_day*24+1)*2)));
 
 dy = gradient(Daily_predicted_0)./gradient(datenum(newT(1:end-(est_day*24+1)*2)));
 dy_tot = [dy_tot, dy'];
 dy_smoothed = smoothdata(dy, 'movmean', est_day*24);
 dy1 = gradient(Daily_predicted_1)./gradient(datenum(newT(1:end-(est_day*24+1)*2)));
 dy1_smoothed = smoothdata(dy1, 'movmean', est_day*24);
 dy2 = gradient(Daily_predicted_2)./gradient(datenum(newT(1:end-(est_day*24+1)*2)));
 dy2_smoothed = smoothdata(dy2, 'movmean', est_day*24);
 
 
%  f = figure;
%  f.Position = [500 300 1200 600];
%  semilogy(time,Daily,'b.');
%  hold on
%  semilogy(newT(1:end-(est_day*24+1)*3), dy(1:end-(est_day*24+1)), 'k');
%  hold on
%  semilogy(newT(end-(est_day*24+1)*3+1:end-(est_day*24+1)*2), dy(end-(est_day*24):end), 'r','LineWidth',2);
%  hold on
%  semilogy(newT(end-(est_day*24+1)*2+1:end-(est_day*24+1)), dy1(end-(est_day*24):end), 'g','LineWidth',2);
%  hold on
%  semilogy(newT(end-(est_day*24+1)+1:end), dy2(end-(est_day*24):end), 'y','LineWidth',2);
%  hold on
%  xline(datetime(tIndex(datetime(tIndex)>=datetime(tIndex(end))-364)), ":")
%  ylabel('Number of cases')
%  xlabel('time (days)')
%  leg = {'Reported', 'Fitted', 'Predicted (Maintained)',...
%      'Predicted (reinforced)','Predicted (relaxed)'};
%  legend(leg{:},'location','eastoutside')
%  set(gcf,'color','w')
%  set(gca,'yscale','lin')
%  title(['Location: ', "Korea"])
%  hold off
%  shg
 


data_cases = [data_cases, round(dy(end-hours(to-from):24:end)'), round(dD(end-hours(to-from):24:end)')];
%% SEIQRDVUP model (V is fully and effectively vaccinated)

 time = tStart1:tEnd;
 
 
 Deaths_time = data_v(indV, [4,8]); % total deaths
 Daily_time = data_v(indV, [4,6]); % new cases
 
 Vaccination1_time = data_v(indV, [4,36]); % people vaccinated
 Vaccination2_time = data_v(indV, [4,37]); % people fully vaccinated
 Vaccination3_time = data_v(indV, [4,38]); % total boosters
 
 Deaths_time((Deaths_time.date<tStart1 | Deaths_time.date>tEnd), :) = [];
 Daily_time((Daily_time.date<tStart1 | Daily_time.date>tEnd), :) = [];
 
 Vaccination1_time((Vaccination1_time.date<tStart1 | Vaccination1_time.date>tEnd), :) = [];
 Vaccination2_time((Vaccination2_time.date<tStart1 | Vaccination2_time.date>tEnd), :) = [];
 Vaccination3_time((Vaccination3_time.date<tStart1 | Vaccination3_time.date>tEnd), :) = [];
 
 Recovered(time<datetime(tStart1) | time>datetime(tEnd)) = [];
 time(time<datetime(tStart1) | time>datetime(tEnd)) = [];
 Deaths = table2array(Deaths_time(:, 2));
 Daily = table2array(Daily_time(:, 2));
 
 Daily(Daily<=0) = NaN;
 nanx = isnan(Daily);
 t = 1:numel(Daily);
 Daily(nanx) = interp1(t(~nanx), Daily(~nanx), t(nanx));
 
 Vaccination1 = table2array(Vaccination1_time(:, 2));
 Vaccination2 = table2array(Vaccination2_time(:, 2));
 Vaccination3 = table2array(Vaccination3_time(:, 2));
 
 Vaccination11 = fillmissing(Vaccination1, "linear", "EndValues", "nearest");
 Vaccination22 = fillmissing(Vaccination2, "linear", "EndValues", "nearest");
 Vaccination33 = fillmissing(Vaccination3, "linear", "EndValues", "nearest");
 
 Vaccination111 = diff(Vaccination11);
 Vaccination222 = diff(Vaccination22);
 Vaccination333 = diff(Vaccination33);
 
 vac_len1 = length(Daily)-length(Vaccination111);
 vac_len2 = length(Daily)-length(Vaccination222);
 vac_len3 = length(Daily)-length(Vaccination333);
 
 Vaccination1111 = [repelem(0,vac_len1), transpose(Vaccination111)];
 Vaccination2222 = [repelem(0,vac_len1), transpose(Vaccination222)];
 Vaccination3333 = [repelem(0,vac_len1), transpose(Vaccination333)];
 
 Confirmed = cumsum(Daily);
 
 [Q,R,D,newT] = getMultipleWaves_2(guess,Npop,time,Confirmed,Recovered,Deaths,...
     Vaccination2222,tStart1,tIndex,pIndex,fit,est_day);
 
 Daily_predicted = Q+R+D;
 
  for kk = 2:length(tIndex)
     
     dup_index0 = find(newT==datetime(tIndex(kk)));
     dup_index = dup_index0(1);
     newT(dup_index) = [];
     Daily_predicted(dup_index) = [];
     D(dup_index) = [];
  end 
 
 % 7 days prediction
 
 
 
 Daily_predicted_0 = Daily_predicted(:,1:(end-(est_day*24+1)*2));
 Daily_predicted_1 = Daily_predicted(:,[1:(end-(est_day*24+1)*3) (end-(est_day*24+1)*2+1):(end-(est_day*24+1))]);
 Daily_predicted_2 = Daily_predicted(:,[1:(end-(est_day*24+1)*3) (end-(est_day*24+1)+1):(end)]);
 
 D_0 = D(:,1:(end-(est_day*24+1)*2));
 dD = gradient(D_0)./gradient(datenum(newT(1:end-(est_day*24+1)*2)));
 
 dy = gradient(Daily_predicted_0)./gradient(datenum(newT(1:end-(est_day*24+1)*2)));  
 dy_tot = [dy_tot, dy'];
 dy_smoothed = smoothdata(dy, 'movmean', est_day*24);
 dy1 = gradient(Daily_predicted_1)./gradient(datenum(newT(1:end-(est_day*24+1)*2)));
 dy1_smoothed = smoothdata(dy1, 'movmean', est_day*24);
 dy2 = gradient(Daily_predicted_2)./gradient(datenum(newT(1:end-(est_day*24+1)*2)));
 dy2_smoothed = smoothdata(dy2, 'movmean', est_day*24);
 
 
%  f = figure;
%  f.Position = [500 300 1200 600];
%  semilogy(time,Daily,'b.');
%  hold on
%  semilogy(newT(1:end-(est_day*24+1)*3), dy(1:end-(est_day*24+1)), 'k');
%  hold on
%  semilogy(newT(end-(est_day*24+1)*3+1:end-(est_day*24+1)*2), dy(end-(est_day*24):end), 'r','LineWidth',2);
%  hold on
%  semilogy(newT(end-(est_day*24+1)*2+1:end-(est_day*24+1)), dy1(end-(est_day*24):end), 'g','LineWidth',2);
%  hold on
%  semilogy(newT(end-(est_day*24+1)+1:end), dy2(end-(est_day*24):end), 'y','LineWidth',2);
%  hold on
%  xline(datetime(tIndex(datetime(tIndex)>=datetime(tIndex(end))-364)), ":")
%  ylabel('Number of cases')
%  xlabel('time (days)')
%  leg = {'Reported', 'Fitted', 'Predicted (Maintained)',...
%      'Predicted (reinforced)','Predicted (relaxed)'};
%  legend(leg{:},'location','eastoutside')
%  set(gcf,'color','w')
%  set(gca,'yscale','lin')
%  title(['Location: ', "Korea"])
%  hold off
%  shg
 

data_cases = [data_cases, round(dy(end-hours(to-from):24:end)'), round(dD(end-hours(to-from):24:end)')];
%% SEIQRDVP (V is fully vaccinated group)

 Deaths_time = data_v(indV, [4,8]);
 Daily_time = data_v(indV, [4,6]);
 
 Vaccination1_time = data_v(indV, [4,36]);
 Vaccination2_time = data_v(indV, [4,37]);
 Vaccination3_time = data_v(indV, [4,38]);
 
 Deaths_time((Deaths_time.date<tStart1 | Deaths_time.date>tEnd), :) = [];
 Daily_time((Daily_time.date<tStart1 | Daily_time.date>tEnd), :) = [];
 
 Vaccination1_time((Vaccination1_time.date<tStart1 | Vaccination1_time.date>tEnd), :) = [];
 Vaccination2_time((Vaccination2_time.date<tStart1 | Vaccination2_time.date>tEnd), :) = [];
 Vaccination3_time((Vaccination3_time.date<tStart1 | Vaccination3_time.date>tEnd), :) = [];
 
 
 
 
 Recovered(time<datetime(tStart1) | time>datetime(tEnd)) = [];
 time(time<datetime(tStart1) | time>datetime(tEnd)) = [];
 Deaths = table2array(Deaths_time(:, 2));
 Daily = table2array(Daily_time(:, 2));
 
 Daily(Daily<=0) = NaN;
 nanx = isnan(Daily);
 t = 1:numel(Daily);
 Daily(nanx) = interp1(t(~nanx), Daily(~nanx), t(nanx));
 
 Vaccination1 = table2array(Vaccination1_time(:, 2));
 Vaccination2 = table2array(Vaccination2_time(:, 2));
 Vaccination3 = table2array(Vaccination3_time(:, 2));
 
 Vaccination11 = fillmissing(Vaccination1, "linear", "EndValues", "nearest");
 Vaccination22 = fillmissing(Vaccination2, "linear", "EndValues", "nearest");
 Vaccination33 = fillmissing(Vaccination3, "linear", "EndValues", "nearest");
 
 Vaccination111 = diff(Vaccination11);
 Vaccination222 = diff(Vaccination22);
 Vaccination333 = diff(Vaccination33);
 
 vac_len1 = length(Daily)-length(Vaccination111);
 vac_len2 = length(Daily)-length(Vaccination222);
 vac_len3 = length(Daily)-length(Vaccination333);
 
 Vaccination1111 = [repelem(0,vac_len1), transpose(Vaccination111)];
 Vaccination2222 = [repelem(0,vac_len1), transpose(Vaccination222)];
 Vaccination3333 = [repelem(0,vac_len1), transpose(Vaccination333)];
 
 Confirmed = cumsum(Daily);
 
 [Q,R,D,newT] = getMultipleWaves_3(guess,Npop,time,Confirmed,Recovered,Deaths,...
     Vaccination2222,tStart1,tIndex,pIndex,fit,est_day);
 
 Daily_predicted = Q+R+D;
 
  for kk = 2:length(tIndex)
     
     dup_index0 = find(newT==datetime(tIndex(kk)));
     dup_index = dup_index0(1);
     newT(dup_index) = [];
     Daily_predicted(dup_index) = [];
     D(dup_index) = [];
  end 
 
 % 7 days prediction
 
 
 
 Daily_predicted_0 = Daily_predicted(:,1:(end-(est_day*24+1)*2));
 Daily_predicted_1 = Daily_predicted(:,[1:(end-(est_day*24+1)*3) (end-(est_day*24+1)*2+1):(end-(est_day*24+1))]);
 Daily_predicted_2 = Daily_predicted(:,[1:(end-(est_day*24+1)*3) (end-(est_day*24+1)+1):(end)]);
 
 D_0 = D(:,1:(end-(est_day*24+1)*2));
 dD = gradient(D_0)./gradient(datenum(newT(1:end-(est_day*24+1)*2)));
 
 dy = gradient(Daily_predicted_0)./gradient(datenum(newT(1:end-(est_day*24+1)*2)));  
 dy_tot = [dy_tot, dy'];
 dy_smoothed = smoothdata(dy, 'movmean', est_day*24);
 dy1 = gradient(Daily_predicted_1)./gradient(datenum(newT(1:end-(est_day*24+1)*2)));
 dy1_smoothed = smoothdata(dy1, 'movmean', est_day*24);
 dy2 = gradient(Daily_predicted_2)./gradient(datenum(newT(1:end-(est_day*24+1)*2)));
 dy2_smoothed = smoothdata(dy2, 'movmean', est_day*24);
 
 
%  f = figure;
%  f.Position = [500 300 1200 600];
%  semilogy(time,Daily,'b.');
%  hold on
%  semilogy(newT(1:end-(est_day*24+1)*3), dy(1:end-(est_day*24+1)), 'k');
%  hold on
%  semilogy(newT(end-(est_day*24+1)*3+1:end-(est_day*24+1)*2), dy(end-(est_day*24):end), 'r','LineWidth',2);
%  hold on
%  semilogy(newT(end-(est_day*24+1)*2+1:end-(est_day*24+1)), dy1(end-(est_day*24):end), 'g','LineWidth',2);
%  hold on
%  semilogy(newT(end-(est_day*24+1)+1:end), dy2(end-(est_day*24):end), 'y','LineWidth',2);
%  hold on
%  xline(datetime(tIndex(datetime(tIndex)>=datetime(tIndex(end))-364)), ":")
%  ylabel('Number of cases')
%  xlabel('time (days)')
%  leg = {'Reported', 'Fitted', 'Predicted (Maintained)',...
%      'Predicted (reinforced)','Predicted (relaxed)'};
%  legend(leg{:},'location','eastoutside')
%  set(gcf,'color','w')
%  set(gca,'yscale','lin')
%  title(['Location: ', "Korea"])
%  hold off
%  shg

data_cases = [data_cases, round(dy(end-hours(to-from):24:end)'), round(dD(end-hours(to-from):24:end)')];
%% SEIQRDV3P
 
 Deaths_time = data_v(indV, [4,8]);
 Daily_time = data_v(indV, [4,6]);
 
 Vaccination1_time = data_v(indV, [4,36]);
 Vaccination2_time = data_v(indV, [4,37]);
 Vaccination3_time = data_v(indV, [4,38]);
 
 Deaths_time((Deaths_time.date<tStart1 | Deaths_time.date>tEnd), :) = [];
 Daily_time((Daily_time.date<tStart1 | Daily_time.date>tEnd), :) = [];
 
 Vaccination1_time((Vaccination1_time.date<tStart1 | Vaccination1_time.date>tEnd), :) = [];
 Vaccination2_time((Vaccination2_time.date<tStart1 | Vaccination2_time.date>tEnd), :) = [];
 Vaccination3_time((Vaccination3_time.date<tStart1 | Vaccination3_time.date>tEnd), :) = [];
 
 
 
 
 Recovered(time<datetime(tStart1) | time>datetime(tEnd)) = [];
 time(time<datetime(tStart1) | time>datetime(tEnd)) = [];
 Deaths = table2array(Deaths_time(:, 2));
 Daily = table2array(Daily_time(:, 2));
 
 Daily(Daily<=0) = NaN;
 nanx = isnan(Daily);
 t = 1:numel(Daily);
 Daily(nanx) = interp1(t(~nanx), Daily(~nanx), t(nanx));
 
 Vaccination1 = table2array(Vaccination1_time(:, 2));
 Vaccination2 = table2array(Vaccination2_time(:, 2));
 Vaccination3 = table2array(Vaccination3_time(:, 2));
 
 Vaccination11 = fillmissing(Vaccination1, "linear", "EndValues", "nearest");
 Vaccination22 = fillmissing(Vaccination2, "linear", "EndValues", "nearest");
 Vaccination33 = fillmissing(Vaccination3, "linear", "EndValues", "nearest");
 
 Vaccination111 = diff(Vaccination11);
 Vaccination222 = diff(Vaccination22);
 Vaccination333 = diff(Vaccination33);
 
 vac_len1 = length(Daily)-length(Vaccination111);
 vac_len2 = length(Daily)-length(Vaccination222);
 vac_len3 = length(Daily)-length(Vaccination333);
 
 Vaccination1111 = [repelem(0,vac_len1), transpose(Vaccination111)];
 Vaccination2222 = [repelem(0,vac_len1), transpose(Vaccination222)];
 Vaccination3333 = [repelem(0,vac_len1), transpose(Vaccination333)];
 
 Confirmed = cumsum(Daily);
 
 [Q,R,D,newT] = getMultipleWaves_5(guess,Npop,time,Confirmed,Recovered,Deaths,...
     Vaccination1111,Vaccination2222,Vaccination3333,tStart1,tIndex,pIndex,fit,est_day);
 
 Daily_predicted = Q+R+D;
 
  for kk = 2:length(tIndex)
     
     dup_index0 = find(newT==datetime(tIndex(kk)));
     dup_index = dup_index0(1);
     newT(dup_index) = [];
     Daily_predicted(dup_index) = [];
     D(dup_index) = [];
  end 
 
 % 7 days prediction
 
 
 
 Daily_predicted_0 = Daily_predicted(:,1:(end-(est_day*24+1)*2));
 Daily_predicted_1 = Daily_predicted(:,[1:(end-(est_day*24+1)*3) (end-(est_day*24+1)*2+1):(end-(est_day*24+1))]);
 Daily_predicted_2 = Daily_predicted(:,[1:(end-(est_day*24+1)*3) (end-(est_day*24+1)+1):(end)]);
 
 D_0 = D(:,1:(end-(est_day*24+1)*2));
 dD = gradient(D_0)./gradient(datenum(newT(1:end-(est_day*24+1)*2)));
 
 dy = gradient(Daily_predicted_0)./gradient(datenum(newT(1:end-(est_day*24+1)*2)));  
 dy_tot = [dy_tot, dy'];
 dy_smoothed = smoothdata(dy, 'movmean', est_day*24);
 dy1 = gradient(Daily_predicted_1)./gradient(datenum(newT(1:end-(est_day*24+1)*2)));
 dy1_smoothed = smoothdata(dy1, 'movmean', est_day*24);
 dy2 = gradient(Daily_predicted_2)./gradient(datenum(newT(1:end-(est_day*24+1)*2)));
 dy2_smoothed = smoothdata(dy2, 'movmean', est_day*24);
 
 
%  f = figure;
%  f.Position = [500 300 1200 600];
%  semilogy(time,Daily,'b.');
%  hold on
%  semilogy(newT(1:end-(est_day*24+1)*3), dy(1:end-(est_day*24+1)), 'k');
%  hold on
%  semilogy(newT(end-(est_day*24+1)*3+1:end-(est_day*24+1)*2), dy(end-(est_day*24):end), 'r','LineWidth',2);
%  hold on
%  semilogy(newT(end-(est_day*24+1)*2+1:end-(est_day*24+1)), dy1(end-(est_day*24):end), 'g','LineWidth',2);
%  hold on
%  semilogy(newT(end-(est_day*24+1)+1:end), dy2(end-(est_day*24):end), 'y','LineWidth',2);
%  hold on
%  xline(datetime(tIndex(datetime(tIndex)>=datetime(tIndex(end))-364)), ":")
%  ylabel('Number of cases')
%  xlabel('time (days)')
%  leg = {'Reported', 'Fitted', 'Predicted (Maintained)',...
%      'Predicted (reinforced)','Predicted (relaxed)'};
%  legend(leg{:},'location','eastoutside')
%  set(gcf,'color','w')
%  set(gca,'yscale','lin')
%  title(['Location: ', "Korea"])
%  hold off
%  shg


data_cases = [data_cases, round(dy(end-hours(to-from):24:end)'), round(dD(end-hours(to-from):24:end)')];
%% SEIQRDV3P + omicron (k=3)
 
 pIndex(14:length(pIndex)) = 3;
 
 Deaths_time = data_v(indV, [4,8]);
 Daily_time = data_v(indV, [4,6]);
 
 Vaccination1_time = data_v(indV, [4,36]);
 Vaccination2_time = data_v(indV, [4,37]);
 Vaccination3_time = data_v(indV, [4,38]);
 
 Deaths_time((Deaths_time.date<tStart1 | Deaths_time.date>tEnd), :) = [];
 Daily_time((Daily_time.date<tStart1 | Daily_time.date>tEnd), :) = [];
 
 Vaccination1_time((Vaccination1_time.date<tStart1 | Vaccination1_time.date>tEnd), :) = [];
 Vaccination2_time((Vaccination2_time.date<tStart1 | Vaccination2_time.date>tEnd), :) = [];
 Vaccination3_time((Vaccination3_time.date<tStart1 | Vaccination3_time.date>tEnd), :) = [];
 
 
 Recovered(time<datetime(tStart1) | time>datetime(tEnd)) = [];
 time(time<datetime(tStart1) | time>datetime(tEnd)) = [];
 Deaths = table2array(Deaths_time(:, 2));
 Daily = table2array(Daily_time(:, 2));
 
 Daily(Daily<=0) = NaN;
 nanx = isnan(Daily);
 t = 1:numel(Daily);
 Daily(nanx) = interp1(t(~nanx), Daily(~nanx), t(nanx));
 
 Vaccination1 = table2array(Vaccination1_time(:, 2));
 Vaccination2 = table2array(Vaccination2_time(:, 2));
 Vaccination3 = table2array(Vaccination3_time(:, 2));
 
 Vaccination11 = fillmissing(Vaccination1, "linear", "EndValues", "nearest");
 Vaccination22 = fillmissing(Vaccination2, "linear", "EndValues", "nearest");
 Vaccination33 = fillmissing(Vaccination3, "linear", "EndValues", "nearest");
 
 Vaccination111 = diff(Vaccination11);
 Vaccination222 = diff(Vaccination22);
 Vaccination333 = diff(Vaccination33);
 
 vac_len1 = length(Daily)-length(Vaccination111);
 vac_len2 = length(Daily)-length(Vaccination222);
 vac_len3 = length(Daily)-length(Vaccination333);
 
 Vaccination1111 = [repelem(0,vac_len1), transpose(Vaccination111)];
 Vaccination2222 = [repelem(0,vac_len1), transpose(Vaccination222)];
 Vaccination3333 = [repelem(0,vac_len1), transpose(Vaccination333)];
 
 Confirmed = cumsum(Daily);
 
 [Q,R,D,newT] = getMultipleWaves_5(guess,Npop,time,Confirmed,Recovered,Deaths,...
     Vaccination1111,Vaccination2222,Vaccination3333,tStart1,tIndex,pIndex,fit,est_day);
 
 Daily_predicted = Q+R+D;
 
  for kk = 2:length(tIndex)
     
     dup_index0 = find(newT==datetime(tIndex(kk)));
     dup_index = dup_index0(1);
     newT(dup_index) = [];
     Daily_predicted(dup_index) = [];
     D(dup_index) = [];
  end 
 
 % 7 days prediction
 
 
 
 Daily_predicted_0 = Daily_predicted(:,1:(end-(est_day*24+1)*2));
 Daily_predicted_1 = Daily_predicted(:,[1:(end-(est_day*24+1)*3) (end-(est_day*24+1)*2+1):(end-(est_day*24+1))]);
 Daily_predicted_2 = Daily_predicted(:,[1:(end-(est_day*24+1)*3) (end-(est_day*24+1)+1):(end)]);
 
 D_0 = D(:,1:(end-(est_day*24+1)*2));
 dD = gradient(D_0)./gradient(datenum(newT(1:end-(est_day*24+1)*2)));
 
 dy = gradient(Daily_predicted_0)./gradient(datenum(newT(1:end-(est_day*24+1)*2)));  
 dy_tot = [dy_tot, dy'];
 dy_smoothed = smoothdata(dy, 'movmean', est_day*24);
 dy1 = gradient(Daily_predicted_1)./gradient(datenum(newT(1:end-(est_day*24+1)*2)));
 dy1_smoothed = smoothdata(dy1, 'movmean', est_day*24);
 dy2 = gradient(Daily_predicted_2)./gradient(datenum(newT(1:end-(est_day*24+1)*2)));
 dy2_smoothed = smoothdata(dy2, 'movmean', est_day*24);
 
 
%  f = figure;
%  f.Position = [500 300 1200 600];
%  semilogy(time,Daily,'b.');
%  hold on
%  semilogy(newT(1:end-(est_day*24+1)*3), dy(1:end-(est_day*24+1)), 'k');
%  hold on
%  semilogy(newT(end-(est_day*24+1)*3+1:end-(est_day*24+1)*2), dy(end-(est_day*24):end), 'r','LineWidth',2);
%  hold on
%  semilogy(newT(end-(est_day*24+1)*2+1:end-(est_day*24+1)), dy1(end-(est_day*24):end), 'g','LineWidth',2);
%  hold on
%  semilogy(newT(end-(est_day*24+1)+1:end), dy2(end-(est_day*24):end), 'y','LineWidth',2);
%  hold on
%  xline(datetime(tIndex(datetime(tIndex)>=datetime(tIndex(end))-364)), ":")
%  ylabel('Number of cases')
%  xlabel('time (days)')
%  leg = {'Reported', 'Fitted', 'Predicted (Maintained)',...
%      'Predicted (reinforced)','Predicted (relaxed)'};
%  legend(leg{:},'location','eastoutside')
%  set(gcf,'color','w')
%  set(gca,'yscale','lin')
%  title(['Location: ', "Korea"])
%  hold off
%  shg


data_cases = [data_cases, round(dy(end-hours(to-from):24:end)'), round(dD(end-hours(to-from):24:end)')];
%% SEIQRDV3P + omicron (k=5)
 
 pIndex(14:length(pIndex)) = 5;
 
 Deaths_time = data_v(indV, [4,8]);
 Daily_time = data_v(indV, [4,6]);
 
 Vaccination1_time = data_v(indV, [4,36]);
 Vaccination2_time = data_v(indV, [4,37]);
 Vaccination3_time = data_v(indV, [4,38]);
 
 Deaths_time((Deaths_time.date<tStart1 | Deaths_time.date>tEnd), :) = [];
 Daily_time((Daily_time.date<tStart1 | Daily_time.date>tEnd), :) = [];
 
 Vaccination1_time((Vaccination1_time.date<tStart1 | Vaccination1_time.date>tEnd), :) = [];
 Vaccination2_time((Vaccination2_time.date<tStart1 | Vaccination2_time.date>tEnd), :) = [];
 Vaccination3_time((Vaccination3_time.date<tStart1 | Vaccination3_time.date>tEnd), :) = [];
 
 
 Recovered(time<datetime(tStart1) | time>datetime(tEnd)) = [];
 time(time<datetime(tStart1) | time>datetime(tEnd)) = [];
 Deaths = table2array(Deaths_time(:, 2));
 Daily = table2array(Daily_time(:, 2));
 
 Daily(Daily<=0) = NaN;
 nanx = isnan(Daily);
 t = 1:numel(Daily);
 Daily(nanx) = interp1(t(~nanx), Daily(~nanx), t(nanx));
 
 Vaccination1 = table2array(Vaccination1_time(:, 2));
 Vaccination2 = table2array(Vaccination2_time(:, 2));
 Vaccination3 = table2array(Vaccination3_time(:, 2));
 
 Vaccination11 = fillmissing(Vaccination1, "linear", "EndValues", "nearest");
 Vaccination22 = fillmissing(Vaccination2, "linear", "EndValues", "nearest");
 Vaccination33 = fillmissing(Vaccination3, "linear", "EndValues", "nearest");
 
 Vaccination111 = diff(Vaccination11);
 Vaccination222 = diff(Vaccination22);
 Vaccination333 = diff(Vaccination33);
 
 vac_len1 = length(Daily)-length(Vaccination111);
 vac_len2 = length(Daily)-length(Vaccination222);
 vac_len3 = length(Daily)-length(Vaccination333);
 
 Vaccination1111 = [repelem(0,vac_len1), transpose(Vaccination111)];
 Vaccination2222 = [repelem(0,vac_len1), transpose(Vaccination222)];
 Vaccination3333 = [repelem(0,vac_len1), transpose(Vaccination333)];
 
 Confirmed = cumsum(Daily);
 
 [Q,R,D,newT] = getMultipleWaves_5(guess,Npop,time,Confirmed,Recovered,Deaths,...
     Vaccination1111,Vaccination2222,Vaccination3333,tStart1,tIndex,pIndex,fit,est_day);
 
 Daily_predicted = Q+R+D;
 
  for kk = 2:length(tIndex)
     
     dup_index0 = find(newT==datetime(tIndex(kk)));
     dup_index = dup_index0(1);
     newT(dup_index) = [];
     Daily_predicted(dup_index) = [];
     D(dup_index) = [];
  end 
 
 % 7 days prediction
 
 
 
 Daily_predicted_0 = Daily_predicted(:,1:(end-(est_day*24+1)*2));
 Daily_predicted_1 = Daily_predicted(:,[1:(end-(est_day*24+1)*3) (end-(est_day*24+1)*2+1):(end-(est_day*24+1))]);
 Daily_predicted_2 = Daily_predicted(:,[1:(end-(est_day*24+1)*3) (end-(est_day*24+1)+1):(end)]);
 
 D_0 = D(:,1:(end-(est_day*24+1)*2));
 dD = gradient(D_0)./gradient(datenum(newT(1:end-(est_day*24+1)*2)));
 
 dy = gradient(Daily_predicted_0)./gradient(datenum(newT(1:end-(est_day*24+1)*2)));  
 dy_tot = [dy_tot, dy'];
 dy_smoothed = smoothdata(dy, 'movmean', est_day*24);
 dy1 = gradient(Daily_predicted_1)./gradient(datenum(newT(1:end-(est_day*24+1)*2)));
 dy1_smoothed = smoothdata(dy1, 'movmean', est_day*24);
 dy2 = gradient(Daily_predicted_2)./gradient(datenum(newT(1:end-(est_day*24+1)*2)));
 dy2_smoothed = smoothdata(dy2, 'movmean', est_day*24);
 
 
%  f = figure;
%  f.Position = [500 300 1200 600];
%  semilogy(time,Daily,'b.');
%  hold on
%  semilogy(newT(1:end-(est_day*24+1)*3), dy(1:end-(est_day*24+1)), 'k');
%  hold on
%  semilogy(newT(end-(est_day*24+1)*3+1:end-(est_day*24+1)*2), dy(end-(est_day*24):end), 'r','LineWidth',2);
%  hold on
%  semilogy(newT(end-(est_day*24+1)*2+1:end-(est_day*24+1)), dy1(end-(est_day*24):end), 'g','LineWidth',2);
%  hold on
%  semilogy(newT(end-(est_day*24+1)+1:end), dy2(end-(est_day*24):end), 'y','LineWidth',2);
%  hold on
%  xline(datetime(tIndex(datetime(tIndex)>=datetime(tIndex(end))-364)), ":")
%  ylabel('Number of cases')
%  xlabel('time (days)')
%  leg = {'Reported', 'Fitted', 'Predicted (Maintained)',...
%      'Predicted (reinforced)','Predicted (relaxed)'};
%  legend(leg{:},'location','eastoutside')
%  set(gcf,'color','w')
%  set(gca,'yscale','lin')
%  title(['Location: ', "Korea"])
%  hold off
%  shg


data_cases = [data_cases, round(dy(end-hours(to-from):24:end)'), round(dD(end-hours(to-from):24:end)')];
%% SEIQRDV3P + omicron (k=7)
 
 pIndex(14:length(pIndex)) = 7;
 
 Deaths_time = data_v(indV, [4,8]);
 Daily_time = data_v(indV, [4,6]);
 
 Vaccination1_time = data_v(indV, [4,36]);
 Vaccination2_time = data_v(indV, [4,37]);
 Vaccination3_time = data_v(indV, [4,38]);
 
 Deaths_time((Deaths_time.date<tStart1 | Deaths_time.date>tEnd), :) = [];
 Daily_time((Daily_time.date<tStart1 | Daily_time.date>tEnd), :) = [];
 
 Vaccination1_time((Vaccination1_time.date<tStart1 | Vaccination1_time.date>tEnd), :) = [];
 Vaccination2_time((Vaccination2_time.date<tStart1 | Vaccination2_time.date>tEnd), :) = [];
 Vaccination3_time((Vaccination3_time.date<tStart1 | Vaccination3_time.date>tEnd), :) = [];
 
 
 Recovered(time<datetime(tStart1) | time>datetime(tEnd)) = [];
 time(time<datetime(tStart1) | time>datetime(tEnd)) = [];
 Deaths = table2array(Deaths_time(:, 2));
 Daily = table2array(Daily_time(:, 2));
 
 Daily(Daily<=0) = NaN;
 nanx = isnan(Daily);
 t = 1:numel(Daily);
 Daily(nanx) = interp1(t(~nanx), Daily(~nanx), t(nanx));
 
 Vaccination1 = table2array(Vaccination1_time(:, 2));
 Vaccination2 = table2array(Vaccination2_time(:, 2));
 Vaccination3 = table2array(Vaccination3_time(:, 2));
 
 Vaccination11 = fillmissing(Vaccination1, "linear", "EndValues", "nearest");
 Vaccination22 = fillmissing(Vaccination2, "linear", "EndValues", "nearest");
 Vaccination33 = fillmissing(Vaccination3, "linear", "EndValues", "nearest");
 
 Vaccination111 = diff(Vaccination11);
 Vaccination222 = diff(Vaccination22);
 Vaccination333 = diff(Vaccination33);
 
 vac_len1 = length(Daily)-length(Vaccination111);
 vac_len2 = length(Daily)-length(Vaccination222);
 vac_len3 = length(Daily)-length(Vaccination333);
 
 Vaccination1111 = [repelem(0,vac_len1), transpose(Vaccination111)];
 Vaccination2222 = [repelem(0,vac_len1), transpose(Vaccination222)];
 Vaccination3333 = [repelem(0,vac_len1), transpose(Vaccination333)];
 
 Confirmed = cumsum(Daily);
 
 [Q,R,D,newT] = getMultipleWaves_5(guess,Npop,time,Confirmed,Recovered,Deaths,...
     Vaccination1111,Vaccination2222,Vaccination3333,tStart1,tIndex,pIndex,fit,est_day);
 
 Daily_predicted = Q+R+D;
 
  for kk = 2:length(tIndex)
     
     dup_index0 = find(newT==datetime(tIndex(kk)));
     dup_index = dup_index0(1);
     newT(dup_index) = [];
     Daily_predicted(dup_index) = [];
     D(dup_index) = [];
  end 
 
 % 7 days prediction
 
 
 
 Daily_predicted_0 = Daily_predicted(:,1:(end-(est_day*24+1)*2));
 Daily_predicted_1 = Daily_predicted(:,[1:(end-(est_day*24+1)*3) (end-(est_day*24+1)*2+1):(end-(est_day*24+1))]);
 Daily_predicted_2 = Daily_predicted(:,[1:(end-(est_day*24+1)*3) (end-(est_day*24+1)+1):(end)]);
 
 D_0 = D(:,1:(end-(est_day*24+1)*2));
 dD = gradient(D_0)./gradient(datenum(newT(1:end-(est_day*24+1)*2)));
 
 dy = gradient(Daily_predicted_0)./gradient(datenum(newT(1:end-(est_day*24+1)*2)));
 dy_tot = [dy_tot, dy'];
 dy_smoothed = smoothdata(dy, 'movmean', est_day*24);
 dy1 = gradient(Daily_predicted_1)./gradient(datenum(newT(1:end-(est_day*24+1)*2)));
 dy1_smoothed = smoothdata(dy1, 'movmean', est_day*24);
 dy2 = gradient(Daily_predicted_2)./gradient(datenum(newT(1:end-(est_day*24+1)*2)));
 dy2_smoothed = smoothdata(dy2, 'movmean', est_day*24);
 
 
%  f = figure;
%  f.Position = [500 300 1200 600];
%  semilogy(time,Daily,'b.');
%  hold on
%  semilogy(newT(1:end-(est_day*24+1)*3), dy(1:end-(est_day*24+1)), 'k');
%  hold on
%  semilogy(newT(end-(est_day*24+1)*3+1:end-(est_day*24+1)*2), dy(end-(est_day*24):end), 'r','LineWidth',2);
%  hold on
%  semilogy(newT(end-(est_day*24+1)*2+1:end-(est_day*24+1)), dy1(end-(est_day*24):end), 'g','LineWidth',2);
%  hold on
%  semilogy(newT(end-(est_day*24+1)+1:end), dy2(end-(est_day*24):end), 'y','LineWidth',2);
%  hold on
%  xline(datetime(tIndex(datetime(tIndex)>=datetime(tIndex(end))-364)), ":")
%  ylabel('Number of cases')
%  xlabel('time (days)')
%  leg = {'Reported', 'Fitted', 'Predicted (Maintained)',...
%      'Predicted (reinforced)','Predicted (relaxed)'};
%  legend(leg{:},'location','eastoutside')
%  set(gcf,'color','w')
%  set(gca,'yscale','lin')
%  title(['Location: ', "Korea"])
%  hold off
%  shg


data_cases = [data_cases, round(dy(end-hours(to-from):24:end)'), round(dD(end-hours(to-from):24:end)')];
%% data export

data_confirmed = table2timetable(array2table(data_cases(:,[1 2 5 7 9 11 13 15 17])), "RowTimes", data_date);
data_deaths = table2timetable(array2table(data_cases(:,[3 4 6 8 10 12 14 16 18])), "RowTimes", data_date);

data_confirmed.Properties.DimensionNames{1} = 'Date'; 
data_deaths.Properties.DimensionNames{1} = 'Date'; 

writetimetable(data_confirmed, strcat("paper_", Iso_code, "_new_cases_SEIR_raw.csv"))
% writetimetable(data_deaths, strcat("\Users\ojhaha\Desktop\설대\랩인턴\COVID-19 수리모델링\code-2\paper\", Iso_code, "_new_deaths_SEIR_raw.csv"))
%% MSE calculation

MSE = [];
RMSE = [];

y = table2array(data_confirmed(end-est_day+1:end, 1));

for jj = 1:8
    y_hat = table2array(data_confirmed(end-est_day+1:end, jj+1));
    mse = mean((y-y_hat).^2);
    MSE(jj) = mse;
    RMSE(jj) = sqrt(mse);
end

% figure
% plot(1:8, RMSE);
% shg

measure = array2table([MSE; RMSE]);

writetable(measure, strcat("measure_", Iso_code, "_new_cases_SEIR_raw.csv"));
%% Figures

f = figure;
f.Position = [500 300 1200 600];
DD = [Daily;Daily_test];
delay = 0; % move x-axis rage
semilogy(data_date(324+delay:end),DD(324+delay:end),'b.');
hold on
semilogy(newT(7753+delay*24:end-(est_day*24+1)*2), dy_tot(end-3288+delay*24:end,1), 'k');
hold on
semilogy(newT(7753+delay*24:end-(est_day*24+1)*2), dy_tot(end-3288+delay*24:end,2), 'r');
hold on
semilogy(newT(7753+delay*24:end-(est_day*24+1)*2), dy_tot(end-3288+delay*24:end,3), 'g');
hold on
semilogy(newT(7753+delay*24:end-(est_day*24+1)*2), dy_tot(end-3288+delay*24:end,4), 'y');
hold on
semilogy(newT(7753+delay*24:end-(est_day*24+1)*2), dy_tot(end-3288+delay*24:end,5), 'b');
hold on
ylabel('Number of cases')
xlabel('time')
leg = {'Reported', 'SEIR', "SEIQR", "SEIQRDVUP", "SEIQRDVP", "SEIQRDV3P"};
legend(leg{:},'location','northwest')
set(gcf,'color','w')
set(gca,'yscale','lin')
hold off
shg

f = figure;
f.Position = [500 300 1200 600];
DD = [Daily;Daily_test];
delay = 0; % move x-axis rage
semilogy(data_date(324+delay:end),DD(324+delay:end),'b.');
hold on
semilogy(newT(7753+delay*24:end-(est_day*24+1)*2), dy_tot(end-3288+delay*24:end,1), 'k');
hold on
semilogy(newT(7753+delay*24:end-(est_day*24+1)*2), dy_tot(end-3288+delay*24:end,5), 'r');
hold on
ylabel('Number of cases')
xlabel('time')
datetick('x','yyyy-mm')
axis auto
leg = {'Reported', 'SEIR', "SEIQRDV3P"};
legend(leg{:},'location','northwest')
set(gcf,'color','w')
set(gca,'yscale','lin')
hold off
shg


f = figure;
f.Position = [500 300 1200 600];
DD = [Daily;Daily_test];
delay = 0; % move x-axis range
semilogy(data_date(324+delay:end),DD(324+delay:end),'b.');
hold on
semilogy(newT(7753+delay*24:end-(est_day*24+1)*2), dy_tot(end-3288+delay*24:end,5), 'k');
hold on
semilogy(newT(7753+delay*24:end-(est_day*24+1)*2), dy_tot(end-3288+delay*24:end,6), 'r');
hold on
semilogy(newT(7753+delay*24:end-(est_day*24+1)*2), dy_tot(end-3288+delay*24:end,7), 'g');
hold on
semilogy(newT(7753+delay*24:end-(est_day*24+1)*2), dy_tot(end-3288+delay*24:end,8), 'y');
hold on
ylabel('Number of cases')
xlabel('time')
leg = {'Reported', "k=1", "k=3", "k=5", "k=7"};
legend(leg{:},'location','northwest')
set(gcf,'color','w')
set(gca,'yscale','lin')
hold off
shg

% f = figure;
% f.Position = [500 300 1200 600];
% DD = [Daily;repelem(0,est_day)'];
% delay = 100;
% semilogy((datetime("2021-12-29"):1:datetime("2022-02-04"))',DD(324+delay:end-7),'b.');
% hold on
% semilogy(newT(7753+delay*24:end-(est_day*24+1)*2), dy_tot(end-3288+delay*24-24*7:end,5), 'k');
% hold on
% semilogy(newT(7753+delay*24:end-(est_day*24+1)*2), dy_tot(end-3288+delay*24-24*7:end,6), 'r');
% hold on
% semilogy(newT(7753+delay*24:end-(est_day*24+1)*2), dy_tot(end-3288+delay*24-24*7:end,7), 'g');
% hold on
% semilogy(newT(7753+delay*24:end-(est_day*24+1)*2), dy_tot(end-3288+delay*24-24*7:end,8), 'y');
% hold on
% ylabel('Number of cases')
% xlabel('time')
% leg = {'Reported', "k=1", "k=3", "k=5", "k=7"};
% legend(leg{:},'location','northwest')
% set(gcf,'color','w')
% set(gca,'yscale','lin')
% hold off
% shg