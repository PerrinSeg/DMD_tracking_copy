function dmd_tracking_fn(filename)

% pause for 6 seconds
% guppy software takes some time to write...
% should use condition to check file integrity...
pause(6);

recordlength = 501*10;
recordlengthV = recordlength;
recordlengthH = recordlength;

recordwidth = 8;
recordwidthV = recordwidth;
recordwidthH = recordwidth;

logV = zeros(recordlengthV,recordwidthV);
logH = zeros(recordlengthH,recordwidthH);

%%

need_overwriteH = 1;
need_overwriteV = 1;

if exist('dmd_tracking_log_hor.csv','file')==2
    log_tempH = csvread('dmd_tracking_log_hor.csv');
    if sum(size(log_tempH)==[recordlengthH,recordwidthH])==2
        need_overwriteH = 0;
        logH = log_tempH;
    end
end

if exist('dmd_tracking_log_vert.csv','file')==2
    log_tempV = csvread('dmd_tracking_log_vert.csv');
    if sum(size(log_tempV)==[recordlengthV,recordwidthV])==2
        need_overwriteV = 0;
        logV = log_tempV;
    end
end

if need_overwriteH==1 
    csvwrite('dmd_tracking_log_hor.csv',zeros(recordlengthH,recordwidthH));
end

if need_overwriteV==1
    csvwrite('dmd_tracking_log_vert.csv',zeros(recordlengthV,recordwidthV));
end

%% shifts for vertical direction

%folderName = 'dmd_stability_exp_running_flowBackOn'; %'FAKE_DMD_cal1';
%filePrefix = '4sdelay-'; 
%listing = dir([folderName '\' filePrefix  '*.fits']);
%filename = 'dmd_guppy\default-0005';

[isWrite shiftsV shiftsH] = dmd_fitting_fn(filename);
% shifts = [latticeDX, latticeDY, dotDX, dotDY]
% projection is the running average

% alex XXXX
shots_to_avg = 3; % 3
% disp(['dmd_tracking_fn - shots_to_avg: ' num2str(shots_to_avg)]);
shots_to_avg_test = 9;

%% change to projected value vert direction
shifts_prev_guessV = mean(logV(recordlengthV-shots_to_avg_test+1:recordlengthV,1:4),1);
shifts_prev_guessV = [shifts_prev_guessV;shifts_prev_guessV;shifts_prev_guessV];
shiftsV = shiftsV - round(shiftsV - shifts_prev_guessV);

%check if the shifts are changing (avoid duplicates from guppy) 
if sum(sum(abs(shiftsV - shifts_prev_guessV))) < 0.001
    display('Shifts are not changing. Check guppy images.')
    isWrite = 0;
end
disp(['isWrite = ' num2str(isWrite)]);

% if max(std(shifts,1,1))> 0.3
%     display('Signal unstable or fitting failed')
% else

% relativeDX = dotDX - latticeDX;

[relxV,ind] = sort(shiftsV(:,3) - shiftsV(:,1));                % if two points are displaced with respect to the third one by 
if find(diff(relxV) > 0.5) == 1; relxV(1) = relxV(1) + 1; end   % approximately one lattice constant, these lines ensure that
if find(diff(relxV) > 0.5) == 2; relxV(3) = relxV(3) - 1; end   % the outlier-point is properly translated, without changing
relxV(ind) = relxV;                                             % the order of the points (Robert)

[relyV,ind] = sort(shiftsV(:,4) - shiftsV(:,2));                % if two points are displaced with respect to the third one by
if find(diff(relyV) > 0.5) == 1; relyV(1) = relyV(1) + 1; end   % approximately one lattice constant, these lines ensure that
if find(diff(relyV) > 0.5) == 2; relyV(3) = relyV(3) - 1; end   % the outlier-point is properly translated, without changing
relyV(ind) = relyV;                                             % the order of the points (Robert)

shifts_allV = [shiftsV, relxV, relyV, zeros(3,2)];

% remove the oldest, append the new value
logV(1:3,:) = [];
logV = [logV; shifts_allV];

% calculate projection, to be write to txt file;
projectxV = mean(logV(recordlengthV-shots_to_avg+1:recordlengthV,5));
projectyV = mean(logV(recordlengthV-shots_to_avg+1:recordlengthV,6));

logV(recordlengthV-2:recordlengthV,7) = [projectxV;projectxV;projectxV];
logV(recordlengthV-2:recordlengthV,8) = [projectyV;projectyV;projectyV];

if isWrite == 1
% write log
csvwrite('dmd_tracking_log_vert.csv',logV);

% write txt
 dlmwrite('dmd_tracking_current_vert.txt',[projectxV,projectyV]);
end
% end

% write dummy to indicate done with file
%dummy = '';
%save([filename '.dummy'],'dummy');


%% shifts for horizontal direction
shifts_prev_guessH = mean(logH(recordlengthH-shots_to_avg_test+1:recordlengthH,1:4),1);
shifts_prev_guessH = [shifts_prev_guessH;shifts_prev_guessH;shifts_prev_guessH];
shiftsH = shiftsH - round(shiftsH - shifts_prev_guessH);

%check if the shifts are changing (avoid duplicates from guppy) 
if sum(sum(abs(shiftsH - shifts_prev_guessH))) < 0.001
    display('Shifts are not changing. Check guppy images.')
    isWrite = 0;
end
disp(['isWrite = ' num2str(isWrite)]);

% if max(std(shifts,1,1))> 0.3
%     display('Signal unstable or fitting failed')
% else

%relativeDX = dotDX - latticeDX;

[relxH,ind] = sort(shiftsH(:,3) - shiftsH(:,1));                % if two points are displaced with respect to the third one by 
if find(diff(relxH) > 0.5) == 1; relxH(1) = relxH(1) + 1; end   % approximately one lattice constant, these lines ensure that
if find(diff(relxH) > 0.5) == 2; relxH(3) = relxH(3) - 1; end   % the outlier-point is properly translated, without changing
relxH(ind) = relxH;                                             % the order of the points (Robert)

[relyH,ind] = sort(shiftsH(:,4) - shiftsH(:,2));                % if two points are displaced with respect to the third one by
if find(diff(relyH) > 0.5) == 1; relyH(1) = relyH(1) + 1; end   % approximately one lattice constant, these lines ensure that
if find(diff(relyH) > 0.5) == 2; relyH(3) = relyH(3) - 1; end   % the outlier-point is properly translated, without changing
relyH(ind) = relyH;                                             % the order of the points (Robert)

shifts_allH = [shiftsH, relxH, relyH, zeros(3,2)];

% remove the oldest, append the new value
logH(1:3,:) = [];
logH = [logH; shifts_allH];

% calculate projection, to be write to txt file;
projectxH = mean(logH(recordlengthH-shots_to_avg+1:recordlengthH,5));
projectyH = mean(logH(recordlengthH-shots_to_avg+1:recordlengthH,6));

logH(recordlengthH-2:recordlengthH,7) = [projectxH;projectxH;projectxH];
logH(recordlengthH-2:recordlengthH,8) = [projectyH;projectyH;projectyH];

if isWrite == 1
% write log
csvwrite('dmd_tracking_log_hor.csv',logH);

% write txt
 dlmwrite('dmd_tracking_current_hor.txt',[projectxH,projectyH]);
 dlmwrite('dmd_tracking_current_combine.txt',[projectxV,projectyV,projectxH,projectyH]);
end
% end

% write dummy to indicate done with file
%dummy = '';
%save([filename '.dummy'],'dummy');


%% plot


start_ind = 5000-100; % to plot less % FEWER!!!!****** moron

latticeDXV = logV(start_ind:end, 1);
latticeDYV = logV(start_ind:end, 2);
dotDXV = logV(start_ind:end, 3);
dotDYV = logV(start_ind:end, 4);
relativeDXV = logV(start_ind:end, 5);
relativeDYV = logV(start_ind:end, 6);

projectedDXV = logV(start_ind:end, 7);
projectedDYV = logV(start_ind:end, 8);



latticeDXH = logH(start_ind:end, 1);
latticeDYH = logH(start_ind:end, 2);
dotDXH = logH(start_ind:end, 3);
dotDYH = logH(start_ind:end, 4);
relativeDXH = logH(start_ind:end, 5);
relativeDYH = logH(start_ind:end, 6);

projectedDXH = logH(start_ind:end, 7);
projectedDYH = logH(start_ind:end, 8);


%% 02/08/2023 Yanfei Li
% warning function that sends email when DMD tracking crashed
% set_alarm(data, check the last n shots, trigger_status)
% when triggered once, one needs to reset trigger in setup_trigger.m

% n = 30; % Until 2023/08/15
n = 10;
trigger = fileread('Safety_trigger_status.txt');
set_alarm(relativeDXH, relativeDXV, relativeDYH, relativeDYV, n, trigger);


%%
%Plotting

range = 1;


% figure(3)
% % title(filename)
% subplot(3,3,1);
% scatter(latticeDXV,latticeDYV)
% axis equal
% axis square
% title ('lattice Vert')
% %xlim ([-0.1 0.1])
% %ylim ([-0.1 0.1])
% % title (['Scatter:' ' xRMS ' num2str(xRMS,3),' lattice sites, yRMS' num2str(xRMS,3),' lattice sites'])
% % xlabel ('xscatter in lattice sites')
% % ylabel ('ylabel in lattice sites')
% 
% subplot(3,3,2);
% scatter(dotDXV,dotDYV)
% title ('dot Vert')
% axis equal
% axis square
% %xlim ([-0.1 0.1])
% %ylim ([-0.1 0.1])
% 
% subplot(3,3,3);
% scatter(relativeDXV,relativeDYV)
% title ('relative Vert')
% axis equal
% axis square
% %xlim ([-0.1 0.1])
% %ylim ([-0.1 0.1])
% 
% subplot(3,3,4)
% plot(latticeDXV,'s')
% title ('latticeDX Vert (all in lattice sites)')
% xlim ([1 recordlengthV-start_ind+1])
% ylim ([latticeDXV(end)-range latticeDXV(end)+range])
% 
% subplot(3,3,7)
% plot(latticeDYV,'s')
% title ('latticeDY Vert')
% xlim ([1 recordlengthV-start_ind+1])
% ylim ([latticeDYV(end)-range latticeDYV(end)+range])
% 
% subplot(3,3,5)
% plot(dotDXV,'s')
% title ('dotDX Vert')
% xlim ([1 recordlengthV-start_ind+1])
% ylim ([dotDXV(end)-range dotDXV(end)+range])
% 
% subplot(3,3,8)
% plot(dotDYV,'s')
% title ('dotDY Vert')
% xlim ([1 recordlengthV-start_ind+1])
% ylim ([dotDYV(end)-range dotDYV(end)+range])
% 
% subplot(3,3,6)
% plot(relativeDXV,'s')
% hold on
% plot(projectedDXV,'-r','LineWidth',2)
% hold off
% title ('relativeDX Vert')
% xlim ([1 recordlengthV-start_ind+1])
% ylim ([projectedDXV(end)-range projectedDXV(end)+range])
% 
% subplot(3,3,9)
% plot(relativeDYV,'s')
% hold on
% plot(projectedDYV,'-r','LineWidth',2)
% hold off
% title ('relativeDY Vert')
% xlim ([1 recordlengthV-start_ind+1])
% ylim ([projectedDYV(end)-range projectedDYV(end)+range])





% figure(4)
% % title(filename)
% subplot(3,3,1);
% scatter(latticeDXH,latticeDYH)
% axis equal
% axis square
% title ('lattice Hor')
% %xlim ([-0.1 0.1])
% %ylim ([-0.1 0.1])
% % title (['Scatter:' ' xRMS ' num2str(xRMS,3),' lattice sites, yRMS' num2str(xRMS,3),' lattice sites'])
% % xlabel ('xscatter in lattice sites')
% % ylabel ('ylabel in lattice sites')
% 
% subplot(3,3,2);
% scatter(dotDXH,dotDYH)
% title ('dot Hor')
% axis equal
% axis square
% %xlim ([-0.1 0.1])
% %ylim ([-0.1 0.1])
% 
% subplot(3,3,3);
% scatter(relativeDXH,relativeDYH)
% title ('relative Hor')
% axis equal
% axis square
% %xlim ([-0.1 0.1])
% %ylim ([-0.1 0.1])
% 
% subplot(3,3,4)
% plot(latticeDXH,'s')
% title ('latticeDX Hor (all in lattice sites)')
% xlim ([1 recordlengthH-start_ind+1])
% ylim ([latticeDXH(end)-range latticeDXH(end)+range])
% 
% subplot(3,3,7)
% plot(latticeDYH,'s')
% title ('latticeDY Hor')
% xlim ([1 recordlengthH-start_ind+1])
% ylim ([latticeDYH(end)-range latticeDYH(end)+range])
% 
% subplot(3,3,5)
% plot(dotDXH,'s')
% title ('dotDX Hor')
% xlim ([1 recordlengthH-start_ind+1])
% ylim ([dotDXH(end)-range dotDXH(end)+range])
% 
% subplot(3,3,8)
% plot(dotDYH,'s')
% title ('dotDY Hor')
% xlim ([1 recordlengthH-start_ind+1])
% ylim ([dotDYH(end)-range dotDYH(end)+range])
% 
% subplot(3,3,6)
% plot(relativeDXH,'s')
% hold on
% plot(projectedDXH,'-r','LineWidth',2)
% hold off
% title ('relativeDX Hor')
% xlim ([1 recordlengthH-start_ind+1])
% ylim ([projectedDXH(end)-range projectedDXH(end)+range])
% 
% subplot(3,3,9)
% plot(relativeDYH,'s')
% hold on
% plot(projectedDYH,'-r','LineWidth',2)
% hold off
% title ('relativeDY Hor')
% xlim ([1 recordlengthH-start_ind+1])
% ylim ([projectedDYH(end)-range projectedDYH(end)+range])

figure(6)

subplot(2,2,1)
plot(relativeDXV,'s')
hold on
plot(projectedDXV,'-r','LineWidth',2)
hold off
title ('relativeDX Vert')
xlim ([1 recordlengthV-start_ind+1])
ylim ([projectedDXV(end)-range projectedDXV(end)+range])

subplot(2,2,3)
plot(relativeDYV,'s')
hold on
plot(projectedDYV,'-r','LineWidth',2)
hold off
title ('relativeDY Vert')
xlim ([1 recordlengthV-start_ind+1])
ylim ([projectedDYV(end)-range projectedDYV(end)+range])

subplot(2,2,2)
plot(relativeDXH,'s')
hold on
plot(projectedDXH,'-r','LineWidth',2)
hold off
title ('relativeDX Hor')
xlim ([1 recordlengthH-start_ind+1])
ylim ([projectedDXH(end)-range projectedDXH(end)+range])

subplot(2,2,4)
plot(relativeDYH,'s')
hold on
plot(projectedDYH,'-r','LineWidth',2)
hold off
title ('relativeDY Hor')
xlim ([1 recordlengthH-start_ind+1])
ylim ([projectedDYH(end)-range projectedDYH(end)+range])

end
    