function dmd_timer_fcn(fileDir, filePrefix)
%(src, eventdata, fileDir, filePrefix)
% This function is meant to be used in conjunction with
% 'atom_statistics_timer.m' When run repeatedly, it will update keep updating
% new incoming files.
%disp('starting atom_timer_fcn.m');

try
 delete('\\rubidium\recycle_bin\greinerlab\Data\DMD_tracking\guppy_save\*.bin');
catch err
end

% persistent dir_length;
persistent inactive;

if isempty(inactive)
    inactive = 0;
end
%%% Check for files %%%
files = dir([fullfile(fileDir, filePrefix) '*dark.bin']);

if isempty(files)
    fprintf('.');
    inactive = inactive +1;
    
    if inactive > 500
        disp('timed out')
        evalin('base', 'stop(a)')
    end
    
else

% if isempty(inactive)
%     inactive = 0;
% end
%     inactive = inactive +1
% if numel(files) > dir_length
%     inactive = 0;
%     dir_length = numel(files);
% elseif numel(files) == dir_length
%     inactive = inactive + 1;
%     fprintf('.')
% elseif numel(files) < dir_length
%     disp('something wrong... fewer files?')
% else
%     disp('something very wrong')
% end

% if inactive > 500
%     disp('timed out')
%     evalin('base', 'stop(a)')
% end


% found_a_file = 0;
loop_index = 1;
% Let's go
% while ~found_a_file 
%     disp('1')
    filename = fullfile(fileDir, files(loop_index).name);
%     disp('2')
    [pathstr, name, ext] = fileparts(filename);
%     current_filePrefix = fullfile(pathstr, name);
%     
    string = filename(1:length(filename)-length('-dark.bin'))
%     dummyfilename = [string '.dummy'];

%     if exist(dummyfilename) ~= 2
        % run script that creates the dummy file
        
        dmd_tracking_fn(string);
%         found_a_file = 1;
        inactive = 0;
%     else
% 	loop_index = loop_index + 1;
%     end

%     if loop_index > numel(files)
%         found_a_file = 1;
%     end

% end

end
%disp('done with atom_timer_fcn.m');
