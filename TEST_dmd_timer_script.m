%clear all
filePrefix = 'default-';
fileDir = 'guppy_save';

set(0,'DefaultFigureColormap',jet);

a = timer('ExecutionMode', 'fixedSpacing', ...
    'Period', 1, ...
    'StartFcn', 'disp(''Lets track the DMD..'')'); %, ...
    %'StopFcn', 'sitter');

% set(a, 'TimerFcn', ...
%     { @dmd_timer_fcn, fileDir, filePrefix});
set(a, 'TimerFcn', 'dmd_timer_fcn(fileDir, filePrefix)');

get(a)
start(a)
