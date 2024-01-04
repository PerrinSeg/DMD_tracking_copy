delete(timerfind);
dmd_timer_script

disp('Babysit dmd tracking...');
bb = timer('ExecutionMode', 'fixedSpacing', ...
    'Period', 30, ...
    'StartFcn', 'sitter', ...
    'StopFcn', 'clear all');

set(bb, 'TimerFcn', 'sitter');

get(bb)
start(bb) 
