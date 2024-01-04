function [ ] = sitter()

disp('sitter called!')
listing = dir('guppy_save/*.bin');

if size(listing,1)> 4
    delete('guppy_save/*.bin');
    %delete(timerfind);
    dmd_timer_script
end

% delete('guppy_save/*.bin');
% dmd_timer_script

end