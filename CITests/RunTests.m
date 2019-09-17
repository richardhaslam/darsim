% RunTests.m
exit_code = 0;

try
    disp('Running tests');    
catch ME
    disp(getReport(ME))
    exit_code = 1;
end

% Ensure that we ALWAYS call exit
exit(exit_code);

