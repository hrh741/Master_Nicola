mkdir(directoryOut);
logFID = fopen([directoryOut,'Log.txt'],'at');
fprintf(logFID,'Simulation from: %s \n',datestr(now));
fprintf(logFID,'Working directory: %s \n\n',pwd);
