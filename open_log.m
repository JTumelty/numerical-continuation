function sd = open_log(sd)
% OPEN_LOG opens the log file and writes the header line.
sd.fileID = fopen(sd.log,'w');
fprintf(sd.fileID,'\t its \t r \t ds \t Newton iteration \t flag \t tol \t residual \t F(u;r) \t iterations needed \n');
end