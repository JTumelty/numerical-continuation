function res = contSH_main()
%CONTSH_MAIN is the main function that performs numerical continuation on
%the Swift--Hohenberg equation.
%
%OUTPUTS: 
%       res:        A structure containing all of the results from the 
%                   continuation 
%
%USAGE:
%   Add the directory containing the source files to the path, navigate to
%   the working directory containing a file all_param.mat and run the
%   command:
%     res = contSH_main()

% Validate properties stored within the structure all_param and set
% unspecified properties to default values.
all_param  = init_param();

% Save all of the parameters in the current working directory
save(strcat(all_param.sd.dir,'/final_param.mat'),'all_param');

% Initialise the results
all_res = init_all_res(all_param);

% Open a log file if necessary
if all_param.sd.log == 1
    all_param.sd  = open_log(all_param.sd);
end

% Perform the numerical continuation
all_res = run_cont(all_res,all_param);

% Save all of final results name
save(strcat(all_param.sd.dir,'/',all_param.sd.dataname),'all_res');

end

