function stitch_together_evolve_langevin_files(outfilename,infileprefix,...
                                               nFiles,idcsFilesRequired)
    %% Load files infileprefix_data_x_of_nFiles*.mat with x in idcsFilesRequired
    %  stitch together mps_series,time_series,noise_series, and save in 
    %  outfilename
    if nargin==3, idcsFilesRequired=1:nFiles; end
    
    mps_series = {};
    time_series = [];
    noise_series = [];
    % noise series will have to take its shape from first file
    
    % get intervals data for mps_series,time_series,noise_series
    first_time_flag = true;
    for idx = idcsFilesRequired
        idxPrefix = sprintf('%s_data_%i_of_%i',infileprefix,idx,nFiles);
        files = dir([idxPrefix,'*.mat']);
        if length(files)==0
            error('No files found with prefix %s',idxPrefix)
        elseif length(files)~=1
            error('More than one file found with prefix %s',idxPrefix)
        end
        
        p=load(sprintf('%s\\%s',files.folder,files.name));
        
        mps_series = {mps_series{:}, p.mps_series{2:end}};
        time_series = [time_series, p.time_series(2:end)];  
        
        if first_time_flag
            noise_series = p.noise_series;
            first_time_flag = false;    
        else
            for n=1:length(noise_series)
                for k=1:length(noise_series{n})
                    noise_series{n}{k} = ...
                        [noise_series{n}{k}, p.noise_series{n}{k}(2:end)];
                end
            end
        end
    end
    
    % get initial workspace data, prepend initial time and mps
    prefix = sprintf('%s_initial_workspace',infileprefix);
    files = dir([prefix,'*.mat']);
    if isempty(files)
        error('No files found with prefix %s',prefix)
    elseif length(files)~=1
        error('More than one file found with prefix %s',prefix)
    end
    
    init_workspace = load(sprintf('%s\\%s',files.folder,files.name));
    mps_series = {init_workspace.mpsIn, mps_series{:}};
    time_series = [0 time_series];
    
    lastwarn('') % clear memory of warnings
    save(outfilename,'mps_series','time_series','noise_series','init_workspace');
    [~,id] = lastwarn;
    if isequal(id, 'MATLAB:save:sizeTooBigForMATFile')
        % file is too big to used compressed filetype
        save(outfilename,'mps_series','time_series','noise_series','init_workspace','-v7.3')
    end
end