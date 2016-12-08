classdef climate_space
    properties
        sitecode;
        swc_val;   % N by 1; SWC bin centers
        T_val;     % M by 1; T bin centers
        swc_depth; % depth of SWC obs for these climate spaces
        year_idx;  % Y by 1; years for each "slice" of year_clim_space
        year_clim_space; % N by M by Y; the climate space data for each year
        alltime_clim_space; % N by M; the climate space across all years
        note; % user-entered text
    end
    
    methods
        
            function obj = climate_space( sitecode, ...
                                          swc_val, ...
                                          T_val, ...
                                          swc_depth, ...
                                          year_idx, ...
                                          year_clim_space, ...
                                          alltime_clim_space, ...
                                          note )
            % class constructor 
            obj.sitecode = sitecode;
            obj.swc_val = swc_val;
            obj.T_val = T_val;
            obj.swc_depth = swc_depth;
            obj.year_idx = year_idx;
            obj.year_clim_space = year_clim_space;
            obj.alltime_clim_space = alltime_clim_space;
            obj.note = note;
            
            end
    end
    end
