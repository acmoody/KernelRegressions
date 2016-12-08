function all_data = assemble_KR_data()
% ASSEMBLE_KR_DATA - 
%   

all_data = cell( 6, 1 );

for sitecode = 1:6

    t0 = now();
    
    data07 = get_kernel_regression_data( sitecode, 2007 );
    data08 = get_kernel_regression_data( sitecode, 2008 );
    data09 = get_kernel_regression_data( sitecode, 2009 );
    data10 = get_kernel_regression_data( sitecode, 2010 );
    data11 = get_kernel_regression_data( sitecode, 2011 );

    data = vertcat( data07, data08, data09, data10, data11 );
    all_data{ sitecode } = data;

    t_secs = ( now() - t0 ) * 86400;
    fprintf( 'done %s (%d)\n', char( UNM_sites( sitecode ) ), t_secs );
end

save( 'kernel_regression_parsed_data.mat', 'all_data' );