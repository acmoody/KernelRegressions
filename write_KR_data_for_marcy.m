function result = write_KR_data_for_marcy()
% WRITE_KR_DATA_FOR_MARCY - 
%   

load( 'tim_out.mat' );
for this_site = 1:6
    this_write = write_single_site( this_site, tim_out{ this_site } );
end

result = 0;

%--------------------------------------------------
function result = write_single_site( sitecode, data )
% WRITE_SINGLE_SITE - 
%   

data( :, 3:4 ) = [];  % this is the mean hour -- meaningless

names = { 'year', 'day', 'air_T', 'precip', 'swc_0_5', 'swc_10_20', ...
          'swc_20_30', 'swc_30plus', 'NEE', 'GPP', 'RE', 'NEE_day', ...
          'NEE_night' };

% write the header
fname = sprintf( '%s_KR_data.txt', char( UNM_sites( sitecode ) ) );
outfile = fopen( fname, 'w' );
fprintf( outfile, '%s\t', names{ : } );
fprintf( outfile, '\n' );
fclose( outfile );

%write the data
dlmwrite( fname, data, '-append', 'delimiter', '\t' );

fprintf( 'wrote %s\n', fname );
result = 0;