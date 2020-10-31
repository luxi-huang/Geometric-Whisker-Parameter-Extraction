%% delet most recent row
% this script is used to delete the most recent rows on SealWhisker.mat 

data_filename = 'SealWhisker.mat';
ma = matfile(data_filename, 'Writable', true);

a = size(ma.row);

% delete row data;
ma.row(a(1), :)= [];

% delete col data 
ma.col(a(1), :)= [];

% delete harbor seal number 
ma.SealNum(a(1), :)= [];

% delete diameter and ratiodata  
ma.D_base(a(1), :)= [];
ma.D_tip(a(1), :)= [];
ma.Ratio_R(a(1), :)= [];

% delete whisker length
ma.length(a(1), :)= [];

% delete whisker data 

ma.whisker_xx(a(1), :)= [];

ma.whisker_yy(a(1), :)= [];


ma.std_upper_concave_xx(a(1), :)= [];

ma.std_upper_concave_yy(a(1), :)= [];

ma.std_upper_convex_xx(a(1), :)= [];

ma.std_upper_convex_yy(a(1), :)= [];


ma.std_lower_concave_xx(a(1), :)= [];

ma.std_lower_concave_yy(a(1), :)= [];

ma.std_lower_convex_xx(a(1), :)= [];

ma.std_lower_convex_yy(a(1), :)= [];

% delete distance from base to tip
ma.dis_tip_base(a(1), :)= [];

% delete ifValid value 
ma.IfValid(a(1)+1, :) = [];
