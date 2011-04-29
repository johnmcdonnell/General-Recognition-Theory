% Plot everything.
% Copyright 2010 John McDonnell

clf
fns = {}; % '1.dat', '7.dat', '11.dat', '13.dat', '17.dat', '23.dat', '27.dat', '29.dat', '33.dat' };
lsout = dir( 'data/*.dat' )
for i = 1:(length(lsout)-1)
    if regexp( lsout(1).name, '\d' )
        fns = [fns lsout(i).name ]
    end
end
n = length( fns )
xfigs = ceil( sqrt(n) )
yfigs = ceil( n / xfigs )

for ( i= 1:length(fns) )
    fns{i}
    subplot( xfigs, yfigs, i)
    plot_fit( load_file(['data/' fns{i}]), 4 )
    title( fns{i} )
end
