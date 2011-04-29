% Based on a demo by Alfonso-Reese, (c) 1995
%% Read in files
format compact;
clc;

[~,files] = unix( 'ls data/<1->.dat' )
[~,count] = sscanf( files, '%s' );

rest=files;

%% Get bics and write to a file
old = []
for i = 1:count
    [file, rest] = strtok( rest);
    file
    thesedata = load_file( file );

    %thesedata = thesedata( 1:20, : ); % TAKING FIRST 20 TRIALS

    bics = []
    for i = 0:4
        [ loglik, bic, params ] = fitsubj( i, thesedata ) 
        bics = [ bics bic ]
    end
    old = [ old; bics]
end

%file = strtok( files )
fid = fopen('fits.dat','w');
rest=files;
for i = 1:count
    [file, rest] = strtok( rest);
    fprintf(fid, strcat( file, ' %0.5f %0.5f %0.5f %0.5f %0.5f\n'), old( i, 1), old(i,2), old(i,3), old(i,4), old(i, 5)  );
end

