function [] = plot_fit( data, fittype )
% Following convention, 'data' is of form [ label, x1, x2 ]
% fittype is the same as in fitsubj

[ loglik, bic, params ] =  fitsubj( fittype, data )

Aindices = data( :, 1) ==1;
Bindices = ~ Aindices;
newplot
hold on
plot( data(Aindices,2), data(Aindices, 3 ), 'ob' );
plot( data(Bindices,2), data(Bindices, 3 ), 'xr' );
switch fittype
    case 2
        interc = - params(3) / params( 2 )
        plot( [interc, interc], [0 1], 'k' )
    case 3
        interc = - params(3) / params( 2 )
        plot( [0 1], [interc, interc], 'k' )
    case 4
        plot2dlinbnd( params(2:length(params)), 'k', [0 1 0 1] )
end
axis([0 1 0 1])
hold off

params
end

