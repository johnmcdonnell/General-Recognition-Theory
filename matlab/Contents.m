% General Recognition Theory Toolbox
%
% Version 2.5  10-March-05  by Leola A. Alfonso-Reese
% Copyright (c) 1995-2005
%
% The GRT Toolbox is a collection of routines and demo scripts for designing and
%   analyzing data based on Ashby's General Recognition Theory (1992) involving
%   two categories.
%
%   Note:  Some analysis routines depend on the Optimization Toolbox for MATLAB.
%
% If you find that this software is helpful for conducting your research
% and you would like to acknowledge use of this toolbox, please cite me
% with the following reference.
%
%    Alfonso-Reese, L.A. (2005) General Recognition Theory Toolbox for MATLAB.
%    Manuscript in preparation.
%
% An example of a citation is as follows: 
%    "Stimuli were generated with the aid of the General Recognition Theory
%    Toolbox (Alfonso-Reese, 2005)."
%
% General Functions
%   bivnormpdf         - compute bivariate normal probability density function
%   bnd2slint          - convert linear bnds to slope & y-intercept pairs (2D)
%   fastnormalcdf      - compute normal cumulative distribution function quickly
%   lindecisbnd        - compute ideal linear decision bound coefficients (nD)
%   linprobcorr        - compute probability correct for optimal linear bound (nD)
%   normalcdf          - compute normal cumulative distribution function
%   normalcdfinv       - compute inverse of normal cumulative distribution
%   normalpdf          - compute normal probability density function
%   quaddecisbnd       - compute ideal quadratic decision bound coefficients (nD)
%
% Stimulus Generation Functions
%   clip_vals          - clip values that are outside of specified range (nD)
%   gensample          - generate normally distributed data points (nD)
%   randrows           - randomize order of rows in a matrix
%   transample         - transform sample params to match desired params (nD)
%
% Subject Simulation Functions
%   lindiscrim1dvals
%   lindiscrim2dvals   - compute linear discriminant values (2D,3D)
%   lindiscrim3dvals
%   sim2dlin           - simulate observer's response data using linear bnd (2D,3D)
%   sim3dlin
%
% Data Analysis Functions
%   angsmat2xyzmat     - convert matrix of angle pairs to x,y,z coordinate matrix 
%   angvec2xymat       - convert vector of angles to x,y coordinate matrix 
%   dprime             - compute d-prime from data (nD)
%   dprimef            - compute d-prime from population parameters (nD)
%   fisherdiscrim1d
%   fisherdiscrim2d    - compute coefficients of Fisher's discriminant fcn (2D,3D)
%   fisherdiscrim3d
%   fit_1dGLC
%   fit_1dnoise
%   fit_2dGLC          - set up and begin search for the GLC (2D,3D)
%   fit_3dGLC
%   fit_2dGQC          - set up and begin search for the GQC (2D,3D)
%   fit_3dGQC
%   getcovs
%   getmeans
%   negloglike_1dGLC
%   negloglike_2dGLC   - compute -loglikelihood of the data for GLC model (2D,3D)
%   negloglike_3dGLC
%   negloglike_2dGQC   - compute -loglikelihood of the data for GQC model (2D,3D)
%   negloglike_3dGQC
%   new2old_2dparams   - convert new format of search parameters to old (2D,3D)
%   new2old_3dparams
%   norm_old_1dparams
%   norm_old_2dparams  - normalize the old parameters:  [noise avec b] (2D,3D)
%   norm_old_3dparams
%   old2new_2dparams   - convert old format of search parameters to new (2D,3D)
%   old2new_3dparams
%   percorr            - determine percent of responses correctly categorized (nD)
%   quadbndpercorr
%   quaddecisbnd
%   reducemat
%   reducevec
%   sphdecisbnd
%   xymat2angvec       - convert matrix of x,y coordinates to vector of angles
%   xyzmat2angmat      - convert matrix of x,y,z coordinates to matrix of angle pairs
%
% Data Graphing Functions
%   colorstr2colorval
%   colorval2colorstr
%   getbndendpts       - return two endpoints of a linear bound (2D)
%   plotbivnormcon     - plot a bivariate normal equivocality contour
%   plot2bivnorms      - plot two bivariate normal probability density functions
%   plot2dlinbnd       - plot a linear boundary on the current figure (2D,3D)
%   plot3dlinbnd
%   plot2dquadbnd      - plot a quadratic boundary on the current figure (2D,3D)
%   plot3dquadbnd
%   plot1dresp
%   plot2dstim
%   plot2dstimcustom
%   plot1dstim
%   plot2dresp         - plot response data points (2D,3D)
%   plot2drespcustom
%   plot3dresp
%   plot3drespcustom
%   plot2dstim         - plot stimulus data points (2D,3D)
%   plot3dstim
%
% Demonstrations
%   For examples of designing experiments, generating stimuli, simulating
%   subject data, analyzing and graphing data, type 'help GRTdemos.'
%
% GRT Toolbox Information
%   GRT_README         - general information about the GRT toolbox

