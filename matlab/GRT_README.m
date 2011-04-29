echo on;
%          INTRODUCTION TO THE GENERAL RECOGNITION THEORY TOOLBOX
%          ======================================================
%
%                             Version 2.5
%
%                       by Leola A. Alfonso-Reese
%
%                          San Diego State University
%                          5500 Campanile Drive
%                          San Diego, CA  92182-4611
%
%                        Copyright (c) 1995-2005
%                     
%
% This is a package of MATLAB programs to generate and analyze data for
% psychology experiments based on the General Recognition Theory.
%
%
%                               Notes
%                               -----
%
% This toolbox is free and the author makes no warranty of any kind with
% regard to the programs or online help contained herein.  The author,
% Leola A. Alfonso-Reese, is not liable for any consequences of installing
% or running these files.
%
% The M-files and any modifications of the M-files contained in this
% toolbox may not be sold, published, or used for commercial gain.
%
% You may copy and modify the GRT Toolbox files as long as the original
% author is acknowledged AND you document any modifications you make
% as your own.  ALSO, IF THIS TOOLBOX IS USED IN ANY PUBLISHED RESEARCH,
% I WOULD APPRECIATE AN ACKNOWLEDGEMENT.  For example, using APA style,
% the toolbox may be referenced as follows:
%
%    Alfonso-Reese, L.A. (2005) General Recognition Theory Toolbox for
%    MATLAB.  Manuscript submitted for publication. 
%
% Address any questions, comments, suggestions, and ESPECIALLY BUG
% REPORTS to:
%
%       Leola A. Alfonso-Reese
%       Department of Psychology
%       San Diego State University
%       5500 Campanile Drive
%       San Diego, CA  92182-4611
%
%       leola@alum.mit.edu
%
% I would enjoy email regarding your reaction to this toolbox.
%
% January 2005
%
%                              Thanks
%                              ------
%
% Thanks to David H. Brainard for helping me get started with MATLAB
% and for teaching me how to speed up my code.
% 
% Thanks to F. Gregory Ashby for many discussions of General Recognition
% Theory.
% 
% Special thanks to Michael Steele Reese for originally suggesting that
% I make GRT Toolbox!
%
%
%                        Installation Notes
%                        ------------------
%
% The toolbox is written to be used with the Macintosh.  Two versions
% are available, one for Mac OS 9.2 running MATLAB 5.2.1, and another
% for Mac OS X running MATLAB 7.0.  Since all the files are simply text
% files, they may also be installed on a UNIX machine or an IBM PC.
% Finally, some data analysis functions depend on the Optimization Toolbox
% which can be purchased from the MathWorks, Inc.
%
% I suggest that the GRT folder, containing many subroutines and a
% GRTdemos subdirectory, be placed with other toolboxes in MATLAB's
% toolbox folder.  Also, the GRT directory and its GRTdemos subdirectory
% must be added to the MATLAB path.
%
%
%                        After Installation
%                        ------------------
%
% For a brief description of all the functions in the GRT Toolbox, type
%     help GRT
%
% For a brief description of the demonstrations, type
%     help GRTdemos
%
% For a description of each GRT function, type 'help' followed by the
% function name, e.g.,
%     help dprime
%
% To see the contents of a function name or demonstration file, type
% 'type' followed by the function or demo name, e.g.,
%     type Comp2dStats
% 
%
%                 Background on General Recognition Theory
%                 ----------------------------------------
%
% General Recognition Theory (GRT; Ashby, 1992) is a multi-dimensional
% version of signal detection theory (Thurstone, 1927; Green & Swets,
% 1966). It was developed by psychologists (Ashby & Townsend, 1986;
% Ashby & Gott, 1988) to study the human categorization process.
%
% In GRT, a stimulus is represented by a point in multi-dimensional
% space. When a subject observes a multi-dimensional stimulus, the
% sensory and perceptual system adds noise to the signal.  This signal
% plus noise is called the perceptual effect.  The perceptual effects
% associated with a particular stimulus are represented by a multi-
% dimensional probability distribution.  A category is represented by
% a mixture of the associated stimulus distributions.  In a categorization
% task, the observer learns to separate the multi-dimensional space into
% regions, one region per category.  Then, when the observer needs to
% respond to a stimulus, he determines the region in which the percept
% falls and emits the category label associated with that region.  Thus,
% the boundaries separating these regions represent the decision rules
% the observer uses to classify stimuli.
%
% Since perceptual effects are not directly observable, it is difficult
% for an experimenter to examine the decision boundaries that subjects
% use to categorize stimuli.  The General Recognition Randomization
% Technique developed by Ashby and Gott (1988) provides a solution to
% this problem.  First, an experimenter randomly generates stimuli by
% adding external noise to each category prototype, where the external
% noise varies according to some distribution.  (The external noise must
% be large enough to make the perceptual noise negligible in comparison.)
% Then, the stimuli are presented to a subject one at a time at high
% contrast and long duration; this further reduces perceptual noise.  For
% each stimulus presentation, the subject's category response is recorded
% and if the task is a supervised one, then the subject receives feedback
% indicating whether his response was correct or not.  After all the
% responses are collected, the experimenter represents each exemplar as
% a point in space, determines the decision boundary that best predicts
% the subjects' responses, and uses this boundary to study the human
% categorization process.
%
%
%            Studying Human Categorization with the GRT Toolbox
%            --------------------------------------------------
%
% This toolbox provides routines to help an experimenter design 
% categorization (or identification) experiments, generate stimuli for
% these experiments, simulate a subject's responses, analyze categorization
% data, and graph results.  The typical user of this toolbox designs
% experiments for two-category tasks, where the categories are specified
% by multi-variate normal distributions.  This version of the toolbox
% provides tools for fitting the General Linear Classifier and the
% General Quadratic Classifier to a data set.
%
% The functions currently provided in this toolbox are
%angsmat2xyzmat.m        lindiscrim1dvals.m      plot2dlinbnd.m
%angvec2xymat.m          lindiscrim2dvals.m      plot2dquadbnd.m
%bivnormpdf.m            lindiscrim3dvals.m      plot2dresp.m
%bnd2slint.m             linprobcorr.m           plot2drespcustom.m
%clip_vals.m             negloglike_1dGLC.m      plot2dstim.m
%colorstr2colorval.m     negloglike_2dGLC.m      plot2dstimcustom.m
%colorval2colorstr.m     negloglike_2dGQC.m      plot3dlinbnd.m
%dprime.m                negloglike_3dGLC.m      plot3dquadbnd.m
%dprimef.m               negloglike_3dGQC.m      plot3dresp.m
%fastnormalcdf.m         new2old_2dparams.m      plot3drespcustom.m
%fisherdiscrim1d.m       new2old_3dparams.m      plot3dstim.m
%fisherdiscrim2d.m       norm_old_1dparams.m     plotbivnormcon.m
%fisherdiscrim3d.m       norm_old_2dparams.m     quadbndpercorr.m
%fit_1dGLC.m             norm_old_3dparams.m     quaddecisbnd.m
%fit_1dnoise.m           normalcdf.m             randrows.m
%fit_2dGLC.m             normalcdfinv.m          reducemat.m
%fit_2dGQC.m             normalpdf.m             reducevec.m
%fit_3dGLC.m             old2new_2dparams.m      sim1dlin.m
%fit_3dGQC.m             old2new_3dparams.m      sim2dlin.m
%gensample.m             percorr.m               sim3dlin.m
%getbndendpts.m          plot1dlinbnd.m          sphdecisbnd.m
%getcovs.m               plot1dresp.m            transample.m
%getmeans.m              plot1dstim.m            xymat2angvec.m
%lindecisbnd.m           plot2bivnorms.m         xyzmat2angmat.m

% The demonstrations currently provided in this toolbox are
%Comp2dStats.m           Find3dGQC.m             Plot2dGRTdata.m
%Design2dGRTexp.m        Gen1dGRTstim.m          Plot3dGRTdata.m
%Find1dGLC.m             Gen2dGRTstim.m          Sim1dSubj.m
%Find2dGLC.m             Gen3dGRTstim.m          Sim2dSubj.m
%Find2dGQC.m             Plot1dGRTdata.m         Sim3dSubj.m
%Find3dGLC.m             Plot2BVNs.m
echo off;

