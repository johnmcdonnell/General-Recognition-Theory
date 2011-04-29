function [ this ] = load_file(fn)

unix( [ 'data/summ_subj.awk ', fn, ' >! this.dat' ] );
% Load raw data
load this.dat;

end
