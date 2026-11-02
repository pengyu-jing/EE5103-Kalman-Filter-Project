% RUN_ALL  EE5103 Computer Control Systems - Kalman filter project (MATLAB version).
%
%   Runs Parts 1-3, prints all numbers of the report to the command window and saves
%   the figures in ./figures and the Part 3 tables in ./results.
%   No toolbox is required (base MATLAB only).
%
%   Parts 1 and 2 use random noise. MATLAB's random-number generator differs from the
%   Python one used for the report, so their statistics differ slightly from the
%   report (the analytic values are identical). Part 3 uses the measured data of the
%   project, so its results are identical to the report (see the self-check).
clear;  close all;  clc;
here   = fileparts(mfilename('fullpath'));
figdir = fullfile(here, 'figures');
addpath(here);

res1 = part1(figdir);
res2 = part2(figdir);
res3 = part3(figdir);
fprintf('\nFigures saved in %s\n', figdir);
