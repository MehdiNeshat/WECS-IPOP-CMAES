%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%     Wave Energy Converters (WEC) position optimization by IPOP-CMA-ES
%%% PhD student : Mehdi Neshat, Optimization and Logistics group, University of Adelaide
%%% neshat.mehdi@gmail.com
%%% mehdi.neshat@adelaide.edu.au
%%% Supervisors : Dr.Alexander and Dr.Wagner
%%%               School of Computer Sceince , University of Adelaide
%%%               https://cs.adelaide.edu.au/~optlog/research/energy.php
%%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%
%%% As the code applies the parallel evaluation of the population model,
%%% the Parallel Computing Toolbox™ is necessary.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Acknowledgement:
%%% The fitness function is programmed and modified by Dr. Nataliia Sergiienko (01/06/2018)
%%% https://www.adelaide.edu.au/directory/nataliia.sergiienko
%%% -----------------------
%%% We would like to express our deep gratitude to Dr.Hansen and Dr.Auger for
%%% publishing the source code of IPOP-CMA-ES.
%%% Auger, A., & Hansen, N. (2005, September). A restart CMA evolution strategy
%%% with increasing population size. In 2005 IEEE congress on evolutionary
%%% computation (Vol. 2, pp. 1769-1776). IEEE.
%%% http://www.cmap.polytechnique.fr/~nikolaus.hansen/cmaes_inmatlab.html#matlab
%%% -------------------------
%%% And also special thanks to
%%%  John D'Errico (2020). fminsearchbnd, fminsearchcon
%%% (https://www.mathworks.com/matlabcentral/fileexchange/8277-fminsearchbnd-fminsearchcon),
%%% MATLAB Central File Exchange. Retrieved April 11, 2020.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Fitness function calculates the total power output of Submerged spheres,
%%% 3-tether arrangement of a wave farm model
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% The results were published in
%%% Neshat, M., Alexander, B., Sergiienko, N., & Wagner, M. (2019).
%%% A new insight into the Position Optimization of Wave Energy Converters
%%% by a Hybrid Local Search. arXiv preprint arXiv:1904.09599.
%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Improvements :
%%% 15/03/2020      Implementation starts
%%%
%%%
%%-------------------------------------------------------------------------

function []=Main()
clc
clear all
close all
warning off
rng('shuffle')% creating random number generator based on the current time

global siteOpts
global buoy
id                = 1 ;       % number of run
Opt.Buoy_Num      = 16;       % WEC number
Opt.WaveModel     = 1;        % 1= Perth, 2=Adelaide, 3=Tasmania, 4=Sydeny
Opt.Safe_Distance = 50;        % the minimum distance between buoys
Opt.PopSize       = 20;
Opt.Area          = round( sqrt(Opt.Buoy_Num *20000));
Opt.LBounds       = 0;       % lower bounds, scalar or Nx1-vector';
Opt.UBounds       = Opt.Area;% upper bounds, scalar or Nx1-vector';
Opt.nVar          = Opt.Buoy_Num*2;
%% IPOP_CMAES options
Opt.Restarts         = 50; % up to computational budget
Opt.TolFun           = 10^3;
Opt.TolHistFun       = 10^3; % stop if back fun-changes smaller TolHistFun';
Opt.TolX             = 0;
Opt.StopOnStagnation = 'on';
Opt.DispModulo       = 5;
Opt.IncPopSize       = 1.5;    % multiplier for population size before each restart';
Opt.ParentNumber     = floor(Opt.PopSize /2);

fitfun               = 'Eval_Power';
xstart               = initialize_array(Opt); %(first array )

insigma              = 0.3* (Opt.UBounds-Opt.LBounds);
%%
[siteOpts,buoy,Opt.siteName]=Compute_Hydrodynamic(Opt);

%% Run the position optimization method
[xmin,fmin, counteval, stopflag, out,bestever,generation] = IPOP_CMAES( fitfun,xstart, insigma, Opt,id)  ;
disp (['xmin=',mat2str(xmin)])
disp (['fmin=',num2str(fmin)])
disp (['bestever.x=',mat2str([bestever.x])])
disp (['bestever.f=',num2str([bestever.f])])
disp (['bestever.eval=',mat2str([bestever.evals])])

save(['IPOP_CMAES_',Opt.siteName,'_NB_',num2str(Opt.Buoy_Num),'_id_',num2str(id),'.mat'],'generation')
end