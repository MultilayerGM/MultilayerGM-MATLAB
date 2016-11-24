% Example use of DirichletDCSBMBenchmark for multiplex and temporal multilayer networks
%
% References:
% 
%       [1] Generative benchmark models for mesoscale structure in multilayer 
%       networks, M. Bazzi, L. G. S. Jeub, A. Arenas, S. D. Howison, M. A. 
%       Porter. arXiv1:608.06196.
%
% Citation: 
%
%       If you use this code, please cite as
%       Lucas G. S. Jeub and Marya Bazzi
%       "A generative model for mesoscale structure in multilayer networks 
%       implemented in MATLAB," [insert website] (2016).



% EXAMPLE 1: Multiplex
% Multiplex example with uniform interlayer dependencies between all pairs
% of layers (interlayer depencency tensor as in FIG 3b of [1] and model parameter
% choices as in FIG 5 and FIG 8 of [1] for specific values of p and mu) 

n_layers=15;
n_nodes=1000;
p = 0.95; mu = 0.1;

L = MultiplexDependencyMatrix(n_layers,p);
[A,S]=DirichletDCSBMBenchmark(n_nodes,n_layers,L,...
'UpdateSteps',200,'theta',1,'communities',10,'q',1,...
'exponent',-2,'kmin',3,'kmax',150,'mu',mu,'maxreject',100);

% EXAMPLE 2: TEMPORAL 
% Temporal example with uniform interlayer dependencies between successive
% layers (interlayer depencency tensor as in FIG 3a of [1] and model parameter
% choices as in FIG 4 and FIG 9 of [1] for specific values of p and mu)

n_layers = 100;
n_nodes = 150;
p = 0.95; mu = 0.1;

L = TemporalDependencyMatrix(n_layers,p);
[A,S]=DirichletDCSBMBenchmark(n_nodes,n_layers,L,...
'UpdateSteps',1,'theta',1,'communities',5,'q',1,...
'exponent',-2,'kmin',3,'kmax',30,'mu',mu,'maxreject',100);


% EXAMPLE 3: TEMPORAL WITH CHANGE POINTS
% Temporal example with nonuniform interlayer dependencies between successive
% layers (interlayer depencency tensor as in FIG 3a of [1] and model parameter
% choices as in FIG 10 of [1] for specific values of p and mu) 

n_layers = 100;
n_nodes = 150;

p = 0.95; p_change = 0; 
p_CP = pvec(pit)*ones(n_layers-1,1);
p_CP(25) = p_change; p_CP(50) = p_change; p_CP(75) = p_change;

L = diag(p_CP,1);

[A,S]=DirichletDCSBMBenchmark(n_nodes,n_layers,L,...
'UpdateSteps',1,'theta',1,'communities',5,'q',1,...
'exponent',-2,'kmin',3,'kmax',30,'mu',mu,'maxreject',100);
        
        
