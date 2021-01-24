% Function for checking the domination between the multi-obj. PSO/GA population
% 
% Input:
% fitness [Np x No] double
%   fitness values for `Np` different particles/individuals of a PSO or GA
%   optimization with `No` different objectives
% 
% Output:
% dom_vector [Np x 1] logical
%   Vector that indicates if each particle is dominated (1) or not
% 
% Extracted from MOPSO.m of https://de.mathworks.com/matlabcentral/
% fileexchange/62074-multi-objective-particle-swarm-optimization-mopso
% Modifications: Use logical data type

% Copyright (c) 2017, Víctor Martínez; All rights reserved.
% License conditions: See MOPSO_license.txt

function dom_vector = pareto_dominance(fitness)
  Np = size(fitness,1);
  if Np == 1, dom_vector = false; return; end % code below does not work for Np=1
  dom_vector = false(Np,1);
  all_perm = nchoosek(1:Np,2);    % Possible permutations
  all_perm = [all_perm; [all_perm(:,2) all_perm(:,1)]]; % also check other way around
  d = dominates(fitness(all_perm(:,1),:),fitness(all_perm(:,2),:)); % check all particles against each other
  dominated_particles = unique(all_perm(d,2));
  dom_vector(dominated_particles) = true;
end

% Function that returns 1 if x dominates y and 0 otherwise
function d = dominates(x,y)
  d = all(x<=y,2) & any(x<y,2);
end