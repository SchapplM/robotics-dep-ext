%ZOrderGet   Get the current z-order position of an existing object of a figure
%
%   [hIndex, hCount] = ZOrderGet(h = gco)
%
%   Set a new z-order position of an existing object of a figure.
%   The lower the z-order, the closer the object to the user eye.
%
%   When adding a new object to a figure, by default it is provided with
%   the first z-order value, so it becomes the most visible one.
%
%   Example:
%      figure(1); clf; hold on;
%      h1 = plot(1, 1, 'o', 'MarkerSize', 10, 'MarkerFaceColor', [1 0 0]);
%      h2 = plot(1, 1, 'o', 'MarkerSize', 20, 'MarkerFaceColor', [0 1 0]);
%      % the red marker (h1) is not visible at this time
%      [i cnt] = ZOrderGet(h1)
%      i =
%           2
%      cnt =
%           2
%
%   Thanks and inspired to AddReorderButtons by Geoffrey K. Adams, November 2007
%
%   Copyright 2011
%
%   v1.0.0 - 28/04/2011
%   Marcello Ferro <marcello.ferro@ilc.cnr.it>
%   http://www.ilc.cnr.it
%
function [hIndex, hCount] = ZOrderGet(h)

% Check params
if(nargin < 1)
    h = gco;
end

% Get the parent of the object
parent = get(h, 'Parent');
if(parent == 0)
    % This means the object is probably the figure.. So just ignore
    hIndex = 0;
    hCount = 0;
    return;
end

% Get all children of the object's parent
children = get(parent, 'Children');

% Find out this object's index number
hIndex = find(h == children);  
hCount = length(children);