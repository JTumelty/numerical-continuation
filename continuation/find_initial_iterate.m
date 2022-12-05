function  nres = find_initial_iterate(nres,fp,cp)
% FIND_INITIAL_ITERATE finds the starting position for each step in the
% continuation process.
    
    if nres.its == 1   
        % Initial iterate is given by the initial conditions.
    elseif nres.its == 2
        
        % Previous iterate gives the initial profile, update bifurcation
        % parameter by adding step size
        nres = update_bifurcation_parameter_its2(nres,cp.mp,cp.bt);
        
    else
        % Computes the tangent to the curve from the previous two iterates
        nres = point_on_tangent(nres,fp.w,cp.mp,cp.bt);
        
        % Update (un,rn,...) with (unt,rnt,...), a point found on the
        % tangent to the bifurcation curve, as the initial iterate
        nres = update_nres_with_initial_iterate(nres,cp.mp,cp.bt);

    end   
end
%-------------------------------------------------------------------------%
function nres = update_bifurcation_parameter_its2(nres,mp,bt)
    % UPDATE_BIFURCATION_PARAMETER_ITS2 updates the bifurcation parameter
    % on the second iteration, depending upon the type of continuation
    
    if bt == 0 && mp == 0
        % Normal continuation: Update r by ds
        nres.rn = nres.rn + nres.ds;    
    elseif (bt == 0 && mp == 1) || (bt == 1 && mp == 0)
        % Non-standard continuation. Vary nu by ds
        nres.nun = nres.nun + nres.ds;
    else
        error('Error. Invalid entry for cp.bt or cp.mp. They should be 0 or 1.')
    end
end
%-------------------------------------------------------------------------%
%% Compute tangent
%-------------------------------------------------------------------------%
function  nres = point_on_tangent(nres,w,mp,bt)
    % Computes the tangent to the curve, using the two previous points,
    % extending the line between them to length ds to give an initial 
    % guess for the next step. The methods depends on the type of continuation  
    
    if mp == 0 && bt == 0
        nres = point_on_tangent_standard(nres,w);
    elseif mp == 1 && bt == 0
        nres = point_on_tangent_mp(nres,w);
    elseif mp == 0 && bt == 1
        nres = point_on_tangent_bt(nres,w);
    else
        error('Error. Check cp.mp and cp.bt.')
    end
end

function nres = point_on_tangent_standard(nres,w)
    % Computes the tangent and a point on it for standard numerical continuation
    
    % Finds the gradients of different variables wrt arclength
    [nres.rdot, nres.udot] = find_gradient_all(nres.sprev,nres.rprev, nres.uprev);
    
    % Computes a normalisation coefficient for the gradient
    nres.lambda= comp_lambda({nres.rdot,w},{nres.udot,1});
    
     % Find a point on the tangent a distance ds from previous iterate for
     % different variables
     [nres.rnt, nres.unt] = find_point(nres.lambda, nres.ds, ...
                                        {nres.rprev,nres.rdot, w},...
                                        {nres.uprev,nres.udot,1});
end

function nres = point_on_tangent_mp(nres,w)
    % Computes the tangent and a point on it when continuing the Maxwell Point
    
    % Finds the gradients of different variables wrt arclength
    [nres.rdot, nres.udot, nres.nudot] = ...
        find_gradient_all(nres.sprev,nres.rprev, nres.uprev, nres.nuprev);
    
    % Computes a normalisation coefficient for the gradient
    nres.lambda= comp_lambda({nres.rdot,w},{nres.udot,1}, {nres.nudot,w});
    
    % Find a point on the tangent a distance ds from previous iterate
    [nres.rnt, nres.unt,nres.nunt] = find_point(nres.lambda, nres.ds, ...
                                               {nres.rprev,nres.rdot, w},...
                                               {nres.uprev,nres.udot,1},...
                                               {nres.nuprev,nres.nudot,w});
end

function nres = point_on_tangent_bt(nres,w)
    % Computes the tangent and a point on it when continuing a bifurcation point
    
    % Finds the gradients of different variables wrt arclength
    [nres.rdot, nres.udot, nres.nudot, nres.vdot] = ...
        find_gradient_all(nres.sprev,nres.rprev, nres.uprev, nres.nuprev, nres.vprev);
    
    % Computes a normalisation coefficient for the gradient
    nres.lambda= comp_lambda({nres.rdot,w}, {nres.nudot,w},...
                                     {nres.vdot,1}, {nres.udot,1});
    
    % Find a point on the tangent a distance ds from previous iterate                        
    [nres.rnt, nres.unt,nres.nunt, nres.vnt] = ...
                                find_point(nres.lambda, nres.ds, ...
                                               {nres.rprev,nres.rdot, w},...
                                               {nres.uprev,nres.udot,1},...
                                               {nres.nuprev,nres.nudot,w},...
                                               {nres.vprev,nres.vdot,1});
end
%-------------------------------------------------------------------------%
%% Update nres with starting point
%-------------------------------------------------------------------------%
function nres = update_nres_with_initial_iterate(nres,mp,bt)
     % Update (un,rn,...) with (unt,rnt,...), a point found on the
     % tangent to the bifurcation curve, as the initial iterate
    
    nres.un = nres.unt;      
    nres.rn = nres.rnt;

    if bt == 1 && mp == 0
        nres.vn = nres.vnt;      
        nres.nun = nres.nunt;
    elseif bt == 0 && mp == 1
        nres.nun = nres.nunt;
    end
end
%-------------------------------------------------------------------------%
%% Utility functions
%-------------------------------------------------------------------------%
function adot = find_gradient(aprev,sprev)
   % Finds the gradient of the line through two previous points
    adot = (aprev(:,2) - aprev(:,1))./(sprev(2)-sprev(1));
end

function varargout = find_gradient_all(sprev,varargin)
    % Finds the gradients for a variable number of input arguments
    varargout = cell(1,nargin-1);
    for i = 1:nargin-1
        varargout{i} = find_gradient(varargin{i},sprev);
    end
end

function lambda = comp_lambda(varargin)
%  Find the normalisation constant of the tangent vector obtained from the
%  two previous points. Each pair in varargin adds e.g. w*|udot|^2 to total
%  norm. 
    diff_norm = 0;
    for i = 1:nargin
        diff_norm = diff_norm + varargin{i}{2}*norm(varargin{i}{1})^2;
    end
    lambda = 1./diff_norm^0.5;
end

function varargout = find_point(lambda, ds, varargin)
     % Find a point on the tangent a distance ds from previous iterate for
     % different variables
     varargout = cell(1,nargin-2);
     for i = 1:nargin-2
         varargout{i} = varargin{i}{1}(:,2) + ...
                                ds*lambda*varargin{i}{2}*varargin{i}{3}^0.5;
     end
end
