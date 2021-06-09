function  nres = find_initial_iterate(nres,fp,cp)
% FIND_INITIAL_ITERATE finds the starting position for each step in the
% continuation process.

    nres.restart = 0;
    nres.its
    
    if nres.its == 1   
        % r and u are given by the initial conditions.
    elseif nres.its == 2
        if cp.bt == 0
            if cp.mp == 1
                nres.nun = nres.nun + nres.ds;
            else
                nres.rn = nres.rn + nres.ds;
            end
        elseif cp.bt == 1
            nres.nun = nres.nun + nres.ds;
        else
            error('Error. Invalid entry for cp.bt. It should be 0 or 1.')
        end
    else             
        tangent();
        nres.un = nres.unt;      
        nres.rn = nres.rnt;
        
        if cp.bt == 1
            nres.vn = nres.vnt;      
            nres.nun = nres.nunt;
        end
        
        if cp.mp == 1
            nres.nun = nres.nunt;
        end
    end   

    nres.itn = 0;
%-------------------------------------------------------------------------%
    function  tangent()
        % Computes the tangent to the curve, using the two previous points,
        % extending the line between them to length ds to give an initial 
        % guess for the next step.   
        nres.rdot    = comp_dot(nres.rprev,nres.sprev);
        nres.udot    = comp_dot(nres.uprev,nres.sprev);
        
        if cp.bt == 0 
            %Compute normalisation coefficient for udot and rdot
            if cp.mp == 1
                nres.nudot = comp_dot(nres.nuprev,nres.sprev);
                nres.lambda= comp_lambda({nres.rdot,fp.w},{nres.udot,1}, {nres.nudot,fp.w});
            else
                nres.lambda= comp_lambda({nres.rdot,fp.w},{nres.udot,1});
            end
        else
            nres.nudot    = comp_dot(nres.nuprev,nres.sprev);
            nres.vdot     = comp_dot(nres.vprev,nres.sprev);
            
            nres.lambda= comp_lambda({nres.rdot,fp.w}, {nres.nudot,fp.w},...
                                     {nres.vdot,1}, {nres.udot,1});
        end
        
        nres.rnt     = nres.rprev(2)   + fp.w^0.5*nres.ds*nres.lambda*nres.rdot;
        nres.unt     = nres.uprev(:,2) + nres.ds*nres.lambda*nres.udot;
        
        if cp.mp == 1
            nres.nunt = nres.nuprev(:,2) + fp.w^0.5*nres.ds*nres.lambda*nres.nudot;
        end
        
        if cp.bt == 1
            nres.nunt  = nres.nuprev(:,2) + fp.w^0.5*nres.ds*nres.lambda*nres.nudot;
            nres.vnt   = nres.vprev(:,2) + nres.ds*nres.lambda*nres.vdot;
        end
    end
end

function adot = comp_dot(aprev,sprev)
% Finds the gradient of a through two previous points
    adot = (aprev(:,2) - aprev(:,1))./(sprev(2)-sprev(1));
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
