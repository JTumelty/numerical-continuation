function kp = init_idrs(nres,kp,ic,fp)
% INIT_IDRS initialises the structures needed for the IDR(s) method if it
% is chosen. 
%
% kp.iLHS.name and kp.iP.name are required field to functions to evaluate 
% LHS and the preconditioner respectively. 
% Other fields in the structure are the parameters needed for the 
% evaluation of these functions. 

    if kp.kmeth == 2
        kp.iLHS     = struct;
        kp.iP       = struct;
        kp.iLHS.name= 'idrs_lhs';
        kp.iP.name  = 'idrs_precon';
        
        kp.iLHS.un  = nres.un;
        kp.iLHS.rn  = nres.rn;
        kp.iLHS.k   = fp.k;
        kp.iLHS.nu  = ic.nu;
        kp.iLHS.N   = fp.N;
        kp.iLHS.w   = fp.w;
        
        kp.iP.rn    = nres.rn;
        kp.iP.k     = fp.k;
        kp.iP.N     = fp.N;   
        
        %  If using dt preconditioner specify dt. 
        if kp.precon == 1
            kp.iP.dt= kp.dt;
        end     
    end
end
