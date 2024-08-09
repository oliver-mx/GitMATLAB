function [c, ceq] = nonlcon(input,option,option_data,option_mesh,option_BVP, minConstraint)
%%  nonlcon  Nonlinear condition in fmincon 
%
%       nonlcon(input,option,option_data,option_mesh,option_BVP, minFW)
%       will be used in fmincon, if a lower bound for FW is needed,
%       when minimizing SEC. Use the same option_data, option_mesh 
%       and option_BVP as used in fmincon. The desired value of freshwater 
%       production rate is minFW.
%
%   Input:
%       input         -   Input for the Data function
%       option        -   'default' or 'FW'
%       option_data   -   Selects Data function used: data(input)
%       option_mesh   -   NMax of bvp5c
%       option_BVP    -   RelTol of bvp5c
%       minConstraint -   lower bound for FW/SEC_net
%
%   Output:
%       c(x) ≤0
%       ceq(x) =0
%
%   Example:
%       @(x)nonlcon(x,'default')
%       or
%       @(x)nonlcon(x,'FW',option_data,option_mesh,option_BVP, minFW)
%

%% nonlcon:
switch (option)
    case 'default'
    c=[];
    case 'FW'
    % output from fun: -FW    
    c= fun_1(input,option_data,'FW',option_mesh,option_BVP) + minConstraint; 
    case 'SEC'
    % output from fun: -SEC_net    
    c= fun_1(input,option_data,'SEC',option_mesh,option_BVP) + minConstraint; 
end
ceq=[];
end