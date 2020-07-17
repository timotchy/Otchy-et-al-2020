function args2timing_vars(allargs,userargs)
% args2timing_vars(allargs,userargs)
%
% Contains default function inputs, takes all required arguments for a
% given function (allargs), and sets values either to user-input arguments
% (in userargs) or defaults. 

conv_crit = 0.01;
loglkflg = 0;
maxiter = 10000;
rotopt = 1;
psiopt = 1;
initN = 100;
loglkN = 5;
crit_loglk = 0.1;
crit_param = 0.01;
waitbaropt = 1;

for argind = 1:2:length(userargs)
    argnm = userargs{argind};
    eval([argnm '= userargs{argind+1};'])
end

for argind = 1:length(allargs)
    argnm = allargs{argind};
    eval(['assignin(''caller'',argnm,' argnm ');'])
end