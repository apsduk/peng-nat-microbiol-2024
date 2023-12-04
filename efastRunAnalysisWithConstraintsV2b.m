function [youtput, X, N, OMi, NS] = efastRunAnalysisWithConstraintsV2(odefun, lb, ub, pmean, pstd, dist, noutputs, K, NR, NS, MI, constrainedkidx)

% Computation of the frequency for the group of interest OMi and the # of sample points NS (here N=NS)
N=NS*K*NR; % wanted no. of sample points
OMi = floor(((N/NR)-1)/(2*MI)/K);
NS = 2*MI*OMi+1;
if(NS*NR < 65)
    fprintf(['Error: sample size must be >= ' ...
        '65 per factor.\n']);
    return;
end

% Pre-allocation of the output matrix Y
% Y will save only the points of interest specified in
% the vector time_points
Y(NS, noutputs, K, NR) = 0;  % pre-allocation

% Loop over k parameters (input factors)
for k = 1:K % i=# of replications (or blocks)
    % Algorithm for selecting the set of frequencies.
    % OMci(i), i=1:k-1, contains the set of frequencies
    % to be used by the complementary group.
    OMci = efastSetFreq(K,OMi/2/MI,k);
    % Loop over the NR search curves.

    % --- Print progress --------------------------------------------
    disp(['... ',num2str(k),' ',datestr(now(),'HH:MM:SS')]);

    for r = 1:NR
        % Setting the vector of frequencies OM
        % for the k parameters
        cj = 1;
        for j=1:K
            if(j==k)
                % For the parameter (factor) of interest
                OM(k) = OMi;
            else
                % For the complementary group.
                OM(j) = OMci(cj);
                cj = cj+1;
            end
        end
        % Setting the relation between the scalar
        % variable S and the coordinates
        % {X(1),X(2),...X(k)} of each sample point.
        FI = rand(1,K)*2*pi; % random phase shift
        S_VEC = pi*(2*(1:NS)-NS-1)/NS;
        OM_VEC = OM(1:K);
        FI_MAT = FI(ones(NS,1),1:K)';
        ANGLE = OM_VEC'*S_VEC+FI_MAT;
        
        X(:,:,k,r) = 0.5+asin(sin(ANGLE'))/pi; % between 0 and 1
        
        % Transform distributions from standard
        % uniform to general.
        X(:,:,k,r) = efastParameterDistWithConstraints( X(:,:,k,r), ub, lb, pmean, pstd, NS, dist, constrainedkidx); %%this is what assigns 'our' values rather than 0:1 dist
                
        % Do the NS model evaluations.
        for s = 1:NS
            
            % ----- Get the parameter set ---------------------------------
            thisX = X(s, :, k, r);
            
            % [thisX(constrainedkidx), sum(thisX(constrainedkidx))]
            
            % ------ Simulate model ---------------------------------------
            performancevector = odefun(thisX);
            
            % ----- Save output -------------------------------------------
            youtput(s, :, k, r) = performancevector;
            
        end %run_num=1:NS
    end % L=1:NR
end % i=1:k
