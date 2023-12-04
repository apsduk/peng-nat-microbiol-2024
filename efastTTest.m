function SAefast_struct = efastttest(Si, rangeSi, Sti, rangeSti, nvar, alpha) % use rangeSi or ramgeSti from efast_sd
%Si=Si(:,time_points);
%Sti=Sti(:,time_points);
%rangeSi=rangeSi(:,time_points,:);
%rangeSti=rangeSti(:,time_points,:);
%Si_struct
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%{
function pairwise_ttest

input: matrix of sensitivity indices Si(# parameters,# search curves)
output: square matrix of p-values, comparing each parameter
%}
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
SAefast_struct=struct;

[k, NR, output] = size(rangeSi); % record k number of parameters, t number of tpoints and NR number of resampled search curves

SAefast_struct.Si=Si;
SAefast_struct.rangeSi=rangeSi;
SAefast_struct.Sti=Sti;
SAefast_struct.rangeSti=rangeSti;
for u = nvar
    %% Compare Si or STi of ith parameter with the dummy uncontected parameter
    for i = 1:k
        %% Si
        [a,p_Si(i,u)] = ttest2(squeeze(rangeSi(i,:,u)), squeeze(rangeSi(k,:,u)), alpha, 'right', 'unequal');
        avg_Si(i,u) = mean(squeeze(rangeSi(i,:,u)));
        std_Si(i,u) = std(squeeze(rangeSi(i,:,u)));
        %% Sti
        [a,p_Sti(i,u)] = ttest2(squeeze(rangeSti(i,:,u)), squeeze(rangeSti(k,:,u)), alpha, 'right', 'unequal');
        avg_Sti(i,u) = mean(squeeze(rangeSti(i,:,u)));
        std_Sti(i,u) = std(squeeze(rangeSti(i,:,u)));
    end % for i
    SAefast_struct.p_Si(:,:,u) = p_Si(:,u);
    SAefast_struct.p_Sti(:,:,u) = p_Sti(:,u);
    SAefast_struct.avg_Si(:,:,u) = avg_Si(:,u);
    SAefast_struct.std_Si(:,:,u) = std_Si(:,u);
    SAefast_struct.avg_Sti(:,:,u) = avg_Sti(:,u);
    SAefast_struct.std_Sti(:,:,u) = std_Sti(:,u);
end
SAefast_struct.p_Si = squeeze(SAefast_struct.p_Si);
SAefast_struct.p_Sti = squeeze(SAefast_struct.p_Sti);

SAefast_struct.avg_Si = squeeze(SAefast_struct.avg_Si);
SAefast_struct.std_Si = squeeze(SAefast_struct.std_Si);
SAefast_struct.avg_Sti = squeeze(SAefast_struct.avg_Sti);
SAefast_struct.std_Sti = squeeze(SAefast_struct.std_Sti);

end
