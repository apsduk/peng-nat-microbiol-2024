function [Si,Sti,rangeSi,rangeSti] = efastsd(Y,OMi,MI,var)
[NS s k NR]=size(Y);

if nargin<4
    display(['ERROR = Choose one or more outputs from var 1 and variable ',num2str(c),' of the model'])
    error('eFAST: the output components for the sensitivity is missing. Not enough input arguments.')
end
for u=var
    ttest_sti=zeros(k,k);
    ttest_si=ttest_sti;
    for i=1:k %loop through parameters
        % Initialize AV,AVi,AVci to zero.
        AV = 0;
        AVi = 0;
        AVci = 0;
        for L=1:NR%length(Y(1,u,i,:))
            Y(:,u,i,L) = (Y(:,u,i,L)-mean(Y(:,u,i,L)))';
            % Fourier coeff. at [1:OMi/2].
            N=length(Y(:,u,i,L));
            NQ = (N-1)/2;
            N0 = NQ+1;
            COMPL = 0;
            Y_VECP = Y(N0+(1:NQ),u,i,L)+Y(N0-(1:NQ),u,i,L);
            Y_VECM = Y(N0+(1:NQ),u,i,L)-Y(N0-(1:NQ),u,i,L);
            for j=1:OMi/2
                ANGLE = j*2*(1:NQ)*pi/N;
                C_VEC = cos(ANGLE);
                S_VEC = sin(ANGLE);
                AC(j) = (Y(N0,u,i,L)+Y_VECP'*C_VEC')/N;
                BC(j) = Y_VECM'*S_VEC'/N;
                COMPL = COMPL+AC(j)^2+BC(j)^2;
            end
            % Computation of V_{(ci)}.
            Vci(L) = 2*COMPL;
            %AVci = AVci+Vci;
            % Fourier coeff. at [P*OMi, for P=1:MI].
            COMPL = 0;
            Y_VECP = Y(N0+(1:NQ),u,i,L)+Y(N0-(1:NQ),u,i,L);
            Y_VECM = Y(N0+(1:NQ),u,i,L)-Y(N0-(1:NQ),u,i,L);
            for j=OMi:OMi:OMi*MI
                ANGLE = j*2*(1:NQ)*pi/N;
                C_VEC = cos(ANGLE');
                S_VEC = sin(ANGLE');
                AC(j) = (Y(N0,u,i,L)+Y_VECP'*C_VEC)/N;
                BC(j) = Y_VECM'*S_VEC/N;
                COMPL = COMPL+AC(j)^2+BC(j)^2;
            end
            % Computation of V_i.
            Vi(L) = 2*COMPL;
            %AVi = AVi+Vi;
            % Computation of the total variance
            % in the time domain.
            V(L) = Y(:,u,i,L)'*Y(:,u,i,L)/N;
        end %L
        % Computation of sensitivity indexes.
        %AV = AV/length(Y(1,1,i,:)); %AV=average total variance
        %AVi = AVi/length(Y(1,1,i,:));
        %AVci = AVci/length(Y(1,1,i,:));
        Si(i,u) = mean(Vi)/mean(V);
        Sti(i,u) = 1-mean(Vci)/mean(V);
        rangeSi(i,:,u) = Vi./V;
        rangeSti(i,:,u) = 1-(Vci./V);
    end %i
end