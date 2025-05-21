clear,clc

A_input = [0	0	0	0	0	0	0	0	3	0
    0	0	0	0	0	0	6	10	0	0
    0	0	0	0	2	0	0	2	0	0
    0	0	0	0	0	0	0	4	0	0
    0	0	2	0	0	10	0	0	4	0
    0	0	0	0	10	0	0	0	0	0
    0	6	0	0	0	0	0	0	0	0
    0	10	2	4	0	0	0	0	0	0
    3	0	0	0	4	0	0	0	0	3
    0	0	0	0	0	0	0	0	3	0];


T = graph(A_input);
subplot(1,2,1)
h1 = plot(T,'EdgeLabel',T.Edges.Weight,'NodeColor',[0.8500 0.3250 0.0980], ...
'EdgeAlpha',0.5,'LineWidth',1,'MarkerSize',7,'EdgeLabelColor',[0 0.4470 0.7410],'NodeFontSize',10);
x = h1.XData;
y = h1.YData;
Omega = EffectiveResistance(A_input)
% Generate a demand matrix
D = EffectiveResitance_withinverseA(A_input)

D2 = distances(T)


[output_Atilde,output_Omega] = IERP2(Omega)

G_out = graph(output_Atilde,'upper');
subplot(1,2,2)
h2 = plot(G_out,'EdgeLabel',G_out.Edges.Weight,'XData', x, 'YData', y,'NodeColor',[0.8500 0.3250 0.0980], ...
'EdgeAlpha',0.5,'LineWidth',1,'MarkerSize',7,'EdgeLabelColor',[0 0.4470 0.7410],'NodeFontSize',10);


sum(sum(abs(output_Omega - Omega)))


function Omega = EffectiveResitance_withinverseA(A)
% the resitance is the inverse of the link weight
    Aforomega = A;
    Aforomega(Aforomega ~= 0) = 1 ./ Aforomega(Aforomega ~= 0);
    Omega = EffectiveResistance(Aforomega);
end



function [output_Atilde,output_Omega] = IERP2(D)
    N = size(D,1);
    Input_Omega = D;
    W_tilde  = Input_Omega;
    W_tilde(W_tilde ~= 0) = 1 ./ W_tilde(W_tilde ~= 0);
    Omega_new = EffectiveResistance(W_tilde);            % Compute the effective resistance Omega
    OmegaDiff = round(D-Omega_new,10);
    diff_change = sum(sum(OmegaDiff));
    val = 1;
    flag = 1;
    A = (W_tilde > 0);
    Gnow = graph(W_tilde,"upper");
    while(flag==1 && val > 0 && all(conncomp(Gnow) == 1))                     % Remove links one by one until we exceed the constraints                            
        previous_change = diff_change;
        % method 1 R = A.*((Omega_new+eye(N)).^-1 - W_tilde)*(D-Omega_new) 
        R = A.*((Omega_new+eye(N)).^-1 - W_tilde).*(D-Omega_new);     % Compute R
        [val,~] = max(max(R));                              % Identify the maximum element
        [row,col] = find(R == val);                        % Identify the link
        A(row(1),col(1)) = 0; A(col(1),row(1)) = 0;        % Remove the link
        W_tilde = A.*D;
        Gnow = graph(W_tilde,"upper");
        W_tilde(W_tilde ~= 0) = 1 ./ W_tilde(W_tilde ~= 0);      % Compute W tilde
        Omega_new = EffectiveResistance(W_tilde);               % Update the shortest path weight matrix
        
        OmegaDiff = round(D-Omega_new,10)                          
        diff_change = sum(sum(OmegaDiff))
        if abs(diff_change) > abs(previous_change)
            flag=0;
        end    
    end
    
    if(row(1) ~= col(1))
        A(row(1),col(1)) = 1; A(col(1),row(1)) = 1;         % Return lastly removed link
        W = A.*D;
    %     W = 2*1/sum(sum(A.*D)).*A.*D;                       % Update the weighted adjacency matrix
    end
    output_Atilde = W;
    output_Omega = EffectiveResitance_withinverseA(W);
end

