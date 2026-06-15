function [output_Atilde, output_Omega, diff_history] = IERP_forconvergencey(D)
    N = size(D,1);
    Input_Omega = D;

    W_tilde  = Input_Omega;
    W_tilde(W_tilde ~= 0) = 1 ./ W_tilde(W_tilde ~= 0);

    Omega_new = EffectiveResistance(W_tilde);  % Compute the effective resistance Omega

    diff_change = sum(sum((abs(D - Omega_new))./(D+(D==0))))/(N*(N-1));

    % Record the initial error
    diff_history = diff_change;

    flag = 1;
    A = (W_tilde > 0);

    % Initialize row and col in case the while loop stops unexpectedly
    row = [];
    col = [];

    while(flag == 1)

        previous_change = diff_change;

        % Compute R
        R = A .* ((Omega_new + eye(N)).^-1 - W_tilde) .* (D - Omega_new);

        % Identify the maximum element
        [val, ~] = max(max(R));

        % Identify the link
        [row, col] = find(R == val);

        % Remove the selected link
        A(row(1), col(1)) = 0;
        A(col(1), row(1)) = 0;

        W_tilde = A .* D;
        W_tilde(W_tilde ~= 0) = 1 ./ W_tilde(W_tilde ~= 0);

        % Update effective resistance
        Omega_new = EffectiveResistance(W_tilde);

        % Compute new error
        diff_change = sum(sum((abs(D - Omega_new))./(D+(D==0))))/(N*(N-1));

        % Record the error after this attempted pruning step
        diff_history(end+1) = diff_change;

        % Stop if the error increases
        if abs(diff_change) > abs(previous_change)
            flag = 0;
        end

    end

    % If the last removed link caused error increase, return it
    if ~isempty(row) && ~isempty(col) && row(1) ~= col(1)

        A(row(1), col(1)) = 1;
        A(col(1), row(1)) = 1;

        W = A .* D;

    else

        W = A .* D;

    end

    W(W ~= 0) = 1 ./ W(W ~= 0);

    output_Atilde = W;
    output_Omega = EffectiveResistance(output_Atilde);

    alpha = alpha_l1_global_para(output_Omega, D);

    output_Omega = alpha * output_Omega;
    output_Atilde = 1/alpha * output_Atilde;

end
