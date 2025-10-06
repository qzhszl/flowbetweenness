function Omega = EffectiveResitance_withinverseA(A)
% the resitance is the inverse of the link weight
    Aforomega = A;
    Aforomega(Aforomega ~= 0) = 1 ./ Aforomega(Aforomega ~= 0);
    Omega = EffectiveResistance(Aforomega);
end
