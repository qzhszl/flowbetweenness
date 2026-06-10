function Omega = EffectiveResistance(Atilde)
m=size(Atilde,1);

DEGREE = sum(Atilde);
Qtilde = diag(DEGREE)-Atilde;

PQ = pinv(Qtilde);

xita = diag(PQ);
u = ones(m,1);

Omega = u*xita.'+xita*u.'-2*PQ;
end
