function test_GL_leastsquares
%TEST_GL_LEASTSQUARES Basic regression tests for the public wrapper.

    Atrue = [
        0, 2, 0, 1
        2, 0, 3, 0
        0, 3, 0, 4
        1, 0, 4, 0
    ];
    D = localEffectiveResistance(Atrue);

    [Alearned, Omega] = GL_leastsquares(D);
    assert(norm(Alearned - Alearned.', 'fro') < 1e-12);
    assert(all(Alearned >= 0, 'all'));
    assert(all(diag(Alearned) == 0));
    assert(norm(Omega - localEffectiveResistance(Alearned), 'fro') < 1e-10);
    assert(norm(Alearned - Atrue, 'fro') / norm(Atrue, 'fro') < 1e-8);
    assert(norm(Omega - D, 'fro') / norm(D, 'fro') < 1e-10);

    Dpartial = D;
    Dpartial(1, 3) = 0;
    Dpartial(3, 1) = 0;
    [Apartial, Omegapartial] = GL_leastsquares(Dpartial);
    constrained = Dpartial > 0;
    relativeConstraintError = norm( ...
        Omegapartial(constrained) - Dpartial(constrained)) / ...
        norm(Dpartial(constrained));
    assert(relativeConstraintError < 1e-5);
    assert(norm(Omegapartial - localEffectiveResistance(Apartial), 'fro') < 1e-10);

    didRejectNegative = false;
    try
        GL_leastsquares([0, -1; -1, 0]);
    catch exception
        didRejectNegative = strcmp(exception.identifier, ...
            'GL_leastsquares:NegativeDemand');
    end
    assert(didRejectNegative);

    fprintf('test_GL_leastsquares: all tests passed.\n');
end

function Omega = localEffectiveResistance(A)
    L = diag(sum(A, 2)) - A;
    Lplus = pinv(L);
    diagonal = diag(Lplus);
    Omega = diagonal + diagonal.' - 2 * Lplus;
    Omega = (Omega + Omega.') / 2;
    Omega(1:size(A, 1)+1:end) = 0;
end
