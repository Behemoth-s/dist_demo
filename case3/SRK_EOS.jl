using PolynomialRoots
using CSV
prop = CSV.read("prop.csv")
function fugacityVapour(T, y...)
    u = 1;
    v = 0;
    omeg = 0.08664;
    psi = 0.42748;
    mi = [0.480 1.574 0.176];
    m = [w.^0  w  -w.^2] * mi';
    alfa = (1 + m .* (1 - (T ./ Tc).^0.5)).^2;
    ai = psi * (R^2) * (Tc.^2) ./ Pc .* alfa;
    bi = omeg * R * Tc ./ Pc;
    Q = ((ai * ai').^0.5) .* (1 - bipk);
    a = y' * Q * y;
    b = y' * bi;
    A = a * p / (R * T)^2;
    B = b * p / (R * T);
    c1 = 1;
    c2 = -(1 + B - u * B);
    c3 = A + v * B^2 - u * B - u * B^2;
    c4 = -(A * B + v * B^2 + v * B^3);
    r = roots([c1,c2,c3,c4])
    z =
    temp = (1 - bipk)' * (z .* ai.^0.5)
end
