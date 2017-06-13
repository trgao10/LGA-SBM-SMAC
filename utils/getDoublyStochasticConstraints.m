function [A,b] = getDoublyStochasticConstraints(d)

T = getTensorTranspose(d,d);
C1 = kron(ones(d,1)',eye(d,d));
C2 = kron(ones(d,1)',eye(d,d)) * T;
A = [C1;C2];
b = ones(2*d,1);

end