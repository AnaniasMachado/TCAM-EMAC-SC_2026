clear;
clc;
t0 = tic;
load('Maragal_2.mat');
A = full(Problem.A);
[m,n] = size(A);
% get linearly independent rows and columns
[~, R_qr, Ecol] = qr(A, 'vector');
tol = 1e-8;
r = sum(abs(diag(R_qr)) > tol);
[~, ~, Erow] = qr(A', 'vector');
C = sort(Ecol(1:r));
R = sort(Erow(1:r));
assert(rank(A(R,C), tol) == r, 'A(R,C) is rank-deficient');
% run local search
[norm,time,swaps,C_out] = LSLAFI_Det_P3 (A,r,n,R,C);
time_tot = toc(t0);
time_tot