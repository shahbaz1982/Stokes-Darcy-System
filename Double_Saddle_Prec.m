% This code is based on the paper entitled "A preconditioning method utilizing an augmented
%Lagrangian approach for the coupled Stokes-Darcy system. by Shahbaz Ahmad
%and Wasif Khan

clc;
clear all;
close all;

load('matrix1.mat'); % Load the coefficient matrix from the file 'matrix1.mat'
load('rhs1.mat'); % Load the RHS load vector from the file 'rhs1.mat'

A11=A(1:153,1:153);
A12T=A(1:153,154:459);
O13=A(1:153,460:504);

A21=A(154:459,1:153);
A22=A(154:459,154:459);
B23T=A(154:459,460:504);

O31=A(460:504,1:153);
B32=A(460:504,154:459);
O33=A(460:504,460:504);

% Generate a random matrix
n = 45;  % Size of the matrix
Q = rand(n);
% Make the matrix symmetric
Q = 0.5*(Q + Q');
% Make the matrix positive definite
Q = Q + n*eye(n); % Adding a multiple of identity matrix to make eigenvalues positive definite

gamma=100;%gamma
AAA=[A11 A12T O13;
    A21 A22+gamma*B23T*inv(Q)*B32 B23T;
    O31 B32 O33];

 %---------------------GMRES Solver---------------------------------------
tol = 1e-8; % tolerance for FGMRES
restart = 10; % restart parameter for FGMRES
[Ar Ac]=size(AAA);
% Initialize FGMRES variables
x0 = zeros(Ar, 1);
maxit = 5; % maximum number of iterations for FGMRES

[x10, flag0, relres0, iter0, resvec0] = gmres(AAA, b, restart, tol, maxit, [], [], x0);

 semilogy(resvec0,'-*')

%------------------Zero Matrices-------------------------------------------
O21=zeros(306,153);
O32=zeros(45,306);
O12T=zeros(153,306);
%-------------------PconD---------------------------------------------------
PCOD=2*[A11 O12T O13;
    O21 A22 B23T;
    O31 B32 O33];
[x1, flag1, relres1, iter1, resvec1] = gmres(AAA, b, restart, tol, maxit, PCOD, [], x0);
hold on
semilogy(resvec1,'-*')
 hold off
 
%-------------------PconT---------------------------------------------------
PCOT=[A11 O12T O13;
    O21 A22 B23T;
    O31 B32 O33];
[x2, flag2, relres2, iter2, resvec2] = gmres(AAA, b, restart, tol, maxit, PCOT, [], x0);
hold on
semilogy(resvec2,'-*')
 hold off
%-------------------PB1---------------------------------------------------
PB1=[A11 A12T O13;
    O21 A22+gamma*B23T*inv(Q)*B32 B23T;
    O31 O32 (-1/gamma)*Q];
[x3, flag3, relres3, iter3, resvec3] = gmres(AAA, b, restart, tol, maxit, PB1, [], x0);
hold on
semilogy(resvec3,'-*')
 hold off
%-------------------PB2---------------------------------------------------
alpha=2*gamma;
PB2=[A11 A12T O13;
    O21 A22+gamma*B23T*inv(Q)*B32 (1-gamma*(1/alpha))*B23T;
    O31 B32 (-1/alpha)*Q];

[x4, flag4, relres4, iter4, resvec4] = gmres(AAA, b, restart, tol, maxit, PB2, [], x0);
hold on
semilogy(resvec4,'-o')
 hold off
%--------------------------------P1----------------------------------
P1=[A11 O12T O13;
    A21 A22+gamma*B23T*inv(Q)*B32-A21*inv(A11)*A12T B23T;
    O31 B32 O33];
[x5, flag5, relres5, iter5, resvec5] = gmres(AAA, b, restart, tol, maxit, P1, [], x0);
hold on
semilogy(resvec5,'-o')
 hold off
%--------------------------------P2----------------------------------
P2=[A11 2*A12T O13;
    A21 A22+gamma*B23T*inv(Q)*B32+A21*inv(A11)*A12T B23T;
    O31 B32 O33];
[x6, flag6, relres6, iter6, resvec6] = gmres(AAA, b, restart, tol, maxit, P2, [], x0);
hold on
semilogy(resvec6,'-o')
 
legend('GMRES','P_{conD}','P_{conT}','P_{B1}','P_{B2}','P_1','P_2')
title('Relative Residual Norms')
hold off

%-----------------------------Eigenvalues-------------------------------------
figure; % Create a new figure
eigenvalues_AAA = eig(full(AAA)); % Compute eigenvalues of the loaded matrix
plot(real(eigenvalues_AAA), imag(eigenvalues_AAA), 'o'); % Plot eigenvalues
title('Eigenvalues of Matrix A');
xlabel('Real Part');
ylabel('Imaginary Part');

figure; % Create a new figure
eigenvalues_PCOD = eig(full(inv(PCOD)*AAA)); % Compute eigenvalues of the loaded matrix
plot(real(eigenvalues_PCOD), imag(eigenvalues_PCOD), 'o'); % Plot eigenvalues
title('Eigenvalues of Matrix P_{conD}^{-1}A');
xlabel('Real Part');
ylabel('Imaginary Part');
ylim([-1 1])

figure; % Create a new figure
eigenvalues_PCOT = eig(full(inv(PCOT)*AAA)); % Compute eigenvalues of the loaded matrix
plot(real(eigenvalues_PCOT), imag(eigenvalues_PCOT), 'o'); % Plot eigenvalues
title('Eigenvalues of Matrix P_{conT}^{-1}A');
xlabel('Real Part');
ylabel('Imaginary Part');
ylim([-1 1])

figure; % Create a new figure
eigenvalues_PB1 = eig(full(inv(PB1)*AAA)); % Compute eigenvalues of the loaded matrix
plot(real(eigenvalues_PB1), imag(eigenvalues_PB1), 'o'); % Plot eigenvalues
title('Eigenvalues of Matrix P_{B1}^{-1}A');
xlabel('Real Part');
ylabel('Imaginary Part');
ylim([-1 1])

figure; % Create a new figure
eigenvalues_PB2 = eig(full(inv(PB2)*AAA)); % Compute eigenvalues of the loaded matrix
plot(real(eigenvalues_PB2), imag(eigenvalues_PB2), 'o'); % Plot eigenvalues
title('Eigenvalues of Matrix P_{B2}^{-1}A');
xlabel('Real Part');
ylabel('Imaginary Part');
ylim([-1 1])

figure; % Create a new figure
eigenvalues_P1 = eig(full(inv(P1)*AAA)); % Compute eigenvalues of the loaded matrix
plot(real(eigenvalues_P1), imag(eigenvalues_P1), 'o'); % Plot eigenvalues
title('Eigenvalues of Matrix P_1^{-1}A');
xlabel('Real Part');
ylabel('Imaginary Part');
ylim([-1 1])

figure; % Create a new figure
eigenvalues_P2 = eig(full(inv(P2)*AAA)); % Compute eigenvalues of the loaded matrix
plot(real(eigenvalues_P2), imag(eigenvalues_P2), 'o'); % Plot eigenvalues
title('Eigenvalues of Matrix P_2^{-1}A');
xlabel('Real Part');
ylabel('Imaginary Part');
ylim([-1 1])
%---------------------------------------------------------------------------