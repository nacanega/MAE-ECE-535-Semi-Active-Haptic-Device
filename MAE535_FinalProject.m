% Nolan Canegallo
% MAE/ECE 535
% Final Project
% Due 28 July 2023
clear; clc; close all;
%% Setup Optimization Problem
d1 = 0.254; % mm
fun = @objectiveFunction2;
nvars = 10;
      %h0 h1 z2 w1 w2 w3 w4 w5 q1 Bmr  
A =   [ 0  0  0  0 -1 -1  0  0  1  0];
b = -d;
Aeq = [];
beq = [];
%     h0 h1 z2 w1 w2 w3 w4 w5 q1 Bmr
lb = [0;0;0;0;0;0;0;0;0;0];
ub = [1;1;1;1;1;1;1;1;1;1];
% nonlcon = @nonlinearConstraint2;
% optOpts = optimoptions('gamultiobj','Display','iter','UseParallel',true,'PlotFcn','gaplotpareto');
% [xList,fvalList,exitflag,output,population,scores] = gamultiobj(fun,nvars,A,b,Aeq,beq,lb,ub,nonlcon);
% save("gaResult.mat","xList","fvalList","exitflag","output","population","scores","-mat")