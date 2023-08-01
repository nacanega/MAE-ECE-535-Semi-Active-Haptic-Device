% Nolan Canegallo
% MAE/ECE 535
% Final Project
% Due 28 July 2023
clear; clc; close all;
%% Setup Optimization Problem
nvars = 10;
      %h0 h1 z2 w1  w2  w3 w4 w5 q1  Ic 
A =   [];
b =   [];
Aeq = [];
beq = [];
%     h0 h1 z2 w1 w2 w3 w4 w5 q1 Ic
lb = [ 0; 0; 0; 0; 0; 0; 0; 0; 0; 0];
ub = [ 1; 1; 1; 1; 1; 1; 1; 1; 1; 1];

% Set nondefault solver options
options = optimoptions("gamultiobj","UseParallel",true,"PopulationSize",200,...
    "MaxGenerations",2000,"PlotFcn",["gaplotrankhist","gaplotpareto",...
    "gaplotparetodistance"]);

% Solve
[solution,objectiveValue] = gamultiobj(@objfunFEMM,nvars,[],[],[],[],lb,ub,...
    [],[],options);

% Clear variables
clearvars options

solution
objectiveValue

[m,n] = size(solution);
ValidResults = struct;
counter = 0;
for i = 1:m
    fprintf('Testing Design %d: ',i)
    [F2,c2,xD2,R2] = objfunFEMM(solution(i,:),false);
    if sum(c2(2:3,:),"all") == 0
        counter = counter + 1;
        saveString = sprintf("result%d",counter); 
        ValidResults.(saveString).objectives = F2;
        ValidResults.(saveString).constraints = c2(1,:);
        ValidResults.(saveString).designString = xD2;
        ValidResults.(saveString).optimizeString = solution(i,:);
        ValidResults.(saveString).extraResistance = R2(1);
        ValidResults.(saveString).Mass = R2(2:3);
        ValidResults.(saveString).Power = R2(4);
        ValidResults.(saveString).Vwaste = R2(5);
        fprintf('Valid\n')
    else
        LC1 = find(c2(2,:));
        LC2 = find(c2(3,:));
        fprintf('Invalid\n');
        fprintf('%d Too Low\n',LC1);
        fprintf('%d Too High\n',LC2);
    end
end

save("results.mat","solution","objectiveValue","ValidResults")