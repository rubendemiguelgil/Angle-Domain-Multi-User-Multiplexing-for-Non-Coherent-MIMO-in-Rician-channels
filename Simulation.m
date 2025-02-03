clear, clc, close all;

%% Parameters
N_users = 2; 
M = 100; % Number of Rx antennas (BS)
L = 100; % Tx length (in bits)

%% Create bit strings
bits = round(rand(L, N_users));

%% Constellation modulation