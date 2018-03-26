function [qd] = TrajGen( u )
%TRAJGEN Summary of this function goes here
%   Detailed explanation goes here

global t0 t1 t2 t3 t4 t5

t0 = 5;
t1 = 10;
t2 = 15;
t3 = 20;
t4 = 25;
t5 = 30;

t = u(4);
A = 45;
w = 5;


if(t<t0)
    q1d = u(1);
    q2d = u(2);
    q3d = u(3);
    q1pd = 0;
    q2pd = 0;
    q3pd = 0;
elseif(t>=t0 && t<t1)
    q1d = 0;
    q2d = 0;
    q3d = 0;
    q1pd = 0;
    q2pd = 0;
    q3pd = 0;
elseif(t>=t1 && t<t2)
    q1d = u(1);
    q2d = u(2);
    q3d = u(3);
    q1pd = 0;
    q2pd = 0;
    q3pd = 0;
elseif(t>=t2 && t<t3)
    q1d = A*sin(w*t) + u(1);
    q2d = A*sin(w*t) + u(2);
    q3d = A*sin(w*t) + u(3);
    q1pd = A*w*cos(w*t);
    q2pd = A*w*cos(w*t);
    q3pd = A*w*cos(w*t);
elseif(t>=t3 && t<t4)
    q1d = A*sin(w*t) + u(1);
    q2d = A*sin(w*t) + u(2);
    q3d = A*sin(w*t) + u(3);
    q1pd = A*w*cos(w*t);
    q2pd = A*w*cos(w*t);
    q3pd = A*w*cos(w*t);
else
    q1d = A*sin(w*t) + u(1);
    q2d = A*sin(w*t) + u(2);
    q3d = A*sin(w*t) + u(3);
    q1pd = A*w*cos(w*t);
    q2pd = A*w*cos(w*t);
    q3pd = A*w*cos(w*t);
end


% Interval between PD, PD+G, and PID+G
% Generate Different trajectories at different times


qd=[q1d;q2d;q3d;q1pd;q2pd;q3pd];


end

