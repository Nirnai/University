function Xef = f(in1)
%F
%    XEF = F(IN1)

%    This function was generated by the Symbolic Math Toolbox version 8.0.
%    21-Feb-2018 18:37:11

% ForwarKinematics x = f(q)
q1 = in1(1,:);
q2 = in1(2,:);
q3 = in1(3,:);
q4 = in1(4,:);
Xef = [sin(q1+pi.*(1.0./4.0)).*3.43492424049175e-1-sin(q1).*(q2+6.9e1./2.0e2)-cos(q1).*(q3+4.269741340859453e-1)-sqrt(2.0).*cos(q4).*cos(q1+pi.*(1.0./4.0)).*(2.1e1./2.0e2);cos(q1+pi.*(1.0./4.0)).*(-3.43492424049175e-1)+cos(q1).*(q2+6.9e1./2.0e2)-sin(q1).*(q3+4.269741340859453e-1)-sqrt(2.0).*cos(q4).*sin(q1+pi.*(1.0./4.0)).*(2.1e1./2.0e2);sqrt(2.0).*sin(q4).*(-2.1e1./2.0e2)+8.697413408594536e-2;angle(cos(q4).*cos(q1+pi.*(1.0./4.0))+cos(q4).*sin(q1+pi.*(1.0./4.0)).*1i);angle(sin(q4).*-1i+sqrt(-sin(q4).^2+1.0));pi.*sign(cos(q4)).*(1.0./2.0)];