%% Author: Nirnai Rao
%% email: rnirnai@gmail.com

% Symbols
syms q1 q2 q3 qd1 qd2 qd3 q1dd q2dd q3dd m1 m2 m3 l1 l2 l3 k1 k2 k3 z0 T1 T2 b1 b2 b3 g real

%% Positions
xcm1 = 0;
ycm1 = 0;
zcm1 = q1;

xcm2 = l2*sin(q2);
ycm2 = 0;
zcm2 = zcm1 + l1 + l2*cos(q2);

xcm3 = xcm2 + l3*sin(q2+q3);
ycm3 = 0;
zcm3 = zcm2 + l3*cos(q2+q3);




%% Velocities

vxcm1 = simplify(diff(xcm1,q1) * qd1);
vycm1 = 0;
vzcm1 = simplify(diff(zcm1,q1) * qd1);
vcm1 = [vxcm1; vycm1; vzcm1];


vxcm2 = simplify(diff(xcm2,q1) * qd1 + diff(xcm2,q2) * qd2);
vycm2 = 0;
vzcm2 = simplify(diff(zcm2,q1) * qd1 + diff(zcm2,q2) * qd2);
vcm2 = [vxcm2; vycm2; vzcm2];

vxcm3 = simplify(diff(xcm3,q1) * qd1 + diff(xcm3,q2) * qd2 + diff(xcm3,q3) * qd3);
vycm3 = 0;
vzcm3 = simplify(diff(zcm3,q1) * qd1 + diff(zcm3,q2) * qd2 + diff(zcm3,q3) * qd3);
vcm3 = [vxcm3; vycm3; vzcm3];

%% Kinetic Energies

K1 = 1/2 * m1 * (vcm1'*vcm1);
K2 = 1/2 * m2 * (vcm2'*vcm2);
K3 = 1/2 * m3 * (vcm3'*vcm3);

KT = simplify(K1+K2+K3);

%% Potential Energies

P1 = -m1*g*q1;
P2 = -m2*g*(q1+l1+l2*cos(q2));
P3 = -m3*g*(q1+l1+l2*cos(q2)+l3*cos(q2+q3));

% P1 = m1 * g * (l3 + l2 + l1 + z0 - q1);
% P2 = m2 * g * (l3 + l2 - l2*cos(q2));
% P3 = m3 * g * (l3 - l3*cos(q2 + q3));

% G = [0;0;-g];
% tcm1_0 = [xcm1;ycm1;zcm1];
% tcm2_0 = [xcm2;ycm2;zcm2];
% tcm3_0 = [xcm3;ycm3;zcm3];
% 
% P1 = m1*G'*tcm1_0;
% P2 = m2*G'*tcm2_0;
% P3 = m3*G'*tcm3_0;


%% Spring Energies
S1 = 1/2*k1*(q1-z0)^2;
S2 = 1/2*k2*(q2-T1)^2;
S3 = 1/2*k3*(q3-T2)^2;

PT = simplify(P1+P2+P3+S1+S2+S3);

%%
% Joint 1
dK_qd1 = simplify(diff(KT,qd1));
dK_qd1_dt = simplify(...
    diff(dK_qd1,q1) * qd1 + diff(dK_qd1,qd1) * q1dd + ...
    diff(dK_qd1,q2) * qd2 + diff(dK_qd1,qd2) * q2dd + ...
    diff(dK_qd1,q3) * qd3 + diff(dK_qd1,qd3) * q3dd);
dK_q1 = simplify(diff(KT,q1));
dP_q1 = simplify(diff(PT,q1));

% Joint2
dK_qd2 = simplify(diff(KT,qd2));
dK_qd2_dt = simplify(...
    diff(dK_qd2,q1) * qd1 + diff(dK_qd2,qd1) * q1dd + ...
    diff(dK_qd2,q2) * qd2 + diff(dK_qd2,qd2) * q2dd + ...
    diff(dK_qd2,q3) * qd3 + diff(dK_qd2,qd3) * q3dd);
dK_q2 = simplify(diff(KT,q2));
dP_q2 = simplify(diff(PT,q2));

% Joint3
dK_qd3 = simplify(diff(KT,qd3));
dK_qd3_dt = simplify(...
    diff(dK_qd3,q1) * qd1 + diff(dK_qd3,qd1) * q1dd + ...
    diff(dK_qd3,q2) * qd2 + diff(dK_qd3,qd2) * q2dd + ...
    diff(dK_qd3,q3) * qd3 + diff(dK_qd3,qd3) * q3dd);
dK_q3 = simplify(diff(KT,q3));
dP_q3 = simplify(diff(PT,q3));

%% Eq Motion

%Eq 1
aux=expand(dK_qd1_dt-dK_q1+dP_q1);
eq1=simplify(aux);
%Eq 2
aux=expand(dK_qd2_dt-dK_q2+dP_q2);
eq2=simplify(aux);
%Eq 3
aux=expand(dK_qd3_dt-dK_q3+dP_q3);
eq3=simplify(aux);

%% Inertia Matrix

M = sym('M', [3,3]);

eq = [eq1,eq2,eq3];
qdd = [q1dd, q2dd, q3dd];

for i = 1:3
    for j = 1:3
        coef = fliplr(coeffs(collect(eq(i),qdd(j)),qdd(j)));
        if(size(coef,2) > 1)
            M(i,j) = coef(1);
            eq(i) = simplify(eq(i) - coef(1) * qdd(j));
        else
            M(i,j) = 0;
        end
    end
end

%% derivative of M
Md = sym('Md',[3,3]);
q = [q1,q2,q3];
for i = 1:3
    for j = 1:3
        Md(i,j) = simplify(diff(M(i,j),q1) * qd1 + diff(M(i,j),q2)*qd2 + diff(M(i,j),q3)*qd3);
    end
end 




%% Centripital and Coriolis Matrix
qd = [qd1,qd2,qd3];
C = sym(zeros(3));
cross = sym(zeros(3,4));
for i = 1:3
    terms = children(eq(i));
    
    % Expand all terms with factor 2 (possible to implement for all factors larger 1)
    for n = 1:size(terms,2)
        factors = children(terms(n));
        for factor = factors
            if(abs(factor) == 2)
                terms(n) = terms(n)/2;
                terms(end+1) = terms(n);
                
            end 
        end
    end
    
    % Analyze terms constant in all possible C's
    for n = 1:size(terms,2)
        [c,t] = coeffs(terms(n), [qd2,qd3],'All');
        % only qd2 terms
        if(size(t,2) == 1 && size(t,1)>1)
            % add coefficient to scond colmn
           C(i,2) = C(i,2) + (c'*(t/t(end-1)));
           terms(n) = 0;
        % only qd3 terms
        elseif(size(t,1) == 1 && size(t,2)>1)
            C(i,3) = C(i,3) + (c * (t/t(end-1))');
            terms(n) = 0;
        % Eliminate Gravity vector terms
        elseif(size(t,1) == 1 && size(t,2) == 1)
            terms(n) = 0;
        end
    end
    
    cross(i,1:size(nonzeros(terms),1)) = nonzeros(terms)'; 
 
end


    
% Analyze cross term (only for pd2*pd3) -> make up different variants
% of C   
% create a vector representing all combinations of extraction
p = size(nonzeros(cross),1);    % number of cross terms
base = 2;                       % only anayzing 2 states qd2 and qd3
c = zeros(base^p,p);
for n = 0:(base^p)-1
    digits = dec2base(n,base,p);
    for k = 1:p
        c(n+1,k) = sym(digits(k));
    end
end

C_saved = sym(zeros(3,3,16));
N_saved = sym(zeros(3,3,16));

% Test all possible C
for n = 1:size(c,1)
    C_test = C;
    [row,col,v] = find(cross);
    cross_matrix = reshape(v, [max(row),max(col)]);
    c_matrix = reshape(c(n,:), [max(row),max(col)]);
    for r = 1:size(c_matrix,1)
        for j = 1:size(c_matrix,2)
            if(c_matrix(r,j) == 0)
                C_test(r,2) = C_test(r,2) + coeffs(cross_matrix(r,j),qd2);
            elseif(c_matrix(r,j) == 1)
                C_test(r,3) = C_test(r,3) + coeffs(cross_matrix(r,j),qd3);
            end
        end
    end
    
    N = simplify(Md - 2*C_test);
    x = sym('x', [3,1],'real');
    test = simplify(x'*N*x);
    if(test == 0)
        C = C_test;
        % Eliminate Terms in original Equations
        eq = simplify((eq' - C*qd')');
        break;
    end
    fprintf('%i \n', n)
end

%% Iterativ Method for C
% C = sym('C', [3,3],'real');
% qd = [qd1, qd2, qd3];
% temp = sym('temp', [1,3],'real');
% for k = 1:3
%     for j = 1:3
%         for i = 1:3
%             temp(i) = (diff(M(k,j),q(i)) + diff(M(k,i),q(j)) - diff(M(i,j),q(k)))*qd(i);
%         end
%         C(k,j) = 1/2 * sum(temp);
%     end
% end
% 
% N = simplify(Md - 2*C);
% x = sym('x', [3,1],'real');
% test = simplify(expand(x'*N*x));
% eq = simplify((eq' - C*qd')');

%% Gravity Vector
% All Mass/Inertia and Centripital/Coriolis effect were already eliminated
G = eq';

%%
u = [q1 q2 q3 qd1 qd2 qd3 l1 l2 l3 m1 m2 m3 k1 k2 k3 z0 T1 T2 b1 b2 b3 g];
matlabFunction(M,Md,C,G, ...
                   'File','dynamic_model', ...
                   'Vars',{u}, ...
                   'Outputs',{'M','Md','C','g'});

