clearvars

syms phi_x phi_y phi_z real;
syms pphi_x pphi_y pphi_z real;

phi = [phi_x, phi_y, phi_z];
pphi = [pphi_x, pphi_y, pphi_z];

T_A = sym(eye(6));

T = [cos(phi_y)*cos(phi_z), -sin(phi_z),        0;
    cos(phi_y)*sin(phi_z), cos(phi_z),          0;
    -sin(phi_y),            0,          1];
    
T_A(4:6,4:6) = T;

fid = fopen('My_TA_UR10ccodeProject.txt','wt');
WriteToFile(fid,T_A,'T_A');



T_Ainv = sym(eye(6));

Tinv = [cos(phi_z)/cos(phi_y),               sin(phi_z)/cos(phi_y), 0;
        -sin(phi_z)      ,               cos(phi_z)       , 0;
        cos(phi_z)*sin(phi_y)/cos(phi_y), sin(phi_z)*sin(phi_y)/cos(phi_y), 1];
    
T_Ainv(4:6,4:6) = Tinv;

WriteToFile(fid,T_Ainv,'T_Ainv');

dT_A = sym(zeros(6,6,3));

for i = 1:3
    dT_A(:,:,i)=  simplify(diff(T_A(:,:),phi(i)));
    
end

pT_A = sym(zeros(6,6));

for i = 1:3
    pT_A=  simplify(pT_A + dT_A(:,:,i)*pphi(i));
    
end

WriteToFile(fid,dT_A,'dT_A');

WriteToFile(fid,pT_A,'pT_A');


syms L1 L2 L3 L4 L5 L6 L7 L8 L9 L10 L11 L12 real;
syms q1 q2 q3 q4 q5 q6 real;
syms q1p q2p q3p q4p q5p q6p t real;
syms lambda real

L = [L1 L2 L3 L4 L5 L6 L7 L8 L9 L10 L11 L12];
q = [q1,q2,q3,q4,q5,q6];
qip = [q1p,q2p,q3p,q4p,q5p,q6p];


Jef=[sin(q1)*(sin(q2)*(L5*cos(q3) + cos(q3)*(L11*cos(q4) + L4*sin(q4)*sin(q5) + L12*sin(q4)*sin(q5)) + sin(q3)*(L4*cos(q4)*sin(q5) - L11*sin(q4) + L12*cos(q4)*sin(q5))) + cos(q2)*(sin(q3)*(L11*cos(q4) + L4*sin(q4)*sin(q5) + L12*sin(q4)*sin(q5)) + L5*sin(q3) - cos(q3)*(L4*cos(q4)*sin(q5) - L11*sin(q4) + L12*cos(q4)*sin(q5))) + L3*sin(q2)) + cos(q1)*(L2 + L4*cos(q5) + L12*cos(q5)),-cos(q1)*(cos(q2)*(L5*cos(q3) + cos(q3)*(L11*cos(q4) + L4*sin(q4)*sin(q5) + L12*sin(q4)*sin(q5)) + sin(q3)*(L4*cos(q4)*sin(q5) - L11*sin(q4) + L12*cos(q4)*sin(q5))) - sin(q2)*(sin(q3)*(L11*cos(q4) + L4*sin(q4)*sin(q5) + L12*sin(q4)*sin(q5)) + L5*sin(q3) - cos(q3)*(L4*cos(q4)*sin(q5) - L11*sin(q4) + L12*cos(q4)*sin(q5))) + L3*cos(q2)),-cos(q1)*(cos(q2)*(L5*cos(q3) + cos(q3)*(L11*cos(q4) + L4*sin(q4)*sin(q5) + L12*sin(q4)*sin(q5)) + sin(q3)*(L4*cos(q4)*sin(q5) - L11*sin(q4) + L12*cos(q4)*sin(q5))) - sin(q2)*(sin(q3)*(L11*cos(q4) + L4*sin(q4)*sin(q5) + L12*sin(q4)*sin(q5)) + L5*sin(q3) - cos(q3)*(L4*cos(q4)*sin(q5) - L11*sin(q4) + L12*cos(q4)*sin(q5)))),-cos(q1)*((L4*cos(q2 + q3 + q4 - q5))/2 - (L12*cos(q2 + q3 + q4 + q5))/2 - (L4*cos(q2 + q3 + q4 + q5))/2 + (L12*cos(q2 + q3 + q4 - q5))/2 + L11*cos(q2 + q3 + q4)),-(L4 + L12)*(sin(q1)*sin(q5) - cos(q1)*cos(q2)*cos(q3)*cos(q4)*cos(q5) + cos(q1)*cos(q2)*cos(q5)*sin(q3)*sin(q4) + cos(q1)*cos(q3)*cos(q5)*sin(q2)*sin(q4) + cos(q1)*cos(q4)*cos(q5)*sin(q2)*sin(q3)),0;
sin(q1)*(L2 + L4*cos(q5) + L12*cos(q5)) - cos(q1)*(sin(q2)*(L5*cos(q3) + cos(q3)*(L11*cos(q4) + L4*sin(q4)*sin(q5) + L12*sin(q4)*sin(q5)) + sin(q3)*(L4*cos(q4)*sin(q5) - L11*sin(q4) + L12*cos(q4)*sin(q5))) + cos(q2)*(sin(q3)*(L11*cos(q4) + L4*sin(q4)*sin(q5) + L12*sin(q4)*sin(q5)) + L5*sin(q3) - cos(q3)*(L4*cos(q4)*sin(q5) - L11*sin(q4) + L12*cos(q4)*sin(q5))) + L3*sin(q2)),-sin(q1)*(cos(q2)*(L5*cos(q3) + cos(q3)*(L11*cos(q4) + L4*sin(q4)*sin(q5) + L12*sin(q4)*sin(q5)) + sin(q3)*(L4*cos(q4)*sin(q5) - L11*sin(q4) + L12*cos(q4)*sin(q5))) - sin(q2)*(sin(q3)*(L11*cos(q4) + L4*sin(q4)*sin(q5) + L12*sin(q4)*sin(q5)) + L5*sin(q3) - cos(q3)*(L4*cos(q4)*sin(q5) - L11*sin(q4) + L12*cos(q4)*sin(q5))) + L3*cos(q2)),-sin(q1)*(cos(q2)*(L5*cos(q3) + cos(q3)*(L11*cos(q4) + L4*sin(q4)*sin(q5) + L12*sin(q4)*sin(q5)) + sin(q3)*(L4*cos(q4)*sin(q5) - L11*sin(q4) + L12*cos(q4)*sin(q5))) - sin(q2)*(sin(q3)*(L11*cos(q4) + L4*sin(q4)*sin(q5) + L12*sin(q4)*sin(q5)) + L5*sin(q3) - cos(q3)*(L4*cos(q4)*sin(q5) - L11*sin(q4) + L12*cos(q4)*sin(q5)))),-sin(q1)*((L4*cos(q2 + q3 + q4 - q5))/2 - (L12*cos(q2 + q3 + q4 + q5))/2 - (L4*cos(q2 + q3 + q4 + q5))/2 + (L12*cos(q2 + q3 + q4 - q5))/2 + L11*cos(q2 + q3 + q4)),-(L4 + L12)*(cos(q2)*cos(q5)*sin(q1)*sin(q3)*sin(q4) - cos(q2)*cos(q3)*cos(q4)*cos(q5)*sin(q1) - cos(q1)*sin(q5) + cos(q3)*cos(q5)*sin(q1)*sin(q2)*sin(q4) + cos(q4)*cos(q5)*sin(q1)*sin(q2)*sin(q3)),0;
0,(L4*sin(q2 + q3 + q4 + q5))/2 + (L12*sin(q2 + q3 + q4 + q5))/2 - L5*sin(q2 + q3) - L3*sin(q2) - (L4*sin(q2 + q3 + q4 - q5))/2 - (L12*sin(q2 + q3 + q4 - q5))/2 - L11*sin(q2 + q3 + q4),(L4*sin(q2 + q3 + q4 + q5))/2 + (L12*sin(q2 + q3 + q4 + q5))/2 - L5*sin(q2 + q3) - (L4*sin(q2 + q3 + q4 - q5))/2 - (L12*sin(q2 + q3 + q4 - q5))/2 - L11*sin(q2 + q3 + q4),(L4*sin(q2 + q3 + q4 + q5))/2 + (L12*sin(q2 + q3 + q4 + q5))/2 - (L4*sin(q2 + q3 + q4 - q5))/2 - (L12*sin(q2 + q3 + q4 - q5))/2 - L11*sin(q2 + q3 + q4),((sin(q2 + q3 + q4 + q5) + sin(q2 + q3 + q4 - q5))*(L4 + L12))/2,0;
0,sin(q1),sin(q1),sin(q1),-sin(q2 + q3 + q4)*cos(q1),cos(q5)*sin(q1) + cos(q2 + q3 + q4)*cos(q1)*sin(q5);
0,-cos(q1),-cos(q1),-cos(q1),-sin(q2 + q3 + q4)*sin(q1),cos(q2 + q3 + q4)*sin(q1)*sin(q5) - cos(q1)*cos(q5);
1,0,0,0,cos(q2 + q3 + q4),sin(q2 + q3 + q4)*sin(q5)];

Ja = T_Ainv*Jef;

Jalambda = Ja*Ja' +sym(eye(6))*lambda^2;

pJa = sym(zeros(6,6));

for i = 1:6
    
    pJa = simplify(diff(Ja,q(i))*qip(i)+ pJa);
end

for j = 1:3
    
    pJa = simplify(diff(Ja,phi(j))*pphi(j)+ pJa);
end

pJalambda = sym(zeros(6,6));

for i = 1:6
    
    pJalambda = simplify(diff(Ja,q(i))*qip(i)+ pJalambda);
end

for j = 1:3
    
    pJalambda = simplify(diff(Ja,phi(j))*pphi(j)+ pJalambda);
end

for i = 1:6
    for j = 1:6
WriteToFile(fid,Ja(i,j),strcat('Ja(',num2str(i),',',num2str(j),')'),'noSpace');
    end
end

fprintf(fid,'\n');
fprintf(fid,'%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%');
fprintf(fid,'\n');
fprintf(fid,'\n');

for i = 1:6
    for j = 1:6
WriteToFile(fid,pJa(i,j),strcat('pJa(',num2str(i),',',num2str(j),')'),'noSpace');
    end
end

fprintf(fid,'\n');
fprintf(fid,'%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%');
fprintf(fid,'\n');
fprintf(fid,'\n');

WriteToFile(fid,Jalambda,'Jalambda');

WriteToFile(fid,pJalambda,'pJalambda');

fclose(fid);
