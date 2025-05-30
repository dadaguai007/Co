function Eo = hybrid_2x4(E1, E2)


C_3dB = 1/sqrt(2).*[1,1i;1i,1];

M =[C_3dB,zeros(2);zeros(2),C_3dB];

Ei = [E1;zeros(size(E1));zeros(size(E1));E2];

U = [1,0,0,0;0,0,1,0;0,1,0,0;0,0,0,1i];

Eo = M*U*M*Ei;



end