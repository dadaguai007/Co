function []=Constellation_plot(input)
DataLen = length(input);
xyz(:,1)=real(input);
xyz(:,2)=imag(input);
z3 = round(rand(1,DataLen)*4+1);
xyz(:,3)=reshape(z3,DataLen,1);
xyz(:,3)=xyz(:,3)-xyz(:,3);
sizeXYZ = size(xyz);
searchR = 0.2;%搜索圆半径
for i=1:sizeXYZ(1)
    index_i = find(xyz(:,1)>(xyz(i,1)-searchR) & xyz(:,1)<xyz(i,1)+searchR & xyz(:,2) > xyz(i,2) -searchR & xyz(:,2)<xyz(i,2)+searchR );
    sizeIndexI = size(index_i);
    xyz(i,3)=sizeIndexI(1); % 找到其周围半径内0.2的点的个数，越大，则颜色的RGB值越大
end
[sortXYZ,sortI]=sort(xyz(:,3));
sz = 5;% 圆圈的大小
scatter(xyz(:,1),xyz(:,2),sz,xyz(:,3),'filled');%% filled表示圆圈为实心
% scatter(xyz(sortI,1),xyz(sortI,2),sz,xyz(sortI,3),'filled');%% filled表示圆圈为实心
hold on ;
    






end