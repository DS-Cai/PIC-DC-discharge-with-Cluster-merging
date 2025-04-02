load("OUTPUT\meanee_xt.dat")
x=N_XT;
y=N_G;
ne_xt_r=reshape(meanee_xt,x,y);
ne_xt_r(find(isnan(ne_xt_r())))=0;
[X,Y]=meshgrid(1:x,1:y);
mesh(X,Y,ne_xt_r.')
m_ne_xt_r=mean(ne_xt_r,2);
NAMEUVP='uvp.dat';
outname=['E:\Matlab code\等离子模拟\1_D discharge\Merging-gai\OUTPUT\',NAMEUVP];
[fileID, message] = fopen(outname,'w');
if fileID < 0
    error( '无法打开我的文件因为: %s' , message);
end
fprintf(fileID,'%s\n','VARIABLES= "X","Y","ne_xt","ni_xt","pot_xt"');
fprintf(fileID,'%s%d%s%d%s\n','ZONE T="Floor", I=', x, ' J=',y,' F=POINT');
for J=1:y
    for I=1:x
        fprintf(fileID,'%d %d %d %d %d\n',I,J,ne_xt((J-1)*x+I),ni_xt((J-1)*x+I),pot_xt((J-1)*x+I));
    end
end
fclose(fileID);