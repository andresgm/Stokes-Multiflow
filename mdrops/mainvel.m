%%
tic
clear;
clc;

%INTEGRATION (se ejecuta después de tener el archivo generado por
%"main_inout"
data=open([cd '/gotaLambda0.2Ca0.35/it1000.mat']); 

CurDir = cd;
    DirRoot = [CurDir '/gotaLambda0.2Ca0.35'];
    NombreSoluSim = ['L' num2str(data.parms.lamda)];
    DirResult = ['/' NombreSoluSim '/'];
    FilePath = [DirRoot '/'];
    Filename=[FilePath 'gotaLambda0.2Ca0.35inout_L' num2str(data.parms.lamda), '.mat'];    
data2=open(Filename);


%especificar el numero de nodos de la malla
nodesmesh=data2.numnodos;

% Llenar el slide    
m=1;
p=1;
slide = 0;
% for i=1:nodesmesh
i = 17;
    for j=1:nodesmesh
        for k=1:nodesmesh
            coordinates(m,1)=data2.coords(i,1);
            coordinates(m,2)=data2.coords(j,2);
            coordinates(m,3)=data2.coords(k,3);
            test(m)=data2.gota(i,j,k);
            m=m+1;
        end

        VelAtNode=integrationf(data,coordinates);   
        
    
        for h=1:m-1
            results(p,1)=coordinates(h,1); 
            results(p,2)=coordinates(h,2);
            results(p,3)=coordinates(h,3);
            results(p,4)=VelAtNode(h,1);     
            results(p,5)=VelAtNode(h,2);
            results(p,6)=VelAtNode(h,3);

            if p~=1 && results(p,1)==results(p-1,1)
                results(p,7)=slide;
            else
            slide=slide+1;
            results(p,7)=slide;
            end


%Para caso de lambda diferente de cero se corrije la velocidad en los
%puntos en el interior de la gota            
             if data.parms.lamda~=0
                 if test(h)==1
                     results(p,4)=0;     
                     results(p,5)=0;
                     results(p,6)=0;
%                      results(p,4)=VelAtNode(h,1)/data.parms.lamda;     
%                      results(p,5)=VelAtNode(h,2)/data.parms.lamda;
%                      results(p,6)=VelAtNode(h,3)/data.parms.lamda;
                 end       
             end
             p=p+1;
        end     
        m=1;    
    end
% end
       

        
%save ('velocidades_Ca0.25_L1_curv1_grav0_bending0_mesh3_65.mat','results')
toc

% graficacion del campo de velocidades 

% nombrearchivo = 'velocidades_Ca0.25_L1_curv1_grav0_bending0_mesh3.mat'
% load(nombrearchivo);
% num del slide (33 de con x1 = 0.0397 no hay slide xon x1 = 0 -seccion
% central de la gota
slidenum = 1;
% numero de posiciones a extraer
pos = 5;
% factor de amplificacion de los vectores
velfact = 1;

slidef = results(results(:,7) == slidenum,:);
slidef = slidef(1:pos:end,:);
quiver3(slidef(:,1),slidef(:,2),slidef(:,3),slidef(:,4),slidef(:,5),slidef(:,6),velfact);
view(90,0); axis auto
% pause
% 
% %% RHEOLOGY
% 
% %alpha=4*pi/(3*5*5*5);
% alpha=1;
% 
% data=open([cd '/Results/L1/Ca0.05.mat']); 
% figure
% trisurf(data.geom.elements,data.geom.nodes(:,1),data.geom.nodes(:,2),data.geom.nodes(:,3),data.geom.curv); colorbar
% data2=open([cd '/Results/L1/Ca0.1.mat']);
% figure
% trisurf(data2.geom.elements,data2.geom.nodes(:,1),data2.geom.nodes(:,2),data2.geom.nodes(:,3),data2.geom.curv); colorbar
% data3=open([cd '/Results/L1/Ca0.15.mat']);
% figure
% trisurf(data3.geom.elements,data3.geom.nodes(:,1),data3.geom.nodes(:,2),data3.geom.nodes(:,3),data3.geom.curv); colorbar
% data4=open([cd '/Results/L1/Ca0.2.mat']);
% figure
% trisurf(data4.geom.elements,data4.geom.nodes(:,1),data4.geom.nodes(:,2),data4.geom.nodes(:,3),data4.geom.curv); colorbar
% data5=open([cd '/Results/L1/Ca0.25.mat']);
% figure
% trisurf(data5.geom.elements,data5.geom.nodes(:,1),data5.geom.nodes(:,2),data5.geom.nodes(:,3),data5.geom.curv); colorbar
% data6=open([cd '/Results/L1/Ca0.3.mat']);
% figure
% trisurf(data6.geom.elements,data6.geom.nodes(:,1),data6.geom.nodes(:,2),data6.geom.nodes(:,3),data6.geom.curv); colorbar
% 
% [deltaf,deltafconst] = DeltaF(data.geom,data.deltafconst);
% [deltaf2,deltafconst2] = DeltaF(data2.geom,data2.deltafconst);
% [deltaf3,deltafconst3] = DeltaF(data3.geom,data3.deltafconst);
% [deltaf4,deltafconst4] = DeltaF(data4.geom,data4.deltafconst);
% [deltaf5,deltafconst5] = DeltaF(data5.geom,data5.deltafconst);
% [deltaf6,deltafconst6] = DeltaF(data6.geom,data6.deltafconst);
% 
% for i=1:data.geom.numnodes
% norma=norm(data.geom.normal(i,:));
% deltafvec(i,1)=(deltaf.mag(i)/norma)*data.geom.normal(i,1);
% deltafvec(i,2)=(deltaf.mag(i)/norma)*data.geom.normal(i,2);
% deltafvec(i,3)=(deltaf.mag(i)/norma)*data.geom.normal(i,3);
% end
% 
% for i=1:data2.geom.numnodes
% norma=norm(data2.geom.normal(i,:));
% deltaf2vec(i,1)=(deltaf2.mag(i)/norma)*data2.geom.normal(i,1);
% deltaf2vec(i,2)=(deltaf2.mag(i)/norma)*data2.geom.normal(i,2);
% deltaf2vec(i,3)=(deltaf2.mag(i)/norma)*data2.geom.normal(i,3);
% end
% 
% for i=1:data3.geom.numnodes
% norma=norm(data3.geom.normal(i,:));
% deltaf3vec(i,1)=(deltaf3.mag(i)/norma)*data3.geom.normal(i,1);
% deltaf3vec(i,2)=(deltaf3.mag(i)/norma)*data3.geom.normal(i,2);
% deltaf3vec(i,3)=(deltaf3.mag(i)/norma)*data3.geom.normal(i,3);
% end
% 
% for i=1:data4.geom.numnodes
% norma=norm(data4.geom.normal(i,:));
% deltaf4vec(i,1)=(deltaf4.mag(i)/norma)*data4.geom.normal(i,1);
% deltaf4vec(i,2)=(deltaf4.mag(i)/norma)*data4.geom.normal(i,2);
% deltaf4vec(i,3)=(deltaf4.mag(i)/norma)*data4.geom.normal(i,3);
% end
% 
% for i=1:data5.geom.numnodes
% norma=norm(data5.geom.normal(i,:));
% deltaf5vec(i,1)=(deltaf5.mag(i)/norma)*data5.geom.normal(i,1);
% deltaf5vec(i,2)=(deltaf5.mag(i)/norma)*data5.geom.normal(i,2);
% deltaf5vec(i,3)=(deltaf5.mag(i)/norma)*data5.geom.normal(i,3);
% end
% 
% for i=1:data6.geom.numnodes
% norma=norm(data6.geom.normal(i,:));
% deltaf6vec(i,1)=(deltaf6.mag(i)/norma)*data6.geom.normal(i,1);
% deltaf6vec(i,2)=(deltaf6.mag(i)/norma)*data6.geom.normal(i,2);
% deltaf6vec(i,3)=(deltaf6.mag(i)/norma)*data6.geom.normal(i,3);
% end
% 
% 
% 
% strtensornode = strtensorreologia(deltafvec,data.geom.nodes,data.VelAtNode,data.geom.normal,data.stokesprobconst.lamda);
% strtensornode2 = strtensorreologia(deltaf2vec,data2.geom.nodes,data2.VelAtNode,data2.geom.normal,data2.stokesprobconst.lamda);
% strtensornode3 = strtensorreologia(deltaf3vec,data3.geom.nodes,data3.VelAtNode,data3.geom.normal,data3.stokesprobconst.lamda);
% strtensornode4 = strtensorreologia(deltaf4vec,data4.geom.nodes,data4.VelAtNode,data4.geom.normal,data4.stokesprobconst.lamda);
% strtensornode5 = strtensorreologia(deltaf5vec,data5.geom.nodes,data5.VelAtNode,data5.geom.normal,data5.stokesprobconst.lamda);
% strtensornode6 = strtensorreologia(deltaf6vec,data6.geom.nodes,data5.VelAtNode,data6.geom.normal,data6.stokesprobconst.lamda);
% 
% rintegral = alpha*trapecioreologia(data.geom.dsi,strtensornode);
% rintegral2 = alpha*trapecioreologia(data2.geom.dsi,strtensornode2);
% rintegral3 = alpha*trapecioreologia(data3.geom.dsi,strtensornode3);
% rintegral4 = alpha*trapecioreologia(data4.geom.dsi,strtensornode4);
% rintegral5 = alpha*trapecioreologia(data5.geom.dsi,strtensornode5);
% rintegral6 = alpha*trapecioreologia(data6.geom.dsi,strtensornode6);
% 
% x(1)=data.stokesprobconst.ca;
% x(2)=data2.stokesprobconst.ca;
% x(3)=data3.stokesprobconst.ca;
% x(4)=data4.stokesprobconst.ca;
% x(5)=data5.stokesprobconst.ca;
% x(6)=data6.stokesprobconst.ca;
% 
% y1(1)=rintegral(1,1)-rintegral(2,2);
% y1(2)=rintegral2(1,1)-rintegral2(2,2);
% y1(3)=rintegral3(1,1)-rintegral3(2,2);
% y1(4)=rintegral4(1,1)-rintegral4(2,2);
% y1(5)=rintegral5(1,1)-rintegral5(2,2);
% y1(6)=rintegral6(1,1)-rintegral6(2,2);
% 
% 
% y2(1)=rintegral(2,2)-rintegral(3,3);
% y2(2)=rintegral2(2,2)-rintegral2(3,3);
% y2(3)=rintegral3(2,2)-rintegral3(3,3);
% y2(4)=rintegral4(2,2)-rintegral4(3,3);
% y2(5)=rintegral5(2,2)-rintegral5(3,3);
% y2(6)=rintegral6(2,2)-rintegral6(3,3);
% 
% y3(1)=rintegral(2,3);
% y3(2)=rintegral2(2,3);
% y3(3)=rintegral3(2,3);
% y3(4)=rintegral4(2,3);
% y3(5)=rintegral5(2,3);
% y3(6)=rintegral6(2,3);
% 
% y4(1)=rintegral(1,1)/3-2;
% y4(2)=rintegral2(1,1)/3-2;
% y4(3)=rintegral3(1,1)/3-2;
% y4(4)=rintegral4(1,1)/3-2;
% y4(5)=rintegral5(1,1)/3-2;
% y4(6)=rintegral6(1,1)/3-2;
% 
% 
% y5(1)=rintegral(2,2)/3-0.09;
% y5(2)=rintegral2(2,2)/3-0.09;
% y5(3)=rintegral3(2,2)/3-0.09;
% y5(4)=rintegral4(2,2)/3-0.09;
% y5(5)=rintegral5(2,2)/3-0.09;
% y5(6)=rintegral6(2,2)/3-0.09;
% 
% y6(1)=rintegral(3,3)/3-0.09;
% y6(2)=rintegral2(3,3)/3-0.09;
% y6(3)=rintegral3(3,3)/3-0.09;
% y6(4)=rintegral4(3,3)/3-0.09;
% y6(5)=rintegral5(3,3)/3-0.09;
% y6(6)=rintegral6(3,3)/3-0.09;
% 
% figure
% 
% Name = ['Tensor de Esfuerzo para \lambda=' num2str(data.stokesprobconst.lamda)];
% plot(x,-y1,'-b+',x,-y2,'-ro',x,y3,'-gs',x,y4,'-c>','LineWidth',2,'MarkerEdgeColor','k','MarkerSize',8)
% xlabel('C','fontsize',30)
% ylabel('\tau_{ij}','fontsize',30)
% set(gca,'fontsize',24)
% title(Name,'fontsize',48);
% grid on;
% leyenda=legend('\tau_{11}-\tau_{22}','\tau_{22}-\tau_{33}','\tau_{12}','(\tau_{ii}/{3})-2',2);
% set(leyenda,'fontsize',18)
% 
% y7(1)=rintegral(1,2)/x(1);
% y7(2)=rintegral2(1,2)/x(2);
% y7(3)=rintegral3(1,2)/x(3);
% y7(4)=rintegral4(1,2)/x(4);
% y7(5)=rintegral5(1,2)/x(5);
% y7(6)=rintegral6(1,2)/x(6);
% 
% %figure
% %plot(x,y7/alpha)
% toc 


