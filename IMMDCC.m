% Co-Clustering via Iterative Multi-mode Discrtization

% Paper: http://dx.doi.org/10.1007/978-3-030-61527-7_7
% YouTube: https://youtu.be/2mLd0n6VBa4
% How to Cite: Fanaee-T H., Thoresen M. (2020) Iterative Multi-mode Discretization: Applications to Co-clustering. In: Appice A., Tsoumakas  G., Manolopoulos Y., Matwin S. (eds) Discovery Science. DS 2020. Lecture Notes in Computer Science, vol 12323. Springer, Cham. https://doi.org/10.1007/978-3-030-61527-7_7

function [ba,noclusters]=IMMDCC(x,th1,th,minClusMem,topHotspotsNo,tscale)

%try
    
xc=x;
[n,m]=size(x);
i=0;

while(true)
    i=i+1;
    nsX{i}=best_scale(xc,1,tscale,tscale);
    [B,I]=min(nsX{i}(:,4));
    obj_X(i)=nsX{i}(I,4); 
    scale=nsX{i}(I,5);
    nsX_scale(i)=scale;
    middle_value=median(1:scale);
    ZX=matrix_disc(xc,scale,1);
    [a,b]=matrix_sep(x,ZX,scale);
    [xc,ZX_record{i}]=matrix_sep(x,ZX,scale);
    xc=xc{middle_value};
    nt=length(find(~isnan(ZX)));
    obj_X_n(i)=length(find(~isnan(ZX_record{i})));
    disp(['Step=',num2str(i),' Sparsity=',num2str(100*nt/(n*m)),'%']);
    if nt==0; break; end
end



xc=x;
i=0;

while(true)
    i=i+1;
    nsY{i}=best_scale(xc,2,tscale,tscale);
    [B,I]=min(nsY{i}(:,4));
    obj_Y(i)=nsY{i}(I,4);
    scale=nsY{i}(I,5);
    nsY_scale(i)=scale;
    middle_value=median(1:scale);
    ZY=matrix_disc(xc,scale,2);
    [xc,ZY_record{i}]=matrix_sep(x,ZY,scale);
    xc=xc{middle_value};
    nt=length(find(~isnan(ZY)));
    obj_Y_n(i)=length(find(~isnan(ZY_record{i})));
    disp(['Step=',num2str(i),' Sparsity=',num2str(100*nt/(n*m)),'%']);
    if nt==0; break; end
end

%**********************************************************************
% Choose best iterations in x,y
%**********************************************************************
clear optm optm2 optm3;
k=0;
for i=length(obj_X):-1:1
    k=k+1;
    for j=length(obj_Y):-1:1
        T=ZX_record{i}+ZY_record{j};
        optm2(i,j)=length(find(~isnan(ZX_record{i})))+length(find(~isnan(ZY_record{j})));
        optm(i,j)=length(find(T==2));
        optm3(i,j)=optm(i,j)/optm2(i,j);
    end
end

[B,I]=min(obj_X);
ZX_final=ZX_record{I};
nsX_scale_sel=I;
[B,J]=min(obj_Y);
ZY_final=ZY_record{J};
nsY_scale_sel=J;


T1=zeros(size(ZX_final));
T1(find(~isnan(ZX_final)))=1;

T2=zeros(size(ZY_final));
T2(find(~isnan(ZY_final)))=1;

T=T1+T2;

idx=find(T==2);
[hotspots(:,1),hotspots(:,2)] = ind2sub([n,m],idx);


GMx=zeros(m,m);
GMy=zeros(n,n);
GM=zeros(n,m);
for i=1:size(hotspots,1)
    hx=hotspots(i,1);
    hy=hotspots(i,2);
    Lx1=find(T1(hx,:));
    Ly2=find(T2(:,hy));
    Lx=find(T(hx,:));
    Ly=find(T(:,hy));
    GM(Ly,Lx)=GM(Ly,Lx)+1;
    GMs{i}=GM;
    GMx(Lx1,Lx1)=GMx(Lx1,Lx1)+1;
    GMy(Ly2,Ly2)=GMy(Ly2,Ly2)+1;
end

for i=1:size(hotspots,1)
    vl=GM(hotspots(i,1),hotspots(i,2));
    if vl==0
        vl=NaN;
    end
    hotspots(i,3)=vl;
end    

[B,I]=sort(hotspots(:,3),'descend');
hotspots=hotspots(I,:);
hotspots=sortrows(hotspots,3,'desc');
[ZH]=zscore(hotspots(:,3));


clear clus;
c=1;
Xmax=max(max(GMx));
Ymax=max(max(GMy));
XYmax=max(max(GM));
clus{c}=hotspots(1,1:2);
topHotspotsNo=min([topHotspotsNo, size(hotspots,1)-1]);
for i=1:topHotspotsNo
	cand=[];
    for j=1:c
            P=hotspots(i+1,3)/XYmax;
            PconnX=GMx(clus{j}(end,2),hotspots(i+1,2))/Xmax;
            PconnY=GMy(clus{j}(end,1),hotspots(i+1,1))/Ymax;
            disp(['i=',num2str(i),' C=',num2str(j),' PX(',num2str(clus{j}(end,2)),',>',num2str(hotspots(i+1,2)),')=',num2str(PconnX),' PY(',num2str(clus{j}(end,1)),',>',num2str(hotspots(i+1,1)),')=',num2str(PconnY)]);
            if PconnX>th && PconnY>th
                cand(j)=PconnX+PconnY;
            end 
    end
    if ~isempty(cand)
            [Pselected,selectedC]=max(cand);
            clus{selectedC}=[clus{selectedC};hotspots(i+1,1:2)];
            disp(['--> i=',num2str(i),' Cluster=',num2str(selectedC)]);
    else
        if P>th1
        	c=c+1;
            clus{c}=hotspots(i+1,1:2);
            disp(['--> i=',num2str(i),' Cluster ',num2str(c), ' added',' P=',num2str(P)]);
        end
    end
end    

for i=1:length(clus)
   cL(i)=size(clus{i},1);
end


[cL,I]=sort(cL,'descend');
k=0;
for i=1:length(clus)
    if size(clus{I(i)},1)>minClusMem
        k=k+1;
        COCLUSTERS{k}=clus{I(i)};
        clus{I(i)};
    end
end


hostspots_selX=length(find(~isnan(ZX_final)))
hostspots_selY=length(find(~isnan(ZY_final)))
ta=0;
xs=zeros(size(x));
for i=1:size(hotspots,1)
    xs(hotspots(i,1:2))=x(hotspots(i,1:2));
    if hotspots(i,1)<=50 && hotspots(i,1)<=50
        ta=ta+1;
    end
end
hotspots
disp(['True Found points=',num2str(ta),'/',num2str(size(hotspots,1)),' ',' P=%',num2str(100*ta/size(hotspots,1))])

if exist('COCLUSTERS')
    noclusters=length(COCLUSTERS);
    for i=1:length(COCLUSTERS)
        ba(i).rows=unique(COCLUSTERS{i}(:,1));
        ba(i).cols=unique(COCLUSTERS{i}(:,2));
    end
else
     noclusters=0
     ba(1).rows=[];
     ba(1).cols=[];
end
