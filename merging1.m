function [x,vx,vy,vz,WEIGHT]=merging1(x,vx,vy,vz,WEIGHT,DX,N_G,constraint_num)
% x=x_e;
% vx=vx_e;
% vy=vy_e;
% vz=vz_e;
% WEIGHT=WEIGHT_E;
for i=1:N_G-1
    [x_i]=find(x>(i-1)*DX & x<i*DX );
    len_x=length(x_i);
    size_min=constraint_num+1;
    size_max=2*size_min-1;
    if len_x>50
        n_clusters=ceil(len_x/size_max); 
        norm_vx=norm(vx(x_i),2);
        norm_vy=norm(vy(x_i),2);
        norm_vz=norm(vy(x_i),2);
        norm_v=(norm_vx+norm_vy+norm_vz)/3;
        v=[vx(x_i)./norm_v,vy(x_i)./norm_v,vz(x_i)./norm_v];
        X = py.numpy.array(v);
        clf =py.k_means_constrained.KMeansConstrained(n_clusters=py.int(n_clusters),...
            size_min=py.int(size_min),size_max=py.int(size_max),random_state=py.int(0));
        clfX=clf.fit(X);
        lables=double(clfX.fit_predict(X));
        cluster_centers_=double(clfX.cluster_centers_);
        for j=1:n_clusters
            clusters_j=x_i(find(lables==j-1));
            cluster_sqr=sum(sum(sqrt((v(find(lables==j-1),:)-cluster_centers_(j,:)).^2)))/3/length(clusters_j)/norm(cluster_centers_(j,1));
            non_W=find(WEIGHT(clusters_j)>1e-5);
            w_i=0; l_nw=length(non_W);
            while length(non_W)>constraint_num
                non_W_k=non_W(1:constraint_num+1);
                wa=zeros(constraint_num+1,1);
                wa=solve_linear_equation2(x(clusters_j(non_W_k)),vx(clusters_j(non_W_k)),...
                    vy(clusters_j(non_W_k)),vz(clusters_j(non_W_k)),WEIGHT(clusters_j(non_W_k)));
                WEIGHT(clusters_j(non_W_k))=wa;
                non_W=find(WEIGHT(clusters_j)>1e-5);
                if w_i>l_nw+2-constraint_num
                    break;
                end
                w_i=w_i+1;
            end
        end
    end
end
non_W=find(WEIGHT<1e-5);
non_W_k=sort(non_W, 'descend');
for k=1:length(non_W_k)
    x(non_W_k(k))=x(end);
    x(end)=[];
    vx(non_W_k(k))=vx(end);
    vx(end)=[];
    vy(non_W_k(k))=vy(end);
    vy(end)=[];
    vz(non_W_k(k))=vz(end);
    vz(end)=[];
    WEIGHT(non_W_k(k))=WEIGHT(end); 
    WEIGHT(end)=[];
end