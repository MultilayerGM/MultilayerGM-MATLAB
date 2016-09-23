function A=DCSBMNetworkGenerator(S,varargin)
% Degree-corrected block-model benchmark
% as suggested in
%   Karrer, B., & Newman, M. E. J. (2011). Stochastic blockmodels and
%   community structure in networks. Physical Review E, 83(1), 016107.
%   http://doi.org/10.1103/PhysRevE.83.016107
%
% Additionally adds rejection of multiedges in the sampling process.
%
% Input:
%
%   S:  multilayer partition
%
% Options:
%
%   parameters for degree_corrected block-model (LFR-like) network
%   generation in any format accepted by OptionSturct (i.e., key-value
%   list, struct, ...)
%
%   exponent: [-3] exponent for the powerlaw degree distribution
%   kmin:  [3] minimum degree
%   kmax:  [50] maximum degree
%   mu: [0.1] mixing parameter
%   maxreject: [100] maximum number of rejections for a single block before
%       bailing out and issuing a warning, leaving the block with less than
%       the desired number of edges (no multi-edges are generated).
%
% Output:
%
%   A:  cell array of adjacency matrices for each layer of the sampled
%       multilayer network


options=OptionStruct('exponent',-3,'kmin',3,'kmax',50,'mu',.1,'maxreject',100);
options.set(varargin);

mu=options.mu;
exponent=options.exponent;
kmin=options.kmin;
kmax=options.kmax;
max_reject=options.maxreject;

[n_nodes,n_layers]=size(S);
A=cell(n_layers,1);
for layer=1:n_layers
    [uc,~,ind]=unique(S(:,layer));
    nc=length(uc);
    G=sparse(1:n_nodes,ind,true,n_nodes,nc);
    group_sizes=sum(G,1);
    k=PowerlawSampler(n_nodes,exponent,kmin,kmax);
    k_g=arrayfun(@(i) sum(k(G(:,i))),1:nc);
    mm=sum(k);
    pos=0;
    neighbours=cell(n_nodes,1);
    for i=1:n_nodes
        neighbours{i}=i;
    end
    for group1=1:nc
        for group2=group1:nc
            % compute expected number of edges based on group
            % strength and mixing parameter
            if group1==group2
                w=((1-mu)*k_g(group1)+mu*(k_g(group1).^2)./mm);
            else
                w=mu*(k_g(group1)*k_g(group2))./mm;
            end
            % actual number of edges is poisson distributed
            m=poissrnd(w);
            
            if group1==group2
                dense=2*m>group_sizes(group1)*(group_sizes(group1)-1);
                m=m/2;
            else
                dense=2*m>group_sizes(group1)*group_sizes(group2);
            end
            
            
            if dense
                disp('dense');
                
                % use binomial sampling if block is dense
                sigma1=k(G(:,group1))./sum(k(G(:,group1)));
                sigma2=k(G(:,group2))./sum(k(G(:,group2)));
                P=sigma1*w*sigma2';
                nodes_g1=find(G(:,group1))';
                nodes_g2=find(G(:,group2))';
                ng1=group_sizes(group1);
                ng2=group_sizes(group2);
                if group1==group2
                    for it=1:ng1
                        nid=find(P(it,it+1:end)>rand(1,ng1-it))+it;
                        neighbours{nodes_g1(it)}=[neighbours{nodes_g1(it)},nodes_g2(nid)];
                        for it2=1:length(nid)
                            neighbours{nodes_g2(nid(it2))}=[neighbours{nodes_g2(nid(it2))},nodes_g1(it)];
                        end
                    end
                else
                    for it=1:ng1
                        nid=find(P(it,:)>rand(1,ng2));
                        neighbours{nodes_g1(it)}=[neighbours{nodes_g1(it)},nodes_g2(nid)];
                        for it2=1:length(nid)
                            neighbours{nodes_g2(nid(it2))}=[neighbours{nodes_g2(nid(it2))},nodes_g1(it)];
                        end
                    end
                end
                
            else
                
                % sample the edges
                sigma1=cumsum(k(G(:,group1)));
                sigma1=sigma1./sigma1(end);
                nodes_g1=find(G(:,group1));
                sigma2=cumsum(k(G(:,group2)));
                sigma2=sigma2./sigma2(end);
                nodes_g2=find(G(:,group2));
                pos=pos+m;
                for e=1:m
                    isneighbour=true;
                    % rejection sampling to avoid multi-edges and
                    % self-loops
                    reject_count=0;
                    while isneighbour&&reject_count<=max_reject
                        decide=rand(2,1);
                        n1=nodes_g1(find(sigma1>decide(1),1));
                        n2=nodes_g2(find(sigma2>decide(2),1));
                        isneighbour=any(n1==neighbours{n2}(:));
                        reject_count=reject_count+1;
                    end
                    if reject_count>max_reject
                        warning('stopping sampling of edges for current block')
                        break;
                    end
                    neighbours{n1}=[neighbours{n1},n2];
                    neighbours{n2}=[neighbours{n2},n1];
                end
            end
        end
    end
    
    % convert neighbour-list to adjacency matrix (note that first
    % neighbour is node itself and should not be included in output)
    indrow=zeros(2*pos,1);
    indcol=zeros(2*pos,1);
    pos=0;
    for i=1:n_nodes
        indrow(pos+1:pos+length(neighbours{i})-1)=i;
        indcol(pos+1:pos+length(neighbours{i})-1)=neighbours{i}(2:end);
        pos=pos+length(neighbours{i})-1;
    end
    indrow=indrow(1:pos);
    indcol=indcol(1:pos);
    A{layer}=sparse(indrow,indcol,1,n_nodes,n_nodes);
end
end





