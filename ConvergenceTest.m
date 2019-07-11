function ConvergenceTest(nodes,layers,types,dependencyMatrix,nullDistribution,varargin)
    % Test convergence of 'PartitionGenerator'.
    %
    % Input:
    %
    %   nodes: Either a scalar specifying the number of physical nodes or
    %       a |V_M|x(1+d) matrix specifying each state node in the format
    %       [node, aspect_1,...,aspect_d] (the latter format allows for missing
    %       nodes in some layers)
    %
    %   layers: Vector giving the number of elements for each aspect in the
    %       format [l_1,...,l_d]. Note that for a multilayer network with a
    %       single aspect this is just a scalar giving the number of layers.
    %
    %   types: 'char' vector specifying the update type for each aspect, 'o' for
    %       ordered and 'r' for random.
    %
    %   dependencyMatrix: A matrix of copying probabilities, corresponding to
    %       the flattened interlayer dependency tensor. This matrix is either
    %       of size lxl for the layer-coupled case or of size (n*l)x(n*l) in
    %       the general case. If state nodes are given explicitly in 'nodes',
    %       'dependencyMatrix' should be of size |V_M|x|V_M|, giving only the
    %       probabilities for state nodes that are actually present in the
    %       network.  The ways the matrix encodes the tensor are as follows:
    %
    %       layer-coupled case:
    %           Mapping for flattening the indices:
    %
    %           a=aspect_1+l_1*(aspect_2-1)+...+l_1*l_2*...*l_(d-1)*(aspect_d-1)
    %
    %           where 'dependencyMatrix(a,b)' is the probability that a node in
    %           layer b copies its community assignment from the same node in
    %           layer a. The rows of 'dependencyMatrix' should sum to a value<1,
    %           where 1-sum(transtionMatrix(:,b)) is the probability of
    %           resampling the community assignment from the specified null
    %           distribution for a node in layer b.
    %
    %       general case:
    %           Mapping for flattening the indeces when state nodes are not
    %           explicitly specified:
    %
    %           i=node+n*(aspect_1-1)+n*l_1*(aspect_2-1)+..._n*l_1*...*l_(d-1)*(aspect_d-1)
    %
    %           When state nodes are given explicitly, the transition matrix
    %           should have the same number of rows as the matrix of state
    %           nodes, where the ith state node corresponds to the ith row of
    %           the transition matrix.
    %
    %           'dependencyMatrix(i,j)' is the probability that
    %           state node j copies its community assignment from state node i.
    %           The rows of 'dependencyMatrix' should sum to a value < 1,
    %           where 1-sum(transtionMatrix(:,j)) is the probability of
    %           resampling the community assignment from the specified null
    %           distribution for node j.
    %
    %   nullDistribution: A function that takes state nodes (i.e., row vectors
    %       of the form [node, aspect_1,...,aspect_d] ) as input and returns a
    %       random community assignment.
    %
    %
    %
    %
    % Options:
    %
    %   Lags: [default: [1,5,10,20]] Lags to plot (note that the unit for the Lags
    %       depends on the number of 'UpdateSteps' specified)
    %
    %   UpdateSteps: [default: 1] number of Gibbs updates between test partitions
    %       (i.e., the expected number of times each state node's community
    %       assignment is updated).
    %
    %   NumberOfSamples: 200
    %
    %   InitialPartition: Optionally specify starting partition (array of
    %       dimension nxl_1x...xl_d).
    %
    %
    %
    % Version: 1.0.1
    % Date: Tue  4 Jul 2017 16:38:06 BST
    % Author: Lucas Jeub
    % Email: ljeub@iu.edu
    %
    %
    % References:
    %
    %       [1] Generative benchmark models for mesoscale structure in multilayer
    %       networks, M. Bazzi, L. G. S. Jeub, A. Arenas, S. D. Howison, M. A.
    %       Porter. arXiv1:608.06196.
    %
    % Citation:
    %
    %       If you use this code, please cite as
    %       Lucas G. S. Jeub and Marya Bazzi
    %       "A generative model for mesoscale structure in multilayer networks
    %       implemented in MATLAB," https://github.com/MultilayerBenchmark/MultilayerBenchmark (2016).

    parseArgs=inputParser();
    addParameter(parseArgs,'InitialPartition',[]);
    addParameter(parseArgs,'UpdateSteps',1);
    addParameter(parseArgs,'NumberOfSamples',200);
    addParameter(parseArgs,'Lags',[1,5,10,20])
    parse(parseArgs,varargin{:});
    options=parseArgs.Results;
    options.isset=@(s) ~isempty(options.(s));

    % determine network shape and set up state nodes if not given explicitly
    l=prod(layers);
    if numel(nodes)==1
        n=nodes;
        nodes=ind2subarray([n,layers],1:n*l);
    else
        n=max(nodes,[],1);
    end

    % Sample partitions
    S=cell(options.NumberOfSamples+1,1);
    if options.isset('InitialPartition')
        S{1}=options.InitialPartition;
    else
        S{1}=zeros([n,layers]);
        for i=1:size(nodes,1)
            S{1}(subarray2ind([n,layers],nodes(i,:)))=nullDistribution(nodes(i,:));
        end
    end
    for i=2:length(S)
        S{i}=PartitionGenerator(nodes,layers,types,dependencyMatrix,nullDistribution,'InitialPartition',S{i-1},'UpdateSteps',options.UpdateSteps);
    end

    for i=1:length(S)
        m{i}=nmi_matrix(S{i});
    end

    lag=options.Lags;
    lag=lag(lag<length(S)/3);
    map=interpolate_map([0,0.5,1],[0.2314,0.298,0.7529;0.3,0.3,0.3;0.71, 0, 0.15],linspace(0,1,length(lag)));
    styles=cycler({'-','--','-','--','-','--'});
    markers=cycler({'o','+','<','d','x','>'});
    sizes=cycler({2,3,3,3,3,3});
    for i=1:length(lag)
        for j=lag(i)+1:length(S)
            nmi_auto{i}(j-lag(i))=nmi(S{j-lag(i)}(:),S{j}(:));
            d_auto{i}(j-lag(i))=mean(abs(m{j-lag(i)}(triu(true(length(m{j})),1))-m{j}(triu(true(length(m{j})),1))));
        end
    end

    f(1)=figure();
    for i=1:length(lag)
        plot(nmi_auto{i},'color',map(i,:),'markerfacecolor',map(i,:),...
            'linestyle',styles(i),'marker',markers(i),'markersize',sizes(i))
        hold on
    end
    legend(arrayfun(@(lag) sprintf('lag=%u',lag),lag,'UniformOutput',false))
    xlabel('step')
    ylabel('mNMI')
    %set(gca,'xlim',[1,200])
    %set(gca,'xtick',[1,50:50:200])
    ylim=get(gca,'ylim');
    ylim(1)=0;
    set(gca,'ylim',ylim);
    %set(gca,'box','off','YGrid','on','YMinorGrid','off');
    %plot(nmi_first,'k')

    f(2)=figure();
    for i=1:length(lag)
        plot(d_auto{i},'color',map(i,:),'markerfacecolor',map(i,:),...
            'linestyle',styles(i),'marker',markers(i),'markersize',sizes(i));
        hold on
    end
    xlabel('step')
    ylabel('$d$')
    %set(gca,'xlim',[1,200])
    ylim=get(gca,'ylim');
    ylim(1)=0;
    set(gca,'ylim',ylim);
    %set(gca,'xtick',[1,50:50:200])
    %set(gca,'box','off','YGrid','on','YMinorGrid','off');
    %legend(arrayfun(@(lag) sprintf('lag=%u',lag),lag,'UniformOutput',false))

    %mark_indices=strjoin(arrayfun(@num2str,[1,10:10:200],'uniformoutput',false),',');
end
