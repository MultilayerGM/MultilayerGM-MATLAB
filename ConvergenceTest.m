function ConvergenceTest(S)
    for i=1:length(S)
        m{i}=nmi_matrix(S{i});
    end

    lag=[1,5,10,20];
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
