function DA_LME_function_scatter_plot(yr_start,reg_sea)

    data = load('All_lme_offsets_and_diurnal_amplitudes.mat');

    num_of_colors        = 13;
    num_of_marker_shapes = 6;
    [out_st,out_col]     = CDF_linest(num_of_colors * num_of_marker_shapes, num_of_marker_shapes);
    
    % *********************************************************************
    % Set parameters
    % *********************************************************************
    N       = 100;              % Number of bootstrapping members
    n_sigma = 2;                % Number of s.d. for error bars
    yr0     = 1879;             % The first year of data
    
    PP.N             = N;
    PP.n_sigma       = n_sigma;
    PP.P.mute_output = 1;
    PP.pic_st        = out_st;
    PP.pic_col       = out_col;
    
    if reg_sea == 1,
        PP.eri_color = 'k';
    elseif ismember(reg_sea,[2 4]),
        PP.eri_color = 'b';
    elseif ismember(reg_sea,[3 5]),
        PP.eri_color = 'r';
    else
        error('Not a valid region-season ID')
    end
       
    % *********************************************************************
    % Prepare for data to be plotted
    % *********************************************************************
    x     = data.da(2:end,yr_start-yr0,reg_sea);
    y     = data.lme(2:end,yr_start-yr0,reg_sea);
    x_std = data.da_std(2:end,yr_start-yr0,reg_sea);
    y_std = data.lme_std(2:end,yr_start-yr0,reg_sea);
    n     = [];
    x_eri = data.da(1,yr_start-yr0,reg_sea);
    y_eri = data.lme(1,yr_start-yr0,reg_sea);
    x_buoy = nan;
    y_buoy = nan;

    % *********************************************************************
    % Generate plot
    % *********************************************************************
    figure((yr_start-yr0-1)*4+reg_sea); clf; hold on;
    DA_LME_plot_square_panels(x,y,x_std,y_std,n,x_eri,y_eri,x_buoy,y_buoy,PP)
    CDF_panel([0 0.6 -0.6 0.6],'','','Diurnal Amplitude (^oC)','Offsets (^oC)','fontsize',18);
    daspect([0.6 1.2 1])  
    set(gcf,'position',[.1 1 7 7],'unit','inches')
    set(gcf,'position',[.1 1 7 7],'unit','inches')


    % *********************************************************************
    % Generate legend
    % *********************************************************************
    group = data.grp(2:end,:); 
    in_name = {};
    for ct = 1:size(group,1);
        name = double(group(ct,:));
        if name(1) > 100
            in_name{ct} = ['      DCK ',num2str(name(3))];
        else
            in_name{ct} = [name(1:2),' DCK ',num2str(name(3))];
        end
    end
    
    figure(1000); clf; hold on;
    CDF_scatter_legend(in_name,out_st,out_col,[1,7,1],'fontsize',15,'mksize',13)
    axis([1 7.3 -14 0])
    set(gcf,'color','w')
    set(gcf,'position',[.1 1 14 4],'unit','inches')
    set(gcf,'position',[.1 1 14 4],'unit','inches')
    
end


function DA_LME_plot_square_panels(x,y,x_std,y_std,n,x_eri,y_eri,x_buoy,y_buoy,PP)
    
    [slope, inter, slope_member, inter_member] =  ...
                                 CDC_yorkfit_bt(y,x,y_std,x_std,0,1,PP.N,PP.P);  
                        
    [slope_member,I] = sort(slope_member);
    inter_member = inter_member(I);
    
    x1_pic = 0;
    x2_pic = 0.6;
    
    % Plot the shading for 2 s.d. 
    x_temp = x1_pic: 0.01:x2_pic;
    clear('yy')
    for ct = 1:numel(x_temp)
        yy(ct,:) = quantile(x_temp(ct) * slope_member + inter_member,[0.025 0.975]);
    end
    patch([x_temp,fliplr(x_temp)],[yy(:,1)' fliplr(yy(:,2)')],[1 1 1]*.6,'facealpha',0.3,'linest','none')

    if isnan(y_buoy),
        plot(x_buoy* [1 1],[-1 1],'-','linewi',3,'markersize',6,'color',PP.eri_color);
    end

    plot([x1_pic x2_pic],[x1_pic x2_pic]*slope+inter,'w-','linewi',4)
    plot([x1_pic x2_pic],[x1_pic x2_pic]*slope+inter,'m-','linewi',2)  
    
    for ct = 1:numel(x)
        plot(x(ct) + [-1 1]*x_std(ct)*PP.n_sigma , y(ct) * [1 1],'-','color',[1 1 1]*.3);
        plot(x(ct) * [1 1] , [-1 1]*y_std(ct)*PP.n_sigma + y(ct),'-','color',[1 1 1]*.3);
        if isempty(n),
            plot(x(ct),y(ct),PP.pic_st(ct),'color',[1 1 1]*.0,'markerfacecolor',...
                PP.pic_col(ct,:),'markersize',15,'linewi',2);
        else
            plot(x(ct),y(ct),PP.pic_st(ct),'color',[1 1 1]*.0,'markerfacecolor',...
                PP.pic_col(ct,:),'markersize',n(ct).^0.3/2,'linewi',2);
        end
    end
    
    plot(x_eri,y_eri,'o','linewi',3,'markersize',6,'color',PP.eri_color);
    plot(x_eri,y_eri,'o','linewi',3,'markersize',15,'color',PP.eri_color);  

    if ~isnan(y_buoy),
        plot(x_buoy,y_buoy,'bo','linewi',3,'markersize',6,'color',[.9 .1 0]);
        plot(x_buoy,y_buoy,'bo','linewi',3,'markersize',15,'color',[.9 .1 0]);
    end  
    
end