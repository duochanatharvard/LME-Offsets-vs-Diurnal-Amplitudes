function DA_LME_function_scatter_plot(method,yr_start,reg_sea,revision)

    if strcmp(method,'bucket')
        if revision == 0
            data = load('All_lme_offsets_and_diurnal_amplitudes.mat');                    
            off_set = revision;
        elseif revision == 1
            data = load('All_lme_offsets_and_diurnal_amplitudes_sens_wmo.mat');  
            off_set = revision;
        elseif revision == 2
            data = load('All_lme_offsets_and_diurnal_amplitudes_sens_not_infer.mat');  
            off_set = revision;
        elseif revision == 3
            data = load('All_lme_offsets_and_diurnal_amplitudes_sens_nation_level.mat');  
            off_set = revision;
        elseif revision == 4
            data = load('All_lme_offsets_and_diurnal_amplitudes.mat');
            data_clim = load('All_lme_offsets_and_diurnal_amplitudes_clim_diurnal.mat');  
            off_set = revision;
            data.da = data.da - data_clim.da;
            data.da_std = sqrt(data.da_std.^2 + data_clim.da_std.^2);
        end
        num_of_colors        = 13;
        num_of_marker_shapes = 6;
    else
        switch method
            case 1
                data = load('All_lme_offsets_and_diurnal_amplitudes_ERI_do_bucket_0.mat');
            case 2
                data = load('All_lme_offsets_and_diurnal_amplitudes_ERI_do_bucket_1.mat');
            case 3
                data = load('All_lme_offsets_and_diurnal_amplitudes_ERIex_do_bucket_0.mat');
            case 4
                data = load('All_lme_offsets_and_diurnal_amplitudes_ERIex_do_bucket_1.mat');
            case 5
                data = load('All_lme_offsets_and_diurnal_amplitudes_ERIex_do_bucket_1_5000_pairs.mat');
        end

        num_of_colors        = 20;
        num_of_marker_shapes = 6;
    end


    [out_st,out_col]     = CDF_linest(num_of_colors * num_of_marker_shapes, num_of_marker_shapes);
    
    % *********************************************************************
    % Set parameters
    % *********************************************************************
    N       = 100;              % Number of bootstrapping members
    n_sigma = 2;                % Number of s.d. for error bars
    if strcmp(method,'bucket')
        yr0     = 1879;             % The first year of data
    else
        yr0     = 1929;             % The first year of data 
    end
    
    if reg_sea == 1
        PP.eri_color = 'k';
    elseif ismember(reg_sea,[2 4])
        PP.eri_color = 'b';
    elseif ismember(reg_sea,[3 5])
        PP.eri_color = 'r';
    else
        error('Not a valid region-season ID')
    end
       
    % *********************************************************************
    % Prepare for data to be plotted
    % *********************************************************************
    if reg_sea ~= 1
        l_legend = ~isnan(nanmean(data.da(2:end,yr_start-yr0,2:5),3));
    else
        l_legend = ~isnan(nanmean(data.da(2:end,[1890:20:1990]-yr0,1),2));
    end
    
    x     = data.da(2:end,yr_start-yr0,reg_sea);
    y     = data.lme(2:end,yr_start-yr0,reg_sea);
    x_std = data.da_std(2:end,yr_start-yr0,reg_sea);
    y_std = data.lme_std(2:end,yr_start-yr0,reg_sea);
    n     = [];
    x_eri = data.da(1,yr_start-yr0,reg_sea);
    y_eri = data.lme(1,yr_start-yr0,reg_sea);
    x_buoy = data.Buoy(yr_start-yr0,reg_sea);
    y_buoy = nan;
    
    l = ~isnan(x);
    x = x(l);
    y = y(l);
    x_std = x_std(l);
    y_std = y_std(l);

    PP.N             = N;
    PP.n_sigma       = n_sigma;
    PP.P.mute_output = 1;
    PP.pic_st        = out_st(l);
    PP.pic_col       = out_col(l,:);
    temp             = data.grp(2:end,:);
    PP.grp           = temp(l,:);
    if strcmp(method,'bucket')
        PP.grp(:,4)           = 0;
    end
   
    % *********************************************************************
    % Generate plot
    % *********************************************************************
    figure((yr_start-yr0-1)*4+reg_sea + off_set); clf; hold on;
    DA_LME_plot_square_panels(x,y,x_std,y_std,n,x_eri,y_eri,x_buoy,y_buoy,PP)
    if strcmp(method,'bucket')
        CDF_panel([0 0.6 -0.6 0.6],'','','Diurnal Amplitude (^oC)','Offsets (^oC)','fontsize',18);
        daspect([0.6 1.2 1])  
    else
        CDF_panel([0 0.4 -0.6 0.8],'','','Diurnal Amplitude (^oC)','Offsets (^oC)','fontsize',18);
        title([num2str(yr_start),' - ',num2str(yr_start+19)])
        daspect([0.4 1.4 1])  
    end
    set(gca,'xtick',[0.1:0.1:0.5],'ytick',[-.4:.2:.4])

    set(gcf,'position',[.1 1 7 7],'unit','inches')
    set(gcf,'position',[.1 1 7 7],'unit','inches')


    % *********************************************************************
    % Generate legend
    % *********************************************************************
    group = data.grp(2:end,:); 
    in_name = {};
    for ct = 1:size(group,1)
        name = double(group(ct,:));
        if numel(name) > 2
            if name(1) > 100
                in_name{ct} = ['      DCK ',num2str(name(3))];
            else
                in_name{ct} = [name(1:2),' DCK ',num2str(name(3))];
            end
        else
            if name(1) > 100
                in_name{ct} = ['DCK ',num2str(name(1))];
            else
                in_name{ct} = [char(name(1:2))];
            end
        end
    end
    
    figure(1000+off_set); clf; hold on;
    if reg_sea == 1 
        if revision ~= 3 % == 3
            CDF_scatter_legend(in_name(l_legend),out_st(l_legend),out_col(l_legend,:),[1,7,1],'fontsize',15,'mksize',13)
            axis([1 7.3 -ceil(nnz(l_legend)/6) 0])
            set(gcf,'color','w')
            set(gcf,'position',[.1 1 14 4],'unit','inches')
            set(gcf,'position',[.1 1 14 4/ceil(66/6)*ceil(size(group,1)/6)],'unit','inches')
        else
            CDF_scatter_legend(in_name(l_legend),out_st(l_legend),out_col(l_legend,:),[1,2,1],'fontsize',15,'mksize',13)
            axis([1 3.3 -ceil(nnz(l_legend)/2) 0])
            set(gcf,'color','w')
            set(gcf,'position',[.1 1 14 4],'unit','inches')
            set(gcf,'position',[.1 1 5 nnz(l_legend)/4],'unit','inches') 
        end
    else
        CDF_scatter_legend(in_name(l_legend),out_st(l_legend),out_col(l_legend,:),[1,1,1],'fontsize',15,'mksize',13)
        axis([1 1.5 -ceil(nnz(l_legend)) 0])
        set(gcf,'color','w')
        set(gcf,'position',[.1 1 14 4],'unit','inches')
        set(gcf,'position',[.1 1 3 nnz(l_legend)/2.8],'unit','inches')
    end
   
    
end


function DA_LME_plot_square_panels(x,y,x_std,y_std,n,x_eri,y_eri,x_buoy,y_buoy,PP)
    
    [slope, inter, slope_member, inter_member] =  ...
                                 CDC_yorkfit_bt(y,x,y_std,x_std,0,1,PP.N,PP.P);  
                        
    [slope_member,I] = sort(slope_member);
    inter_member = inter_member(I);
    
    x1_pic = 0;
    x2_pic = 0.6;
    
    if 0    % To plot the residual offsets
        y = y - slope * x - inter;  
        slope_member = slope_member - slope;
        slope = 0;
        inter_member = inter_member - inter;
        inter = 0;
    end
    
    % Plot the shading for 2 s.d. 
    x_temp = x1_pic: 0.01:x2_pic;
    clear('yy')
    for ct = 1:numel(x_temp)
        yy(ct,:) = quantile(x_temp(ct) * slope_member + inter_member,[0.025 0.975]);
    end
    patch([x_temp,fliplr(x_temp)],[yy(:,1)' fliplr(yy(:,2)')],[1 1 1]*.6,'facealpha',0.3,'linest','none')

    if isnan(y_buoy)
        plot(x_buoy* [1 1],[-1 1],'-','linewi',3,'markersize',6,'color',PP.eri_color);
    end

    plot([x1_pic x2_pic],[x1_pic x2_pic]*slope+inter,'w-','linewi',4)
    plot([x1_pic x2_pic],[x1_pic x2_pic]*slope+inter,'m-','linewi',2)  
    
    for ct = 1:numel(x)
        plot(x(ct) + [-1 1]*x_std(ct)*PP.n_sigma , y(ct) * [1 1],'-','color',[1 1 1]*.3);
        plot(x(ct) * [1 1] , [-1 1]*y_std(ct)*PP.n_sigma + y(ct),'-','color',[1 1 1]*.3);
        if isempty(n)
            if ismember(PP.grp(ct,4),[0 1])
                edge_color = [0 0 0];
            else
                edge_color = [0 0 1];
            end
            plot(x(ct),y(ct),PP.pic_st(ct),'color',edge_color,'markerfacecolor',...
                PP.pic_col(ct,:),'markersize',15,'linewi',2);
        else
            if ismember(PP.grp(ct,4),[0 1])
                edge_color = [0 0 0];
            else
                edge_color = [0 0 1];
            end
            plot(x(ct),y(ct),PP.pic_st(ct),'color',edge_color,'markerfacecolor',...
                PP.pic_col(ct,:),'markersize',n(ct).^0.3/2,'linewi',2);
        end
    end
    
    % plot(x_eri,y_eri,'o','linewi',3,'markersize',6,'color',PP.eri_color);
    % plot(x_eri,y_eri,'o','linewi',3,'markersize',15,'color',PP.eri_color);  

    if ~isnan(y_buoy)
        plot(x_buoy,y_buoy,'bo','linewi',3,'markersize',6,'color',[.9 .1 0]);
        plot(x_buoy,y_buoy,'bo','linewi',3,'markersize',15,'color',[.9 .1 0]);
    end  
    
end