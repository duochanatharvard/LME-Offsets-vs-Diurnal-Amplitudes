% diurnal = DIURNAL_function_compute_morak_DA(mon, lat, hr, u, C, do_trd)
% 
% Compute expected buoy diurnal cycles from MB16.
% 
% Please make sure that inputs have the same dimensions
% do_trd: whether to compute a trend and form an 'open' cycle
% 
% A linear trend representing seasonal cycles should not be included 
% if input diurnal cycles has two nights 
% 
% Reference:
%    Morak?Bozzo, S., Merchant, C. J., Kent, E. C., Berry, D. I., & Carella, G. (2016). 
%    Climatological diurnal variability in sea surface temperature characterized 
%    from drifting buoy data. Geoscience Data Journal, 3(1), 20-28.

function diurnal = DIURNAL_function_compute_morak_DA(mon, lat, hr, u, C, do_trd)

    % *********************************************************************
    % Convert month into seasons
    % *********************************************************************
    sea = fix(mon/3) + 1;
    sea(sea == 5) = 1;
    
    % *********************************************************************
    % Convert latitude into bands
    % *********************************************************************
    lat_id = discretize(lat,[-90:10:90]);
    
    % *********************************************************************
    % Load fitting parameters
    % *********************************************************************
    [fit_coeffs,avg_u] = DIURNAL_function_morak_DA_coefficient;
    
    % *********************************************************************
    % Set parameters
    % *********************************************************************
    w = 2*pi/24;
    id = sub2ind([5,18],sea,lat_id);
    clear('ct')
    for ct = 1:14
        temp = squeeze(fit_coeffs(ct,:,:));
        c{ct} = temp(id);
    end
    u_avg = avg_u(id);
    
    % *********************************************************************
    % Compute for diurnal cycles
    % *********************************************************************
    if do_trd == 1,
        d1 = c{1} + c{2}.*u + c{3}.*C + c{4}.*hr + c{5}.*u.*hr + c{6}.*C.*hr;
    else
        d1 = c{1} + c{2}.*u + c{3}.*C;
    end
    d2 = c{7}.*sin(w*hr) + c{8}.*C.*sin(w*hr) + c{9}.*cos(w*hr) + c{10}.*C.*cos(w*hr);
    d3 = c{11}.*sin(2*w*hr) + c{12}.*C.*sin(2*w*hr) + c{13}.*cos(2*w*hr) + c{14}.*C.*cos(2*w*hr);
    diurnal  = d1 + (d2 + d3) .* exp(-u./u_avg);

end

% *************************************************************************
% FOR DEBUG ... 
% *************************************************************************
% clear;
% file = 'diurnal_sst_anomalies_from_drifting_buoys_v1.1.nc';
% sstano_fit = ncread(file,'sstano_fit');
% uu = ncread(file,'u');
% cc = ncread(file,'C');
% 
% % fit_coeffs = ncread(file,'fit_coeffs');
% % avg_u = squeeze(nanmean(nanmean(uu,3),1));
% [fit_coeffs,avg_u] = DIURNAL_function_morak_DA_coefficient;
% 
% cat = 1;
% lat = 5;
% sea = 1;
% b = sstano_fit(:,sea,cat,lat);
% c = squeeze(fit_coeffs(:,sea,lat));
% u  = uu(:,sea,cat,lat);
% C  = cc(:,sea,cat,lat);
% 
% u_avg = avg_u(sea,lat);
% w  = 2*pi/24;
% 
% t = [0.5:1:23.5]';
% d1 = c(1) + c(2)*u + c(3)*C + c(4)*t + c(5)*u.*t + c(6)*C.*t;
% d2 = c(7)*sin(w*t) + c(8)*C.*sin(w*t) + c(9)*cos(w*t) + c(10)*C.*cos(w*t);
% d3 = c(11)*sin(2*w*t) + c(12)*C.*sin(2*w*t) + c(13)*cos(2*w*t) + c(14)*C.*cos(2*w*t);
% d  = d1 + (d2 + d3) .* exp(-u./u_avg);
% clf; hold on; plot(b); plot(d)