clear;

dir_save  = 'Users/duochan/Data/LME_intercomparison/Step_04_LME_output/';
file_load = [dir_save,'LME_Ship_vs_Ship_all_measurements_1850_2014_Full_SST_Global.mat'];

load(file_load,'out_rnd')
load(file_load,'out')

% ************************************************************************
% Set Parameters  
% ************************************************************************
yint = 1;
do_sort = 1;
ysrt    = 1850;
alpha   = 0.05;
do_NpD  = 1;
figure(3);

clear('test_mean','test_random')
NY = size(out.bias_decade,1);
for i =  1:5
    out.bias_decade_annual(i:5:NY*5,:) = out.bias_decade;
    out_rnd.bias_decade_rnd_annual(i:5:NY*5,:,:) = out_rnd.bias_decade_rnd;
end
out.bias_decade = out.bias_decade_annual;
out_rnd.bias_decade_rnd = out_rnd.bias_decade_rnd_annual;
out.bias_decade(end+1:165,:) = nan;
out_rnd.bias_decade_rnd(end+1:165,:,:) = nan;

N_group = size(out.unique_grp,1);
N_rnd   = size(out_rnd.bias_fixed_random,1);

fixed_mean   = repmat(out.bias_fixed,1,165)';
test = fixed_mean + out.bias_decade;

% *******************************************************
% Find nations that are to be marked by '*' or '**'    **
% *******************************************************
if yint == 1
    temp = out_rnd.bias_fixed_random;
    temp  = temp - repmat(nanmean(temp,1),size(temp,1),1);
    out.bias_fixed_std = sqrt(nansum(temp.^2,1) / (size(temp,1) - 1))';
end
Z_stats = (1 - normcdf(abs(out.bias_fixed ./ out.bias_fixed_std)))*2;
Sig_Nat_90 = Z_stats < alpha;
Sig_Nat_BF = Z_stats < alpha/46;
List_Nat = out.unique_grp;

% ***************************************************************************
% Add a section that reads the number of measurements from individual groups
% ***************************************************************************
dir_home = LME_OI('home');
file_stats = [dir_home,'Stats_All_ships.mat'];
load(file_stats,'unique_grp','Stats_glb');
Stats_glb = squeeze(nansum(Stats_glb,2));

% *******************************************************
% Prepare Data to be plotted                           **
% *******************************************************
if do_sort == 1
    [~,I] = sort(out.bias_fixed);
    test = test(:,I);
    Stats_glb = Stats_glb(:,I);
    Sig_Nat_90 = Sig_Nat_90(I);
    Sig_Nat_BF = Sig_Nat_BF(I);
    List_Nat   = List_Nat(I,:);
    
    l = nansum(~isnan(test),1) > 1;
    test = test(:,l);
    Stats_glb = Stats_glb(:,l);
    Sig_Nat_90 = Sig_Nat_90(l);
    Sig_Nat_BF = Sig_Nat_BF(l);
    List_Nat   = List_Nat(l,:);
end

% *******************************************************
% Find the color scheme                                **
% *******************************************************
col = colormap_CD([ .5 .67; .05 0.93],[.9 .2],[0 0],10);
cc  = discretize(test',[-inf -0.45:0.05:0.45 inf]);
wid = (log10(Stats_glb')+1)/15;

% ********************************************
% A list of decks that appear from 1908-1941
% ********************************************
l = any(~isnan(out.bias_decade([1908:1941]-1849,:)),1);
deck_bf = out.unique_grp(l,:);

%% *******************************************************
% Generate Figures                                     **
% *******************************************************
figure(3);clf;hold on;
num_row = 20;
num_col = 6;
pic = test';

for ct = 1:num_col
    sty_in{ct} = [1 num_row-2 ct ct];
end
sty_in{num_col+1} = [num_row num_row 1 num_col-2];
sty_in{num_col+2} = [num_row num_row num_col-1 num_col];
out_ly = CDF_layout([num_row,num_col],sty_in);

% for ct = 1:num_col
%     subplot(num_row,num_col,out_ly{ct}), hold on
%     h = patch([1908 1941 1941 1908],[-160 -160 0 0],[1 1 1]*.7);
%     set(h,'EdgeColor','none','facealpha',0.5);
% end

clear('y_label_text')
for i = 1:size(List_Nat,1)
    [p1,p2] = ind2sub([ceil(size(test,2)/num_col), num_col],i);
    
    subplot(num_row,num_col,out_ly{p2}),
    for j = 1:size(pic,2)
        if ~isnan(cc(i,j))
            patch(ysrt + [0 1 1 0]*yint + yint * (j-1),[-1 -1 1 1]*wid(i,j) - p1,...
                col(cc(i,j),:),'linest','none');
        end
    end
    
    if Sig_Nat_BF(i) && ismember(List_Nat(i,:),deck_bf,'rows')
        surfix = '** ';
    elseif Sig_Nat_90(i)
        surfix = '* ';
    else
        surfix = '';
    end
    
    if do_NpD == 1
        clear('surfix_col')
        switch List_Nat(i,end)
            case 0
                surfix_col = '\color[rgb]{0,0,1}';
            case 1
                surfix_col = '\color[rgb]{1,0,0}';
            case 13
                surfix_col = '\color[rgb]{.5,.5,1}';
            case 14
                surfix_col = '\color[rgb]{1,.5,.5}';
            case -1
                surfix_col = '\color[rgb]{0.5,0.5,0.5}';
            case 3
                surfix_col = '\color[rgb]{1,.6,0}';
        end
    end

    if all(List_Nat(i,1) > 100)
        y_label_text{i,p2} = [surfix_col, surfix,'---- D ',num2str(List_Nat(i,1))];
    else
        y_label_text{i,p2} = [surfix_col, surfix,char(List_Nat(i,1:2)),' D ',num2str(List_Nat(i,3))];
    end
end
%%
NN = ceil(size(test,2)/num_col);
if yint == 1
    for i = 1:num_col
        subplot(num_row,num_col,out_ly{i}),
        lgd_text = flipud(y_label_text((i-1)*NN+1:min(i*NN,size(List_Nat,1)),i));
        CDF_panel([1850 2015 -ceil(size(test,2)/num_col)-1 -0],[],{},'Year','');
        set(gca,'ytick',[-numel(lgd_text):1:-1],'yticklabel',lgd_text,'fontsize',10);
        set(gca,'xtick',1850:10:2010,'xticklabel',{'1850','','','','','1900','','','','','1950','','','','','2000',''});
        set(gca,'fontsize',10,'fontweight','bold');   
    end
end

figure(3);  subplot(num_row,num_col,out_ly{end-1}),
for i = 1:size(col,1)
    patch([0 1 1 0]+i - 1,[0 0 1 1],col(i,:),'linewi',1);
end
set(gca,'xtick',0:2:20,'xticklabel',[-0.5:0.1:0.5]);
set(gca,'ytick',[]);
set(gca,'fontsize',16,'fontweight','normal')
xlabel('SST offsets between groups (^oC)');
set(gcf,'color','w')

set(gcf,'position',[1 1 15 8]*1.2,'unit','inches');
set(gcf,'position',[-1 7.5 15 8]*1.7,'unit','inches');

%%
dir_save = '/Users/duochan/Dropbox/Research/01_SST_Bucket_Intercomparison/02_Manuscript/04_ERI_Hull_2020/';
dir_save = '/Users/duochan/Desktop/';
file_save = [dir_save,'LME_all_Ships.png'];
CDF_save(3,'png',400,file_save);