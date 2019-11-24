% RGB_out = CDF_colormap(col,brt,gry,num,random)
% generate the colorbar given certain parameters

function RGB_out = CDF_colormap(col,brt,gry,num,random)

    if(size(col,2)==1)
        col = [col col];
    end
    
    if(size(brt,2)==1)
        brt = [brt brt];
    end
    
    if(size(gry,2)==1)
        gry = [gry gry];
    end
    
    RGB = nan(num(1),3,size(col,1));
    
    for ll = 1:size(col,1)
        
        for i = 1:num(1)
            
            if(size(col,2)==1)
                Hue_temp = col(ll,1)+randn(1)*random/10;
            else
                Hue_temp = col(ll,1)+(i-1)*(col(ll,2)-col(ll,1))/(num-1)+randn(1)*random/10;
            end
            
            Str_temp = gry(1,1)+(i-1)*(gry(1,2)-gry(1,1))/(num-1)+randn(1)*random/10;
            
            Brt_temp = brt(1,1)+(i-1)*(brt(1,2)-brt(1,1))/(num-1)+randn(1)*random/10;
            
            RGB(i,:,ll)=RGB_generator_RGB([Hue_temp Str_temp Brt_temp],1);
        end
        
    end
    
    if(size(col,1) == 1)
        RGB_out = RGB;
    else
        RGB_out = [flipud(RGB(:,:,2)); RGB(:,:,1)];      
    end
    
    colormap(gca,RGB_out);

end