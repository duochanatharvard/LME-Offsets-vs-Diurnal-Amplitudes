 function [output] = RGB_generator_RGB(input,mode)

    if  nargin < 2
        mode = 1;
    end

    Hue_temp = input(:,1);
    Str_temp = input(:,2);
    Brt_temp = input(:,3);

    Hue_temp = rem(Hue_temp + 1000,1);

    Hue_RGB = nan(size(input));
    Str_RGB = nan(size(input));
    if mode == 1;
        Brt_RGB = nan(size(input));
    else
        Gry_RGB = nan(size(input));
    end

    logic = Hue_temp>=0 & Hue_temp<1/6;
    Hue_RGB(logic,:)=repmat([1 0 0],nnz(logic),1) + repmat((Hue_temp(logic,:) - 0/6),1,3) .* 6 .* repmat([0 1 0],nnz(logic),1);

    logic = Hue_temp>=1/6 & Hue_temp<2/6;
    Hue_RGB(logic,:)=repmat([1 1 0],nnz(logic),1) - repmat((Hue_temp(logic,:) - 1/6),1,3) .* 6 .* repmat([1 0 0],nnz(logic),1);

    logic = Hue_temp>=2/6 & Hue_temp<3/6;
    Hue_RGB(logic,:)=repmat([0 1 0],nnz(logic),1) + repmat((Hue_temp(logic,:) - 2/6),1,3) .* 6 .* repmat([0 0 1],nnz(logic),1);

    logic = Hue_temp>=3/6 & Hue_temp<4/6;
    Hue_RGB(logic,:)=repmat([0 1 1],nnz(logic),1) - repmat((Hue_temp(logic,:) - 3/6),1,3) .* 6 .* repmat([0 1 0],nnz(logic),1);

    logic = Hue_temp>=4/6 & Hue_temp<5/6;
    Hue_RGB(logic,:)=repmat([0 0 1],nnz(logic),1) + repmat((Hue_temp(logic,:) - 4/6),1,3) .* 6 .* repmat([1 0 0],nnz(logic),1);

    logic = Hue_temp>=5/6 & Hue_temp<=6/6;
    Hue_RGB(logic,:)=repmat([1 0 1],nnz(logic),1) - repmat((Hue_temp(logic,:) - 5/6),1,3) .* 6 .* repmat([0 0 1],nnz(logic),1);

    Str_RGB=(Hue_RGB - .5) .* repmat(Str_temp,1,3) + .5;
    
    if(mode == 1)
        
        logic = Brt_temp > 0.5;
        Brt_RGB(logic,:) = Str_RGB(logic,:) + (1 - Str_RGB(logic,:)) .* repmat(((Brt_temp(logic,:)-0.5)./0.5),1,3);
        
        logic = Brt_temp <= 0.5;
        Brt_RGB(logic,:) = Str_RGB(logic,:) - (Str_RGB(logic,:)) .* repmat(((0.5 - Brt_temp(logic,:))./0.5),1,3);
        
        output = Brt_RGB;
        
    else
    
        Gry_temp = Brt_temp;
        
        xishu1 = 0.2989;
        xishu2 = 0.5870;
        xishu3 = 0.1140;
        
        Gry_input = sqrt(xishu1 * Str_RGB(:,1) + xishu2 * Str_RGB(:,2) + xishu3 * Str_RGB(:,3));
        
        logic = Gry_input > Gry_temp;
        
        Gry_insuf = Gry_input(logic);
        Gry_RGB(logic,:) = Str_RGB(logic,:) .* repmat(Gry_temp(logic)./Gry_insuf,1,3);
        
        logic = Gry_input < Gry_temp;
        Gry_insuf = Gry_input(logic);
        kk = repmat((Gry_temp(logic)-1)./(Gry_insuf-1),1,3);
        Gry_RGB(logic,:) = (Str_RGB(logic,:)-1).*kk + 1;
        
        logic = Gry_input == Gry_temp;
        Gry_RGB(logic,:) = Str_RGB(logic,:);

        output = Gry_RGB;
    end
    
end