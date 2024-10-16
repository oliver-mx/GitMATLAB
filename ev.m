function ev(a,I)
    % a is output1
    % I is vector with indices
    if nargin == 1
        I=1:1:length(a);
    end
    for index=1:length(I)
        i=I(index);
        if I(index)> length(a)
        else
            if i == 1;str1='SEC_net = ';str2=' [kWh/m^3]';end
            if i == 2;str1='FW      =  ';str2=' [m^3/h]';end
            if i == 3;str1='Rev     =  ';str2=' [$/h]';end
            if i == 4;str1='REC_RO  =  ';str2=' [%]';end
            if i == 5;str1='REC_PRO =  ';str2=' [%]';end
                if isnan(a(4))
                    if i == 6;str1='PD_net  =  ';str2=' [W^2/m^2]';end
                    if i == 7;str1='P_d_out =   ';str2=' [bar]';end
                else
                    if i == 6;str1='RO_draw_in    = ';str2=' [m^3/s]';end
                    if i == 7;str1='Permeate_out  = ';str2=' [m^3/s]';end
                end
            if i == 8;str1='Wastewater_in = ';str2=' [m^3/s]';end
            if i == 9;str1='C_permeate    = ';str2=' [ppm]';end
            if i == 10;str1='C_brine       = ';str2=' [%]';end
            if i == 11; str1='C_dilluted    = ';str2=' [%]';end
            if i == 12;str1='W_net = ';str2=' [kWh]';end
            if i == 13;str1='W_p1  = ';str2=' [kWh]';end
            if i == 14;str1='W_p2  = ';str2=' [kWh]';end
            if i == 15;str1='W_p3  = ';str2=' [kWh]';end
            if i == 16;str1='W_p4  = ';str2=' [kWh]';end
            if i == 17;str1='W_t   = ';str2=' [kWh]';end
            if i == 18;str1='mix_M1  =  ';str2=' ';end
            if i == 19;str1='F_1-F_2 =  ';str2=' [kg m^2/s^3]';end % [kg m^2/s^3] = [W]
            if i == 20;str1='mix_M2  =  ';str2=' ';end
            if i == 21;str1='f_1-f_2 =  ';str2=' [kg m^2/s^3]';end
            disp([str1, num2str(a(i)), str2])
        end
    end
end