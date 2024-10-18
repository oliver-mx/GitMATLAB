function ev(a,I)
    % a is output1
    % I is vector with indices
if nargin == 0
clc;
disp('Output argument overview:')
disp('-------------------------------------')
disp('1  -- SEC_net    | 15 -- J_d^RO(0)')
disp('2  -- FW         | 16 -- C_d^RO(0)')
disp('3  -- Rev        | 17 -- J_d^RO(L)')
disp('4  -- REC_RO     | 18 -- C_d^RO(L)')
disp('5  -- REC_PRO    | 19 -- J_f^RO(L)')
disp('6  -- PD_net     | 20 -- C_f^RO(L)')
disp('-------------------------------------')
disp('7  -- mix_M1     | 21 -- Q_d^PRO(in)')
disp('8  -- mix_M2     | 22 -- c_d^PRO(in)')
disp('9  -- W_net      | 23 -- Q_d^PRO(out)')
disp('10 -- W_p1       | 24 -- c_d^PRO(out)')
disp('11 -- W_p2       | 25 -- Q_f^PRO(0)')
disp('12 -- W_p3       | 26 -- c_f^PRO(0)')
disp('13 -- W_p4       | 27 -- Q_f^PRO(L)')
disp('14 -- W_t        | 28 -- c_f^PRO(L)')
disp('-------------------------------------')
else
    if nargin == 1
        I=1:1:length(a);
    end
    for index=1:length(I)
        i=I(index);
        if I(index)> length(a)
        else
            %--------------------------------------------------------------
            if i == 1;str1='SEC_net = ';str2=' [kWh/m^3]';end
            if i == 2;str1='FW      =  ';str2=' [m^3/h]';end
            if i == 3;str1='Rev     =  ';str2=' [$/h]';end
            if i == 4;str1='REC_RO  =  ';str2=' [%]';end
            if i == 5;str1='REC_PRO =  ';str2=' [%]';end
            %--------------------------------------------------------------
            if i == 6;str1='PD_net  =  ';str2=' [W/m^2]';end
            %--------------------------------------------------------------
            if i == 7;str1='mix_M1  =  ';str2=' ';end
            if i == 8;str1='mix_M2  =  ';str2=' ';end
            %--------------------------------------------------------------
            if i == 9;str1='W_net = ';str2=' [kWh]';end
            if i == 10;str1='W_p1  = ';str2=' [kWh]';end
            if i == 11;str1='W_p2  = ';str2=' [kWh]';end
            if i == 12;str1='W_p3  = ';str2=' [kWh]';end
            if i == 13;str1='W_p4  = ';str2=' [kWh]';end
            if i == 14;str1='W_t   = ';str2=' [kWh]';end
            if i == 15;str1='J_d^RO(0) = ';str2=' [kg/sm]';end
            if i == 16;str1='C_d^RO(0) = ';str2=' [%]';end
            if i == 17;str1='J_d^RO(L) = ';str2=' [kg/sm]';end
            if i == 18;str1='C_d^RO(L) = ';str2=' [%]';end
            if i == 19;str1='J_f^RO(L) = ';str2=' [kg/sm]';end
            if i == 20;str1='C_f^RO(L) = ';str2=' [ppm]';end
            %--------------------------------------------------------------
            if i == 21;str1='Q_d^PRO(in)  = ';str2=' [kg/sm]';end
            if i == 22;str1='c_d^PRO(in)  = ';str2=' [%]';end
            if i == 23;str1='Q_d^PRO(out) = ';str2=' [kg/sm]';end
            if i == 24;str1='c_d^PRO(out) = ';str2=' [%]';end
            if i == 25;str1='Q_f^PRO(0)   = ';str2=' [kg/sm]';end
            if i == 26;str1='c_f^PRO(0)   = ';str2=' [ppm]';end
            if i == 27;str1='Q_f^PRO(L)   = ';str2=' [kg/sm]';end
            if i == 28;str1='c_f^PRO(L)   = ';str2=' [ppm]';end
            %--------------------------------------------------------------
            disp(['% ',str1, num2str(a(i)), str2])
        end
    end
end
end