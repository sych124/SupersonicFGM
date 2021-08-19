function [E, alph, rho, cond, poiss] = MatTemp(Mat)
if strcmp(Mat,'SiN')
    E = [0, 348.43e9 , -3.070e-4 , 2.160e-7 ,-8.946e-11]; % Si3N4
    alph = [0, 5.8723e-6, 9.095e-4, 0, 0] ;
    cond = [0, 13.723, -1.032e-3, 5.466e-7, -7.876e-11];
    rho = [0, 2370, 0 ,0, 0];
    poiss = [0, 0.24, 0,0,0];
end
if strcmp(Mat,'SUS')
    E = [0, 201.04e9, 3.079e-4, -6.534e-7, 0]; % SUS304
    alph = [0, 12.33e-6, 8.086e-4, 0, 0 ];
    cond = [0, 15.379, -1.264e-3, 2.092e-6, -7.223e-10];
    rho = [0, 8166, 0, 0 ,0];  
    poiss = [0, 0.3262, -2.002e-4, 3.797e-7,0];
end
if strcmp(Mat,'AlO')
    E = [0,   349.55e9, -3.853e-4, 4.027e-7, -1.673e-10];
    alph = [0, 6.8269e-6, 1.838e-4, 0, 0 ] ;
    cond = [-1123.6, -14.087, -6.227e-3, 0, 0];
    rho = [0, 3800, 0 ,0, 0];
    poiss = [0, 0.26, 0,0,0]; % Al2O3 properties
end
if strcmp(Mat,'Al')
    E = [0, 20.179e9, 1.73e-2, -3.4698e-5, 1.4808e-8];
    alph = [0, 23e-6, 0, 0, 0 ]; 
    cond = [0, 204, 0, 0, 0 ]; 
    rho = [0, 2700, 0, 0 ,0];  
    poiss = [0, 0.3, 0,0,0];% Al properties
end
if strcmp(Mat,'Zr')
    E = [0, 132.2e9,-3.805e-4, -6.127e-14, 0]; %Check!!!!!
    alph = [0, 0, 0, 0, 0 ];     
    cond = [0, 1.710, 1.228e-4, 6.884e-8, 0 ]; %Check!!!!    
    rho = [0,0, 0, 0 ,0];  
    poiss = [0, 0.2882, 1.133e-4, 0, 0];% ZrO2 properties
end
if strcmp(Mat,'Ti')
    E = [0, 122.7e9, -4.605e-4, 0, 0]; %Check!
    alph = [0, 7.5788e-6, 6.638e-4, -3.147e-6, 0 ];          
    cond = [0, 1, 1.704e-2, 0, 0 ];    
    rho = [0, 0, 0, 0 ,0];  
    poiss = [0,	0.2884,	1.121e-4, 0, 0];% TiAlV properties
end

if strcmp(Mat,'SiC')
    E = [0, 122.7e9, -4.605e-4, 0, 0]; %Check!
    alph = [0, 7.5788e-6, 6.638e-4, -3.147e-6, 0 ];          
    cond = [0, 1, 1.704e-2, 0, 0 ];    
    rho = [0, 0, 0, 0 ,0];  
    poiss = [0,	0.2884,	1.121e-4, 0, 0];% TiAlV properties
end

end