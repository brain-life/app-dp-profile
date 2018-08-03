function [tract_number] = Get_tract_number(tract_name)
switch tract_name
    case 'Thal_Rad_L'
        tract_number = 1;
    case 'Thal_Rad_R'
        tract_number = 2;  
    case 'CST_L'
        tract_number = 3;
    case 'CST_R'
        tract_number = 4;   
    case 'Cing_L'
        tract_number = 5;
    case 'Cing_R'
        tract_number = 6;     
    case 'Hipp_L'
        tract_number = 7;
    case 'Hipp_R'
        tract_number = 8;    
    case 'Call_Maj'
        tract_number = 9;
    case 'Call_Min'
        tract_number = 10;  
    case 'IFOF_L'
        tract_number = 11;
    case 'IFOF_R'
        tract_number = 12;  
    case 'ILF_L'
        tract_number = 13;
    case 'ILF_R'
        tract_number = 14;    
    case 'SLF_L'
        tract_number = 15;
    case 'SLF_R'
        tract_number = 16;  
    case 'Unc_L'
        tract_number = 17;
    case 'Unc_R'
        tract_number = 18;  
    case 'ARC_L'
        tract_number = 19;
    case 'ARC_R'
        tract_number = 20;           
    otherwise
      disp(['ERROR: ',tract_name, ' not found.'])
end    
end