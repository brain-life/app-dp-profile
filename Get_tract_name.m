function [tract_name] = Get_tract_name(tract_number)
switch tract_number
    case 1
        tract_name = 'Thal_Rad_L';
    case 2
        tract_name = 'Thal_Rad_R';  
    case 3
        tract_name = 'CST_L';
    case 4
        tract_name = 'CST_R';   
    case 5
        tract_name = 'Cing_L';
    case 6
        tract_name = 'Cing_R';     
    case 7
        tract_name = 'Hipp_L';
    case 8
        tract_name = 'Hipp_R';    
    case 9
        tract_name = 'Call_Maj';
    case 10
        tract_name = 'Call_Min';  
    case 11
        tract_name = 'IFOF_L';
    case 12
        tract_name = 'IFOF_R';  
    case 13
        tract_name = 'ILF_L';
    case 14
        tract_name = 'ILF_R';    
    case 15
        tract_name = 'SLF_L';
    case 16
        tract_name = 'SLF_R';  
    case 17
        tract_name = 'Unc_L';
    case 18
        tract_name = 'Unc_R';  
    case 19
        tract_name = 'ARC_L';
    case 20
        tract_name = 'ARC_R';           
    otherwise
      disp(['ERROR: ',tract_number, ' not found.'])
end    
end