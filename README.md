# FIS_Analysis_3D
Matlab code for computing Facial Improvement Score (FIS)

Function that generates random groups (with different ammount of treated subjects) and perform the Facial Improvement Score (FIS) on them.                              

As output a FIS histogram is calculated and saved as figure.

% Input:                                                        
%   res_path: Result path where the histigrams will be stored.  
%   type_sample: 0 for humans , 1 for mice.                     
%   age_range: [min max] vector defining the range of ages.     
%   nameOutput: name for the stored histogram.                  
%   treat_level: Treatment dosis level; 1 for low, 2 for high.  
%                                                               
% Notes:                                                        
%   All humans have an predefined treatment dosis level of 1.   
%   All mice have an predefined age of 1.                       
%                                                               
% Example:                                                      
%   FIS_analysis_3D('results\',1,[0 2],'mice_low',1)            
%                                                               
%%%%%%%%%%%%%%%

This code is available as Supplementary Material to the article 
"Epigallocatechin-3-Gallate Stimulates Face Remodeling in Down Syndrome" 
authored by: John M. Starbuck, Sergi Llambrich, Rubèn González, Julia Albaigès, 
Anna Sarlé, Jens Wouters, Alejandro González, Xavier Sevillano, James Sharpe, 
Rafael De La Torre, Mara Dierssen, Greetje Vande Velde, and Neus Martínez-Abadías
