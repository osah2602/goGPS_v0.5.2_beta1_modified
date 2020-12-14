
%**************************************************************************+
%*****DESCRIPTION:                                                         + 
%                 "readGPTgrid" is a routine that searches and reads ...  +
%                 the various Global Pressure and Temperature(GPT) model...+ 
%                 series(GPT2,GPT2w,GPT3) grid files in the current ...    + 
%                 Directory/folder and  and over gives the respective ...  +
%                 parameters in a cell.A situation where Grid files are not+   
%                 found in the current folder,it searches recursively ...  +  
%                 through all Sub-Directories/Folders of the current.      +

%USAGE:                                                                    +
%      [gpt_grid] = readGPTgrid(GRIDfile,GPTmodel,grid_res)               +
%=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-+
%INPUT:                                                                    +
%1.    GRIDfile : GPT model grid file in single quote (e.g:'gpt2_1w.grd',  +
%                 'gpt2_5w.grd','gpt2_5.grd','gpt3_1.grd','gpt3_5.grd')    +
%      NOTE:                                                               +             
%            Gridfiles can also be extracted or imported into matlab 'mat' +
%            format/file with the extension '.mat'(e.g:'gpt2_1w.mat',      +
%                 'gpt2_5w.mat','gpt2_5.mat','gpt3_1.mat','gpt3_5.mat')    +

%2.   GPTmodel : GPT model type. indicate by ('GPT2' or 'GPT2w' or 'GPT3') +
%3.   grid_res : Grid resolution. 1 for 1x1 and 5 for 5x5 resolutions resp.+

%OUTPUT:                                                                   +                                                                   +
%1.     gpt_grid  -->    Gridfile parameters in a cell                     +
%=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-+
%=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-+
%WRITTEN BY: OSAH SAMUEL, MSC GEOMATIC ENGINEERING                         +
%     Email:osahsamuel@yahoo.ca/osahsamuel@gmail.com                       +
%     Phone:+233(0)246137410 / +233(0)509438484                            +
%=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-+
%**************************************************************************+
%**************************************************************************+

function [gpt_grid] = readGPTgrid(GRIDfile,GPTmodel,grid_res)

%***********SEARCH FOR GRIDFILE 
%Call the "searchGRID.m" function
[gpt_fn, gptfile]= searchGRID(GRIDfile);

%CHECK IF EITHER OF searchGRID OUTPUT IS EMPTY([])
if any([isempty(gpt_fn) isempty(gptfile)])
    
   %*****CHECK WHICH FILE TYPE WAS SEARCHED FOR
   %1 '.MAT FILE '
   if regexpi(GRIDfile, '\w*.mat$')
       
     %******CHECK MODEL TYPE
     if strcmpi(GPTmodel,'GPT2')
         
        %*****SEARCH FOR '.grd' file
        %Call the "searchGRID.m" function
        [gpt_fn, gptfile]= searchGRID('gpt2_5.grd');
        
     elseif strcmpi(GPTmodel,'GPT2w')
            %*****SEARCH FOR '.grd' file
            
            %CHECK GRID RESOLUTION 
            if grid_res == 1
               %Call the "searchGRID.m" function
               [gpt_fn, gptfile]= searchGRID('gpt2_1w.grd');
            
               if any([isempty(gpt_fn) isempty(gptfile)])
                  %Call the "searchGRID.m" function
                  [gpt_fn, gptfile]= searchGRID('gpt2_1wA.grd');
               end
               
            elseif grid_res == 5
                   %Call the "searchGRID.m" function
                   [gpt_fn, gptfile]= searchGRID('gpt2_5w.grd');                  
            end
            
     elseif  strcmpi(GPTmodel,'GPT3')
             %*******SEARCH FOR '.grd' file
            
            %CHECK GRID RESOLUTION 
            if grid_res == 1
               %Call the "searchGRID.m" function
               [gpt_fn, gptfile]= searchGRID('gpt3_1.grd');
                                      
            elseif grid_res == 5
                   %Call the "searchGRID.m" function
                   [gpt_fn, gptfile]= searchGRID('gpt3_5.grd');                  
            end
            
     end %//if strcmpi(GPTmodel,'GPT2')
     
   elseif regexpi(GRIDfile, '\w*.grd$')
            
          %******CHECK MODEL TYPE
          if strcmpi(GPTmodel,'GPT2')
         
             %*****SEARCH FOR '.grd' file
             %Call the "searchGRID.m" function
             [gpt_fn, gptfile]= searchGRID('gpt2_5.mat');
        
          elseif  strcmpi(GPTmodel,'GPT2w')
                  %*****SEARCH FOR '.grd' file
            
                  %CHECK GRID RESOLUTION 
                  if grid_res == 1
                     %Call the "searchGRID.m" function
                     [gpt_fn, gptfile]= searchGRID('gpt2_1w.mat');
            
                     if any([isempty(gpt_fn) isempty(gptfile)])
                        %Call the "searchGRID.m" function
                        [gpt_fn, gptfile]= searchGRID('gpt2_1wA.mat');
                     end 
               
                  elseif grid_res == 5
                         %Call the "searchGRID.m" function
                         [gpt_fn, gptfile]= searchGRID('gpt2_5w.mat');                  
                  end 
            
          elseif  strcmpi(GPTmodel,'GPT3')
                  %*******SEARCH FOR '.grd' file
            
                  %CHECK GRID RESOLUTION 
                  if grid_res == 1
                     %Call the "searchGRID.m" function
                     [gpt_fn, gptfile]= searchGRID('gpt3_1.mat');
                                      
                  elseif grid_res == 5
                         %Call the "searchGRID.m" function
                         [gpt_fn, gptfile]= searchGRID('gpt3_5.mat');                  
                  end 
            
          end  %//if strcmpi(GPTmodel,'GPT2')
          
   end %//if regexpi(GRIDfile, '\w*.mat$')
   
end %//if any([isempty(gpt_fn) isempty(gptfile)])
 
%IF THE SPECIFIED GRID RESOLUTION IS NOT AVAILABLE,USE THE ONE AVAILABLE IN
%=================================THE SYSTEM
if any([~exist('gptfile','var') ~exist('gpt_fn','var')]) 
    %GET USER INPUT GRID RESOLUTION
    if grid_res == 5 %IF RESOLUTION IS 5, CHECK 1
        
       if strcmpi(GPTmodel,'GPT2w')
          %Call the "searchGRID.m" function
          [gpt_fn, gptfile]= searchGRID('gpt2_1w.mat');
          grid_res = 1;%NEW GRID RESOLUTION
          
          if any([isempty(gptfile) isempty(gpt_fn)])%IF MAT FILE ISN'T AVAILABLE,USE GRID FILE
             %Call the "searchGRID.m" function
             [gpt_fn, gptfile]= searchGRID('gpt2_1w.grd');  
          end
          
       elseif strcmpi(GPTmodel,'GPT3')
              %Call the "searchGRID.m" function
              [gpt_fn, gptfile]= searchGRID('gpt3_1.mat');
              grid_res = 1;%NEW GRID RESOLUTION
          
              if any([isempty(gptfile) isempty(gpt_fn)])%IF MAT FILE ISN'T AVAILABLE,USE GRID FILE
                 %Call the "searchGRID.m" function
                 [gpt_fn, gptfile]= searchGRID('gpt3_1.grd');  
              end 
       end 
       
    elseif grid_res == 1 %IF RESOLUTION IS 1, CHECK 5
        
           if strcmpi(GPTmodel,'GPT2w')%if model is GPT2w
              %Call the "searchGRID.m" function
              [gpt_fn, gptfile]= searchGRID('gpt2_5w.mat');
              grid_res = 5;%NEW GRID RESOLUTION
          
              if any([isempty(gptfile) isempty(gpt_fn)])%IF MAT FILE ISN'T AVAILABLE,USE GRID FILE
                 %Call the "searchGRID.m" function
                 [gpt_fn, gptfile]= searchGRID('gpt2_5w.grd');  
              end 
          
           elseif strcmpi(GPTmodel,'GPT3')%if model is GPT3
                  %Call the "searchGRID.m" function
                  [gpt_fn, gptfile]= searchGRID('gpt3_5.mat');
                  grid_res = 5;%NEW GRID RESOLUTION
          
                 if any([isempty(gptfile) isempty(gpt_fn)])%IF MAT FILE ISN'T AVAILABLE,USE GRID FILE
                    %Call the "searchGRID.m" function
                    [gpt_fn, gptfile]= searchGRID('gpt3_5.grd');  
                 end  
          
           end  
        
    end
else
    if any([isempty(gptfile) isempty(gpt_fn)])
       %GET USER INPUT GRID RESOLUTION
       if grid_res == 5 %IF RESOLUTION IS 5, CHECK 1
        
          if strcmpi(GPTmodel,'GPT2w')
             %Call the "searchGRID.m" function
             [gpt_fn, gptfile]= searchGRID('gpt2_1w.mat');
             grid_res = 1;%NEW GRID RESOLUTION
          
             if any([isempty(gptfile) isempty(gpt_fn)])%IF MAT FILE ISN'T AVAILABLE,USE GRID FILE
                %Call the "searchGRID.m" function
                [gpt_fn, gptfile]= searchGRID('gpt2_1w.grd');  
             end  
          
          elseif strcmpi(GPTmodel,'GPT3')
                %Call the "searchGRID.m" function
                [gpt_fn, gptfile]= searchGRID('gpt3_1.mat');
                grid_res = 1;%NEW GRID RESOLUTION
          
                if any([isempty(gptfile) isempty(gpt_fn)])%IF MAT FILE ISN'T AVAILABLE,USE GRID FILE
                   %Call the "searchGRID.m" function
                  [gpt_fn, gptfile]= searchGRID('gpt3_1.grd');  
                end 
          
          end 
       
       elseif grid_res == 1 %IF RESOLUTION IS 1, CHECK 5
        
           if strcmpi(GPTmodel,'GPT2w')%if model is GPT2w
               
              %Call the "searchGRID.m" function
              [gpt_fn, gptfile]= searchGRID('gpt2_5w.mat');
              grid_res = 5;%NEW GRID RESOLUTION
          
              if any([isempty(gptfile) isempty(gpt_fn)])%IF MAT FILE ISN'T AVAILABLE,USE GRID FILE
                 %Call the "searchGRID.m" function
                [gpt_fn, gptfile]= searchGRID('gpt2_5w.grd');  
              end 
          
           elseif strcmpi(GPTmodel,'GPT3')
                  %Call the "searchGRID.m" function
                  [gpt_fn, gptfile]= searchGRID('gpt3_5.mat');
                  grid_res = 5;%NEW GRID RESOLUTION
          
                 if any([isempty(gptfile) isempty(gpt_fn)])%IF MAT FILE ISN'T AVAILABLE,USE GRID FILE
                    %Call the "searchGRID.m" function
                    [gpt_fn, gptfile]= searchGRID('gpt3_5.grd');  
                 end   
          
           end   
        
       end 
       
    end
    
end

%IF GRID FILES ARE STILL NOT FOUND, RETURN A MESSAGE
%=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-+
if any([isempty(gptfile) isempty(gpt_fn)])         
   beep %Give a beep sound
   errmsg{1}=sprintf('Grid file was not found in Directory.\n');                 
   errmsg{2}='Please ensure that grid file is in Directory & Try Again.';    
   errordlg(errmsg,'Grid file input Error','modal') 
   
    %Return Empty matrix([]) as output
    gpt_grid=[]; 
    
    return 
    
end    
%=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-+

%****************************READ GRID FILE 

%CONVERT CELL ARRAY TO CHARACTER ARRAY
if iscell(gpt_fn)
   gpt_fn  = gpt_fn{1};
   gptfile = gptfile{1};
end
   
%******CHECK EXTENSION TYPE('.mat or '.grd')
if regexpi(gpt_fn, '\w*.mat$')
    
   %***Import the MAT file
   newData1 = load('-mat', gptfile);

   %Create new variables in the base workspace from those fields.
   vars = fieldnames(newData1);
   for i = 1:length(vars)
       assignin('base', vars{i}, newData1.(vars{i}));
   end 
   dataout = zeros(1,length(newData1));
   field = fieldnames(newData1);

   vec=getfield(newData1,field{:});
   
   %***************READ DATA
   %****CHECK GPT MODEL
  if strcmpi(GPTmodel,'GPT2') 
     %*****READ MEAN AND AMPLITUDE VALUES
      p_grid  = vec(:,3:7);          % pressure in Pascal
      T_grid  = vec(:,8:12);         % temperature in Kelvin
      Q_grid  = vec(:,13:17)./1000;  % specific humidity in kg/kg
      dT_grid = vec(:,18:22)./1000;  % temperature lapse rate in Kelvin/m
      u_grid  = vec(:,23);           % geoid undulation in m
      Hs_grid = vec(:,24);           % orthometric grid height in m  
      ah_grid = vec(:, 25:29)/1000;  % hydrostatic mapping function coefficient, dimensionless
      aw_grid = vec(:, 30:34)/1000;  % wet mapping function coefficient, dimensionless
      
      %COMBINE ALL DATA TO ONE CELL GRID
      gpt_grid = {p_grid, T_grid, Q_grid, dT_grid, u_grid, Hs_grid,ah_grid,aw_grid};  
      
  elseif strcmpi(GPTmodel,'GPT2w')
      
         if any([isequal(grid_res,1), isequal(grid_res,5)])
             
            %*****READ MEAN AND AMPLITUDE VALUES
            p_grid  = vec(:,3:7);          % pressure in Pascal
            T_grid  = vec(:,8:12);         % temperature in Kelvin
            Q_grid  = vec(:,13:17)./1000;  % specific humidity in kg/kg
            dT_grid = vec(:,18:22)./1000;  % temperature lapse rate in Kelvin/m
            u_grid  = vec(:,23);           % geoid undulation in m
            Hs_grid = vec(:,24);           % orthometric grid height in m 
            ah_grid = vec(:,25:29)./1000;  % hydrostatic mapping function coefficient, dimensionless
            aw_grid = vec(:,30:34)./1000;  % wet mapping function coefficient, dimensionless
            la_grid = vec(:,35:39);    	   % wet mapping function coefficient, dimensionless
            Tm_grid = vec(:,40:44);    % mean temperature in Kelvin    

            %COMBINE ALL DATA TO ONE CELL GRID
            gpt_grid = {p_grid, T_grid, Q_grid, dT_grid, u_grid, Hs_grid, ah_grid, aw_grid, la_grid, Tm_grid };    
            
         end 
         
  elseif strcmpi(GPTmodel,'GPT3')
      
         if any([isequal(grid_res,1), isequal(grid_res,5)])
 
            p_grid    = vec(:,3:7);          % pressure in Pascal
            T_grid    = vec(:,8:12);         % temperature in Kelvin
            Q_grid    = vec(:,13:17)/1000;   % specific humidity in kg/kg
            dT_grid   = vec(:,18:22)/1000;   % temperature lapse rate in Kelvin/m
            u_grid    = vec(:,23);           % geoid undulation in m
            Hs_grid   = vec(:,24);           % orthometric grid height in m
            ah_grid   = vec(:,25:29)/1000;   % hydrostatic mapping function coefficient, dimensionless
            aw_grid   = vec(:,30:34)/1000;   % wet mapping function coefficient, dimensionless
            la_grid   = vec(:,35:39);    	   % water vapor decrease factor, dimensionless
            Tm_grid   = vec(:,40:44);        % mean temperature in Kelvin
            Gn_h_grid = vec(:,45:49)/100000;    % hydrostatic north gradient in m
            Ge_h_grid = vec(:,50:54)/100000;    % hydrostatic east gradient in m
            Gn_w_grid = vec(:,55:59)/100000;    % wet north gradient in m
            Ge_w_grid = vec(:,60:64)/100000;    % wet east gradient in m
    
            %COMBINE ALL DATA TO ONE CELL GRID
            gpt_grid = {p_grid, T_grid, Q_grid, dT_grid, u_grid, Hs_grid, ah_grid, aw_grid, la_grid, Tm_grid , Gn_h_grid , Ge_h_grid , Gn_w_grid , Ge_w_grid};    
   
         end
         
  end %//if strcmpi(GPTmodel,'GPT2')
  
elseif regexpi(gpt_fn, '\w*.grd$')
           
       %****CHECK GPT MODEL
       if strcmpi(GPTmodel,'GPT2') %if GPT2 MODEL
           
          %READ GRIDFILE
          fid = fopen(gptfile,'r');
          vec = textscan( fid, '%f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f', 'HeaderLines', 1 , 'CollectOutput', true );
          vec = vec{1};
          fclose (fid);
          
          %*****READ MEAN AND AMPLITUDE VALUES
          p_grid  = vec(:,3:7);          % pressure in Pascal
          T_grid  = vec(:,8:12);         % temperature in Kelvin
          Q_grid  = vec(:,13:17)./1000;  % specific humidity in kg/kg
          dT_grid = vec(:,18:22)./1000;  % temperature lapse rate in Kelvin/m
          u_grid  = vec(:,23);           % geoid undulation in m
          Hs_grid = vec(:,24);           % orthometric grid height in m  
          ah_grid = vec(:, 25:29)/1000;  % hydrostatic mapping function coefficient, dimensionless
          aw_grid = vec(:, 30:34)/1000;  % wet mapping function coefficient, dimensionless
          
          %COMBINE ALL DATA TO ONE CELL GRID
          gpt_grid = {p_grid, T_grid, Q_grid, dT_grid, u_grid, Hs_grid,ah_grid,aw_grid};
        
          
       elseif  strcmpi(GPTmodel,'GPT2w') %if GPT2w MODEL
          
               %CHECK GRID RESOLUTION
               if any([isequal(grid_res,1), isequal(grid_res,5)]) %if either 1x1 or 5x5 grid file
          
                  %READ GRIDFILE
                  fid = fopen(gptfile,'r');
                  vec = textscan( fid, '%f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f', 'HeaderLines', 1 , 'CollectOutput', true );
                  vec = vec{1};
                  fclose (fid);
                  
                 %*****READ MEAN AND AMPLITUDE VALUES
                 p_grid  = vec(:,3:7);          % pressure in Pascal
                 T_grid  = vec(:,8:12);         % temperature in Kelvin
                 Q_grid  = vec(:,13:17)./1000;  % specific humidity in kg/kg
                 dT_grid = vec(:,18:22)./1000;  % temperature lapse rate in Kelvin/m
                 u_grid  = vec(:,23);           % geoid undulation in m
                 Hs_grid = vec(:,24);           % orthometric grid height in m 
                 ah_grid = vec(:,25:29)./1000;  % hydrostatic mapping function coefficient, dimensionless
                 aw_grid = vec(:,30:34)./1000;  % wet mapping function coefficient, dimensionless
                 la_grid = vec(:,35:39);    	   % wet mapping function coefficient, dimensionless
                 Tm_grid = vec(:,40:44);    % mean temperature in Kelvin    

                 %COMBINE ALL DATA TO ONE CELL GRID
                 gpt_grid = {p_grid, T_grid, Q_grid, dT_grid, u_grid, Hs_grid, ah_grid, aw_grid, la_grid, Tm_grid };    
             
               end
               
       elseif strcmpi(GPTmodel,'GPT3')
      
              if any([isequal(grid_res,1), isequal(grid_res,5)])
                   
                  %READ GRIDFILE
                  fid = fopen(gptfile,'r');
                  vec = textscan( fid, '%f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f', 'HeaderLines', 1 , 'CollectOutput', true );
                  vec = vec{1};
                  fclose (fid);
            
                  %*****READ MEAN AND AMPLITUDE VALUES     
                  p_grid    = vec(:,3:7);          % pressure in Pascal
                  T_grid    = vec(:,8:12);         % temperature in Kelvin
                  Q_grid    = vec(:,13:17)/1000;   % specific humidity in kg/kg
                  dT_grid   = vec(:,18:22)/1000;   % temperature lapse rate in Kelvin/m
                  u_grid    = vec(:,23);           % geoid undulation in m
                  Hs_grid   = vec(:,24);           % orthometric grid height in m
                  ah_grid   = vec(:,25:29)/1000;   % hydrostatic mapping function coefficient, dimensionless
                  aw_grid   = vec(:,30:34)/1000;   % wet mapping function coefficient, dimensionless
                  la_grid   = vec(:,35:39);    	   % water vapor decrease factor, dimensionless
                  Tm_grid   = vec(:,40:44);        % mean temperature in Kelvin
                  Gn_h_grid = vec(:,45:49)/100000;    % hydrostatic north gradient in m
                  Ge_h_grid = vec(:,50:54)/100000;    % hydrostatic east gradient in m
                  Gn_w_grid = vec(:,55:59)/100000;    % wet north gradient in m
                  Ge_w_grid = vec(:,60:64)/100000;    % wet east gradient in m
    
                  %COMBINE ALL DATA TO ONE CELL GRID
                  gpt_grid = {p_grid, T_grid, Q_grid, dT_grid, u_grid, Hs_grid, ah_grid, aw_grid, la_grid, Tm_grid , Gn_h_grid , Ge_h_grid , Gn_w_grid , Ge_w_grid};    
                 
              end   
               
                  
       end %//if strcmpi(GPTmodel,'GPT2')
       
end %//if regexpi(gpt_fn, '\w*.mat$')
  
%SAVE OUTPUT
setappdata(0,'gpt_grid',gpt_grid )

%=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-+

%A.******SUBROUTINE TO CHECK IF GRID FILE IS IN CURRENT DIRECTORY

function [FileName, FilePath]= searchGRID(varargin)

file=varargin{1};%get input file

if ischar(file)%Check if file is character
  
  %GET LIST OF FILE(S) FROM DIRECTORY
  listFile = dir(file);%List of file(s) 
  
  if isempty(listFile) %if file is not in Current Directory 
      
     %**CHECK GPT grids from goGPS data FOLDER & SUB-FOLDERS  
     folder ='../data/TropoGRIDS/VMF files/GPT grids';
     
     %****FIRST CHECK IF FOLDER IS EMPTY OR NOT
     SizeListing=size(dir(folder)); %Size of listFiles

     %********CHECK IF FOLDER IS EMPTY OR NOT
     if any([SizeListing(1,1)==2,SizeListing(1,1)==0])% size of 2 means folder is empty & 0 means folder doesn't exist
          
        %IF FOLDER IS EMPTY,CHECK MAIN DIRECTORY 
        %---------------------------------------------------------
        %IMPORT GPT FILES FROM CURRENT DIRECTORY/FOLDER
        folder =pwd; %Get the folder from the Parent Directory containing GPT GRID FILES
        
        %=====================EXTRACT FOLDER CONTENT
        %***GET THE LIST OF FOLDER CONTENTS:
        %===(listFiles is a struct array)===  
        listFiles = dir(folder);%Get List of Folder/current directory Contents
       
     else
         %=====================EXTRACT FOLDER CONTENT
         %***GET THE LIST OF FOLDER CONTENTS:
         %===(listFiles is a struct array)===  
         listFiles = dir(folder);%Get List of Folder/current directory Contents        
     end
     
     dirIndex = [listFiles.isdir];%#Find the index for directories(zeros(0's) indicates files,& ones(1s) for folders or pointers)
     
     %GET LIST OF FILES
     fileList = {listFiles(~dirIndex).name}'; %Get list of files & remove pointers('.','..')
     
     %GET LIST OF SUB-FOLDERS/DIRECTORIES
     subDirs={listFiles(dirIndex).name};%Get list of Sub-folders/directories within parent folder 
                                        %NOTE:filenames={listFiles.name} gives the the list files & subfolders 
                                        %subdirs={listFiles(dirIndex).name} gives list of folders & poniter('.','..') OR                                             
                                        %subdirs1 = filenames([listFiles.isdir])

     %GET INDEX OF SUB-DIRECTORIES
     validIndex = ~ismember(subDirs,{'.','..'});%Find index of subdirectories that are not '.' or '..'                                    
     
     if ~isempty(fileList)% If fileList is not empty   
         
       %SORT FILES & REMOVE DUPLICATE FILES
       fileList=unique(fileList); 
       
       %***CONSTRUCT FULL PATH OF FILES
       try
          Filepath = cellfun(@(x) fullfile(folder,x),fileList,'UniformOutput',false);%Prepend path to files
         
       catch
            Filepath=fullfile(strcat(folder,'\',fileList));%Prepend path to files      
       end
                     
       %COMPARE INPUT GRIDFILE(file) WITH WHAT IS IN DIRECTORY(fileList) AND  
       %EXTRACT CORRESPONDING ROW(S) OF DATA (Cfile) IN EACH FILE
       %WHERE:id & iu are the INDEX OF  fileList in directory(d) & User
       %Input file(u) respectively-----------------------------------------
       [Cfile,id,iu] = intersect(fileList,file);%Grid  filenames (e.g.'gpt2_5w.mat')
       
       CfilePATH = Filepath(id);%Grid filpath(e.g.C:\Users\'gpt2_5w.mat')
       
       %ASSIGNMENT
       FileName = Cfile;
       FilePath = CfilePATH;
       
       %CONVERT CELL ARRAYS TO CHARACTER ARRAYS
       if iscell(FilePath)
          FilePath=char(FilePath);
       end       
       if iscell(FileName)
          FileName=char(FileName);
       end   
        
     else %IF THERE ARE NO FILEs IN THE MAIN FOLDER,GIVE AN INDICATION        
         empty_fileList = 1; %#ok<*NASGU> %flag to indicate there are no files in main folder 
         
     end %//if ~isempty(fileList)% If fileList is not empty 
         
       %WHAT IF THERE WERE FILES BUT WERE NOT GPT GRID FILES
       if exist('Cfile','var')
          
          if isempty(Cfile) 
             empty_fileList = 1; %flag to indicate there is No corresponding file in main folder 
       
          end 
          
       else 
           empty_fileList = 1; %flag to indicate there are no files in main folder 
           
       end 
       
      %***OPEN ALL SUB-FOLDERS & GET GPT GRID FILES IF THERE WERE NO MATCH
      %IN THE MAIN FOLDER-------------------------------------------------
        
      if exist('empty_fileList','var')
           
        p=1;%Loop index       
                
       for iDir = find(validIndex) % Loop over valid sub-directories  
           
          %***GET THE LIST OF SUB-FOLDER CONTENTS:
          %   ===(sublist is a struct array)===    
          subList = dir(strcat(folder,'\',subDirs{iDir})); 
                                                                                            
          %SIZE OF SUB-FOLDER CONTENT
          Sizesublist=size(subList);

          %**IF  SUB-FOLDER NOT EMPTY 
          if Sizesublist(1,1)~=2 % size of 2 means folder is empty 
                        
            %************CHECK FOR ONLY FILES IN THE SUB-FOLDER
            subList([subList.isdir]) = [];%remove non-directories /current and parent directory.
            subIndex = [subList.isdir];%#Find the index for directories  
                      
            %GET LIST OF FILES
            fileLists = {subList(~subIndex).name}';%Get list of files & remove pointers('.','..')                                            
              
            %CONSTRUCT THE FULLFILE PATH USING THE SUB-DIRECTORY & THE FILENAME
            try
               Filepath=fullfile(strcat(folder,'\',subDirs{iDir},'\',fileLists));%Prepend path to files           
            catch  
                 Filepath= cellfun(@(x) fullfile(strcat(folder,'\',subDirs{iDir}),x),...  
                                                     fileLists,'UniformOutput',false); %Prepend path to files                     
            end %try
                              
          end %if Sizesublist(1,1)~=2
          
          Filelists{p,1}= fileLists; %#ok<*AGROW> %File Lists
          Filepaths{p,1}=Filepath; %File Paths
           
          p=p+1;%Update index
                    
         %COMBINE ALL FILES FROM VARIOUS SUB-FOLDERS  
         filelists = vertcat(Filelists{:}) ;   
         filepaths=vertcat(Filepaths{:}) ; 
                
       end  %for iDir = find(validIndex)
       
       %SORT FILES & REMOVE DUPLICATE FILES
       [fileList,i_fl] = unique(filelists);
        filePaths=(filepaths(i_fl));
        
       %COMPARE INPUT FILE WITH LIST OF FILEs IN SUB-DIRECTORIEs
       try
          %COMPARE INPUT GRIDFILE(file) WITH WHAT IS IN DIRECTORY(fileList) AND  
          %EXTRACT CORRESPONDING ROW(S) OF DATA (Cfile) IN EACH FILE
          %WHERE:id & iu are the INDEX OF  fileList in directory(d) & User
          %Input file(u) respectively-----------------------------------------
          [Cfile,id,iu] = intersect(fileList,file);%Grid  filenames (e.g.'gpt2_5w.mat')
       
          CfilePATH = filePaths(id);%Grid filpath(e.g.C:\Users\'gpt2_5w.mat')
       
          %ASSIGNMENT
          FileName = Cfile;
          FilePath = CfilePATH; 
          
          %CONVERT CELL ARRAY TO CHARACTER
          if iscell(FilePath)
             FilePath=char(FilePath);
          end   
                    
          if iscell(FileName)
             FileName=char(FileName);
          end     
          
       catch
             for k =1:length(fileList)%Loop over list of files
                        
                 if strcmpi(file,fileList(k))%Compare files
                    FileName=fileList(k);
                    FilePath=filePaths(k);
              
                    %CONVERT CELL ARRAY TO CHARACTER
                    if iscell(FilePath)
                       FilePath=char(FilePath);
                    end  
                    
                      if iscell(FileName)
                         FileName=char(FileName);
                      end   
            
                 end  
                        
             end   %for k =1:length(filelists)
       
       end 
         
         if ~exist('FileName','var') | ~exist('FilePath','var')
            FileName=[];
            FilePath=[];
         end 
       
       end %//if exist('empty_fileList','var')
                     
  else
      FileName=file;
      FilePath=file;
      
      %CONVERT CELL ARRAY TO CHARACTER
      if iscell(FilePath)
         FilePath=char(FilePath);
      end        
      if iscell(FileName)
         FileName=char(FileName);
      end    
          
  end %if isempty(listFile)
  
end %if ischar(file)
%******************************************END OF searcGRID.m 
%=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-