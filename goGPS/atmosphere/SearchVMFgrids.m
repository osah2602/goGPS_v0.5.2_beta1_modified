%***************************************************************************
%***************************************************************************
%DESCRIPTION:                                                              *                                              
%            The "SearchVMFgrids" subroutine searches Vienna Mapping...    +
%             Function Grid (VMF)files in a given Directory/Folder         +
%            as indicated in the sub-routine. Situation where VMFG files   + 
%            are not found in the default directory/folder,it searches     +   
%            recursively through all Sub-Directories/Folders of the...     +  
%            given Directory. VMFG files are extracted by looping ...      + 
%            through all the listed files in the provided folder or ...    +
%            sub-folders.Finally,if VMFG files are still not found in the..+ 
%            default directory/folder and its sub-folders,the search for   +
%            VMFG files is extended to the current directory.              +
%=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
%%USAGE:                                                                   +
%      [VMF_grid_found,VMFgrids] = SearchVMFgrids()                        +
%=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
%NOTE:                                                                     *
%     THE FF FOLDERS/DIRECTORIES ARE SEARCHED FOR VMF GRID FILES           *
%1.   'VMF grids' folder; main folder for VMF GRIDs                        *
%2.   'VMF files' folder; main folder for VMF GRIDs                        *
%3.   'TropoGRIDS' folder; main folder for Tropospheric & Tm GRIDs         *
%4.   'data' folder; folder for goGPS data                                 *
%5.   'pwd'; Current folder/Directory. In goGPS, it is the 'goGPS' folder  *
%=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
%--------------------------------------------------------------------------
%OUTPUTs                                                                   *
%1.VMF_grid_found : Flag to indicate presence of VMF grid file.1 to ...*
%                   indicate there is file in the folder & 0 to mean ...   *
%                   absence of file                                        *
%2.VMFgrids     : Structure array of VMF grid files which includes:      *
%1.    h0files  : Zero(0) Hour VMFG File Name(s) eg:'VMFG_20180101.H00'    * 
%2.    h6files  : Six(6) Hour VMFG File  Name(s) eg:'VMFG_20180101.H06'    *  
%3.    h12files : Twelve(12) Hour VMFG File Name(s) eg:'VMFG_20180101.H12' *  
%4.    h18files : Eighteen(18) Hour VMFG File Name(s)eg:'VMFG_20180101.H18'*
%5.    H0files  : Zero(0) Hour VMFG File Path(s)                           * 
%                 eg:'C:\Users\...data\VMFG_20180101.H00'                  * 
%6.    H6files  : Six(6) Hour VMFG File  Path(s)                           *  
%7.    H12files : Twelve(12) Hour VMFG File Path(s)                        *  
%8.    H18files : Eighteen(18) Hour VMFG File Path(s)                      * 
%9.    Orography: Orography file                                           *
%=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-+
%=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-+
% WRITTEN BY: 
%            OSAH SAMUEL, MSC GEOMATIC ENGINEERING (PhD STUDENT)           +
%            Email:osahsamuel@yahoo.ca/osahsamuel@gmail.com                +
%            Phone:+233(0)246137410 / +233(0)509438484                     +
%=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-+
%**************************************************************************+

function [VMF_grid_found,VMFgrids] = SearchVMFgrids()

%         ------------------------------------------------
%*********SEARCH FOR VMF GRID FILES FROM VARIOUS DIRECTORY
%         ------------------------------------------------

%1.SEARCH VMF GRIDs from 'VMF grids' folder (A SUB-FOLDER]
[VMF_grid_found,VMFgrids] = SearchVMFgrid('../data/TropoGRIDS/VMF files/VMF grids');

if isequal(VMF_grid_found,0)%if VMF grids are not in 'VMF grids' folder, try the 'VMF files' folder
    
   %2.SEARCH VMF GRIDs from 'VMF files' folder
   [VMF_grid_found,VMFgrids] = SearchVMFgrid('../data/TropoGRIDS/VMF files');

   if isequal(VMF_grid_found,0)%if VMF grids are not in 'VMF files' folder, try the 'TropoGRIDS' folder
       
      %3.SEARCH VMF GRIDs from 'TropoGRIDS' folder
      [VMF_grid_found,VMFgrids] = SearchVMFgrid('../data/TropoGRIDS');

      if isequal(VMF_grid_found,0)%if VMF grids are not in 'TropoGRIDS' folder, try the 'data' folder
    
         %4.SEARCH VMF GRIDs from 'data' folder 
         [VMF_grid_found,VMFgrids] = SearchVMFgrid('../data');%searching the 'data' folder
   
         if isequal(VMF_grid_found,0) %if VMF grids are not in the 'data' folder, try the 'current' directory
       
            %5.SEARCH VMF GRIDs from 'current' folder  
            [VMF_grid_found,VMFgrids] = SearchVMFgrid(pwd);%if VMF grids are not in the 'current' directory, then issue a message
      
            if isequal(VMF_grid_found,0)%if VMF grids are not in the 'current' directory, then issue a message
          
                beep %Give a beep sound
                errmsg2{1}=sprintf('No VMF grid file(s) found in goGPS directory.\n');
                errmsg2{2}=sprintf('Please Provide VMF grid file(s) & Try Again.\n');
                errmsg2{3}=sprintf('VMF grid(s) can as well be downloaded from : http://vmf.geo.tuwien.ac.at/trop_products/GRID/\n');
                warndlg(errmsg2,'VMF grid file Error','modal')
                
                 return
                    
            end   
                               
         end 
                                
      end    
      
   end   
   
end  
%===================================END OF Search FOR FOLDER WITH VMF grids
%=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
      
function [VMF_grid_found,VMFgrids] = SearchVMFgrid(directory)

%ASSIGNMENT
folder = directory;  %folder/Directory suppose to be containing VMFG FILES in goGPS

%****FIRST CHECK IF FOLDER IS EMPTY OR NOT
SizeListing=size(dir(folder)); %Size of folder content

%********CHECK IF FOLDER IS EMPTY OR NOT
if any([SizeListing(1,1)==2,SizeListing(1,1)==0])% size of 2 means folder is empty & 0 means folder doesn't exist
                  
   VMF_grid_found = 0; %Flag to indicate the absence of VMF grid file in directory
   
   %SET VARIOUS FILE NAMES & PATHS TO EMPTY([])
   h0files  = []; h6files  = []; h12files = []; h18files = [];
   H0files  = []; H6files  = []; H12files = []; H18files = [];
   orogFILE = [];
   orogPATH = [];
   VMFgrids = [];
   
   return             
      
else  
     %***GET THE LIST OF FOLDER CONTENTS:
     %===(listFiles is a struct array)=== 
     listFiles = dir(folder);%Get List of Folder/current directory Contents
                     
end   %//if any([SizeListing(1,1)==2,SizeListing(1,1)==0])
                          
if exist('listFiles','var') %If Folder is not Empty
                  
   %**CHECK IF FOLDER HAS FILES & SUBFOLDERS
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
          
       %*****GET INDEXEs FOR THE VARIOUS GRID FILES(i.e.Integer Val /[])
       H0Index = regexpi(fileList, '\w*.H00$'); %Index for 0 hour files
       H6Index = regexpi(fileList, '\w*.H06$'); %Index for 6 hour files
       H12Index = regexpi(fileList, '\w*.H12$');%Index for 12 hour files
       H18Index = regexpi(fileList, '\w*.H18$');%Index for 18 hour files
          
       %*********LOOK FOR GRID FILEs WITH .txt EXTENSION
       %NB:some browser during file download will save grid file as text
       %   file with the extension '.txt'.
       H0_1Index = regexpi(fileList, '\w*.H00.txt$'); %Index for 0 hour files
       H6_1Index = regexpi(fileList, '\w*.H06.txt$'); %Index for 6 hour files
       H12_1Index = regexpi(fileList, '\w*.H12.txt$');%Index for 12 hour files
       H18_1Index = regexpi(fileList, '\w*.H18.txt$');%Index for 18 hour files
          
       %************GET THE GRID FILES
       h0Files  = fileList(~cellfun(@isempty, H0Index));%  0 hour files(eg:'VMFG_20180101.H00')
       h6Files  = fileList(~cellfun(@isempty, H6Index)); % 6 hour files(eg:'VMFG_20180101.H06')
       h12Files = fileList(~cellfun(@isempty, H12Index));%12 hour files(eg:'VMFG_20180101.H12')
       h18Files = fileList(~cellfun(@isempty, H18Index));%18 hour files(eg:'VMFG_20180101.H18')
          
       %*********THOSE WITH '.txt' EXTENSION
       h0_1Files  = fileList(~cellfun(@isempty, H0_1Index));%  0 hour files(eg:'VMFG_20180101.H00.txt')
       h6_1Files  = fileList(~cellfun(@isempty, H6_1Index)); % 6 hour files(eg:'VMFG_20180101.H06.txt')
       h12_1Files = fileList(~cellfun(@isempty, H12_1Index));%12 hour files(eg:'VMFG_20180101.H12.txt')
       h18_1Files = fileList(~cellfun(@isempty, H18_1Index));%18 hour files(eg:'VMFG_20180101.H18.txt')
 
       %****************WORK ON VARIOUS GRID FILEs
          
       %IF Any of the Files is Found(i.e. not empty([])),save its  fullfile path 
       %============                       +                        ============
       %1.******0 HOUR GRID FILE
       if ~isempty(h0Files) & ~isempty(h0_1Files)% If h0files & h0_1files are not empty
                                          
           %CONCATENATE h0files & h0_1files(i.e. Fname0)
           h0files = [h0Files;h0_1Files];
             
           %***CONSTRUCT THE FULL PATH       
           H0files = cellfun(@(x) fullfile(folder,x),h0Files,'UniformOutput',false);%Prepend path to files  
             
       elseif ~isempty(h0Files) & isempty(h0_1Files)% If h0files is not empty & h0_1files is empty
                              
              h0files = h0Files;
              %***CONSTRUCT THE FULL PATH        
              H0files = cellfun(@(x) fullfile(folder,x),h0Files,'UniformOutput',false);%Prepend path to files
             
       elseif  isempty(h0Files) & ~isempty(h0_1Files)% If h0files is empty & h0_1files is not empty
                                 
               h0files = h0_1Files;
               %***CONSTRUCT THE FULL  PATH      
               H0files = cellfun(@(x) fullfile(folder,x),h0_1Files,'UniformOutput',false);%Prepend path to files
             
       else   %if they are both empty
              
              h0files = [];
              H0files = [];%Assign empty([]) matrix
              
       end    %if ~isempty(h0files)
           
       %2.******6 HOUR GRID FILE          
       if ~isempty(h6Files) & ~isempty(h6_1Files)% If h6files & h6_1files are not empty
                                        
          %CONCATENATE h6files & h6_1files(i.e. Fname6)
                         
          h6files = [h6Files;h6_1Files];
             
          %********CONSTRUCT THE FULL PATH 
          H6files = cellfun(@(x) fullfile(folder,x),h6Files,'UniformOutput',false);%Prepend path to files 
          
       elseif ~isempty(h6Files) & isempty(h6_1Files)% If h6files is not empty & h6_1files is empty
                             
              h6files = h6Files;
              %***CONSTRUCT THE FULL PATH        
              H6files = cellfun(@(x) fullfile(folder,x),h6Files,'UniformOutput',false);%Prepend path to files
             
       elseif  isempty(h6Files) & ~isempty(h6_1Files)% If h6files is empty & h6_1files is not empty
                             
               h6files = h6_1Files;
               %***CONSTRUCT THE FULL PATH       
               H6files = cellfun(@(x) fullfile(folder,x),h6_1Files,'UniformOutput',false);%Prepend path to files
             
       else %if they are both empty
                          
            h6files = [];
            H6files = [];%Assig empty([]) matrix
              
       end  %if ~isempty(h6files)
          
       %3.******12 HOUR GRID FILE                    
       if ~isempty(h12Files) & ~isempty(h12_1Files)% If h12files & h12_1files are not empty
                                     
          %CONCATENATE h12files & h12_1files(i.e. Fname12)
          h12files = [h12Files;h12_1Files];
             
          %********CONSTRUCT THE FULL PATH
          H12files = cellfun(@(x) fullfile(folder,x),h12Files,'UniformOutput',false);%Prepend path to files  
            
       elseif ~isempty(h12Files) & isempty(h12_1Files)% If h12files is not empty & h12_1files is empty
                               
              h12files = h12Files;
              %***CONSTRUCT THE FULL PATH        
              H12files = cellfun(@(x) fullfile(folder,x),h12Files,'UniformOutput',false);%Prepend path to files
             
       elseif  isempty(h12Files) & ~isempty(h12_1Files)% If h12files is empty & h12_1files is not empty
                             
               h12files = h12_1Files;
               %***CONSTRUCT THE FULL PATH       
               H12files = cellfun(@(x) fullfile(folder,x),h12_1Files,'UniformOutput',false);%Prepend path to files
             
       else %if they are both empty  
                          
            h12files = [];
            H12files = [];%Assig empty([]) matrix
              
       end   %if ~isempty(h12files)
          
       %4.******18 HOUR GRID FILE           
       if ~isempty(h18Files) & ~isempty(h18_1Files)% If h18files & h18_1files are not empty
                          
          %CONCATENATE h18files & h18_1files(i.e. Fname18)
          h18files = [h18Files;h18_1Files];
             
          %********CONSTRUCT THE FULL PATH
          H18files = cellfun(@(x) fullfile(folder,x),h18Files,'UniformOutput',false);%Prepend path to files  
            
       elseif  ~isempty(h18Files) & isempty(h18_1Files)% If h18files is not empty & h18_1files is empty
                               
               h18files = h18Files;
               %***CONSTRUCT THE FULL PATH        
               H18files = cellfun(@(x) fullfile(folder,x),h18Files,'UniformOutput',false);%Prepend path to files
             
       elseif  isempty(h18Files) & ~isempty(h18_1Files)% If h18files is empty & h18_1files is not empty
                            
               h18files = h18_1Files;
               %***CONSTRUCT THE FULL PATH       
               H18files = cellfun(@(x) fullfile(folder,x),h18_1Files,'UniformOutput',false);%Prepend path to files
             
       else  %if they are both empty  
                          
             h18files = [];
             H18files = [];%Assig empty([]) matrix
              
       end   %if ~isempty(h18files)
                      
                      
      %***************GET OROGRAPHY FILE AS WELL 
      if any([~isempty(H0files),~isempty(H6files),~isempty(H12files),~isempty(H18files)]) | all ([~isempty(H0files),~isempty(H6files),~isempty(H12files),~isempty(H18files)])
                     
         %*****CONSTRUCT FILEs PATH
         Filepath=fullfile(strcat(folder,'\',fileList));%Prepend path to files
                     
         u=1;%Loop index for orography file
          
         for r=1:length(fileList) 
                                       
             if strncmpi(fileList(r),'orography_ell',13)| strncmpi(fileList(r),'orography',9)| strncmpi(fileList(r),'orog',4)
                
                orogFILE=fileList(r,1);%Orography file
                                 
                %********CONSTRUCT THE FULL PATH  
                try
                    orogPATH(u,1)=Filepath(r,1);%Orography file  path
                                 
                catch     
                      orogPATH = cellfun(@(x) fullfile(folder,x),orogFILE,'UniformOutput',false);%Prepend path to files  
                end     
                 
                u=u+1;%Update index
                                                  
             end  %//if strncmpi(fileList(r),'orography_ell',13)|   
             
         end  %//for r=1:length(fileList)     
                      
      else 
           %SET OROGRAPHY FILE NAME & PATH TO EMPTY([])
           orogFILE = [];
           orogPATH = [];
              
      end   %if any([~isempty(H0files),~isempty(H6files),~isempty(H12files)
                                                       
    else
        %SET VMF GRID FILE & OROGRAPHY FILE NAMES & PATHS TO EMPTY([])  
        h0files  = []; h6files  = []; h12files = []; h18files = [];
        H0files  = []; H6files  = []; H12files = []; H18files = [];
        orogFILE = [];
        orogPATH = [];
           
    end  %if ~isempty(fileList)
    %----------------------------------------------------------------------
    %----------------------------------------------------------------------
    
    %IF NO VMFG FILES EXIST IN THE CURRENT FOLDER, CHECK SUB-FOLDERS/              
    %SUB-DIRECTORIES FOR VMF GRID FILEs
    %                              +                        ==========
     if isempty(fileList)
              
        checkSubfolders=1; %Flag to search Sub-Directories/folders
           
     elseif all([isempty(h0Files),isempty(h6Files),isempty(h12Files),isempty(h18Files)])
                 
            checkSubfolders=1;%Flag to search Sub-Directories/folders 
                         
     else
          if any([~isempty(h0files),~isempty(h6files),~isempty(h12files),~isempty(h18files)])
                  
             checkSubfolders=0;%Flag not to search Sub-Directories/folders
             
             VMF_grid_found = 1; %Flag to indicate the presence of VMF grid file in directory
             
          end
                  
                                                                  
     end    %\\if isempty(fileList)
                
     
     if checkSubfolders==1 %if this condition is true,then either files in the current
                           %Directory were not VMF grid files or the there were no files at all.        
                           %System has to look through Sub-Folders or Sub-directories for VMF Grid Files                                                              
                                                                                                                                                                                                                                                                                  
        %***OPEN ALL SUB-FOLDERS & GET VMF GRID FILES
                
        p=1;%Loop index
        q=1;
                
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
                   
              end  %try
                       
              Filelists{p,1}= fileLists;
           
              Filepaths{p,1}=Filepath;
           
              p=p+1;%Update index
                       
              %COMBINE ALL RINEX FILES FROM VARIOUS SUB-FOLDERS  
              filelists = vertcat(Filelists{:});     
              filepaths=vertcat(Filepaths{:}) ;                     
                      
              if ~isempty(fileLists)   % If fileLists is not empty
            
                 j=1;k=1;l=1;m=1; %Loop index for the VMF grid file types
                 jj=1;kk=1;ll=1;mm=1; %Loop index for the VMF grid file types
                 u=1;%Loop index for orography file
                        
                 for i = 1:length(filelists)%****LOOP OVER THE FILES & CHECK FOR VMF GRID FILES
                                            %Identify the needed files and extract them separately
                                            %into different file names(H0files,H6files,H12files & H18files )
                                            %Note that dir also lists the directories,so you have to check for them. 
                                                         
                      %***GET NAME & PATH OF FILES FROM THE LIST OF FILES FOR EACH LOOP               
                           
                      FileName=filelists{i,1};
                      filepath=filepaths{i,1};
                                          
                      %****MATCH THE END OF EACH FILE IN THE 'FileName'
                      %For 0 Hour file look for 'H00'
                      %For 6 Hour file look for 'H06' ...
                      %For 12 Hour file look for 'H12' And...
                      %For 18 Hour file look for 'H18' 
                 
                      %1.******0 HOUR GRID FILE 
                      if regexpi(FileName,'\w*.H00$') 
  
                         %When Tested & the file is a 0h file,save...
                         %its fullfile path into the H0files cell array
                         %============            +        ============
                         h0Files{j,1}=FileName;
                         H0filess{j,1}=filepath;%#ok<*AGROW> %0h file in cells                             
    
                         j=j+1;%Update index
                                    
                      %2.******6 HOUR GRID FILE 
                   
                      elseif  regexpi(FileName,'\w*.H06$' ) %Else test the file if it is a 6 hour file
      
                              %When Tested & the file is a 6h file,save...
                              %its fullfile path into the H6files cell array
                              %============            +        ============
                              h6Files{k,1}=FileName;
                              H6filess{k,1}=filepath;%6 Hour file in cells
                                   
                              k=k+1; %update index
                 
                      %3.******12 HOUR GRID FILE          
                      elseif  regexpi(FileName,'\w*.H12$')
      
                              %When Tested & the file is a 12h file,save...
                              %its fullfile path into the H12files cell array
                              %============            +        ============ 
                              h12Files{l,1}=FileName;
                              H12filess{l,1}=filepath;%12 Hour file in cells
                                 
                              l=l+1;  %update index  
                                 
                      %4.******18 HOUR GRID FILE          
                      elseif regexpi(FileName,'\w*.H18$')
      
                             %When Tested & the file is a 18h file,save...
                             %its fullfile path into the H18files cell array
                             %============            +        ============  
                             h18Files{m,1}=FileName;
                             H18filess{m,1}=filepath;%18 Hour file in cells
                                 
                             m=m+1;  %update index 
                                 
                      %*********THOSE WITH '.txt' EXTENSION       
                      %1.******0 HOUR GRID FILE 
                      elseif  regexpi(FileName,'\w*.H00.txt$') 
                              h0_1Files{jj,1}=FileName;                         
                              H0_1filess{jj,1}=filepath;%0h file in cells                             
    
                              jj=jj+1;%Update index
                                  
                      %2.******6 HOUR GRID FILE 
                      elseif  regexpi(FileName,'\w*.H06.txt$') 
                              h6_1Files{kk,1}=FileName;
                              H6_1filess{kk,1}=filepath;%6h file in cells                             
    
                              kk=kk+1;%Update index  
                                   
                      %3.******12 HOUR GRID FILE 
                      elseif  regexpi(FileName,'\w*.H12.txt$') 
                              h12_1Files{ll,1}=FileName;
                              H12_1filess{ll,1}=filepath;%12h file in cells                             
    
                              ll=ll+1;%Update index          
                                        
                      %4.******18 HOUR GRID FILE 
                      elseif  regexpi(FileName,'\w*.H18.txt$') 
                              h18_1Files{mm,1}=FileName;
                              H18_1filess{mm,1}=filepath;%18h file in cells                             
    
                              mm=mm+1;%Update index 
                              
                                                         
                      end   %//if regexpi(FileName,'\w*.H00$')
                           
                      %IMPORT OROGRAPHY FILE AS WELL
                                                                                                    
                      if strncmpi(FileName,'orography_ell',13)| strncmpi(FileName,'orography',9)| strncmpi(FileName,'orog',4)
                
                         orogFILE{u,1}=FileName;%#ok<*NASGU> %Orography file 
                                  
                         %********CONSTRUCT THE FULL PATH                                           
                         orogPATH{u,1}=filepath;%Orography file  path
                                                          
                         u=u+1;%Update index
                                                                                                
                      end       
                      
                 end %//for i = 1:length(fileLists) 
                                                                                                               
              end    %//if ~isempty(fileLists)
                      
            else                          
                 %WRITE ERROR MESSAGE
                 erMsg{q,1}=sprintf('The Folder %s and  It''s Sub-Folder, %s , are empty / Contain No VMF Grid files .\n',folder,subDirs{iDir});
                         
                 q=q+1; %update index
                        
            end   %if Sizesublist(1,1)~=2
                   
        end    %for iDir = find(validIndex)
        
      if all([~exist('H0Files','var'),~exist('H6Files','var'),~exist('H12Files','var'),~exist('H18Files','var'),...
              ~exist('h0_1Files','var'),~exist('h6_1Files','var'),~exist('h12_1Files','var'),~exist('h18_1Files','var')]) 
         
         VMF_grid_found = 0;%Flag to indicate the absence of VMF grid file in directory
         
         VMFgrids = [];
         
         return
      
      else 
          VMF_grid_found = 1; %Flag to indicate the presence of VMF grid file in directory
          
          %*********************COMBINE VMF GRID FILEs
          %********0 HOUR FILEs
          if exist('H0filess','var') & exist('H0_1filess','var')
             h0files = [h0Files;h0_1Files];
             H0files = [H0filess;H0_1filess];
               
          elseif exist('H0filess','var') & ~exist('H0_1filess','var')
                 h0files = h0Files;
                 H0files = H0filess;
                  
          elseif ~exist('H0filess','var') & exist('H0_1filess','var')
                 h0files = h0_1Files;
                 H0files = H0_1filess;
                                         
          else  
              h0files = [];
              H0files = [];
                                                                                                
          end %\\if exist('H0filess','var')
                        
          %********6 HOUR FILEs
          if exist('H6filess','var') & exist('H6_1filess','var')
             h6files = [h6Files;h6_1Files];                                         
             H6files = [H6filess;H6_1filess];
                 
          elseif exist('H6filess','var') & ~exist('H6_1filess','var')
                 h6files = h6Files;  
                 H6files = H6filess;
                 
          elseif ~exist('H6filess','var') & exist('H6_1filess','var')
                 h6files = h6_1Files;  
                 H6files = H6_1filess;
          else     
              h6files = [];  
              H6files = [];                   
                                                                                     
          end   %\\if exist('H6filess','var')
            
          %********12 HOUR FILEs
          if exist('H12filess','var') & exist('H12_1filess','var')
             h12files = [h12Files;h12_1Files];                                  
             H12files = [H12filess;H12_1filess];
                
          elseif exist('H12filess','var') & ~exist('H12_1filess','var')
                 h12files = h12Files;                              
                 H12files = H12filess;
                                      
          elseif ~exist('H12filess','var') & exist('H12_1filess','var') 
                 h12files = h12_1Files;
                 H12files = H12_1filess;  
                  
          else   
              h12files = [];
              H12files = [];
               
          end %//if exist('H12filess','var') & exist('H12_1filess','var')
                 
          %********18 HOUR FILEs
          if exist('H18filess','var') & exist('H18_1filess','var')
             h18files = [h18Files;h18_1Files];                                     
             H18files = [H18filess;H18_1filess];
                
          elseif exist('H18filess','var') & ~exist('H18_1filess','var')
                 h18files = h18Files;                            
                 H18files = H18filess;
                                      
          elseif ~exist('H18filess','var') & exist('H18_1filess','var') 
                 h18files = h18_1Files;
                 H18files = H18_1filess;   
                                      
          else  
              h18files = [];
              H18files = []; 
                                              
          end  %\\if exist('H18filess','var')
          
         
      end %//if all([~exist('H0Files','var'),~exist('H6Files','var'),~exist('H12Files','var'),~exist('H18Files','var'),...
          %          ~exist('h0_1Files','var'),~exist('h6_1Files','var'),~exist('h12_1Files','var'),~exist('h18_1Files','var')]) 
              
           
     end %//if checkSubfolders==1
     
end %//if exist('listFiles','var')

%********SORT VMF grid FILES & REMOVE DUPLICATE FILES
        
%***********0 HOUR FILES(0h)
if ~isempty(H0files)
   
   [h0files,i_h0] = unique(h0files);%GET 0 HOUR FILES with no repetitions & INDEX(i_h0)
   H0files =(H0files(i_h0));%FILES WITH PATH          
end 
                
%***********6 HOUR FILES(6h)
if ~isempty(H6files)
   [h6files,i_h6] = unique(h6files);%GET 6 HOUR FILES with no repetitions & INDEX(i_h6)
    H6files =(H6files(i_h6));%FILES WITH PATH           
end 
        
%***********12 HOUR FILES(12h)
if ~isempty(H12files)
   [h12files,i_h12] = unique(h12files);%GET 12 HOUR FILES with no repetitions & INDEX(i_h12)
   H12files =(H12files(i_h12));%FILES WITH PATH             
end 
        
%***********18 HOUR FILES(18h)
if ~isempty(H18files)
   [h18files,i_h18] = unique(h18files);%GET 18 HOUR FILES with no repetitions & INDEX(i_h18)
   H18files =(H18files(i_h18));%FILES WITH PATH
end        
                      
%***********OROGRAPHY FILES(orography_ell)
if exist('orogFILE','var') 
    
   if ~isempty(orogPATH)
       
      [orofile,i_oro] = unique(orogFILE);%GET orography file with no repetitions & INDEX(i_oro)
      oropath =(orogPATH(i_oro));%FILES WITH PATH 
       
   else
       orofile = []; %SET orography FILE TO EMPTY([]) IF IT DOES NOT EXIST IN FOLDER   
       oropath = [];%SET orography FILE PATH TO EMPTY([])
   end
   
end     

%***********CHECK IF ALL GRID FILES ARE EMPTY([]) 
if all([isempty(H0files),isempty(H6files),isempty(H12files),isempty(H18files)])
   
   VMFgrids       = [];
   VMF_grid_found = 0;%Flag to indicate the absence of VMF grid file in directory
   
else
    
     %CREATE A STRUCTURE ARRAY(STRUCT) OF ALL THE GRID FILES 
     %***00 HOUR
     h00.name = h0files;%FILE NAME[e.g.:'VMFG_20180101.H00']
     h00.path = H0files;%FILE PATH['../data/TropoGRIDS/VMF files/VMF grids/'VMFG_20180101.H00'

     %****6 HOUR
     h06.name=h6files;%FILE NAME[e.g.:'VMFG_20180101.H06']
     h06.path=H6files;%FILE PATH['../data/TropoGRIDS/VMF files/VMF grids/'VMFG_20180101.H06'

     %****12 HOUR
     h12.name = h12files;%FILE NAME[e.g.:'VMFG_20180101.H12']
     h12.path = H12files;%FILE PATH['../data/TropoGRIDS/VMF files/VMF grids/'VMFG_20180101.H12'

     %****18 HOUR
     h18.name = h18files;%FILE NAME[e.g.:'VMFG_20180101.H18']
     h18.path = H18files;%FILE PATH['../data/TropoGRIDS/VMF files/VMF grids/'VMFG_20180101.H18'

     %OROGRAPHY FILE
     Oro.name = orofile;%FILE NAME[e.g.:'orography_ell']
     Oro.path = oropath;%FILE PATH['../data/TropoGRIDS/VMF files/VMF grids/'orography_ell'

     %COMBINE ALL INTO ONE STRUCTURE(STRUCT)-VMFgrids
     VMFgrids.H00 = h00;
     VMFgrids.H06 = h06;
     VMFgrids.H12 = h12;
     VMFgrids.H18 = h18;
     VMFgrids.Orography = Oro;
end 

%=========================================END OF SearchVMFgrid.m
%=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=- 
  