%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% Library load %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

DB = importdata('materialDB_CMOS.csv',',');

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%% Here is a patch to be able to load the tables in Matlab AND Octave %%%%%%
% Matlab see the header in multiple cells while Octave see the header in one cell only
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

M='Si';
if length(DB.textdata(1,:))==1  %% Octave data load

    DB.textdata{1,1}=[DB.textdata{1,1} ',']; % patch, add a comma "," at the end
    idxM=strfind(DB.textdata{1,1},',');
    idx=strfind(DB.textdata{1,1},[',' M ',']);
    idxM=find(idxM==idx);
    
    Si = DB.data(:,idxM)';
else  %% Matlab data load

    for i=1:length(DB.textdata(1,:))
      idx=strcmp(DB.textdata{1,i},M);
      if idx==1
        Si=DB.data(:,i-1)';
        % break % removing the break makes it slower but more compatible between Matlab and Octave
      end
    end 
    
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

M='Oxide';
if length(DB.textdata(1,:))==1  %% Octave data load

    DB.textdata{1,1}=[DB.textdata{1,1} ',']; % patch, add a comma "," at the end
    idxM=strfind(DB.textdata{1,1},',');
    idx=strfind(DB.textdata{1,1},[',' M ',']);
    idxM=find(idxM==idx);
    
    Oxide = DB.data(:,idxM)';
else  %% Matlab data load

    for i=1:length(DB.textdata(1,:))
      idx=strcmp(DB.textdata{1,i},M);
      if idx==1
        Oxide=DB.data(:,i-1)';
        % break % removing the break makes it slower but more compatible between Matlab and Octave
      end
    end 
    
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

M='Poly';
if length(DB.textdata(1,:))==1  %% Octave data load

    DB.textdata{1,1}=[DB.textdata{1,1} ',']; % patch, add a comma "," at the end
    idxM=strfind(DB.textdata{1,1},',');
    idx=strfind(DB.textdata{1,1},[',' M ',']);
    idxM=find(idxM==idx);
    
    Poly = DB.data(:,idxM)';
else  %% Matlab data load

    for i=1:length(DB.textdata(1,:))
      idx=strcmp(DB.textdata{1,i},M);
      if idx==1
        Poly=DB.data(:,i-1)';
        % break % removing the break makes it slower but more compatible between Matlab and Octave
      end
    end 
    
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%