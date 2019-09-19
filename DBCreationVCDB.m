%% THIS is for database creation using FPSetExtractionAlt2.m 
%TSec es el tamaño del buffer en segundos
%MaxNumBuffer indica cuantos minutos máximos del video son procesados
warning('off','all');
TSec=10; %buffer of 4 seconds
nscales=4;

nslices=24;
nbB=nslices*2;
% VideoNum=0;
% iFpA=0;
% iFpbV=0;
% iFpV=zeros(nscales,1); %number of total scales for ORB descriptors
% Scales=[.8 1 1.2 1.4]; %(scale factor 1.2 )
% FpA=zeros(1,2+nbB);
% FpbV=zeros(1,902);
% FpV=zeros(nscales,1,34); % 8 scales,1 video, 1 keframe, 34= 32 bytes +2 
% LUT_DB=zeros(nscales,28,1,2); %4 escalas, 28 videos, 1 segments para empezar, 2: indice inicial y final
% VName=cell(500,1);
% 
% % times
% tAFp=0;
% tKf=0;
% tbVFp=0;
% tVFp=0;
% tA=0;
% tbV=0;
% tV=zeros(8,1);
% Nquery=0;


FpathO='D:\Testing time\Video Dataset Real-world\AVI Sorted';
listingD=dir(FpathO);
numofDirs=length(listingD);  %-2 '.' & '..'
indexPerType=zeros(numofDirs-3,2);% diferent types, iniIdx, endIdx ;total segments in scale(2)=no scale
idx1Type=0;

for iDir=2:numofDirs-3 %1 
    %el nombre del sud directorio
	FpathsubD=[FpathO '\' listingD(iDir+2).name '\'];
    listingSubD=dir(FpathsubD);
	numofSubDirs=length(listingSubD);  
    for iSubDir=1:numofSubDirs-2
        FpathD2=[FpathsubD listingSubD(iSubDir+2).name '\'];
        disp(FpathD2);
        listingF=dir(FpathD2);
        numofFiles=length(listingF);
        idx1Type=idx1Type+1; %indice para guardar el inicio del directorio por tipos 
        for iVideo=1:numofFiles-2
            %el nombre del archivo
            wasSegment=false;
            VideoNum=VideoNum+1;
            AuxName=listingF(iVideo+2).name;
            FName=[FpathD2 AuxName];
            VName{VideoNum}=AuxName(1:length(AuxName)-4);
            disp(['Video: ' num2str(VideoNum) ' --> ' FName]);

            %% Read video handler
            videoFReader = vision.VideoFileReader;
            videoFReader.Filename = FName;
            videoFReader.ImageColorSpace = 'Intensity';
            videoFReader.VideoOutputDataType = 'uint8';
            videoFReader.AudioOutputPort = true;

            VideoInfo=info(videoFReader); 
            FRVal= VideoInfo.VideoFrameRate;
            Vsize=VideoInfo.VideoSize;
            if VideoInfo.Audio
                 SR=VideoInfo.AudioSampleRate;
                 ch=VideoInfo.AudioNumChannels;
            end

            bufferLength= round(FRVal*TSec); %segment of TSec seconds
            BufferA=[];
            BufferV=zeros(Vsize(2),Vsize(1),bufferLength,'uint8');

            FNum=0;
            bufferNum=0;
            NumMaxBuffer=180/TSec; %3 minute
            countTSb=0;
            %countTSv=zeros(8,1);

           %% Video segmentation
            while ~isDone(videoFReader) && (bufferNum<=NumMaxBuffer)
              [V,A] = step(videoFReader);
              if ~isempty(V)
                  FNum=FNum+1;
                  if FNum<= bufferLength 
                    BufferV(:,:,FNum)= V; 
                    BufferA=[BufferA;  A]; 
                  else
                    wasSegment=true;
                    bufferNum=bufferNum+1;
                    Nquery=Nquery+1;
                    disp(['Segment number: ' num2str(bufferNum)]);

                    %% Fingerprints set
                    if VideoInfo.Audio
                        [bV,bA,dORB,t1,t2,t3,t4]=FPSetExtractionVCDB(BufferV, FRVal, BufferA,SR,ch,true);
                    else
                        [bV,bA,dORB,t1,t2,t3,t4]=FPSetExtractionVCDB(BufferV, FRVal, [],0,0,true); %No audio
                    end

                    %for timing report
                    tAFp=tAFp+t1;
                    tKf=tKf+t2;
                    tbVFp=tbVFp+t3;
                    tVFp=tVFp+t4;

                    %% Reset counters for video and audio buffers
                    BufferA=[];
                    BufferV=zeros(Vsize(2),Vsize(1),bufferLength,'uint8');
                    FNum=0;

                    %% save into DBs
                    iFpA=iFpA+1;
                    FpA(iFpA,1)=VideoNum;
                    FpA(iFpA,2)=bufferNum; %timestamp increments every Tseconds
                    if VideoInfo.Audio
                        FpA(iFpA,3:2+nbB)=bA;
                    else
                        FpA(iFpA,3:2+nbB)=-1;
                    end

                    iFpbV=iFpbV+1;
                    FpbV(iFpbV,1)=VideoNum;
                    FpbV(iFpbV,2)=bufferNum;
                    FpbV(iFpbV,3:902)=bV;

                    %Using arrays instead of cells

                    for sdo=1:nscales %sdo index for number of scales do=descriptorsORB
                        ido=1;  %number of keyframes
                        if ~isempty(dORB{ido,sdo})
                            [nd,nB]= size(dORB{ido,sdo});
                            LUT_DB(sdo,VideoNum,bufferNum,1)=iFpV(sdo)+1;
                            for jdo=1:nd   %number of descriptors 
                                iFpV(sdo)=iFpV(sdo)+1;
                                FpV(sdo,iFpV(sdo),1)=VideoNum;
                                FpV(sdo,iFpV(sdo),2)=bufferNum;
                                FpV(sdo,iFpV(sdo),3:2+nB)=dORB{ido,sdo}(jdo,:);
                            end
                            LUT_DB(sdo,VideoNum,bufferNum,2)=iFpV(sdo);
                        end
                    end
                  end

              end
            end
            release(videoFReader);
            %save indexes for LUT
            save('DB_VCDB410sec.mat');

        end %For all the videos inside each directory
        %save('DB_VCDB.mat', 'LUT_DB', 'FpA','FpbV','FpV', 'iFpA','iFpbV', 'iFpV','tAFp','tKf','tbVFp','tVFp');
        indexPerType(iSubDir,1)=idx1Type;
        indexPerType(iSubDir,2)=VideoNum;
        idx1Type=VideoNum;
    end %for all subdirectories
        
%To generate lut using video numbers and ts for all directories

tAFp=tAFp/Nquery;
tKf=tKf/Nquery;
tbVFp=tbVFp/Nquery;
tVFp=tVFp/Nquery;

% tA=tA/Nquery; 
% tbV=tbV/Nquery; 
% tV=tV/Nquery; 
save('DB_VCDB4sec.mat');

AnnotationIdx;
%save('DB_VCDB.mat', 'LUT_DB', 'FpA','FpbV','FpV', 'iFpA','iFpbV', 'iFpV','tAFp','tKf','tbVFp','tVFp');