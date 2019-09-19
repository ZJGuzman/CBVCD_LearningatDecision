
%This script search the closest fingerprints pertaining to a video segment
%Previously the fingerprints of "original" video segments were extracted
%   and storaged in a LUT simulating a database (using DBCreationVCDB.m).  
%I used cross validation. Each video in VCDB dataset was considered original while the rest of the videos 
%   in the dataset are the queries to find possible copies. 
%The searching strategy is using K-NN based on distances or similarities of the video fingerprints 
%   (K is at most set to 10 in this example) 
%The decision strategy is based on and adaptation of the traditional
%   Q-Learning algorithm, which generates and update the Q-matrix for each query
%

%% Metrics functions
warning('off', 'all');
%diary('SearchPool4sec.txt');
P=@(tp,fp)(tp/(tp+fp));
R=@(tp,fn)(tp/(tp+fn));
Fb=@(b,P,R)((1+b^2)*((P*R)/((b^2*P)+R)));

%% The directory for audio files
TVideos=VideoNum-1;
[TSegments,~]=size(FpA);
Knn=max(10,round(TVideos/10));
MaxNumKF=max(max(FpbV(:,2)),max(FpA(:,2)));% asignarlo a 100 por ejemplo 

Gamma=1;%0.8;
ThresH=0.1;
%This was for empirical assessment of different thresholds:
%ThresH=[0.03 0.1 0.2 0.3 0.4 0.5 0.6 0.7 0.8 0.9];
%ConfMatrix(1:11)=struct('TP',0,'TN',0,'FP',0,'FN',0, 'TR',0); %TR is the number of groundtruth
ConfMatrix=struct('TP',0,'TN',0,'FP',0,'FN',0); %TR is the number of groundtruth
Metrics_perQuery=struct('P',0,'R',0,'F1',0);

NqueryComparisons=0;

%time 
tA=0;
tbV=0;
tV=zeros(8,1);
tDesicion=0;

wasfp=false; 
nbB=nslices*2;  %same as in FPSetExtractionVCDB.m
nbB1=30*30;
nscales=8;
%% Para cada fingerprint (segmento en cada video de tamaño TSec)

Reward(1:TVideos,1:MaxNumKF)=-1;
selScale=ones(nscales,1);
bufferNum=0;
NumSearch=0;

%This was used for assessment by batches
FpA=FpA_Saved;%(1:3759,:);
FpbV=FpbV_Saved;%(1:3759,:);
FpA_Saved=FpA;
FpbV_Saved=FpbV;
FpV_Saved=FpV;

idxVideo1=0;
for idxSeg=1:TSegments
   
    disp(['Number of segment in DB: ' num2str(idxSeg) ' out of ' num2str(TSegments)]);
    bufferNum=bufferNum+1;
	numQVideo=FpbV_Saved(idxSeg,1);
    numQSeg=FpbV_Saved(idxSeg,2);
    
    %reset tables for q-learning decision
    if numQVideo~=idxVideo1 %reset Reward for each video
        clear('Reward');
        Reward(1:TVideos,1:MaxNumKF)=-1;
        selScale=ones(nscales,1);
        idxVideo1=FpbV_Saved(idxSeg,1);
        %copy the databases without the current video 
		FpbV=kron((FpbV_Saved(:,1)~=idxVideo1),ones(1,902)).*FpbV_Saved;  
        FpA=kron((FpA_Saved(:,1)~=idxVideo1),ones(1,50)).*FpA_Saved;     
        for sc=1:nscales
			FpV1(:,:)=FpV_Saved(sc,:,:);
            FpV1=kron((FpV1(:,1)~=idxVideo1),ones(1,34)).*FpV1;
            FpV(sc,:,:)=FpV1;
		
        end

        bufferNum=1;
        NumSearch=0;
        BestScale=0;
		
    end
    %uses auxiliar arrays to save the best results 
    Raux=Reward;
    Qaux=zeros(TVideos,MaxNumKF);
    ChosenScaleAux=zeros(TVideos,MaxNumKF);
    %FlagMultiMode=false;
		
    %% get query fingerprint and Search them
    wasAudio=(sum(FpA_Saved(idxSeg,3:nbB+2))>0); % if there is not fingerprint:  -1's
    if wasAudio
       disp(['Video: ' num2str(numQVideo) ' Segment number : ' num2str(numQSeg) ]);
       fpAudio=FpA_Saved(idxSeg,3:nbB+2);
       
       % searching audio fingerprint
        tic;
        [Aidx,Ad] = knnsearch(FpA(:,3:nbB+2),fpAudio,'K',Knn,'IncludeTies',true,'NSMethod','exhaustive','Distance','hamming');
        tA=tA+toc;
        %save detected audio keyframes
        [~,knna]=size(Aidx{1,1}); %cuantos son los K-más cercanos
        for iknna=1:knna
            if Ad{1,1}(1,iknna)<1 
                VideoIdx=FpA(Aidx{1,1}(1,iknna),1);
                Videots=FpA(Aidx{1,1}(1,iknna),2);
                Qaux(VideoIdx,Videots)=1-Ad{1,1}(1,iknna); % score in [0:1]; greater is better
                if Raux(VideoIdx,Videots)<0
                    Raux(VideoIdx,Videots)=0; %first reward, =0 if it was found
                end

            end %If match audio
         end %For Knn of A

    end
    fpbVisual=FpbV_Saved(idxSeg,3:nbB1+2);
    %searching global visual
    tic;
    [bVidx,bVd] = knnsearch(FpbV(:,3:nbB1+2),fpbVisual,'K',Knn,'IncludeTies',true,'NSMethod','exhaustive','Distance','hamming');
    tbV=tbV+toc; 
    %save detected visual keyframes
    [~,knnbV]=size(bVidx{1,1});  
    for iknnbV=1:knnbV
        if bVd{1,1}(1,iknnbV)<1 
            VideoIdx=FpbV(bVidx{1,1}(1,iknnbV),1);
            Videots=FpbV(bVidx{1,1}(1,iknnbV),2);
            if VideoIdx>0 && Videots>0
                Qaux(VideoIdx,Videots)=Qaux(VideoIdx,Videots)+ 1-bVd{1,1}(1,iknnbV); % score in [0:2]; greater is better
                if Raux(VideoIdx,Videots)<0
                    Raux(VideoIdx,Videots)=0; %first reward, =0 if it was found
                end
            end
        end %If match audio
     end
    
    %find the Nearest Neighbour with orb
    idxVisual1=LUT_DB(4,numQVideo,numQSeg,1);
    idxVisual2=LUT_DB(4,numQVideo,numQSeg,2);
    clear('fpVisual');
    %Look in the next scale if orb was not find
    if idxVisual1==0 || idxVisual2==0
        idxVisual1=LUT_DB(5,numQVideo,numQSeg,1);
        idxVisual2=LUT_DB(5,numQVideo,numQSeg,2);
        if idxVisual1==0 || idxVisual2==0
            therewasORB=false;
        else
            therewasORB=true;
            if idxVisual1==idxVisual2 %for the assignment of one descriptor
                fpVisual(1,1:32)=FpV_Saved(5,idxVisual1:idxVisual2,3:34); %4=no scale
            else %save the descriptors for the scale and keyframe from the reference database
                fpVisual(:,1:32)=FpV_Saved(5,idxVisual1:idxVisual2,3:34); %4=no scale
            end
            
        end
    else
        therewasORB=true;
        if idxVisual1==idxVisual2 %for the assignment of one descriptor
            fpVisual(1,1:32)=FpV_Saved(4,idxVisual1:idxVisual2,3:34); 
        else %save the descriptors for the scale and keyframe from the reference database
            fpVisual(:,1:32)=FpV_Saved(4,idxVisual1:idxVisual2,3:34); 
        end
        
    end
    
    if therewasORB
        %fpVisual=FpV(4,idxVisual1:idxVisual2,3:34); 
        detScale=0;
        numDesc=idxVisual2-idxVisual1+1;
               
       %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%   
        if numDesc>0   %if there are orb descriptors
           for iv=1:TVideos % for all  the videos
               if max(Raux(iv,:))>=0 % for only previously selected videos
                    ScoreByTs=zeros(2, MaxNumKF);
                    for idxTs=1:MaxNumKF
                        if Raux(iv,idxTs)>=0  %if this frame was previously selected
                            ScoreByScale=zeros(nscales,2);
                            for idxScale=1:nscales
                                if selScale(idxScale)>=1 %if the scale is approved 
                                    idx1=LUT_DB(idxScale,iv,idxTs,1); %initial index
                                    idx2=LUT_DB(idxScale,iv,idxTs,2); %ending index
                                    if idx1~=0 && idx2~=0 %if there are descriptors in this keyframe
                                        clear('ORBofVideoX');
                                        if idx1==idx2 %for the assignment of one descriptor
                                            ORBofVideoX(1,1:34)=FpV(idxScale,idx1:idx2,:);
                                        else %save the descriptors for the scale and keyframe from the reference database
                                            ORBofVideoX(:,1:34)=FpV(idxScale,idx1:idx2,:);
                                        end
                                        [Vidx,Vd] = knnsearch(ORBofVideoX(:,3:34),fpVisual,'K',1,'IncludeTies',false,'NSMethod','exhaustive','Distance','hamming');
                                        numDescRef=idx2-idx1+1;
                                        [C,ia,~]=unique((Vd<1).*Vidx); %C is the vector of non repited descriptors
                                        numDescDet=length(C);
                                        if numDesc<numDescRef
                                            ScoreByScale(idxScale,1)=numDescDet/numDescRef; %distance (normalized)
                                        else
                                            ScoreByScale(idxScale,1)=numDescDet/numDesc; %distance (normalized)
                                        end
                                        ScoreByScale(idxScale,2)=1- mean(Vd(ia)); %score in [0:1] greater is better

                                    end %if there are orb descriptors in this keyframe
                                end %if idxScale was previously selected or is the first keyframe )
                            end %for all the scales in the same keyframe
                            [maxScore, chosenScale]=max(ScoreByScale(:,1).*ScoreByScale(:,2));
                            ScoreByTs(1,idxTs)=maxScore;
                            ScoreByTs(2,idxTs)=chosenScale;
                        end
                    end %for all the keyframes
                    Qaux(iv,:)=Qaux(iv,:)+ ScoreByTs(1,:); % score in [0:2]; greater is better
                    ChosenScaleAux(iv,:)=ScoreByTs(2,:); %all the scales per keyframe
               end  %end if this video was selected
            end   % end for all videos
        else
            detScale=-1;
        end %end if there were orb descriptors
    else
        detScale=-1;
    end
    %% Results analysis and saving
    %update Q
	tic;
    Q=Raux+Gamma*Qaux;
    %To find the best detection scores and to generates the confusion
    %matrix
    
	%Primero genera una matriz lógica con 1's para los valores de la tabla de anotación que deberían ser detectados
    Copies=AnnTable{numQVideo,numQSeg,1};
	numCopies=length(Copies);
    sg1=AnnTable{numQVideo,numQSeg,2}; %the beginning of the coincidence segment
    sg2=AnnTable{numQVideo,numQSeg,3}; %the end of the coincidence segment	
	AnnotationMatrix=false(TVideos,MaxNumKF);
	for c=1:numCopies
		AnnotationMatrix(Copies(c),sg1(c):sg2(c))=true;
	end
	
	%tDesicion=toc+tDesicion;  %partí el tiempo porque en el sig segmento ejecuto la búsqueda 11 veces para probar los %umbrales
    %choose best score in Q
	
	%tic;
	%for th=1:11
		%Q=(Q>=ThresH(th)).*Q;
        Q=(Q>=ThresH).*Q;
		[Valmax, idxmax]=max(Q(:));
		if Valmax>0
			Module=mod(idxmax,TVideos);
			Division=floor(idxmax/TVideos);
			if Module==0
				detVideo=TVideos;
				detTs=Division;
			else
				detVideo=Module;
				detTs=Division+1; 
			end
			if detScale>=0
				detScale=ChosenScaleAux(detVideo,detTs);
				selScale(detScale)=selScale(detScale)+1;
			end
		
			if AnnotationMatrix(detVideo,detTs)
				%ConfMatrix(th).TP=ConfMatrix(th).TP+1;
                ConfMatrix.TP=ConfMatrix.TP+1;
				disp(['TP >> Best similarity: Query ' num2str(numQVideo) ' Segment ' num2str(bufferNum) ...
				' Detected Video ' num2str(detVideo) ' Segment ' num2str(detTs) ' Scale ' num2str(detScale)]);
            else
                %ConfMatrix(th).FP=ConfMatrix(th).FP+1;
                ConfMatrix.FP=ConfMatrix.FP+1;
				disp(['FP << Best similarity: Query ' num2str(numQVideo) ' Segment ' num2str(bufferNum) ...
				' Detected Video ' num2str(detVideo) ' Segment ' num2str(detTs) ' Scale ' num2str(detScale)]);
			end
		else %if there was no detection over the threshold
            detVideo=0;
            detTs=0;
			if max(AnnotationMatrix(:))>0 %#ok<ALIGN>
                disp('FN >> No similarity');
                %ConfMatrix(th).FN=ConfMatrix(th).FN+1;
                ConfMatrix.FN=ConfMatrix.FN+1;
            else
                disp('TN >> No reference in database');
                %ConfMatrix(th).TN=ConfMatrix(th).TN+1;
                ConfMatrix.TN=ConfMatrix.TN+1;
            end
		end
	%end
	
    %update Reward
    if detVideo>0 && detTs>0
        if numQSeg==1 
            Reward(detVideo,detTs)=Q(detVideo,detTs)/2;
        else
            Reward(detVideo,detTs)=(Reward(detVideo,detTs)+Q(detVideo,detTs))/2;
        end
   
        Reward(detVideo,detTs+1:MaxNumKF)=Reward(detVideo,detTs);
    end
%     for ii=1:TVideos %este ciclo replica el mayor valor de la tabla Reward a los siguientes segmentos del mismo video detectado previamente
%         [value,idxM]=max(Reward(ii,:));
%         Reward(ii,idxM:MaxNumKF)=value;
%     end

%     %normalize Reward and Q
%     Reward=Reward/max(Reward(:));
%     Q=Q/max(Q(:));

	tDesicion=toc+tDesicion; 
	%tDesicion=(toc/11)+tDesicion; 
    %save('SearchDBFps4sec.mat');
end

%% calculate the detection metrics

    ConfMatrix.TP=ConfMatrix.TP+TVideos;
	Metrics_perQuery.P=P(ConfMatrix.TP,ConfMatrix.FP);
	Metrics_perQuery.R=P(ConfMatrix.TP,ConfMatrix.FN);
	Metrics_perQuery.F1=Fb(1,Metrics_perQuery.P, Metrics_perQuery.R);

% for t=1:11
% 	Metrics_perQuery(t).P=P(ConfMatrix(t).TP,ConfMatrix(t).FP);
% 	Metrics_perQuery(t).R=P(ConfMatrix(t).TP,ConfMatrix(t).FN);
% 	Metrics_perQuery(t).F1=Fb(1,Metrics_perQuery(t).P, Metrics_perQuery(t).R);
% end

tA=tA/TSegments; 
tbV=tbV/TSegments; 
tV=tV/TSegments; 
tDesicion=tDesicion/TSegments;

disp(['Average execution time for searching Audio fingerprint: ' num2str(tA)]);
disp(['Average execution time for searching Visual (GLOBAL) fingerprint: ' num2str(tbV)]);
disp(['Average execution time for decision: ' num2str(tDesicion)]);

%save('SearchDBFps4sec.mat');
save('ConfMatrix_4sTh0_3.mat', 'ConfMatrix','Metrics_perQuery');
diary off;

