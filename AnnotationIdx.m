%TSec=10;
TVideos=VideoNum-1;
numSegAnn=1+60/TSec;
%TSegments=max(FpA(:,2));%+5; %+5 para anchar el arreglo y no falte espacio
AnnTable=cell(TVideos,numSegAnn,3);% Numero de video 28,numero de keyframe máximo 43, 
%encontrado en video numero y keyframe inicial y final  
FpathAnnotation='D:\Testing time\Video Dataset Real-world\annotation\';
listingA=dir(FpathAnnotation);
numofDirs=length(listingA);
[numV, ~]=size(VName); %(1:VideoNum) tiene los nombres de los videos numerados
for iDir=1:numofDirs-2
    %el nombre del directorio
    FAname=[FpathAnnotation listingA(iDir+2).name];
    fid=fopen(FAname, 'r');
	Annotation=textscan(fid, '%s %s %s %s %s %s', 'delimiter', ',');
    [nV,~]=size(Annotation{1,1});
    for iA=1:nV
        nlength=length(Annotation{1,1}{iA});
        NameVideoQuery=Annotation{1,1}{iA}(1:nlength-4);
        nlength=length(Annotation{1,2}{iA});
        NameVideoRef=Annotation{1,2}{iA}(1:nlength-4);
		
		[~,IdxV1]=max(strcmp(VName, NameVideoQuery)); 
        [~,IdxV2]=max(strcmp(VName, NameVideoRef)); %obtiene el numero de video de la columna 1 y 2
		T=Annotation{1,3}{iA};	
		Ki1=hms2kf(T,TSec);	%Initial time for query video
        T=Annotation{1,4}{iA};	
		Ke1=hms2kf(T,TSec);	%ending time for query video
		T=Annotation{1,5}{iA};	
		Ki2=hms2kf(T,TSec);	%Initial time for reference video
        T=Annotation{1,6}{iA};	
		Ke2=hms2kf(T,TSec);	%ending time for reference video	
		
        %para guardar las correspondencias de los videos el video 1 con el
        %video2 y a la vez el video 2 con el video 1
        if Ki1<=numSegAnn && Ke1<=numSegAnn && Ki2<=numSegAnn && Ke2<=numSegAnn
            for k=Ki1:Ke1
                AnnTable{IdxV1,k,1}=[AnnTable{IdxV1,k,1} IdxV2];
                AnnTable{IdxV1,k,2}=[AnnTable{IdxV1,k,2} Ki2];
                AnnTable{IdxV1,k,3}=[AnnTable{IdxV1,k,3} Ke2];
            end
            for k=Ki2:Ke2
                AnnTable{IdxV2,k,1}=[AnnTable{IdxV2,k,1} IdxV1];
                AnnTable{IdxV2,k,2}=[AnnTable{IdxV2,k,2} Ki1];
                AnnTable{IdxV2,k,3}=[AnnTable{IdxV2,k,3} Ke1];
            end
        end
        
				
		
    end
    fclose(fid);
    
   
end

        
                
                
                
