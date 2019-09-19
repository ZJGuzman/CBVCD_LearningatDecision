function [bF,bA,FpORB,t1,t2,t3,t4]=FPSetExtractionVCDB(f,FRVal, a,SR,ch,FlagScales)
        nslices=24;
        nSR=8000;
        %% THIS extracts features from frames 'f' and audio 'a' 
        %       and generates their fingerprints 
        if isempty(a)
            bA=[];
            t1=0;
        else
            tic;
            bA=BinSpectrogram(AudioStandardize(a));
            t1=toc;
        end
        
        tic;
        kF=frameSelection(f);
       
        t2=toc;
         %imshow(kF);
        tic;
        %bF=BinHashI(kF);
        bF=ImageToVector(kF);
        t3=toc;
        
        tic;
        FpORB= GetORB(kF); %FAST+PATCH+SCALING+ORB
        t4=toc;
        
               
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    function a=AudioStandardize(a)
        %% THIS converts audio 'a' into mono channel and downsample it to 8KHz 
        %nSR=8000;
        
        %Mono channel conversion
        if ch>1
            a = sum(a,2);
            a = a/ch; 
        end
        % Resampling
        if (nSR ~= SR)
            a = resample (a,nSR,SR,100); 
            %resample applies an antialiasing FIR lowpass filter to Amono and compensates for the delay introduced by the filter.
            %uses an antialiasing filter of order 2 × 100 × max(SR,Fs).
        end
    end


    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    function [BSa]=BinSpectrogram(a)
        %% THIS generates and binarize the spectrogram of 'a' using OUALI2016
        %nslices=24;
        fm=1.6;
        
        Freqmin=300;    %minimum frequency required
        Freqmax=3000;   %maximum frequency required
        Freqbins=257;   %number of total frequency bins
        Fs=nSR;       %Fs is the default sample frequency 
        wlen=Fs*0.096;  %window length measured in seconds
        Noverlap=wlen-(Fs*0.048); %overlaping factor = 50% of the window length
        inci=(Freqmax-Freqmin+2)/Freqbins; %to generate only 257 frequencies
        j=1;
        
        freqV=zeros(Freqbins,1); %the frequencies vector for filtering
        for i=Freqmin:inci:Freqmax
            freqV(j)=i*2/Fs; %normalized frequency rad/sample
            j=j+1;
        end

        S=abs(spectrogram(a,blackman(wlen),Noverlap, freqV,'yaxis'));
        
        %normalize S
        s=S-min(S(:));
        s=s./max(S(:));
        
        [w,h]=size(s);
        finc=round(w/nslices);  %size of interval x
        tinc=round(h/nslices);  %size of interval y
        FpGH=zeros(1,nslices);
        FpGV=zeros(1,nslices);
        
        GMHash = im2bw(s,fm*mean(s(:)));
       
        for sliceidx=1:nslices
           FpGH(sliceidx)=sum(sum(GMHash((sliceidx-1)*finc+1:min(sliceidx*finc,w),:)));
           FpGV(sliceidx)=sum(sum(GMHash(:,(sliceidx-1)*tinc+1:min(sliceidx*tinc,h))));
        end

        BSa=[FpGH FpGV];
    
    end
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    function ImgV=ImageToVector(k)
        k1=imresize(k,[30 30]);
        ImgV=k1(:)';
    end
%    function bI=BinHashI(k)
%         %% THIS converts an image s into a binary string of nslices*2 bits 
%         %based on the global mean as binarization threshold
%         %nslices=24;
%         [h,w]=size(k);
%         finc=round(w/nslices);  %size of interval x
%         tinc=round(h/nslices);  %size of interval y
%         FpGH=zeros(1,nslices);
%         FpGV=zeros(1,nslices);
%         GM1=mean(k(:));
%         GMHash=(k>GM1);
%         for sliceidx=1:nslices
%            FpGH(sliceidx)=sum(sum(GMHash(:,(sliceidx-1)*finc+1:min(sliceidx*finc,w))));
%            FpGV(sliceidx)=sum(sum(GMHash((sliceidx-1)*tinc+1:min(sliceidx*tinc,h),:)));
%         end
%         bI=[FpGH FpGV];
%         
%         
%     end

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    function f1=frameSelection(f)
        %% THIS selects and folds the representative frames from the frames set 'f' 
        %cada 2 segundos  un shot, entonces promediar los shots
       [~,w, nf]=size(f);
        Shot=round(FRVal*2);
        halfShot=round(FRVal/2);
        iniShot=halfShot+Shot;
        midw=w/2;
        wei=0.5;
        Kf1=f(:,:,halfShot);
        for j=iniShot:Shot:nf
            if mean(Kf1(:))<1
                Kf=f(:,:,j);
            else
                Kf=wei*f(:,:,j)+(1-wei)*Kf1;
            end
            Kf1=Kf;
        end
%          for j=2:nf
%             Kf=wei*f(:,:,j-1)+(1-wei)*f(:,:,j);
%             f(:,:,j)=Kf;
%         end
        f1(:,:)=0.5*Kf(:,1:midw)+0.5*(fliplr(Kf(:,midw+1:w)));
        
    end  

    function FpORB= GetORB(f)
        %% THIS extracts FAST corners, divides the image into 16 block, 
        % in each block, if there are corners, then extract ORB at
        % differente scales
        [h,w]=size(f);
        patchSizeX=floor(w/4);
        patchSizeY=floor(h/4);
        Scales=[.57 .69 .8 1 1.2 1.4 1.7 2]; %(scale factor 1.2 )
        FpORB=cell(1,8);
        if FlagScales
            iniS=1;
            endS=8;
        else    %No scaling
            iniS=4;
            endS=4;
        end

        %figure; imshow(f);
        Corners=cv.FAST(f(:,:)); 
        numCorners=length(Corners);
        %hold on;
        for iy=1:patchSizeY:h%-patchSizeY
            for ix=1:patchSizeX:w%-patchSizeX
                clear('patchImage');
                iy1=iy+patchSizeY-1;
                ix1=ix+patchSizeX-1;
                patchImage(:,:)=f(iy:iy1, ix:ix1);

                flagCorners=false;
                c=1;
                while (~flagCorners && c<=numCorners)
                     if inpolygon(Corners(c).pt(1),Corners(c).pt(2),[ix ix1],[iy iy1])
                         flagCorners=true;
                         %plot(Corners(c).pt(1),Corners(c).pt(2), 'g+');
                     else
                         c=c+1;
                     end
                end
                if flagCorners

                    for iscales=iniS:endS
                        %scales
                        I=imresize(patchImage,Scales(iscales));
                        [~,ORBdescriptor]=cv.ORB(I);
                        if ~isempty(ORBdescriptor)
                            FpORB{iscales}=[FpORB{iscales}; ORBdescriptor];
                        end
                    end
                    %rectangle('Position',[ix iy patchSizeX patchSizeY],'EdgeColor', [1,0,0]);
                    %plot(Corners(c).pt(1),Corners(c).pt(2), 'g+');
                end
            end

        end
        
    end
end

