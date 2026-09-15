clear all
close all
files={  'mouse260--sal-';     'mouse260-rev--sal-'; 
     'mouse261--sal-';     'mouse261-rev--sal-';
     'mouse262--sal-';  'mouse262-rev--sal-';
     'mouse263--sal-';     'mouse263-rev--sal-';
     'mouse264--sal-';     'mouse264-rev--sal-';  
     'mouse256-rev--sal-';
    'mouse257--sal-';'mouse257-rev--sal-';
    'mouse258--sal-';   'mouse258-rev--sal-'; 
    'mouse259--sal-';
    'mouse295--sal-';    'mouse295-rev--sal-';
    };%  

   mss=[1 1 2 2 3 3 4 4 5 5 6 7 7 8 8 9 10 10 ];
   sxx=[1 1 1 1 1 1 1 1 1 1 1 2 2 1 1 2 1  1  ]; %dosal=1; clear dodcz;

%  files={  'mouse260--dcz-';     'mouse260-rev--dcz-'; 
%       'mouse261--dcz-';     'mouse261-rev--dcz-';
%       'mouse262--dcz-';
%       'mouse263--dcz-';     'mouse263-rev--dcz-';
%       'mouse264--dcz-';     'mouse264-rev--dcz-';  
%       'mouse256-rev--dcz-';
%      'mouse257--dcz-';
%         'mouse258-rev--dcz-'; 
%      'mouse259-rev--dcz-';
%      'mouse295--dcz-';    'mouse295-rev--dcz-';
%     };% 
% mss=[1 1 2 2 3 4 4 5 5 6 7 8 9 10 10 ];
% sxx=[1 1 1 1 1 1 1 1 1 1 2 1 2 1  1  ]; %dodcz=1; clear dosal; 

cellcount=0; cellscontrol=[]; cellsOB=[]; cellscontrolSTD=[]; cellsOBSTD=[]; mouseid=[];sexid=[]; FOVid=[];
basedefMIN=1:14; basedefMEAN=10:14; stimdefON=15:28; stimdefOFF=31:39;%these are for the standardizing step, which determines whether cells are "responsive"
exclloco=1; limtrials=6; rejectingdoubles=1;

for fil=1:size(files,1)

    load(strcat('CA_workspace_SORTED2_',files{fil}))
    load(strcat('reject_cells_',files{fil}));
    for rejectingdoubles=rejectingdoubles;
        if rejectingdoubles==1;
            data1=data1(reject_cells==0,:);
            dfofanalyse=dfofanalyse(reject_cells==0,:);
            dfofanalyse2=dfofanalyse2(reject_cells==0,:);
            vec=vectrials{1};
            vecdat=vec{1,1};
            vecdat=vecdat(reject_cells==0,:,:);
            vec{1,1}=vecdat;
            vectrials{1}=vec;

            vec=vectrials{2};
            vecdat=vec{1,1};
            vecdat=vecdat(reject_cells==0,:,:);
            vec{1,1}=vecdat;
            vectrials{2}=vec;
        end
    end

    samp=1/framerate;
    timeaxis=-(samp*(framerate/2))+samp:samp:(samp*(framerate*1))+samp;

    %do the control run
    vec=vectrials{1};
    vecdat=vec{1,1};
    vecstim=vec{1,2};
    ords=vec{1,3};
    locs=vec{1,4};

    %limits to the 5th or 4th in sequence and eliminates repeats in the control run
    vecstim(1)=0;
    vecstim1=vecstim;
    for s=2:size(vecstim,2);
        if and(ords(s)>1,vecstim(s-1)==vecstim(s))
            vecstim1(s)=0;
        end
        if ords(s)<4
            vecstim1(s)=0;
        end
    end
    vecstim=vecstim1;
    
    %eliminates locomotion, if wanted
    for exclloco=exclloco;
        if exclloco==1;
            for s=1:size(vecstim,2);
                if locs(1,s)>0;
                    vecstim(1,s)=0;
                end
            end
        end
    end

   for limtrials=limtrials; 
       if limtrials>0;
           counts=zeros(1,8);
           for tr=1:size(vecstim,2);
               tr1=(size(vecstim,2)+1)-tr;
               aa=vecstim(1,tr1); 
               if aa>0;
               counts(1,aa)=counts(1,aa)+1;
               if counts(1,aa)>limtrials;
                   vecstim(1,tr1)=0;
               end
               end
           end
       end
   end

    numtrials(fil,1)=sum(vecstim==1);numtrials(fil,2)=sum(vecstim==2);numtrials(fil,3)=sum(vecstim==5);
  
    dfofanalyse1=dfofanalyse; 

    for cell=1:size(vecdat,1);
        celltmp=[]; celltrialtmp=[];
        for ori=1:3;
            celltmp(:,ori)=squeeze(mean(vecdat(cell,:,vecstim==trtype(ori)),3));
        end
        ssss=vecdat(cell,:,:); ssss=ssss(:); sss=ssss(ssss>-0); ss=sss(sss<prctile(sss,50));
        for ori=1:3;
            cellscontrolSTD(cellcount+cell,:,ori)=celltmp(:,ori)./std(ss(:));   
        end

        for ori=1:3;
            cellscontrol(cellcount+cell,:,ori)=celltmp(:,ori)-mean(celltmp(basedefMEAN,ori));
        end

        %for trialwise stats
        for ori=1:3;
            celltrialON{cellcount+cell,ori}=squeeze(mean(vecdat(cell,stimdefON,vecstim==trtype(ori)),2))./std(ss(:));
            celltrialOFF{cellcount+cell,ori}=squeeze(mean(vecdat(cell,stimdefOFF,vecstim==trtype(ori)),2))./std(ss(:));
        end
    end
    %make "inc" variable to figure out which cells show a true response
    for cella=1:cell;
        inc_cntON(cellcount+cella,:)=squeeze(mean(cellscontrolSTD(cellcount+cella,stimdefON,:),2)-mean(cellscontrolSTD(cellcount+cella,basedefMEAN,:),2));
        inc_cntOFF(cellcount+cella,:)=squeeze(mean(cellscontrolSTD(cellcount+cella,stimdefOFF,:),2)-mean(cellscontrolSTD(cellcount+cella,basedefMEAN,:),2));
      
    end

    %do the oddball run
    vec=vectrials{2};
    vecdat=vec{1,1};
    vecstim=vec{1,2};
    ords=vec{1,3};
    locs=vec{1,4}; 
    
    %part one
    for limtrials=limtrials; %dont need to do this for the 6 and8, but only 4
       if limtrials>0;
           %first, limit 4s to the sequences right before the 6 and 8
           vecstim1=vecstim;
           for tr=1:size(vecstim,2)-5;
                if vecstim(tr)==4
                    if vecstim(tr+5)==4;
                        vecstim1(tr)=0;
                    end
                end
           end
       end
    end
    vecstim=vecstim1;
     %eliminates locomotion, if wanted
    for exclloco=exclloco;
        if exclloco==1;
            for s=1:size(vecstim,2);
                if locs(1,s)>0;
                    vecstim(1,s)=0;
                end
            end
        end
    end
    %part 2 limtrials (after locomotion exclude. need this
    for limtrials=limtrials; %dont need to do this for the 6 and8, but only 4
       if limtrials>0;
           %vecstim=vecstim1;
           %then, limit to last 6 occurrances
           counts=zeros(1,8);
           for tr=1:size(vecstim,2);
               tr1=(size(vecstim,2)+1)-tr;
               aa=vecstim(1,tr1); 
               if aa>0
               counts(1,aa)=counts(1,aa)+1;
               if counts(1,aa)>6;
                   vecstim(1,tr1)=0;
               end
               end
           end
       end
   end
    trs=[4 6 8];
    numtrials(fil,4)=sum(vecstim==4);numtrials(fil,5)=sum(vecstim==6);numtrials(fil,6)=sum(vecstim==8);


   

    for cell=1:size(vecdat,1);
        celltmp=[];
        for ori=1
            celltmp(:,1)=squeeze(mean(vecdat(cell,:,and(ords>2,vecstim==2)),3));
        end

        ssss=vecdat(cell,:,:); ssss=ssss(:); sss=ssss(ssss>-0); ss=sss(sss<prctile(sss,50));
        for ori=1;
            cellsOBSTD(cellcount+cell,:,ori)=celltmp(:,ori)./std(ss(:));
        end

        for ori=1
            cellsOB(cellcount+cell,:,1)=celltmp(:,1)-mean(celltmp(basedefMEAN,ori));
        end

    end
    for cell=1:size(vecdat,1);
        celltmp=[];
        for ori=1:3
            celltmp(:,ori)=squeeze(mean(vecdat(cell,:,vecstim==trs(ori)),3));
        end

         ssss=vecdat(cell,:,:); sss=ssss(ssss>-0); ss=sss(sss<prctile(sss,50));
        for ori=1:3;
            cellsOBSTD(cellcount+cell,:,ori+1)=celltmp(:,ori)./std(ss(:));
        end

        for ori=1:3
            cellsOB(cellcount+cell,:,ori+1)=celltmp(:,ori)-mean(celltmp(basedefMEAN,ori));
        end

        %for trialwise stats
        for ori=1:3;
            celltrialON{cellcount+cell,ori+3}=squeeze(mean(vecdat(cell,stimdefON,vecstim==trs(ori)),2))./std(ss(:));
            celltrialOFF{cellcount+cell,ori+3}=squeeze(mean(vecdat(cell,stimdefOFF,vecstim==trs(ori)),2))./std(ss(:));
        end

    end

    %make "inc" variable to figure out which cells show a true response
    for cella=1:cell;
        inc_OBON(cellcount+cella,:)=squeeze(mean(cellsOBSTD(cellcount+cella,stimdefON,:),2)-mean(cellsOBSTD(cellcount+cella,basedefMEAN,:),2));
        inc_OBOFF(cellcount+cella,:)=squeeze(mean(cellsOBSTD(cellcount+cella,stimdefOFF,:),2)-mean(cellsOBSTD(cellcount+cella,basedefMEAN,:),2));
    
        
    end

    cellcount=cell+cellcount;
    mouseid=vertcat(mouseid,zeros(cell,1)+mss(fil));
    sexid=vertcat(sexid,zeros(cell,1)+sxx(fil));
    FOVid=vertcat(FOVid,zeros(cell,1)+fil);

end

%deterimining which cells to include 
stdcut1=1.97; %p<.05 %this added subcriterion adds only 3 total neurons
stdcut2=2.6; %p<.01 %making them both .05 has same statistical results, but waters down most effects.
doLMEstats=1; outliercut=23.23; %5 stds
incRDNT=zeros(1,cellcount); incDEV1=zeros(1,cellcount); incDEV2=zeros(1,cellcount);
scale_to_cntMAX=0; mousewise=0;%only for plotting purposes.
do_onsets=1;
do_offsets=0; viewtruerdnt=0;

for cell=1:cellcount;
    if do_onsets==1;
        if mouseid(cell)<6; %gcamp8m doesnt have as large transients, so teh cut is slightly more liberal
            if and(or(and(inc_cntON(cell,1)>stdcut1, inc_OBON(cell,3)>stdcut1),or(inc_cntON(cell,1)>stdcut2, inc_OBON(cell,3)>stdcut2)),and(abs(inc_cntON(cell,1))<outliercut,inc_OBON(cell,3)<outliercut))
                incRDNT(1,cell)=1;
            end
            if and(or(and(inc_cntON(cell,2)>stdcut1, inc_OBON(cell,2)>stdcut1),or(inc_cntON(cell,2)>stdcut2, inc_OBON(cell,2)>stdcut2)),and(abs(inc_cntON(cell,2))<outliercut,inc_OBON(cell,2)<outliercut))
                incDEV1(1,cell)=1;
            end
            if and(or(and(inc_cntON(cell,3)>stdcut1, inc_OBON(cell,4)>stdcut1),or(inc_cntON(cell,3)>stdcut2, inc_OBON(cell,4)>stdcut2)),and(abs(inc_cntON(cell,3))<outliercut,inc_OBON(cell,4)<outliercut))
                incDEV2(1,cell)=1;
            end
            
        else
             if and(or(and(inc_cntON(cell,1)>stdcut1, inc_OBON(cell,3)>stdcut1),or(inc_cntON(cell,1)>stdcut1, inc_OBON(cell,3)>stdcut1)),and(abs(inc_cntON(cell,1))<outliercut,inc_OBON(cell,3)<outliercut))
                incRDNT(1,cell)=1;
            end
            if and(or(and(inc_cntON(cell,2)>stdcut1, inc_OBON(cell,2)>stdcut1),or(inc_cntON(cell,2)>stdcut1, inc_OBON(cell,2)>stdcut1)),and(abs(inc_cntON(cell,2))<outliercut,inc_OBON(cell,2)<outliercut))
                incDEV1(1,cell)=1;
            end
            if and(or(and(inc_cntON(cell,3)>stdcut1, inc_OBON(cell,4)>stdcut1),or(inc_cntON(cell,3)>stdcut1, inc_OBON(cell,4)>stdcut1)),and(abs(inc_cntON(cell,3))<outliercut,inc_OBON(cell,4)<outliercut))
                incDEV2(1,cell)=1;
            end
           
        end
    end
     if do_offsets==1;
        if mouseid(cell)<6; 
            if and(or(and(inc_cntOFF(cell,1)>stdcut1, inc_OBOFF(cell,3)>stdcut1),or(inc_cntOFF(cell,1)>stdcut2, inc_OBOFF(cell,3)>stdcut2)),and(inc_cntOFF(cell,1)<outliercut,inc_OBOFF(cell,1)<outliercut))
                incRDNT(1,cell)=1;
            end
            if and(or(and(inc_cntOFF(cell,2)>stdcut1, inc_OBOFF(cell,2)>stdcut1),or(inc_cntOFF(cell,2)>stdcut2, inc_OBOFF(cell,2)>stdcut2)),and(inc_cntOFF(cell,2)<outliercut,inc_OBOFF(cell,2)<outliercut))
                incDEV1(1,cell)=1;
            end
            if and(or(and(inc_cntOFF(cell,3)>stdcut1, inc_OBOFF(cell,4)>stdcut1),or(inc_cntOFF(cell,3)>stdcut2, inc_OBOFF(cell,4)>stdcut2)),and(inc_cntOFF(cell,3)<outliercut,inc_OBOFF(cell,4)<outliercut))
                incDEV2(1,cell)=1;
            end
         
        else
             if and(or(and(inc_cntOFF(cell,1)>stdcut1, inc_OBOFF(cell,3)>stdcut1),or(inc_cntOFF(cell,1)>stdcut1, inc_OBOFF(cell,3)>stdcut1)),and(inc_cntOFF(cell,1)<outliercut,inc_OBOFF(cell,1)<outliercut))
                incRDNT(1,cell)=1;
            end
            if and(or(and(inc_cntOFF(cell,2)>stdcut1, inc_OBOFF(cell,2)>stdcut1),or(inc_cntOFF(cell,2)>stdcut1, inc_OBOFF(cell,2)>stdcut1)),and(inc_cntOFF(cell,2)<outliercut,inc_OBOFF(cell,2)<outliercut))
                incDEV1(1,cell)=1;
            end
            if and(or(and(inc_cntOFF(cell,3)>stdcut1, inc_OBOFF(cell,4)>stdcut1),or(inc_cntOFF(cell,3)>stdcut1, inc_OBOFF(cell,4)>stdcut1)),and(inc_cntOFF(cell,3)<outliercut,inc_OBOFF(cell,4)<outliercut))
                incDEV2(1,cell)=1;
            end
           
        end
    end
end

disp(strcat('proportion of RDNT responsive:',num2str(sum(incRDNT)/cellcount)))
disp(sum(incRDNT))
disp(strcat('proportion of DEV1 responsive:',num2str(sum(incDEV1)/cellcount)))
disp(sum(incDEV1))
disp(strcat('proportion of DEV2 responsive:',num2str(sum(incDEV2)/cellcount)))
disp(sum(incDEV2))

stimON=15:28; stimOFF=31:39; %early on is 15:20; lateon is 21:25; fullon is 15:28
for usestandardized=1;
    if usestandardized==1;
        cellsOB=cellsOBSTD;
        cellscontrol=cellscontrolSTD;
    end
end

for localDDanalysis=1;
    OB=cellsOB(incDEV1==1,:,2);
    CT=cellscontrol(incDEV1==1,:,2);
    for cell=1:size(OB,1) %BASE CORRECTION
        
            OB(cell,:)=OB(cell,:)-mean(OB(cell,basedefMEAN));%
            CT(cell,:)=CT(cell,:)-mean(CT(cell,basedefMEAN));%
       
    end
   
            stdzvals=ones(size(CT,1),1);
       
    clear OBtmp CTtmp
    for cell=1:size(OB,1); OBtmp(cell,:)=OB(cell,:)./stdzvals(cell,1);CTtmp(cell,:)=CT(cell,:)./stdzvals(cell,1);end

    OBmn=mean(OBtmp,1); CTmn=mean(CTtmp,1);
    stdev=std(OBtmp-CTtmp,0,1)./sqrt(size(OB,1)-1);
    stdev1=std(OBtmp)./sqrt(size(OB,1)-1);  stdev2=std(CTtmp)./sqrt(size(CT,1)-1);
    figure; shadedErrorBar(timeaxis,CTmn,stdev2,'lineProps','k'); hold on;
    shadedErrorBar(timeaxis,OBmn,stdev1,'lineProps','r');hold on;
     xlim([-.2 1]); title ('local deviant, global standard', 'FontSize', 16, 'FontWeight', 'bold');
    ylabel('estFiringRate','FontSize',16,'FontWeight','bold'); 
    xlabel('time (sec)','FontSize',16,'FontWeight','bold'); set(gca,'FontSize',16,'FontWeight','bold');
    set(gcf,'Color','w');
     
    figure(3); subplot(1,2,1); scatter(mean(CT(:,stimON),2),mean(OB(:,stimON),2),'fill'); title('localDEV globalRDNT ONSET responses'); hold on; 
   clear x;  tmp1=mean(OB(:,stimON),2); tmp2=mean(CT(:,stimON),2); aa=horzcat(tmp1,tmp2); x(1)=min(aa(:))*1.02; x(2)=max(aa(:))*1.02;
    y = x; plot(x,y);[h,p,ci,stats]=ttest(tmp1,tmp2); disp(strcat('local DD global rdnt. t=',num2str(stats.tstat),'. p=',num2str(p)))
     ylabel(strcat('deviant: p=',num2str(p)),'FontSize',16,'FontWeight','bold'); 
    xlabel('control','FontSize',16,'FontWeight','bold');
    figure(4); subplot(1,4,1); errorbar_groups([mean(tmp2) mean(tmp1)]',[std(tmp2)./sqrt(size(tmp1,1)-1) std(tmp1)./sqrt(size(tmp1,1)-1)]', 'bar_colors',[.5 .5 .5; 1 0 0]);
        ylabel('eFR','FontSize',16,'FontWeight','bold'); xlabel('onsets LDGr','FontSize',16,'FontWeight','bold'); set(gca,'FontSize',16,'FontWeight','bold');
        set(gcf,'Color','w');
        props1=((tmp1-tmp2))./((abs(tmp1)+abs(tmp2))./2); 
        figure(4); subplot(1,4,2); errorbar_groups([mean(props1)]',[std(props1)./sqrt(size(tmp1,1)-1)]', 'bar_colors',[1 0 0]); 
       

     for doLMEstats=doLMEstats;
         if doLMEstats==1;
             tmpmouse = mouseid(incDEV1==1);
             tmpsex   = sexid(incDEV1==1);
             tmpFOV   = FOVid(incDEV1==1);
             n = size(tmp1, 1);

             CellID = (1:n)';   % unique ID per neuron, since tmp1/tmp2 are row-matched

             control_vals = tmp1; deviant_vals = tmp2;
             mice = tmpmouse; sexes = tmpsex; fovs = tmpFOV;

             Value = [control_vals; deviant_vals];
             Condition = [repmat({'control'}, n, 1); repmat({'deviant'}, n, 1)];
             Mouse = [mice; mice];
             Sex = [sexes; sexes];
             FOV = [fovs; fovs];
             Cell = [CellID; CellID];   % same cell ID repeated for both conditions

             T = table(Value, categorical(Condition), categorical(Mouse), categorical(Sex), ...
                 categorical(FOV), categorical(Cell), ...
                 'VariableNames', {'Value', 'Condition', 'Mouse', 'Sex', 'FOVid', 'CellID'});

             lme = fitlme(T, 'Value ~ Condition + (1|Mouse) + (1|Mouse:FOVid) + (1|Sex) + (1|CellID)' );
                anova_results = anova(lme, 'DFMethod', 'satterthwaite');
             disp(strcat('local DD global rdnt. LME F=', num2str(anova_results.FStat(2)), '. p=',  num2str(anova_results.pValue(2)),' df=',num2str(anova_results.DF2(2))))
            try dosal=dosal; T_sal_LD=T; save T_sal_LD T_sal_LD; end
            try dodcz=dodcz; T_dcz_LD=T; save T_dcz_LD T_dcz_LD; end
         end
     end
     eff=(mean(tmp1)-mean(tmp2))./std(vertcat(tmp1,tmp2)); disp(strcat('effsizeON=',num2str(eff)));

   figure(3); subplot(1,2,2); scatter(mean(CT(:,stimOFF),2),mean(OB(:,stimOFF),2),'fill'); title('localDEV globalRDNT OffSET responses'); hold on; 
    tmp1=mean(OB(:,stimOFF),2); tmp2=mean(CT(:,stimOFF),2); %aa=horzcat(tmp1,tmp2); x(1)=min(aa(:))*1.02; x(2)=max(aa(:))*1.02;
    y = x; plot(x,y);[h,p,ci,stats]=ttest(tmp1,tmp2);disp(strcat('OFFlocal DD global rdnt. t=',num2str(stats.tstat),'. p=',num2str(p)))
     ylabel(strcat('deviant: p=',num2str(p)),'FontSize',16,'FontWeight','bold'); 
    xlabel('control','FontSize',16,'FontWeight','bold'); 
     figure(4); subplot(1,4,3); errorbar_groups([mean(tmp2) mean(tmp1)]',[std(tmp2)./sqrt(size(tmp1,1)-1) std(tmp1)./sqrt(size(tmp1,1)-1)]', 'bar_colors',[.5 .5 .5; 1 0 0]);
        ylabel('eFR','FontSize',16,'FontWeight','bold'); xlabel('offsets LDGr','FontSize',16,'FontWeight','bold'); set(gca,'FontSize',16,'FontWeight','bold');
        set(gcf,'Color','w');
        props2=((tmp1-tmp2))./((abs(tmp1)+abs(tmp2))./2);
        figure(4); subplot(1,4,4); errorbar_groups([mean(props2)]',[std(props2)./sqrt(size(tmp1,1)-1)]', 'bar_colors',[1 0 0]);
        
        for doLMEstats=doLMEstats;
            if doLMEstats==1;
                tmpmouse = mouseid(incDEV1==1);
                tmpsex   = sexid(incDEV1==1);
                tmpFOV   = FOVid(incDEV1==1);
                n = size(tmp1, 1);

                CellID = (1:n)';   % unique ID per neuron, since tmp1/tmp2 are row-matched

                control_vals = tmp1; deviant_vals = tmp2;
                mice = tmpmouse; sexes = tmpsex; fovs = tmpFOV;

                Value = [control_vals; deviant_vals];
                Condition = [repmat({'control'}, n, 1); repmat({'deviant'}, n, 1)];
                Mouse = [mice; mice];
                Sex = [sexes; sexes];
                FOV = [fovs; fovs];
                Cell = [CellID; CellID];   % same cell ID repeated for both conditions

                T = table(Value, categorical(Condition), categorical(Mouse), categorical(Sex), ...
                    categorical(FOV), categorical(Cell), ...
                    'VariableNames', {'Value', 'Condition', 'Mouse', 'Sex', 'FOVid', 'CellID'});

                lme = fitlme(T, 'Value ~ Condition + (1|Mouse) + (1|Mouse:FOVid) + (1|Sex) + (1|CellID)');

                anova_results = anova(lme, 'DFMethod', 'satterthwaite');
             %disp(strcat('OFFlocal DD global rdnt. LME F=', num2str(anova_results.FStat(2)), '. p=',  num2str(anova_results.pValue(2)),' df=',num2str(anova_results.DF2(2))))
            end
        end
     eff=(mean(tmp1)-mean(tmp2))./std(vertcat(tmp1,tmp2)); disp(strcat('effsizeOFF=',num2str(eff)));



   
    titles={'control run';'resp to B in AAAAB when trained on AAAAB'}; 
    figure;
    for rasterpoot=1;
        tmp1=mean(OB(:,stimON),2); tmp2=mean(CT(:,stimON),2);
        
        [ix,b]=sort(tmp1+tmp2,'descend'); 
       
        for z1=1;
            subplot(1,2,z1); imagesc(timeaxis,1:size(CT,1),CT(b,:)); colormap gray; set(gca,'FontSize',14,'FontWeight','bold');
           xlim([-.12 1.03]);
            title(titles{z1});

            if z1==1; xlabel('sec','FontSize',14,'FontWeight','bold');
                ylabel('neurons','FontSize',14,'FontWeight','bold');
            end
        end
        for z1=2;
            subplot(1,2,z1); imagesc(timeaxis,1:size(OB,1),OB(b,:)); colormap gray; set(gca,'FontSize',14,'FontWeight','bold');
           xlim([-.12 1.03]);
            title(titles{z1});

            if z1==1; xlabel('sec','FontSize',14,'FontWeight','bold');
                ylabel('neurons','FontSize',14,'FontWeight','bold');
            end
        end
        for settingcolorscale=1;
            
            tmpmax=0; tmpmin=0;
            for z1=1:2;
                subplot(1,2,z1);
                cb=caxis; if cb(1)<tmpmin; tmpmin=cb(1); end
                if cb(2)>tmpmax; tmpmax=cb(2); end
              % tmpmin=0;
            end
            for z1=1:2;
                subplot(1,2,z1);
                cb=[tmpmin tmpmax];
                caxis(cb*.2); %caxis([0.02 .25]);
            end
        end
    end
        
end

for globalDDanalysis=1;
    OB=cellsOB(incDEV2==1,:,4);
    CT=cellscontrol(incDEV2==1,:,3);
     for cell=1:size(OB,1) %BASE CORRECTION
        
            OB(cell,:)=OB(cell,:)-mean(OB(cell,basedefMEAN));%
            CT(cell,:)=CT(cell,:)-mean(CT(cell,basedefMEAN));%
        
     end
     
            stdzvals=ones(size(CT,1),1);
       
    clear OBtmp CTtmp
    for cell=1:size(OB,1); OBtmp(cell,:)=OB(cell,:)./stdzvals(cell,1);CTtmp(cell,:)=CT(cell,:)./stdzvals(cell,1);end

    OBmn=mean(OBtmp,1); CTmn=mean(CTtmp,1);
    stdev=std(OBtmp-CTtmp,0,1)./sqrt(size(OB,1)-1);
    stdev1=std(OBtmp)./sqrt(size(OB,1)-1);  stdev2=std(CTtmp)./sqrt(size(CT,1)-1);
    figure; shadedErrorBar(timeaxis,CTmn,stdev2,'lineProps','k'); hold on;
    shadedErrorBar(timeaxis,OBmn,stdev1,'lineProps','r');hold on;
     xlim([-.2 1]); title ('local deviant, global deviant', 'FontSize', 16, 'FontWeight', 'bold');
    ylabel('estFiringRate','FontSize',16,'FontWeight','bold'); 
    xlabel('time (sec)','FontSize',16,'FontWeight','bold'); set(gca,'FontSize',16,'FontWeight','bold');
    set(gcf,'Color','w');
    figure(7); subplot(1,2,1); scatter(mean(CT(:,stimON),2),mean(OB(:,stimON),2),'fill'); title('localDEV globalDEV ONSET responses'); hold on; 
   clear x;  tmp1=mean(OB(:,stimON),2); tmp2=mean(CT(:,stimON),2); aa=horzcat(tmp1,tmp2); x(1)=min(aa(:))*1.02; x(2)=max(aa(:))*1.02;
    y = x; plot(x,y); [h,p,ci,stats]=ttest(tmp1,tmp2);disp(strcat('local DD global DD. t=',num2str(stats.tstat),'. p=',num2str(p)))
    ylabel(strcat('deviant: p=',num2str(p)),'FontSize',16,'FontWeight','bold'); 
    xlabel('control','FontSize',16,'FontWeight','bold');
     figure(8); subplot(1,4,1);errorbar_groups([mean(tmp2) mean(tmp1)]',[std(tmp2)./sqrt(size(tmp1,1)-1) std(tmp1)./sqrt(size(tmp1,1)-1)]', 'bar_colors',[.5 .5 .5; 1 0 0]);
        ylabel('eFR','FontSize',16,'FontWeight','bold'); xlabel('onsets LDGD','FontSize',16,'FontWeight','bold'); set(gca,'FontSize',16,'FontWeight','bold');
        set(gcf,'Color','w');
        props3=((tmp1-tmp2))./((abs(tmp1)+abs(tmp2))./2);
        figure(8); subplot(1,4,2); errorbar_groups([mean(props3)]',[std(props3)./sqrt(size(tmp1,1)-1)]', 'bar_colors',[1 0 0]);
       
     for doLMEstats=doLMEstats;
         if doLMEstats==1;
             tmpmouse = mouseid(incDEV2==1);
             tmpsex   = sexid(incDEV2==1);
             tmpFOV   = FOVid(incDEV2==1);
             n = size(tmp1, 1);

             CellID = (1:n)';   % unique ID per neuron, since tmp1/tmp2 are row-matched

             control_vals = tmp1; deviant_vals = tmp2;
             mice = tmpmouse; sexes = tmpsex; fovs = tmpFOV;

             Value = [control_vals; deviant_vals];
             Condition = [repmat({'control'}, n, 1); repmat({'deviant'}, n, 1)];
             Mouse = [mice; mice];
             Sex = [sexes; sexes];
             FOV = [fovs; fovs];
             Cell = [CellID; CellID];   % same cell ID repeated for both conditions

             T = table(Value, categorical(Condition), categorical(Mouse), categorical(Sex), ...
                 categorical(FOV), categorical(Cell), ...
                 'VariableNames', {'Value', 'Condition', 'Mouse', 'Sex', 'FOVid', 'CellID'});

             lme = fitlme(T, 'Value ~ Condition + (1|Mouse) + (1|Mouse:FOVid) + (1|Sex) + (1|CellID)');
               
                anova_results = anova(lme, 'DFMethod', 'satterthwaite');
             disp(strcat('local DD global DD. LME F=', num2str(anova_results.FStat(2)), '. p=',  num2str(anova_results.pValue(2)),' df=',num2str(anova_results.DF2(2))))
            try dosal=dosal; T_sal_GD=T; save T_sal_GD T_sal_GD; end
            try dodcz=dodcz; T_dcz_GD=T; save T_dcz_GD T_dcz_GD; end
         end
     end
     eff=(mean(tmp1)-mean(tmp2))./std(vertcat(tmp1,tmp2)); disp(strcat('effsizeON=',num2str(eff)));


    figure(7); subplot(1,2,2); scatter(mean(CT(:,stimOFF),2),mean(OB(:,stimOFF),2),'fill'); title('localDEV globalDEV OffSET responses'); hold on; 
    tmp1=mean(OB(:,stimOFF),2); tmp2=mean(CT(:,stimOFF),2);% aa=horzcat(tmp1,tmp2); x(1)=min(aa(:))*1.02; x(2)=max(aa(:))*1.02;
    y = x; plot(x,y);[h,p,ci,stats]=ttest(tmp1,tmp2);;disp(strcat('OFFlocal DD global DD. t=',num2str(stats.tstat),'. p=',num2str(p)))
     ylabel(strcat('deviant: p=',num2str(p)),'FontSize',16,'FontWeight','bold'); 
    xlabel('control','FontSize',16,'FontWeight','bold');
     figure(8); subplot(1,4,3); errorbar_groups([mean(tmp2) mean(tmp1)]',[std(tmp2)./sqrt(size(tmp1,1)-1) std(tmp1)./sqrt(size(tmp1,1)-1)]', 'bar_colors',[.5 .5 .5; 1 0 0]);
        ylabel('eFR','FontSize',16,'FontWeight','bold'); xlabel('offsets LDGD','FontSize',16,'FontWeight','bold'); set(gca,'FontSize',16,'FontWeight','bold');
        set(gcf,'Color','w');
         props4=((tmp1-tmp2))./((abs(tmp1)+abs(tmp2))./2);
        figure(8); subplot(1,4,4); errorbar_groups([mean(props4)]',[std(props4)./sqrt(size(tmp1,1)-1)]', 'bar_colors',[1 0 0]);
      
      for doLMEstats=doLMEstats;
         if doLMEstats==1;
             tmpmouse = mouseid(incDEV2==1);
             tmpsex   = sexid(incDEV2==1);
             tmpFOV   = FOVid(incDEV2==1);
             n = size(tmp1, 1);

             CellID = (1:n)';   % unique ID per neuron, since tmp1/tmp2 are row-matched

             control_vals = tmp1; deviant_vals = tmp2;
             mice = tmpmouse; sexes = tmpsex; fovs = tmpFOV;

             Value = [control_vals; deviant_vals];
             Condition = [repmat({'control'}, n, 1); repmat({'deviant'}, n, 1)];
             Mouse = [mice; mice];
             Sex = [sexes; sexes];
             FOV = [fovs; fovs];
             Cell = [CellID; CellID];   % same cell ID repeated for both conditions

             T = table(Value, categorical(Condition), categorical(Mouse), categorical(Sex), ...
                 categorical(FOV), categorical(Cell), ...
                 'VariableNames', {'Value', 'Condition', 'Mouse', 'Sex', 'FOVid', 'CellID'});

             lme = fitlme(T, 'Value ~ Condition + (1|Mouse) + (1|Mouse:FOVid) + (1|Sex) + (1|CellID)');

                anova_results = anova(lme, 'DFMethod', 'satterthwaite');
       end
     end
    eff=(mean(tmp1)-mean(tmp2))./std(vertcat(tmp1,tmp2)); disp(strcat('effsizeOFF=',num2str(eff)));

    titles={'control run';'resp to C in AAAAC when trained on AAAAB'}; 
    figure;
    for rasterpoot=1;
        tmp1=mean(OB(:,stimON),2); tmp2=mean(CT(:,stimON),2);
        
         [ix,b]=sort(tmp2+tmp1,'descend'); %comment this out if wanting to sort according to first raster
        for z1=1;
            subplot(1,2,z1); imagesc(timeaxis,1:size(CT,1),CT(b,:)); colormap gray; set(gca,'FontSize',14,'FontWeight','bold');
           xlim([-.12 1.03]);
            title(titles{z1});

            if z1==1; xlabel('sec','FontSize',14,'FontWeight','bold');
                ylabel('neurons','FontSize',14,'FontWeight','bold');
            end
        end
        for z1=2;
            subplot(1,2,z1); imagesc(timeaxis,1:size(OB,1),OB(b,:)); colormap gray; set(gca,'FontSize',14,'FontWeight','bold');
           xlim([-.12 1.03]);
            title(titles{z1});

            if z1==1; xlabel('sec','FontSize',14,'FontWeight','bold');
                ylabel('neurons','FontSize',14,'FontWeight','bold');
            end
        end
        for settingcolorscale=1;
            
            tmpmax=0; tmpmin=0;
            for z1=1:2;
                subplot(1,2,z1);
                 if cb(2)>tmpmax; tmpmax=cb(2); end
            end
            for z1=1:2;
                subplot(1,2,z1);
                
                caxis(cb*.2); 
            end
        end
    end
    
end

for redundnatanalysis=1;
    if viewtruerdnt==1; 
    OB=cellsOB(incRDNT==1,:,1);
    CT=cellscontrol(incRDNT==1,:,1);
    FRs=cellsOB(incRDNT==1,:,1);
    else
    OB=cellsOB(incRDNT==1,:,3);
    CT=cellscontrol(incRDNT==1,:,1);
    FRs=cellsOB(incRDNT==1,:,1);
    end

     for cell=1:size(OB,1) %BASE CORRECTION
        
            OB(cell,:)=OB(cell,:)-mean(OB(cell,basedefMEAN));%
            CT(cell,:)=CT(cell,:)-mean(CT(cell,basedefMEAN));%
         
     end
    
            stdzvals=ones(size(CT,1),1);
       
    clear OBtmp CTtmp
    for cell=1:size(OB,1); OBtmp(cell,:)=OB(cell,:)./stdzvals(cell,1);CTtmp(cell,:)=CT(cell,:)./stdzvals(cell,1);end

    OBmn=mean(OBtmp,1); CTmn=mean(CTtmp,1);
    stdev=std(OBtmp-CTtmp,0,1)./sqrt(size(OB,1)-1);
    stdev1=std(OBtmp)./sqrt(size(OB,1)-1);  stdev2=std(CTtmp)./sqrt(size(CT,1)-1);
    figure; shadedErrorBar(timeaxis,CTmn,stdev2,'lineProps','k'); hold on;
    shadedErrorBar(timeaxis,OBmn,stdev1,'lineProps','r');hold on;
    %shadedErrorBar(timeaxis,FRsmn,stdev,'lineProps','b');hold on;
     
    xlim([-.2 1]); title ('local standard, global deviant', 'FontSize', 16, 'FontWeight', 'bold');
    ylabel('estFiringRate','FontSize',16,'FontWeight','bold'); 
    xlabel('time (sec)','FontSize',16,'FontWeight','bold'); set(gca,'FontSize',16,'FontWeight','bold');
    set(gcf,'Color','w');

    figure(11); subplot(1,2,1); scatter(mean(CT(:,stimON),2),mean(OB(:,stimON),2),'fill'); title('localRDNT globalDEV ONSET responses'); hold on; 
   clear x;  tmp1=mean(OB(:,stimON),2); tmp2=mean(CT(:,stimON),2); aa=horzcat(tmp1,tmp2); x(1)=min(aa(:))*1.02; x(2)=max(aa(:))*1.02;
    y = x; plot(x,y);[h,p,ci,stats]=ttest(tmp1,tmp2);disp(strcat('local rdnt global DD. t=',num2str(stats.tstat),'. p=',num2str(p)))
ylabel(strcat('deviant: p=',num2str(p)),'FontSize',16,'FontWeight','bold'); 
    xlabel('control','FontSize',16,'FontWeight','bold');
     figure(12);subplot(1,4,1);  errorbar_groups([mean(tmp2) mean(tmp1)]',[std(tmp1)./sqrt(size(tmp2,1)-1) std(tmp1)./sqrt(size(tmp1,1)-1)]', 'bar_colors',[.5 .5 .5; 1 0 0]);
        ylabel('eFR','FontSize',16,'FontWeight','bold'); xlabel('onsets LrGD','FontSize',16,'FontWeight','bold'); set(gca,'FontSize',16,'FontWeight','bold');
        set(gcf,'Color','w');
        props5=((tmp1-tmp2))./((abs(tmp1)+abs(tmp2))./2);
        figure(12); subplot(1,4,2); errorbar_groups([mean(props5)]',[std(props5)./sqrt(size(tmp1,1)-1)]', 'bar_colors',[1 0 0]);
        ylabel('eFR','FontSize',16,'FontWeight','bold'); xlabel('offsets LDGr','FontSize',16,'FontWeight','bold'); set(gca,'FontSize',16,'FontWeight','bold');
        set(gcf,'Color','w');
       
          for doLMEstats=doLMEstats;
            if doLMEstats==1;
                tmpmouse = mouseid(incRDNT==1);
                tmpsex   = sexid(incRDNT==1);
                tmpFOV   = FOVid(incRDNT==1);
                n = size(tmp1, 1);

                CellID = (1:n)';   % unique ID per neuron, since tmp1/tmp2 are row-matched

                control_vals = tmp1; deviant_vals = tmp2;
                mice = tmpmouse; sexes = tmpsex; fovs = tmpFOV;

                Value = [control_vals; deviant_vals];
                Condition = [repmat({'control'}, n, 1); repmat({'deviant'}, n, 1)];
                Mouse = [mice; mice];
                Sex = [sexes; sexes];
                FOV = [fovs; fovs];
                Cell = [CellID; CellID];   % same cell ID repeated for both conditions

                T = table(Value, categorical(Condition), categorical(Mouse), categorical(Sex), ...
                    categorical(FOV), categorical(Cell), ...
                    'VariableNames', {'Value', 'Condition', 'Mouse', 'Sex', 'FOVid', 'CellID'});

                lme = fitlme(T, 'Value ~ Condition + (1|Mouse) + (1|Mouse:FOVid) + (1|Sex) + (1|CellID)');
               
                anova_results = anova(lme, 'DFMethod', 'satterthwaite');
             disp(strcat('local rdnt global DD. LME F=', num2str(anova_results.FStat(2)), '. p=',  num2str(anova_results.pValue(2)),' df=',num2str(anova_results.DF2(2))))

            end
        end
    eff=(mean(tmp1)-mean(tmp2))./std(vertcat(tmp1,tmp2)); disp(strcat('effsizeON=',num2str(eff)));


    figure(11); subplot(1,2,2); scatter(mean(CT(:,stimOFF),2),mean(OB(:,stimOFF),2),'fill'); title('localRDNT globalDEV OffSET responses'); hold on; 
      tmp1=mean(OB(:,stimOFF),2); tmp2=mean(CT(:,stimOFF),2); %aa=horzcat(tmp1,tmp2); x(1)=min(aa(:))*1.02; x(2)=max(aa(:))*1.02;
    y = x; plot(x,y);[h,p,ci,stats]=ttest(tmp1,tmp2);;disp(strcat('OFFlocal rdnt global DD. t=',num2str(stats.tstat),'. p=',num2str(p)))
     ylabel(strcat('deviant: p=',num2str(p)),'FontSize',16,'FontWeight','bold'); 
    xlabel('control','FontSize',16,'FontWeight','bold'); 
    figure(12);subplot(1,4,3); errorbar_groups([mean(tmp2) mean(tmp1)]',[std(tmp2)./sqrt(size(tmp1,1)-1) std(tmp1)./sqrt(size(tmp1,1)-1)]', 'bar_colors',[.5 .5 .5; 1 0 0]);
        ylabel('eFR','FontSize',16,'FontWeight','bold'); xlabel('offsets LDGr','FontSize',16,'FontWeight','bold'); set(gca,'FontSize',16,'FontWeight','bold');
        set(gcf,'Color','w');
        props6=((tmp1-tmp2))./((abs(tmp1)+abs(tmp2))./2);
   figure(12); subplot(1,4,4); errorbar_groups([mean(props6)]', [std(props6)./sqrt(size(tmp1,1)-1)]', 'bar_colors',[1 0 0]);
       ylabel('eFR','FontSize',16,'FontWeight','bold'); xlabel('offsets LDGr','FontSize',16,'FontWeight','bold'); set(gca,'FontSize',16,'FontWeight','bold');
       set(gcf,'Color','w');
     
         for doLMEstats=doLMEstats;
            if doLMEstats==1;
                tmpmouse = mouseid(incRDNT==1);
                tmpsex   = sexid(incRDNT==1);
                tmpFOV   = FOVid(incRDNT==1);
                n = size(tmp1, 1);

                CellID = (1:n)';   % unique ID per neuron, since tmp1/tmp2 are row-matched

                control_vals = tmp1; deviant_vals = tmp2;
                mice = tmpmouse; sexes = tmpsex; fovs = tmpFOV;

                Value = [control_vals; deviant_vals];
                Condition = [repmat({'control'}, n, 1); repmat({'deviant'}, n, 1)];
                Mouse = [mice; mice];
                Sex = [sexes; sexes];
                FOV = [fovs; fovs];
                Cell = [CellID; CellID];   % same cell ID repeated for both conditions

                T = table(Value, categorical(Condition), categorical(Mouse), categorical(Sex), ...
                    categorical(FOV), categorical(Cell), ...
                    'VariableNames', {'Value', 'Condition', 'Mouse', 'Sex', 'FOVid', 'CellID'});

                lme = fitlme(T, 'Value ~ Condition + (1|Mouse) + (1|Mouse:FOVid) + (1|Sex) + (1|CellID)');

                anova_results = anova(lme, 'DFMethod', 'satterthwaite');
            

            end
        end
    eff=(mean(tmp1)-mean(tmp2))./std(vertcat(tmp1,tmp2)); disp(strcat('effsizeOFF=',num2str(eff)));





    titles={'control run';'resp to A in AAAAA when trained on AAAAB'}; 
    figure;
    for rasterpoot=1;
        tmp1=mean(OB(:,stimON),2); tmp2=mean(CT(:,stimON),2);
       [ix,b]=sort(tmp2+tmp1,'descend'); %coment out if wanting to sort byb first raster
        for z1=1;
            subplot(1,2,z1); imagesc(timeaxis,1:size(CT,1),CT(b,:)); colormap gray; set(gca,'FontSize',14,'FontWeight','bold');
           xlim([-.12 1.03]);
            title(titles{z1});

            if z1==1; xlabel('sec','FontSize',14,'FontWeight','bold');
                ylabel('neurons','FontSize',14,'FontWeight','bold');
            end
        end
        for z1=2;
            subplot(1,2,z1); imagesc(timeaxis,1:size(OB,1),OB(b,:)); colormap gray; set(gca,'FontSize',14,'FontWeight','bold');
           xlim([-.12 1.03]);
            title(titles{z1});

            if z1==1; xlabel('sec','FontSize',14,'FontWeight','bold');
                ylabel('neurons','FontSize',14,'FontWeight','bold');
            end
        end
        for settingcolorscale=1;
            
            tmpmax=0; tmpmin=0;
            for z1=1:2;
                subplot(1,2,z1);
                if cb(2)>tmpmax; tmpmax=cb(2); end
            end
            for z1=1:2;
                subplot(1,2,z1);
                 caxis(cb*.2); 
            end
        end
    end
    
        
end
totalresponsivecells=287; %figure this out separately. 

for determine_propsBOTH=1;
    criter=.025;
   if and(do_onsets==1,do_offsets==1)
        
        for redunndantanalysis=1;
            type=1;
            cnts=0;
            for cell=1:size(incRDNT,2);
                if incRDNT(1,cell)==1;
                    cnts=1+cnts;
                    tmpcnts=celltrialON{cell,type};
                    tmpob=celltrialON{cell,type+3};
                    [h,p1,ci,stats]=ttest2(tmpcnts,tmpob,"Tail","left"); %left means second is greate
                    tmpcnts=celltrialOFF{cell,type};
                    tmpob=celltrialOFF{cell,type+3};
                    [h,p2,ci,stats]=ttest2(tmpcnts,tmpob,"Tail","left"); %left means second is greater
                    DDrdnt(cnts,1)=or(p1<criter,p2<criter)*1;
                    % tmpcnts1=celltrialON{cell,type};
                    % tmpob1=celltrialON{cell,type+3};
                    % tmpcnts2=celltrialOFF{cell,type};
                    % tmpob2=celltrialOFF{cell,type+3};
                    % [clusters, p_values, t_sums, permutation_distribution] = permutest(vertcat(tmpcnts1,tmpcnts2)',vertcat(tmpob1,tmpob2)',false);
                    % DDrdnt(cnts,1)=(p_values<criter)*1;
                end
            end
        end
         disp(strcat('proportion of rdntDD cells:',num2str(mean(DDrdnt))))  
         tot=totalresponsivecells; resp=cnts; dev=sum(DDrdnt);
        figure; 
        ax = gca(); 
        pieData = [(resp-dev)/tot (tot-(resp))/tot dev/tot]; 
        h = pie(ax, pieData); 
        % Define 3 colors, one for each of the 3 wedges
        newColors = [...
        0.5,       0.5,       .5;
        1, 1, 1;   
        1,       0, 0];  
        ax.Colormap = newColors; 
        title(strcat('proportion of dev detecting cells: AAAA-A', num2str(dev/tot)))

        for localDDanalysis=1;
            type=2;
            cnts=0;
            for cell=1:size(incDEV1,2);
                if incDEV1(1,cell)==1;
                    cnts=1+cnts;
                     tmpcnts=celltrialON{cell,type};
                    tmpob=celltrialON{cell,type+3};
                    [h,p1,ci,stats]=ttest2(tmpcnts,tmpob,"Tail","left"); %left means second is greate
                    tmpcnts=celltrialOFF{cell,type};
                    tmpob=celltrialOFF{cell,type+3};
                    [h,p2,ci,stats]=ttest2(tmpcnts,tmpob,"Tail","left"); %left means second is greater
                    DDdev1(cnts,1)=or(p1<criter,p2<criter)*1;
                    % tmpcnts1=celltrialON{cell,type};
                    % tmpob1=celltrialON{cell,type+3};
                    % tmpcnts2=celltrialOFF{cell,type};
                    % tmpob2=celltrialOFF{cell,type+3};
                    % [clusters, p_values, t_sums, permutation_distribution] = permutest(vertcat(tmpcnts1,tmpcnts2)',vertcat(tmpob1,tmpob2)',false);
                    % DDdev1(cnts,1)=(p_values<criter)*1;
                end
            end
        end
        disp(strcat('proportion of local DD cells:',num2str(mean(DDdev1))))
            tot=totalresponsivecells; resp=cnts; dev=sum(DDdev1);
        figure; 
        ax = gca(); 
        pieData = [(resp-dev)/tot (tot-(resp))/tot dev/tot]; 
        h = pie(ax, pieData); 
        % Define 3 colors, one for each of the 3 wedges
        ax.Colormap = newColors; 
        title(strcat('proportion of dev detecting cells: AAAA-B', num2str(dev/tot)))

        for globalDDanalysis=1;
            type=3;
            cnts=0;
            for cell=1:size(incDEV2,2);
                if incDEV2(1,cell)==1;
                    cnts=1+cnts;
                     tmpcnts=celltrialON{cell,type};
                    tmpob=celltrialON{cell,type+3};
                    [h,p1,ci,stats]=ttest2(tmpcnts,tmpob,"Tail","left"); %left means second is greate
                    tmpcnts=celltrialOFF{cell,type};
                    tmpob=celltrialOFF{cell,type+3};
                    [h,p2,ci,stats]=ttest2(tmpcnts,tmpob,"Tail","left"); %left means second is greater
                    DDdev2(cnts,1)=or(p1<criter,p2<criter)*1;
                    % tmpcnts1=celltrialON{cell,type};
                    % tmpob1=celltrialON{cell,type+3};
                    % tmpcnts2=celltrialOFF{cell,type};
                    % tmpob2=celltrialOFF{cell,type+3};
                    % [clusters, p_values, t_sums, permutation_distribution] = permutest(vertcat(tmpcnts1,tmpcnts2)',vertcat(tmpob1,tmpob2)',false);
                    % DDdev2(cnts,1)=(p_values<criter)*1;
                end
            end
        end
        disp(strcat('proportion of glob DD cells:',num2str(mean(DDdev2))))
        tot=totalresponsivecells; resp=cnts; dev=sum(DDdev2);
        figure; 
        ax = gca(); 
        pieData = [(resp-dev)/tot (tot-(resp))/tot dev/tot]; 
        h = pie(ax, pieData); 
        ax.Colormap = newColors; 
        title(strcat('proportion of dev detecting cells: AAAA-C', num2str(dev/tot)))


    end
end


