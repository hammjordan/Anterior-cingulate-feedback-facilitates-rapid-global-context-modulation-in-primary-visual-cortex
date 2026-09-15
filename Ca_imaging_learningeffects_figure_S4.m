clear all
close all
files={'013024_m1_pyrs_global_002_pre_train';'101523_m1_pyrs_global_001_pre_train';'013124_m1_pyrs_global_002_pre_train';
    '012924_m1_pyrs_global_001_pre_train';'112923_m1_pyrs_global_001_pre_train'
    '013024_m1_pyrs_global_002_pre_trainREV';'101523_m1_pyrs_global_001_pre_trainREV';'013124_m1_pyrs_global_002_pre_trainREV';
    '012924_m1_pyrs_global_001_pre_trainREV';'112923_m1_pyrs_global_001_pre_trainREV'};

filesCNT={'013024_m1_pyrs_global_002_pre';'101523_m1_pyrs_global_rev_001_pre';'013124_m1_pyrs_global_002_pre';
    '012924_m1_pyrs_global_001_pre';'112923_m1_pyrs_global_001_pre';
    '013024_m1_pyrs_global_rev_002_pre';'101523_m1_pyrs_global_rev_001_pre';'013024_m1_pyrs_global_rev_002_pre';
    '012924_m1_pyrs_global_rev_001_pre';'112923_m1_pyrs_global_001_pre'};
revs=[0 0 0 0 0 1 1 1 1 1];
%during oddball phase, the 4's are the local oddball (e.g. 2),

cellcount=0; cellscontrol=[]; cellsOB=[]; cellscontrolSTD=[]; cellsOBSTD=[]; mouseid=[];fovid=[];dataall=[];dataallSTD=[];trialAVSall=[];trialAVSallSTD=[];
basedefMEAN=10:14; stimdefON=15:28; stimdefOFF=31:39;%these are for the standardizing step, which determines whether cells are "responsive"
stdnorm=1;
rejectingdoubles=1;
cut1=1.97; cut2=2.6; bins=5;


%mainly for selecting out responsive neurons
for fil=1:size(files,1)

    load(strcat('CA_workspace_SORTED_train',files{fil}))
    if revs(fil)==0;
        load(strcat('reject_cells_',files{fil}(1:end-6)));
    else
        load(strcat('reject_cells_',files{fil}(1:end-9)));
    end


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

        end
    end

    samp=1/framerate;
    timeaxis=-(samp*(framerate/2))+samp:samp:(samp*(framerate*1))+samp;

    vec=vectrials{1};
    vecdat=vec{1,1};
    vecstim=vec{1,2}; vecstim1=vecstim;
    ords=vec{1,3};
    locs=vec{1,4};
    trialindex=1:size(ords,2);


    numtrials(fil,1)=sum(vecstim==2);numtrials(fil,2)=sum(vecstim==4);numtrials(fil,3)=sum(vecstim==8);

    %for normalizing all cell activity by its stdev
    dfofanalyse1=dfofanalyse; vecdatSTD=vecdat;
    for cell=1:size(vecdat,1);
        ssss=dfofanalyse1(cell,:); sss=ssss(ssss>-0); ss=sss(sss<prctile(sss,50));
        vecdatSTD(cell,:,:)=vecdat(cell,:,:)./std(ss);
    end


    dfofanalyse1=dfofanalyse;


    %for excluding the first redundant (only for seleicting responsive
    %neurons)
    vecstim=horzcat(vecstim,zeros(1,4));
    for tt=1:size(vecstim,2)-4
        if and(vecstim(1,tt)==2,vecstim(1,tt+4)==4)
            vecstim(1,tt)=0;
        end
    end
    vecstim=vecstim(1,1:end-4);

    %for baseadj
    for cell=1:size(vecdat,1);
        for tr=1:size(vecdat,3);
            vecdat(cell,:,tr)=vecdat(cell,:,tr)-mean(vecdat(cell,basedefMEAN,tr));
            vecdatSTD(cell,:,tr)=vecdatSTD(cell,:,tr)-mean(vecdatSTD(cell,basedefMEAN,tr));
        end
    end

    dattmp=[];stimIDs=[];dattmpstd=[];stimIDsON=[];trAVs=[];trAVsSTD=[];
    for tr=1:max(trialindex(:));
        dattmp=horzcat(dattmp,vecdat(:,11:31,tr));
        dattmpstd=horzcat(dattmpstd,vecdatSTD(:,11:31,tr));
        stimIDs=horzcat(stimIDs,zeros(1,21)+vecstim(tr));
        tmpa=zeros(1,21)+vecstim(tr);tmpa(1,1:6)=zeros(1,6);%tmpa(1,18:21)=zeros(1,4);
        stimIDsON=horzcat(stimIDsON,tmpa);
        trAVsSTD=horzcat(trAVsSTD,mean(vecdatSTD(:,stimdefON,tr),2));
        trAVs=horzcat(trAVs,mean(vecdat(:,stimdefON,tr),2));
    end
    dataall=vertcat(dataall,dattmp);
    dataallSTD=vertcat(dataallSTD,dattmpstd);
    trialAVSall=vertcat(trialAVSall,trAVs);
    trialAVSallSTD=vertcat(trialAVSallSTD,trAVsSTD);

    for cell=1:size(dattmpstd,1);
        stepa=floor(size(stimIDsON,2)./bins);
        for binna=1:bins;
            stimtmp=stimIDsON-stimIDsON;
            stimtmp(1,((binna-1)*stepa)+1:(binna*stepa))=stimIDsON(1,((binna-1)*stepa)+1:(binna*stepa));
            incsRDNT(cellcount+cell,binna)=mean(dattmpstd(cell,stimtmp==2));
            incsOB(cellcount+cell,binna)=mean(dattmpstd(cell,stimtmp==4));

        end
        stimtmp=stimIDsON(1,:);
        incsRDNT(cellcount+cell,bins+1)=mean(dattmpstd(cell,stimtmp==2));
        incsOB(cellcount+cell,bins+1)=mean(dattmpstd(cell,stimtmp==4));
        locIDs(cellcount+cell,:)=locs;

    end
    mouseid((cellcount+1):(cell+cellcount),1)=zeros(cell,1)+fil;
    fovid((cellcount+1):(cell+cellcount),1)=zeros(cell,1)+fil;
    cellcount=cell+cellcount;


end
if stdnorm==1;
    dataall=dataallSTD;
    trialAVSall=trialAVSallSTD;
end
incsRDNTtmp=[]; incsOBtmp=[];
for c=1:size(incsRDNT,1);
    if or(or(incsRDNT(c,1)>cut2,incsRDNT(c,end-1)>cut2),and(incsRDNT(c,1)>cut1,incsRDNT(c,end-1)>cut1))
        incsRDNTtmp(c,1)=1;
    else incsRDNTtmp(c,1)=0;
    end
    if or(or(incsOB(c,1)>cut2,incsOB(c,end-1)>cut2),and(incsOB(c,1)>cut1,incsOB(c,end-1)>cut1))
        incsOBtmp(c,1)=1;
    else incsOBtmp(c,1)=0;
    end

end
cut=.99; incsRDNT=incsRDNTtmp; incsOB=incsOBtmp;
vecstim=vecstim1;






%%%main analysis
cellcount1=0; cellcount2=0;  cellsOB=[]; cellsOBSTD=[]; dataall=[];dataallSTD=[];trialAVSall=[];trialAVSallSTD=[];
cellsCNT=[]; cellsCNTSTD=[]; incsCNT=[]; locoOB=[]; locoCNT=[];
for fil=1:size(files,1)

    load(strcat('CA_workspace_SORTED_train',files{fil}))
    if revs(fil)==0;
        load(strcat('reject_cells_',files{fil}(1:end-6)));
    else
        load(strcat('reject_cells_',files{fil}(1:end-9)));
    end

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


        end
    end

    samp=1/framerate;
    timeaxis=-(samp*(framerate/2))+samp:samp:(samp*(framerate*1))+samp;
    vec=vectrials{1};
    vecdat=vec{1,1};

    bingroups=round(0:((size(vecdat,3))/bins):(size(vecdat,3)))
    for bin=1:bins;
        vec=vectrials{1};
        vecdat=vec{1,1};
        vecstim=vec{1,2};
        ords=vec{1,3};
        locs=vec{1,4};

        vecdat=vecdat(:,:,(bingroups(bin)+1):(bingroups(bin+1)));
        vecstim=vecstim(1,(bingroups(bin)+1):(bingroups(bin+1)));
        ords=ords(1,(bingroups(bin)+1):(bingroups(bin+1)));
        locs=locs(1,(bingroups(bin)+1):(bingroups(bin+1)));


        locoOB(fil,bin)=mean(locs(:));



        dfofanalyse1=dfofanalyse;



        for cell=1:size(vecdat,1);
            celltmp=[]; celltrialtmp=[];
            for ori=1:5
                celltmp(:,ori)=squeeze(mean(vecdat(cell,:,ords==ori),3));
            end
            ssss=dfofanalyse1(cell,:); sss=ssss(ssss>-0); ss=sss(sss<prctile(sss,50));
            for ori=1:5
                cellsOBSTD(cellcount1+cell,:,ori,bin)=(celltmp(:,ori)-mean(celltmp(basedefMEAN,ori)))./std(ss(:));
            end

            for ori=1:5
                cellsOB(cellcount1+cell,:,ori,bin)=celltmp(:,ori)-mean(celltmp(basedefMEAN,ori));
            end


        end

    end
    cellcount1=cell+cellcount1;


    %do the control run
    load(strcat('CA_workspace_SORTED_',filesCNT{fil}));
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


    for bin=1:bins;
        vec=vectrials{1};
        vecdat=vec{1,1};
        vecstim=vec{1,2};
        ords=vec{1,3};
        locs=vec{1,4};

        vecdat=vecdat(:,:,(bingroups(bin)+1):(bingroups(bin+1)));
        vecstim=vecstim(1,(bingroups(bin)+1):(bingroups(bin+1)));
        ords=ords(1,(bingroups(bin)+1):(bingroups(bin+1)));
        locs=locs(1,(bingroups(bin)+1):(bingroups(bin+1)));



        locoCNT(fil,bin)=mean(locs(:));



        dfofanalyse1=dfofanalyse;


        %for limiting to the right orientations
        for s=1:size(vecstim,2);
            if vecstim(1,s)>2;
                %ords(1,s)=0;
            end
        end
        for cell=1:size(vecdat,1);
            celltmp=[]; celltrialtmp=[];
            for ori=1:5
                celltmp(:,ori)=squeeze(mean(vecdat(cell,:,ords==ori),3));
            end
            ssss=dfofanalyse1(cell,:); sss=ssss(ssss>-0); ss=sss(sss<prctile(sss,50));
            for ori=1:5
                cellsCNTSTD(cellcount2+cell,:,ori,bin)=(celltmp(:,ori)-mean(celltmp(basedefMEAN,ori)))./std(ss(:));
            end

            for ori=1:5
                cellsCNT(cellcount2+cell,:,ori,bin)=celltmp(:,ori)-mean(celltmp(basedefMEAN,ori));
            end


        end
        %make "inc" variable to figure out which cells show a true response
        for cella=1:size(vecdat,1)
            inc_cntON(cellcount2+cella,:,bin)=squeeze(mean(cellsCNTSTD(cellcount2+cella,stimdefON,:,bin),2)-mean(cellsCNTSTD(cellcount2+cella,basedefMEAN,:,bin),2));

        end


    end
    cellcount2=cell+cellcount2;




end

for stdnorm=stdnorm;
    if stdnorm==1;
        cellsOB=cellsOBSTD;cellsCNT=cellsCNTSTD;
    end
end
for z=1:size(mouseid,1);
    if mouseid(z)>5;
        if mouseid(z)~=8;
            mouseid(z)=mouseid(z)-5;
        end
    end
end

figure;
%plot controls. the "cut" is different here than in the next plot, bc i
%already sorted out the include/exclude for teh OB session in the first
%part of the script. it's cut2.
for bin=1:bins;
    clear CT
    OB=cellsCNT(inc_cntON(:,5,bin)>cut2,:,5,bin);
    m5ct{bin}=mouseid(inc_cntON(:,5,bin)>cut2);
    fovct5{bin}=fovid(inc_cntON(:,5,bin)>cut2);
    for stim=1:4;
        CT{stim}=cellsCNT(inc_cntON(:,stim,bin)>cut2,:,stim,bin);
        mmCT{stim,bin}=mouseid(inc_cntON(:,stim,bin)>cut2);
        fvCT{stim,bin}=fovid(inc_cntON(:,stim,bin)>cut2);
    end

    OB1=OB(:,11:41,1); CTmn=[];stdevCT=[];
    CT5stat{bin}=mean(OB1(:,5:18),2);
    for stim=1:4
        CC=CT{stim};
        CCmn=mean(CC(:,11:41),1);stdevCC=std(CC(:,11:41),0,1)./sqrt(size(CC,1)-1);
        CTmn=horzcat(CTmn,CCmn);
        stdevCT=horzcat(stdevCT,stdevCC);
        CT1stat{stim,bin}=mean(CC(:,5:18),2);
    end
    OB=OB1;
    OBmn=mean(OB,1);
    stdevOB=std(OB,0,1)./sqrt(size(OB,1)-1);
    mn=horzcat(CTmn,OBmn);
    st=horzcat(stdevCT,stdevOB);
    subplot(1,bins,bin); shadedErrorBar(1:size(mn,2),mn,st,'lineProps','k'); hold on;
    axis;
    axes(bin+5,:)=ans(1,3:4);
    subplot(1,bins,bin); shadedErrorBar((size(mn,2)-31):size(mn,2),mn(1,(size(mn,2)-31):size(mn,2)),st(1,(size(mn,2)-31):size(mn,2)),'lineProps','k'); hold on;
    if bin==1;
        title('first trials during training');
    elseif bin==bins;
        title('last trials');
    end
    make_eps_saveable


end
%plot OBs
for bin=1:bins;
    OB=cellsOB(incsOB>cut,:,5,bin);
    m5ob=mouseid(incsOB>cut);
    fovOB5=fovid(incsOB>cut);
    CT=cellsOB(incsRDNT>cut,:,1:4,bin);
    mmob=mouseid(incsRDNT>cut);
    fvob=fovid(incsRDNT>cut);

    OB1=OB(:,11:41,1); CT1=[];
    OBstat{bin}=mean(OB1(:,5:18),2);
    for stim=1:4;
        CT1=horzcat(CT1,CT(:,11:41,stim));
        RDstat{stim,bin}=mean(CT(:,5:18),2);
    end

    OB=OB1; CT=CT1;
    OBmn=mean(OB,1); CTmn=mean(CT,1);
    mn=horzcat(CTmn,OBmn);
    stdevOB=std(OB,0,1)./sqrt(size(OB,1)-1);
    stdevCT=std(CT,0,1)./sqrt(size(CT,1)-1);
    st=horzcat(stdevCT,stdevOB);
    subplot(1,bins,bin); shadedErrorBar(1:size(mn,2),mn,st,'lineProps','b'); hold on;
    axis;
    axes(bin,:)=ans(1,3:4);
    subplot(1,bins,bin); shadedErrorBar((size(mn,2)-31):size(mn,2),mn(1,(size(mn,2)-31):size(mn,2)),st(1,(size(mn,2)-31):size(mn,2)),'lineProps','r'); hold on;
    make_eps_saveable

end

for bin=1:bins;

    subplot(1,bins,bin); ylim([min(axes(:,1)) max(axes(:,2))]);

end

disp(strcat('proportion of cells=',num2str(mean(incsOB>cut))));


%do stats


doLMEstats=1;
for doLMEstats=doLMEstats;
    if doLMEstats==1;
        for bin=1:5;
            for stim=1:4;
                a=RDstat{stim,bin};
                b=CT1stat{stim,bin};
                tmpmouseA=mmob;tmpmouseB=mmCT{stim,bin};
                fovA=fvob;fovB=fvCT{stim,bin};
                values = [a; b];
                condition = [repmat("A", length(a), 1);
                    repmat("B", length(b), 1)];
                mouse = [tmpmouseA; tmpmouseB];
                fovs=[fovA; fovB];

                % 2. Put into a table
                tbl = table(values, condition, mouse, fovs);
                tbl.mouse = categorical(tbl.mouse);
                tbl.condition = categorical(tbl.condition);

                % 3. Fit linear mixed-effects model:
                % Fixed effect: condition (A vs B)
                % Random effect: mouse (random intercept)
                lme = fitlme(tbl, 'values ~ condition + (1|mouse)+(1|fovs)');

                % 4. Display results
                % disp(lme);
                X=anova(lme, 'DFMethod', 'satterthwaite');
                tvals(bin,stim)=X{2,2};pvals(bin,stim)=X{2,5};
                dfvals(bin,stim)=X{2,4};
            end
            for stim=5;
                a=OBstat{1,bin};
                b=CT5stat{1,bin};
                tmpmouseA=m5ob;tmpmouseB=m5ct{bin};
                fovA=fovOB5;fovB=fovct5{bin};
                values = [a; b];
                condition = [repmat("A", length(a), 1);
                    repmat("B", length(b), 1)];
                mouse = [tmpmouseA; tmpmouseB];
                fovs=[fovA; fovB];


                % 2. Put into a table
                tbl = table(values, condition, mouse, fovs);
                tbl.mouse = categorical(tbl.mouse);
                tbl.condition = categorical(tbl.condition);

                % 3. Fit linear mixed-effects model:
                % Fixed effect: condition (A vs B)
                % Random effect: mouse (random intercept)
                lme = fitlme(tbl, 'values ~ condition + (1|mouse)+(1|fovs)');


                X=anova(lme);
                tvals(bin,stim)=X{2,2};pvals(bin,stim)=X{2,5};
                dfvals(bin,stim)=X{2,4};
            end
        end
    end
end

disp(pvals(1,4));
disp(pvals(1,5));
disp(pvals(5,4));
disp(pvals(5,5));


for locomotion_analysis=1;
    % figure; plot(mean(locoOB(revs==0,:)),'r'); hold on;
    % plot(mean(locoOB(revs==1,:)),'--r'); hold on;
    % plot(mean(locoCNT(revs==0,:)),'k'); hold on;
    % plot(mean(locoCNT(revs==1,:)),'--k'); hold on;
    clear X
    mm=[1 2 3 4 5 1 2 6 4 5];
    for z=1:6;
        locoOB1(z,:)=mean(locoOB(mm==z,:),1);
        locoCNT1(z,:)=mean(locoCNT(mm==z,:),1);
    end

    % figure; plot(mean(locoOB1,1),'r'); hold on;
    % plot(mean(locoCNT1,1),'k'); hold on;
    %  anova_rm({locoOB1 locoCNT1})

    % locoOB1 and locoCNT1 are 6x5 (mice x timebins)

    % Compute summary stats
    meanOB1 = mean(locoOB1,1);
    meanCNT1 = mean(locoCNT1,1);
    semOB1 = std(locoOB1,[],1)./sqrt(size(locoOB1,1));
    semCNT1 = std(locoCNT1,[],1)./sqrt(size(locoCNT1,1));

    timeBins = 1:5;
    jit = 0.1; % jitter spacing for visibility

    % Colors
    colOB  = [0.85 0.1 0.1];   % red
    colCNT = [0.2 0.2 0.2];    % dark gray
    spagColOB  = [0.85 0.1 0.1 0.2]; % faint spaghetti
    spagColCNT = [0.2 0.2 0.2 0.2];

    figure; hold on;

    % --- Spaghetti lines for individual mice ---
    for m = 1:size(locoOB1,1)
        plot(timeBins - jit, locoOB1(m,:), '-', 'Color', spagColOB);
        plot(timeBins + jit, locoCNT1(m,:), '-', 'Color', spagColCNT);
    end

    % --- Individual mouse dots ---
    scatter(timeBins - jit, locoOB1, 45, colOB,  'filled', 'MarkerFaceAlpha', 0.8);
    scatter(timeBins + jit, locoCNT1, 45, colCNT, 'filled', 'MarkerFaceAlpha', 0.8);

    % --- Mean ± SEM (larger marker) ---
    errorbar(timeBins - jit, meanOB1, semOB1, 'o', 'Color','k', ...
        'MarkerFaceColor',colOB, 'MarkerEdgeColor','k', ...
        'LineStyle','none', 'LineWidth',1.5, 'MarkerSize',10);

    errorbar(timeBins + jit, meanCNT1, semCNT1, 'o', 'Color','k', ...
        'MarkerFaceColor',colCNT, 'MarkerEdgeColor','k', ...
        'LineStyle','none', 'LineWidth',1.5, 'MarkerSize',10);

    % --- Formatting ---
    xlim([0.5 5.5]);
    xlabel('Time bin');
    ylabel('Locomotion (rate)');
    legend({'OB1 mice','CNT1 mice','OB1 mean ± SEM','CNT1 mean ± SEM'}, 'Location','best');
    title('Locomotion Over Time: OB1 vs CNT1');
    set(gca,'FontSize',12,'Box','off');
end


