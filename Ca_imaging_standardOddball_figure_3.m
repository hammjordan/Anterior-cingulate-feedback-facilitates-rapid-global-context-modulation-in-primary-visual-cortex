clcl
for file_load=1;


files={'012924_m1_pyrs_oddball_001_pre';'013024_m1_pyrs_oddball_002_pre';'013124_m1_pyrs_oddball_002_pre';};
doingmatch=1;onestim=0; %if 0, then plot responses to both stimuli (e.g. 45 and 135 deg). if 1, then plot only to first orientation. if 2, then plot to second
realign=[1 1 1]; %the oddball run got out of alignment. 

%%% 
%note, for matching to global local run, run with onestim=1 first, and then
%with onestim=2. stack those matrices and they will align with the global
%lcoal data from tehse same mice. also, in global local data, you have rev
%and regular. you use the reverse to map to onestim=1, and regular goes to
%onestim=2.
end

limtrials1=6; %limit the number of trials analysed for each stimulus type(control/redundandt/deviant). IF set high, then it normalizes to the smallest number of trials in any condition (cont, rdnt, dev)
baseBEGIN=13; baseEND=16; %what you call the baseline
BEGINo=18; ENDo=31; %these are the analysis windows. note: stimulus starts at 18 and finishes at 31
BEGINoff=32; ENDoff=40; 
dobeforeafter=0; %if this is "1", then it compares the first trials to the last for the test. if it's zero, it does even/odd
exclloco=1; rejectingdoubles=1;
%%%%%%%this compiles all the files into three matrices
for buildingmatrices=1;
    effectsallBUILD=[]; %even
    effectsallTEST=[]; %odd
    singtrialsall=[];
    count=0; counttest=0;
    effectsall=[];  %this is the main varable
    effectsallSTD=[];  %this is for limiting trials
    actcont=[];
    mousevec=[];
    countall=0;
end
for file=1:size(files,1);
   
    load(strcat('CA_workspace_SORTED_',files{file}),'vectrials','framerate','phzstd','phases','dfofanalyse','stimvisMAIN','stimvistypeMAIN');
      
    for eliminating_repeats=1;
        for stim=1:2;
            timestamps=vectrials{stim,3};
            trialtypes=vectrials{stim,2};
            timesince=horzcat(10000,diff(timestamps));
            for tr=1:size(trialtypes,2);
                if trialtypes(1,tr)==phzstd(stim);
                    if timesince(1,tr)<45;
                        trialtypes(1,tr)=99;
                    end
                end
            end
            vectrials{stim,2}=trialtypes;
        end
    end

    for exclloco=exclloco;
        if exclloco==1;
            for stim=1:2;
                locos=vectrials{stim,4};
                trialtypes=vectrials{stim,2};
                for tr=1:size(trialtypes,2)
                    if locos(1,tr)==1;
                        trialtypes(1,tr)=99;
                    end
                end
                vectrials{stim,2}=trialtypes;
            end
        end
    end
    %figure; plot(stimvistypeMAIN); title(num2str(file));
    disp(files{file});
    for remove_first_trial_of_each_type=1;
        if remove_first_trial_of_each_type==1; 
            if limtrials1>0;
                limtrials=limtrials1;
                for ss1=1:2;
                    vecstim=vectrials{ss1,2};
                    tmps1=[];
                    for tr1=1:size(vecstim,2);
                        
                        if and(sum(vecstim(tr1)==tmps1)==0,vecstim(tr1)~=77)
                            tmps1=horzcat(tmps1,vecstim(tr1));
                            vecstim(tr1)=99;
                        end
                        
                    end
                    for additionally_remove_first_few_controls=1; %counting backward from the end
                        tmps1=[]; gh=size(vecstim,2)+1;
                        disp(sum(vecstim==(20+(ss1-1)))); 
                        for tr1=1:size(vecstim,2);
                            if and(vecstim(gh-tr1)>15,vecstim(gh-tr1)<22)
                                if sum(tmps1==1)<(limtrials)
                                    tmps1=horzcat(tmps1,1);
                                else
                                    vecstim(gh-tr1)=99;
                                end
                            end
                        end
                    end
                    vectrials{ss1,2}=vecstim;
                end
            
               
            end
        end
    end
    
    load(strcat('reject_cells_',files{file}));
    for rejectingdoubles=rejectingdoubles;
        if rejectingdoubles==1;
            vecdat=vectrials{1,1};
            vecdat=vecdat(reject_cells==0,:,:);
            vectrials{1,1}=vecdat;

            vecdat=vectrials{2,1};
            vecdat=vecdat(reject_cells==0,:,:);
            vectrials{2,1}=vecdat;
        end
    end


    clear effects effectsSTD
    for analysis=1
        
        samp=1/framerate;
        timeaxis=-((samp*(framerate/2))+(samp*2)):samp:((samp*framerate)-(samp*2)); 
        
        for stim=1:2;
            limtrials=limtrials1;
            vecdat=vectrials{stim,1};
            mousetmp=zeros(size(vecdat,1),1)+file; 
            vecstim=vectrials{stim,2};
            if realign(file)==1;
                for tr=1:size(vecdat,3);                 
                    if vecstim(tr)==15;
                         vecdat(:,3:43,tr)=vecdat(:,1:41,tr);
                     end
                end
            end
            if limtrials1>0; 
            for checkingfortoo_few_trials=1;
                ss=[1 2 3 4 15 phzstd(stim)];
                for gss=1:6; 
                    sss(gss)=sum(vecstim==ss(gss));
                end
                if min(sss)<limtrials1;
                   limtrials=min(sss); 
                end               
                disp(limtrials); 
            end
            
            effects=[];
            for z1=1:8
                tmp1=vecdat(:,:,vecstim==z1);
                if size(tmp1,3)>(limtrials-1)
                    effects(:,:,z1+1)=mean(tmp1(:,:,1:limtrials),3);
                else
                    effects(:,:,z1+1)=mean(tmp1(:,:,1:end),3);
                end
            end
            tmp1=vecdat(:,:,vecstim==15);
            if size(tmp1,3)>(limtrials-1)
                effects(:,:,10)=mean(tmp1(:,:,1:limtrials),3);
            else
                effects(:,:,10)=mean(tmp1(:,:,1:end),3);
            end
            tmp1=vecdat(:,:,vecstim==phzstd(stim));
            
                if size(tmp1,3)>(limtrials-1)
                    effects(:,:,1)=mean(tmp1(:,:,1:limtrials),3);
                     else
                    effects(:,:,1)=mean(tmp1(:,:,1:end),3);
                end
                
            else
                 vecdat=vectrials{stim,1};
                vecstim=vectrials{stim,2};
                effects=[];
                stimax=[phzstd(stim) 1 2 3 4 5 6 7 8 15];
                for z1=1:10
                    effects(:,:,z1)=mean(vecdat(:,:,vecstim==stimax(z1)),3);
                end
            end
                
            
            for c1=1:size(effects,1) %baseline correction
                basos=[1 2 2 2 2 2 2 2 2 3; 1 3 3 3 3 3 3 3 3 2];
                
                for z1=1:size(effects,3)
                    mnn=mean(effects(c1,baseBEGIN:baseEND,z1),2); 
                    tmpbasos=dfofanalyse(c1,phases==basos(stim,z1));
                    tmpbasos=tmpbasos(tmpbasos>0); tmpbasos=tmpbasos(tmpbasos<prctile(tmpbasos,50));
                    stt=std(tmpbasos); 
                    
                    if stt<.001; stt=std(effects(c1,ENDo:end,z1)); end
                    if stt<.001; stt=10000; end
                    for t1=1:size(effects,2)
                        effectsSTD(c1,t1,z1)=(effects(c1,t1,z1)-mnn)./stt;
                        effects(c1,t1,z1)=(effects(c1,t1,z1)-mnn);
                    end
                    mnns(file,stim,c1,z1)=mnn;
                    stts(file,stim,c1,z1)=stt;
                end
                
            end
            
                
            
            if onestim==0;
                effectsall((1+countall):(countall+size(effects,1)),1:size(effects,2),2:10,stim)=effects(:,1:end,2:10);
                effectsall((1+countall):(countall+size(effects,1)),1:size(effects,2),1,stim)=effects(:,1:end,1);
                
                mousevec((1+countall):(countall+size(effects,1)),stim)=mousetmp;
            elseif onestim==1;
                if stim==1;
                    effectsall((1+countall):(countall+size(effects,1)),1:size(effects,2),2:10,stim)=effects(:,1:end,2:10);
                    effectsall((1+countall):(countall+size(effects,1)),1:size(effects,2),1,stim)=effects(:,1:end,1);
                    
                mousevec((1+countall):(countall+size(effects,1)),stim)=mousetmp;
                else
                    effectsall((1+countall):(countall+size(effects,1)),1:size(effects,2),2:10,stim)=effects(:,1:end,2:10)-effects(:,1:end,2:10);
                    effectsall((1+countall):(countall+size(effects,1)),1:size(effects,2),1,stim)=effects(:,1:end,1)-effects(:,1:end,1);
                   
                mousevec((1+countall):(countall+size(effects,1)),stim)=mousetmp-mousetmp;
                end
            elseif onestim==2;
                if stim==2;
                    effectsall((1+countall):(countall+size(effects,1)),1:size(effects,2),2:10,stim)=effects(:,1:end,2:10);
                    effectsall((1+countall):(countall+size(effects,1)),1:size(effects,2),1,stim)=effects(:,1:end,1);
                mousevec((1+countall):(countall+size(effects,1)),stim)=mousetmp;
                else
                    effectsall((1+countall):(countall+size(effects,1)),1:size(effects,2),2:10,stim)=effects(:,1:end,2:10)-effects(:,1:end,2:10);
                    effectsall((1+countall):(countall+size(effects,1)),1:size(effects,2),1,stim)=effects(:,1:end,1)-effects(:,1:end,1);
                mousevec((1+countall):(countall+size(effects,1)),stim)=mousetmp-mousetmp;
                end
            end

            if onestim==0;
                effectsallSTD((1+countall):(countall+size(effects,1)),1:size(effects,2),2:10,stim)=effectsSTD(:,1:end,2:10);
                effectsallSTD((1+countall):(countall+size(effects,1)),1:size(effects,2),1,stim)=effectsSTD(:,1:end,1);
                
            elseif onestim==1;
                if stim==1;
                    effectsallSTD((1+countall):(countall+size(effects,1)),1:size(effects,2),2:10,stim)=effectsSTD(:,1:end,2:10);
                    effectsallSTD((1+countall):(countall+size(effects,1)),1:size(effects,2),1,stim)=effectsSTD(:,1:end,1);
                    
                else
                    effectsallSTD((1+countall):(countall+size(effects,1)),1:size(effects,2),2:10,stim)=effectsSTD(:,1:end,2:10)-effectsSTD(:,1:end,2:10);
                    effectsallSTD((1+countall):(countall+size(effects,1)),1:size(effects,2),1,stim)=effectsSTD(:,1:end,1)-effectsSTD(:,1:end,1);
                    
                end
            elseif onestim==2;
                if stim==2;
                    effectsallSTD((1+countall):(countall+size(effects,1)),1:size(effects,2),2:10,stim)=effectsSTD(:,1:end,2:10);
                    effectsallSTD((1+countall):(countall+size(effects,1)),1:size(effects,2),1,stim)=effectsSTD(:,1:end,1);
                  
                else
                    effectsallSTD((1+countall):(countall+size(effects,1)),1:size(effects,2),2:10,stim)=effectsSTD(:,1:end,2:10)-effectsSTD(:,1:end,2:10);
                    effectsallSTD((1+countall):(countall+size(effects,1)),1:size(effects,2),1,stim)=effectsSTD(:,1:end,1)-effectsSTD(:,1:end,1);
                    
                end
            end
            
        end
        
        countall=countall+size(effects,1);
    end
    
    for build=1
        
        for stim=1:2
            limtrials=limtrials1;
            if limtrials1>0;
                vecdat=vectrials{stim,1};
                vecstim=vectrials{stim,2};
                for checkingfortoo_few_trials=1;
                    ss=[1 2 3 4 15 phzstd(stim)];
                    for gss=1:6;
                        sss(gss)=sum(vecstim==ss(gss));
                    end
                    if min(sss)<limtrials1;
                        limtrials=min(sss);
                    end
                  
                end
                
                effects=[];
                for z1=1:8
                    tmp1=vecdat(:,:,vecstim==z1);
                    ex=[];
                    for t1=1:size(tmp1,3)
                        if rem(t1,2)==1;
                            ex(t1)=1;
                        else
                            ex(t1)=0;
                        end
                    end
                    if dobeforeafter==1; pt=ceil(limtrials/2); clear ex; ex(1,1:pt)=ones(1,pt); try ex(1,(pt+1):limtrials)=zeros(1,pt); catch ex(1,(pt+1):limtrials)=zeros(1,pt-1); end; end
                    
                    if size(tmp1,3)>(limtrials-1)
                        atmp=(tmp1(:,:,1:limtrials));
                        ex=ex(1:limtrials);
                        effects(:,:,z1+1)=mean(atmp(:,:,ex==1),3);
                    else
                        atmp=(tmp1(:,:,1:end));
                        ex=ex(1:size(atmp,3));
                        effects(:,:,z1+1)=mean(atmp(:,:,ex==1),3);
                    end
                end
                
                tmp1=vecdat(:,:,vecstim==15);
                ex=[];
                for t1=1:size(tmp1,3)
                    if rem(t1,2)==1
                        ex(t1)=1;
                    else
                        ex(t1)=0;
                    end
                end
                if dobeforeafter==1; pt=ceil(limtrials/2); clear ex; ex(1,1:pt)=ones(1,pt); try ex(1,(pt+1):limtrials)=zeros(1,pt); catch ex(1,(pt+1):limtrials)=zeros(1,pt-1); end; end
                
                if size(tmp1,3)>(limtrials-1)
                    atmp=(tmp1(:,:,1:limtrials));
                    ex=ex(1:limtrials);
                    effects(:,:,10)=mean(atmp(:,:,ex==1),3);
                else
                    atmp=(tmp1(:,:,1:end));
                    ex=ex(1:size(atmp,3));
                    effects(:,:,10)=mean(atmp(:,:,ex==1),3);
                end
                tmp1=vecdat(:,:,vecstim==phzstd(stim));
                ex=[];
                for t1=1:size(tmp1,3)
                    if rem(t1,2)==1;
                        ex(t1)=1;
                    else
                        ex(t1)=0;
                    end
                end
                if dobeforeafter==1; pt=ceil(limtrials/2); clear ex; ex(1,1:pt)=ones(1,pt); try ex(1,(pt+1):limtrials)=zeros(1,pt); catch ex(1,(pt+1):limtrials)=zeros(1,pt-1); end; end
                
                atmp=(tmp1(:,:,1:limtrials));
                ex=ex(1:limtrials);
                effects(:,:,1)=mean(atmp(:,:,ex==1),3);
            else
                vecdat=vectrials{stim,1};
                vecstim=vectrials{stim,2};
                effects=[];
                stimax=[phzstd(stim) 1 2 3 4 5 6 7 8 15];
                for z1=1:10;
                    tmp1=vecdat(:,:,vecstim==stimax(z1));
                    ex=[];
                    for t1=1:size(tmp1,3)
                        if rem(t1,2)==1;
                            ex(t1)=1;
                        else
                            ex(t1)=0;
                        end
                    end
                    tmp1=tmp1(:,:,ex==1);
                    effects(:,:,z1)=mean(tmp1,3);
                end
            end
                
                
            for c1=1:size(effects,1) %baseline correction
                
                
                for z1=1:size(effects,3)
               
                    mnn=mnns(file,stim,c1,z1);
                    stt=stts(file,stim,c1,z1);
                    for t1=1:size(effects,2)
                        effects(c1,t1,z1)=(effects(c1,t1,z1)-mnn)./stt;
                    end
                end
                
            end
            
            if onestim==0;
                effectsallBUILD((1+count):(count+size(effects,1)),1:size(effects,2),2:10,stim)=effects(:,1:end,2:10);
                effectsallBUILD((1+count):(count+size(effects,1)),1:size(effects,2),1,stim)=effects(:,1:end,1);
            elseif onestim==1;
                if stim==1;
                    effectsallBUILD((1+count):(count+size(effects,1)),1:size(effects,2),2:10,stim)=effects(:,1:end,2:10);
                    effectsallBUILD((1+count):(count+size(effects,1)),1:size(effects,2),1,stim)=effects(:,1:end,1);
                else
                    effectsallBUILD((1+count):(count+size(effects,1)),1:size(effects,2),2:10,stim)=effects(:,1:end,2:10)-effects(:,1:end,2:10);
                    effectsallBUILD((1+count):(count+size(effects,1)),1:size(effects,2),1,stim)=effects(:,1:end,1)-effects(:,1:end,1);
                end
            elseif onestim==2;
                if stim==2;
                    effectsallBUILD((1+count):(count+size(effects,1)),1:size(effects,2),2:10,stim)=effects(:,1:end,2:10);
                    effectsallBUILD((1+count):(count+size(effects,1)),1:size(effects,2),1,stim)=effects(:,1:end,1);
                else
                    effectsallBUILD((1+count):(count+size(effects,1)),1:size(effects,2),2:10,stim)=effects(:,1:end,2:10)-effects(:,1:end,2:10);
                    effectsallBUILD((1+count):(count+size(effects,1)),1:size(effects,2),1,stim)=effects(:,1:end,1)-effects(:,1:end,1);
                end
            end
            
        end
        
        count=count+size(effects,1);
    end
    
    
    for test=1
        samp=1/framerate;
        for stim=1:2
            limtrials=limtrials1;
            if limtrials1>0;
                vecdat=vectrials{stim,1};
                vecstim=vectrials{stim,2};
                for checkingfortoo_few_trials=1;
                    ss=[1 2 3 4 15 phzstd(stim)];
                    for gss=1:6;
                        sss(gss)=sum(vecstim==ss(gss));
                    end
                    if min(sss)<limtrials1;
                        limtrials=min(sss);
                    end
                   
                end
                
                effects=[];
                for z1=1:8
                    tmp1=vecdat(:,:,vecstim==z1);
                    ex=[];
                    for t1=1:size(tmp1,3)
                        if rem(t1,2)==1;
                            ex(t1)=0;
                        else
                            ex(t1)=1;
                        end
                    end
                    if dobeforeafter==1; pt=ceil(limtrials/2); clear ex; ex(1,1:pt)=zeros(1,pt); try ex(1,(pt+1):limtrials)=ones(1,pt); catch ex(1,(pt+1):limtrials)=ones(1,pt-1); end; end
                    
                    if size(tmp1,3)>(limtrials-1)
                        atmp=(tmp1(:,:,1:limtrials));
                        ex=ex(1:limtrials);
                        effects(:,:,z1+1)=mean(atmp(:,:,ex==1),3);
                    else
                        atmp=(tmp1(:,:,1:end));
                        ex=ex(1:size(atmp,3));
                        effects(:,:,z1+1)=mean(atmp(:,:,ex==1),3);
                    end
                end
                
                tmp1=vecdat(:,:,vecstim==15);
                ex=[];
                for t1=1:size(tmp1,3)
                    if rem(t1,2)==1;
                        ex(t1)=0;
                    else
                        ex(t1)=1;
                    end
                end
                if dobeforeafter==1; pt=ceil(limtrials/2); clear ex; ex(1,1:pt)=zeros(1,pt); try ex(1,(pt+1):limtrials)=ones(1,pt); catch ex(1,(pt+1):limtrials)=ones(1,pt-1); end; end
                
                if size(tmp1,3)>(limtrials-1)
                    atmp=(tmp1(:,:,1:limtrials));
                    ex=ex(1:limtrials);
                    effects(:,:,10)=mean(atmp(:,:,ex==1),3);
                else
                    atmp=(tmp1(:,:,1:end));
                    ex=ex(1:size(atmp,3));
                    effects(:,:,10)=mean(atmp(:,:,ex==1),3);
                end
                tmp1=vecdat(:,:,vecstim==phzstd(stim));
                ex=[];
                for t1=1:size(tmp1,3)
                    if rem(t1,2)==1;
                        ex(t1)=0;
                    else
                        ex(t1)=1;
                    end
                end
                if dobeforeafter==1; pt=ceil(limtrials/2); clear ex; ex(1,1:pt)=zeros(1,pt); try ex(1,(pt+1):limtrials)=ones(1,pt); catch ex(1,(pt+1):limtrials)=ones(1,pt-1); end; end
                
                atmp=(tmp1(:,:,1:limtrials));
                ex=ex(1:limtrials);
                effects(:,:,1)=mean(atmp(:,:,ex==1),3);
            else
                vecdat=vectrials{stim,1};
                vecstim=vectrials{stim,2};
                effects=[];
                stimax=[phzstd(stim) 1 2 3 4 5 6 7 8 15];
                for z1=1:10;
                    tmp1=vecdat(:,:,vecstim==stimax(z1));
                    ex=[];
                    for t1=1:size(tmp1,3)
                        if rem(t1,2)==1;
                            ex(t1)=0;
                        else
                            ex(t1)=1;
                        end
                    end
                    tmp1=tmp1(:,:,ex==1);
                    effects(:,:,z1)=mean(tmp1,3);
                end
            end
            
            for c1=1:size(effects,1) %baseline correction
                
                
                for z1=1:size(effects,3)
                  mnn=mnns(file,stim,c1,z1);
                    stt=stts(file,stim,c1,z1);
                    for t1=1:size(effects,2)
                        effects(c1,t1,z1)=(effects(c1,t1,z1)-mnn)./stt;
                    end
                end
                
            end
            
            
            if onestim==0;
                effectsallTEST((1+counttest):(counttest+size(effects,1)),1:size(effects,2),2:10,stim)=effects(:,1:end,2:10);
                effectsallTEST((1+counttest):(counttest+size(effects,1)),1:size(effects,2),1,stim)=effects(:,1:end,1);
            elseif onestim==1;
                if stim==1;
                    effectsallTEST((1+counttest):(counttest+size(effects,1)),1:size(effects,2),2:10,stim)=effects(:,1:end,2:10);
                    effectsallTEST((1+counttest):(counttest+size(effects,1)),1:size(effects,2),1,stim)=effects(:,1:end,1);
                else
                    effectsallTEST((1+counttest):(counttest+size(effects,1)),1:size(effects,2),2:10,stim)=effects(:,1:end,2:10)-effects(:,1:end,2:10);
                    effectsallTEST((1+counttest):(counttest+size(effects,1)),1:size(effects,2),1,stim)=effects(:,1:end,1)-effects(:,1:end,1);
                end
            elseif onestim==2;
                if stim==2;
                    effectsallTEST((1+counttest):(counttest+size(effects,1)),1:size(effects,2),2:10,stim)=effects(:,1:end,2:10);
                    effectsallTEST((1+counttest):(counttest+size(effects,1)),1:size(effects,2),1,stim)=effects(:,1:end,1);
                else
                    effectsallTEST((1+counttest):(counttest+size(effects,1)),1:size(effects,2),2:10,stim)=effects(:,1:end,2:10)-effects(:,1:end,2:10);
                    effectsallTEST((1+counttest):(counttest+size(effects,1)),1:size(effects,2),1,stim)=effects(:,1:end,1)-effects(:,1:end,1);
                end
            end
            
        end
        
        counttest=counttest+size(effects,1);
    end
    
    
end
%close all
titles={'Control';' ';' ';' ';'redundants';' ';' ';' ';' ';'DEVIANT'};
redundantnumberA=3; % this is the redundant in the chain that you will plot/analyse; if it's 0, then average over all redundants
include_cutoff=2.6; %number of standard deviations above baseline for inclusion
exclude_high=100; %get rid of the cell if any responses are this big
exclude_low=-20; %get rid of the cell if any responses are this negative

%if you want to just plot averages accross all trials, make stringent test
%%% equal to 0. if you want to do a hard test of clusters,make it 1... use
%%% 1 only if you are actually testing teh reality of a subcluster
if doingmatch==1; include_cutoff=-100; end;
for process_step_2=1;
        %%%%%%%%%figures out relative responses to different stimuli
        simpleeffects=[]; simpleeffectsPRE=[];simpleeffectsPOST=[];c1=1; c2=1; cols=[]; newindices=[];notindices=[];
        for stim=1:2
            for c=1:size(effectsallSTD,1)
                if redundantnumberA==0; 
                	r=mean(mean(effectsallSTD(c,BEGINo:ENDo,2:9,stim)));
                else
                    r=mean(effectsallSTD(c,BEGINo:ENDo,redundantnumberA+1,stim));
                end
                cn=mean(effectsallSTD(c,BEGINo:ENDo,1,stim));
                d=mean(effectsallSTD(c,BEGINo:ENDo,10,stim)); 
                if redundantnumberA==0; 
                	rff=mean(mean(effectsallSTD(c,BEGINoff:ENDoff,2:9,stim)));
                else
                    rff=mean(effectsallSTD(c,BEGINoff:ENDoff,redundantnumberA+1,stim));
                end
                cnff=mean(effectsallSTD(c,BEGINoff:ENDoff,1,stim)); %cn=cnff; r=rff;
                dff=mean(effectsallSTD(c,BEGINoff:ENDoff,10,stim)); %d=dff;
                rff=r; dff=d; cnff=cn; %forget the offsets 
                 if redundantnumberA==0; 
                	rP=mean(mean(effectsallBUILD(c,BEGINo:ENDo,2:9,stim)));
                else
                    rP=mean(effectsallBUILD(c,BEGINo:ENDo,redundantnumberA+1,stim));
                end
                cnP=mean(effectsallBUILD(c,BEGINo:ENDo,1,stim));
                dP=mean(effectsallBUILD(c,BEGINo:ENDo,10,stim));
                if redundantnumberA==0; 
                	rT=mean(mean(effectsallTEST(c,BEGINo:ENDo,2:9,stim)));
                else
                    rT=mean(effectsallTEST(c,BEGINo:ENDo,redundantnumberA+1,stim));
                end
                cnT=mean(effectsallTEST(c,BEGINo:ENDo,1,stim));
                dT=mean(effectsallTEST(c,BEGINo:ENDo,10,stim));
                if and(and(or(or(or(r>include_cutoff,cn>include_cutoff),d>include_cutoff),or(or(rff>include_cutoff,cnff>include_cutoff),dff>include_cutoff)),...
                        and(and(r<exclude_high,cn<exclude_high),d<exclude_high)),and(and(r>exclude_low,cn>exclude_low),d>exclude_low))
                    simpleeffects(c1,1)=(cn);
                    simpleeffects(c1,2)=(r);
                    simpleeffects(c1,3)=(d);
                    simpleeffectsPRE(c1,1)=(cnP);
                    simpleeffectsPRE(c1,2)=(rP);
                    simpleeffectsPRE(c1,3)=(dP);
                    simpleeffectsPOST(c1,1)=(cnT);
                    simpleeffectsPOST(c1,2)=(rT);
                    simpleeffectsPOST(c1,3)=(dT);
                    newindices(c1,1)=c;
                    newindices(c1,2)=stim;
                    c1=1+c1;
                else
                    notindices(c2,1)=c;
                    notindices(c2,2)=stim;
                    c2=1+c2;
                    
                end
                
            end
        end
        for makingcolors=1;
            for c=1:(c1-1)
                s=simpleeffects(c,1)./max(simpleeffects(c,:));
                cols(c,2)=s;
                s=simpleeffects(c,3)./max(simpleeffects(c,:));
                cols(c,1)=s;
                s=simpleeffects(c,2)./max(simpleeffects(c,:));
                cols(c,3)=s;
            end
            for c=1:(c1-1);
                a=cols(c,:)-min(cols(c,:));
                cols(c,:)=cols(c,:)./max(cols(c,:));
            end
        end
       
 
end

for countcellstotal=1; % for only focusing on one response per cell
    numcells=size(unique(newindices(:,1)),1);
    disp(strcat('all cells=',num2str(size(effectsall,1)),'; responsive=',num2str(numcells),'; percent=',num2str(numcells/size(effectsall,1))));
    ex1=zeros(size(simpleeffects,1),1);
    for c=1:size(effectsall,1);
        if sum(newindices(:,1)==c)>1;
            clear tmp1 tmp2; cnt1=1;
            for c1=1:size(newindices,1);
                if newindices(c1,1)==c;
                    tmp1(cnt1,1)=mean(simpleeffects(c1,[1 2 3]));
                    tmp2(cnt1,1)=c1;
                    cnt1=1+cnt1;
                end
            end
            if tmp1(1,1)>tmp1(2,1);
                ex1(tmp2(2,1),1)=1;
            else
                ex1(tmp2(1,1),1)=1;
            end
        end
    end
    if onestim>0
        if onestim==2;
            ex1=vertcat(ones(sum(newindices(:,2)==1),1),zeros(sum(newindices(:,2)==2),1));
        elseif onestim==1;
            
            ex1=vertcat(zeros(sum(newindices(:,2)==1),1),ones(sum(newindices(:,2)==2),1));;
        end
    end
    for onlylookatonePERcell=1;
        if onlylookatonePERcell==1;
            simpleeffectsPRE=simpleeffectsPRE(ex1==0,:);
            simpleeffects=simpleeffects(ex1==0,:);
            simpleeffectsPOST=simpleeffectsPOST(ex1==0,:);
            newindices=newindices(ex1==0,:);
        end
    end
    for makingcolors=1;
        c1=size(simpleeffects,1)+1; clear cols
        for c=1:(c1-1)
            s=simpleeffects(c,1)./max(simpleeffects(c,:));
            cols(c,2)=s;
            s=simpleeffects(c,3)./max(simpleeffects(c,:));
            cols(c,1)=s;
            s=simpleeffects(c,2)./max(simpleeffects(c,:));
            cols(c,3)=s;

        end
        for c=1:(c1-1);
            %a=cols(c,:)-min(cols(c,:));
            cols(c,:)=cols(c,:)./max(cols(c,:));
            cols(c,2)=0;%cols(c,3)=0;%cols(c,1)=0;
        end
    end

end

effectsall=effectsallSTD;
redundantnumber=3; %if you want to plot (but not "select by") a redudnant later in the stream, use this line
for plottingeverything=1; 
%%%%for selecting all or subset of cells for subesequent plots
%%%for selecting all or subset of cells for subesequent plots
plotanti=0; %if =1, then plot everything outside of the selected box
selections=1; %if the cells you want to select are scattered more than a single box can capture, make this 2 or 3 and select multiple regions.

for plotbOLTH=1;

    pind=[];
    for selecto=1:selections;
        figure; scatter(simpleeffects(:,1),simpleeffects(:,3),50,cols,'fill'); hold on
        ylabel('resp to dev','FontSize',12,'FontWeight','bold'); xlabel('resp to ctrl','FontSize',12,'FontWeight','bold');
        clear x;  tmp1=simpleeffects(:,1); tmp2=simpleeffects(:,3); aa=horzcat(tmp1,tmp2); x(1)=min(aa(:))*1.02; x(2)=max(aa(:))*1.02;
        y = x; plot(x,y);
       
        if selecto>1;
            hold on;
            scatter([x11(1) x11(1) x11(2) x11(2)],[y11(1) y11(2) y11(1) y11(2)]);
        end
        figure; scatter((.5*simpleeffects(:,1)+.5*simpleeffects(:,3))-simpleeffects(:,2),simpleeffects(:,3)-simpleeffects(:,1),50,cols,'fill');

        ylabel('Deviance Detection','FontSize',12,'FontWeight','bold'); xlabel('Stim specific adaptation','FontSize',12,'FontWeight','bold');
     
       x=[-1000 1000]; y=[1000 -1000];

        for c=1:size(simpleeffects,1);
            if and(((.5*simpleeffects(c,1)+.5*simpleeffects(c,3))-simpleeffects(c,2))>x(1),((.5*simpleeffects(c,1)+.5*simpleeffects(c,3))-simpleeffects(c,2))<x(2));
                if and((simpleeffects(c,3)-simpleeffects(c,1))<y(1),(simpleeffects(c,3)-simpleeffects(c,1))>y(2))
                    if sum(c==pind)==0;
                        pind=vertcat(pind,c);
                    end
                end
            end
        end
    end
    effectsallBOTH=effectsall; %(effectsallBUILD+effectsallTEST)/2;
    tmpeff=[]; mousey=[];
    for c=1:numcells;
        tmpeff(c,:,:)=effectsallBOTH(newindices(pind(c),1),:,:,(newindices(pind(c),2)));
        %THIS part is for developing the within-obs errorbars
        for t=1:size(effectsallBOTH,2)
            if redundantnumber>0;
                mn=mean(mean(effectsallBOTH(newindices(pind(c),1),t,[redundantnumber+1 1 10],newindices(pind(c),2)),3),4);
                a1(c,t)=mean(effectsallBOTH(newindices(pind(c),1),t,redundantnumber+1,newindices(pind(c),2)),4)-mn;
            else
                
                a1tmp=mean(mean(effectsallBOTH(newindices(pind(c),1),t,2:9,newindices(pind(c),2)),3),4);
                a2tmp=mean(effectsallBOTH(newindices(pind(c),1),t,1,newindices(pind(c),2)),4);
                a3tmp=mean(effectsallBOTH(newindices(pind(c),1),t,10,newindices(pind(c),2)),4);
                mn=(a1tmp+a2tmp+a3tmp)/3;
                a1(c,t)=mean(mean(effectsallBOTH(newindices(pind(c),1),t,2:9,newindices(pind(c),2)),3),4)-mn;
                mousey(c)=mousevec(newindices(pind(c),1),1);
            end
            a2(c,t)=mean(effectsallBOTH(newindices(pind(c),1),t,1,newindices(pind(c),2)),4)-mn;
            a3(c,t)=mean(effectsallBOTH(newindices(pind(c),1),t,10,newindices(pind(c),2)),4)-mn;
                mousey(c,1)=mousevec(newindices(pind(c),1),1);
        end
    end
   

    x=timeaxis;
    y1=mean(mean(mean(tmpeff(:,:,redundantnumber+1),3),4),1);
    if redundantnumber==0; y1=mean(mean(mean(tmpeff(:,:,2:9),3),4),1); end;
    y1a=y1;%
    err1=std(a1,0,1)./sqrt(numcells);
    
    y2=mean(mean(mean(tmpeff(:,:,1),3),4),1);
    y2a=y2;%
    err2=std(a2,0,1)./sqrt(numcells);
    
    y3=mean(mean(mean(tmpeff(:,:,10),3),4),1); 
    y3a=y3;%
    err3=std(a3,0,1)./sqrt(numcells);
 
    if numcells>1;
        figure; shadedErrorBar(x,y1a,err1,'lineProps','b'); hold on;shadedErrorBar(x,y2a,err2,'lineProps','k');hold on;shadedErrorBar(x,y3a,err3,'lineProps','r');
        xlim([-.2 1]); title ('Average of cell responses', 'FontSize', 16, 'FontWeight', 'bold');ylabel('z-scores','FontSize',16,'FontWeight','bold'); xlabel('time (sec)','FontSize',16,'FontWeight','bold'); set(gca,'FontSize',16,'FontWeight','bold');
        set(gcf,'Color','w');
        
        yy=ylim;
        figure; errorbar_groups([mean(mean(y2a(:,BEGINo:ENDo))) mean(mean(y1a(:,BEGINo:ENDo))) mean(mean(y3a(:,BEGINo:ENDo)))]',[mean(mean(err2(:,BEGINo:ENDo))) mean(mean(err1(:,BEGINo:ENDo))) mean(mean(err3(:,BEGINo:ENDo)))]', 'bar_colors',[.5 .5 .5; 0 0 1; 1 0 0]);
        ylabel('z-scores','FontSize',16,'FontWeight','bold'); xlabel('time (sec)','FontSize',16,'FontWeight','bold'); set(gca,'FontSize',16,'FontWeight','bold');
        set(gcf,'Color','w');
        for statstat=1;
            r1=mean(mean(tmpeff(:,BEGINo:ENDo,redundantnumber+1),3),2)-mean(mean(tmpeff(:,baseBEGIN:baseEND,redundantnumber+1),3),2);
            if redundantnumber==0;
                 r1=mean(mean(tmpeff(:,BEGINo:ENDo,2:9),3),2)-mean(mean(tmpeff(:,baseBEGIN:baseEND,2:9),3),2);
            end
            c1=mean(mean(tmpeff(:,BEGINo:ENDo,1),3),2)-mean(mean(tmpeff(:,baseBEGIN:baseEND,1),3),2);
            d1=mean(mean(tmpeff(:,BEGINo:ENDo,10),3),2)-mean(mean(tmpeff(:,baseBEGIN:baseEND,10),3),2);
            
            [h,p1,ci,stats1]=ttest(r1,c1);
            [h,p2,ci,stats2]=ttest(d1,c1);

            for dolme=1;


                % Stack the dependent variable
                Y = [
                    c1; ...
                    d1; ];
                n=size(d1,1);

                % Fixed-effect factors
                Context = categorical([ ...
                    repmat("control", n, 1); ...
                    repmat("deviant", n, 1); ]);
                cellID = categorical(horzcat(1:n,1:n)');
                mouse1=categorical(vertcat(mousey,mousey));

                % Create table
                tbl = table(Y, Context, cellID, mouse1);

                % Fit linear mixed effects model
                lme = fitlme(tbl, ...
                    'Y ~ Context + (1|mouse1) + (1|cellID)');

                % View full model summary
                fe = lme.Coefficients;
                disp(fe)
                anova_results = anova(lme, 'DFMethod', 'satterthwaite');
            end
            disp(anova_results)
        end
        title(strcat('T_S_S_A','=',num2str(round(stats1.tstat*100)/100),', p=',num2str(round(p1*1000)/1000),'; T_D_D','=',num2str(round(stats2.tstat*100)/100),', p=',num2str(round(p2*1000)/1000)));
    else
        figure; plot(x,y1,'b','LineWidth',3); hold on;plot(x,y2,'k','LineWidth',3); hold on;plot(x,y3,'r','LineWidth',3);yy=ylim;
        xlim([-.2 1]); title ('Average of cell responses', 'FontSize', 16, 'FontWeight', 'bold');
    end
    
    for manyplotthing=1;
        figure;
        for cnd=1:10;
            y1=mean(tmpeff(:,:,cnd),1);
            y1=y1-mean(y1(1,baseBEGIN:baseEND));
            if cnd==1; err1=std(a2,0,1)./sqrt(numcells);
                err=std(a2,0,1)./sqrt(numcells);
                
            elseif cnd==10;
                err=std(a3,0,1)./sqrt(numcells);
            else
                err=std(a1,0,1)./sqrt(numcells);
            end
            subplot(1,10,cnd); shadedErrorBar(x,y1,err,'lineProps','k'); set(gcf,'Color','w'); xlim([-.2 .8]);
            if cnd==1;
                ylabel('z-scores','FontSize',16,'FontWeight','bold'); xlabel('time (sec)','FontSize',16,'FontWeight','bold'); set(gca,'FontSize',16,'FontWeight','bold');
                title(titles{cnd}, 'FontSize', 18, 'FontWeight', 'bold');
            else
                title(titles{cnd}, 'FontSize', 18, 'FontWeight', 'bold'); yticks([]);set(gca,'FontSize',16,'FontWeight','bold');
            end
        end
        for fixaxes=1
            
            ymin1=1; ymax1=-1;
            for cnd=1:10
                subplot(1,10,cnd); yy=ylim;
                if yy(1)<ymin1; ymin1=yy(1); end
                if yy(2)>ymax1; ymax1=yy(2); end
                
            end
            for cnd=1:10
                subplot(1,10,cnd); ylim([ymin1 ymax1]);
            end
            
        end
    end
    
    for manyplotthing=2;
        
        if redundantnumber==0;
            figure;
            ccnds={1;2:9;10};
            cccf=[1 5 10];
            
            colars={'k' 'b' 'r'};
            for cnd=1:3;
                y1=mean(mean(tmpeff(:,:,ccnds{cnd}),3),1);
                y1=y1-mean(y1(1,baseBEGIN:baseEND));
                if cnd==1; err1=std(a2,0,1)./sqrt(numcells/2);
                    err=std(a2,0,1)./sqrt(numcells/2);
                    
                elseif cnd==3;
                    err=std(a3,0,1)./sqrt(numcells/2);
                else
                    err=std(a1,0,1)./sqrt(numcells/2);
                end
                subplot(1,3,cnd); shadedErrorBar(x,y1,err,'lineProps',colars{cnd}); set(gcf,'Color','w'); xlim([-.2 .8]);
                if cnd==1;
                    ylabel('z-scores','FontSize',16,'FontWeight','bold'); xlabel('time (sec)','FontSize',16,'FontWeight','bold'); set(gca,'FontSize',16,'FontWeight','bold');
                    title(titles{cccf(cnd)}, 'FontSize', 18, 'FontWeight', 'bold');
                    cont1=y1./max(y1); cont1err=err; y11=y1;
                else
                    if cnd==2;
                        title(titles{5}, 'FontSize', 18, 'FontWeight', 'bold'); yticks([]);set(gca,'FontSize',16,'FontWeight','bold');
                        rdnt1=y1./max(y11); rdnt1err=err;
                    else
                        title(titles{cccf(cnd)}, 'FontSize', 18, 'FontWeight', 'bold'); yticks([]);set(gca,'FontSize',16,'FontWeight','bold');
                        dev1=y1./max(y11); dev1err=err;
                        
                    end
                end
            end
            for fixaxes=1
                
                ymin1=1; ymax1=-1;
                for cnd=1:3
                    subplot(1,3,cnd); yy=ylim;
                    if yy(1)<ymin1; ymin1=yy(1); end
                    if yy(2)>ymax1; ymax1=yy(2); end
                    
                end
                for cnd=1:3
                    subplot(1,3,cnd); ylim([ymin1 ymax1]);
                end
                
            end
        else
            figure;
            ccnds=[1 redundantnumber+1 10];
            
            colars={'k' 'b' 'r'};
            for cnd=1:3;
                y1=mean(mean(tmpeff(:,:,ccnds(cnd)),3),1);
                y1=y1-mean(y1(1,baseBEGIN:baseEND));
                if cnd==1; err1=std(a1,0,1)./sqrt(numcells/2);
                    err=std(a1,0,1)./sqrt(numcells/2);
                    
                elseif cnd==3;
                    err=std(a3,0,1)./sqrt(numcells/2);
                else
                    err=std(a2,0,1)./sqrt(numcells/2);
                end
                subplot(1,3,cnd); shadedErrorBar(x,y1,err,'lineProps',colars{cnd}); set(gcf,'Color','w'); xlim([-.2 .8]);
                if cnd==1;
                    ylabel('z-scores','FontSize',16,'FontWeight','bold'); xlabel('time (sec)','FontSize',16,'FontWeight','bold'); set(gca,'FontSize',16,'FontWeight','bold');
                    title(titles{ccnds(cnd)}, 'FontSize', 18, 'FontWeight', 'bold');
                    cont1=y1./max(y1); cont1err=err; y11=y1;
                else
                    if cnd==2;
                        title(titles{5}, 'FontSize', 18, 'FontWeight', 'bold'); yticks([]);set(gca,'FontSize',16,'FontWeight','bold');
                        rdnt1=y1./max(y11); rdnt1err=err;
                    else
                        title(titles{ccnds(cnd)}, 'FontSize', 18, 'FontWeight', 'bold'); yticks([]);set(gca,'FontSize',16,'FontWeight','bold');
                        dev1=y1./max(y11); dev1err=err;
                    end
                end
            end
            for fixaxes=1
                
                ymin1=1; ymax1=-1;
                for cnd=1:3
                    subplot(1,3,cnd); yy=ylim;
                    if yy(1)<ymin1; ymin1=yy(1); end
                    if yy(2)>ymax1; ymax1=yy(2); end
                    
                end
                for cnd=1:3
                    subplot(1,3,cnd); ylim([ymin1 ymax1]);
                end
                
            end
        end
    end
     
    
 end



cols3=cols-cols;
try 
   DDs1=simpleeffectsPRE(:,3)-(mean(simpleeffectsPRE(:,1),2));DDs2=simpleeffectsPOST(:,3)-(mean(simpleeffectsPOST(:,1),2));
   DDs3=simpleeffectsPRE(:,3)-(mean(simpleeffectsPRE(:,2),2));DDs4=simpleeffectsPOST(:,3)-(mean(simpleeffectsPOST(:,2),2));
    sm=(DDs1>1.97)+(DDs2>1.97);%+(DDs3>1.67)+(DDs4>1.67);%
    figure; 
    ax = gca(); 
pieData = [mean(sm==2) mean(sm~=2)]; 
h = pie(ax, pieData); 
newColors = [...
    1,       0, 0;   
    .5,       .5,       0.5];  
    ax.Colormap = newColors; 
    colDD=cols-cols;
    for c=1:size(sm,1)
        if sm(c,1)==4;
            colDD(c,1)=1;
        end
    end
    title(num2str(mean(sm~=2)));
%     for determiningcombinedPvals=1;
%         pp=1;
%         for p1=.01:.01:1;
%             pvals(pp)=chi2pdf((-2)*((log(p1)+log(p1))),4);pp=1+pp;
%         end
%         figure; plot(.01:.01:1,pvals);
%     end
end
%axis(sctax)


   sortby=[1 redundantnumber+1 10];
   BOOSTscalefactor=80; %between 0 and 100. higher is brighter
   figure; 
   
   for rasters=1
       
       effects1=tmpeff;
       for c=1:size(tmpeff,1);
           for cnd=1:10;
               m=mean(effects1(c,baseBEGIN:baseEND,cnd));
               effects1(c,:,cnd)=effects1(c,:,cnd)-m;
           end
       end
       ak1=mean(effects1(:,:,sortby),3); [ix,b]=sort(mean(ak1(:,BEGINo:ENDo),2),'descend');
       effects2=effects1(b,:,:);
       
       for z1=1:size(effects2,3)
           subplot(1,size(effects2,3),z1); imagesc(timeaxis,1:size(effects2,1),effects2(:,:,z1)); colormap('Gray'); set(gca,'FontSize',14,'FontWeight','bold');
           xlim([-.2 1]);
           title(titles{z1});
           
           if z1==1; xlabel('sec','FontSize',14,'FontWeight','bold');
               ylabel('neurons','FontSize',14,'FontWeight','bold');
           end
       end
       for settingcolorscale=1;
           scalefactor=(BOOSTscalefactor/100); if scalefactor==1; scalefactor=.99; end
           tmpmax=0; tmpmin=0;
           for z1=1:size(effects2,3)
               subplot(1,size(effects2,3),z1);
               cb=caxis; if cb(1)<tmpmin; tmpmin=cb(1); end
               if cb(2)>tmpmax; tmpmax=cb(2); end
           end
           for z1=1:size(effects2,3)
               subplot(1,size(effects2,3),z1);
               cb=[tmpmin tmpmax];
               caxis(cb*(1-scalefactor))
           end
       end
   end


end
