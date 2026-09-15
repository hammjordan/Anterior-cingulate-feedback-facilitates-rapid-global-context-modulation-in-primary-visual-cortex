%clcl 
clear all
limtrials1=6; 
exclloco=1; 
for fileload=1;
       
      
     files={'1_sal_mouse137';'2_sal_mouse138';'1_sal_mouse139';'2_sal_mouse140';'1_sal_mouse145';'2_sal_mouse146';'1_sal_mouse152';'2_sal_mouse154';'1_sal_mouse162';'2_sal_mouse163';'1_sal_mouse164';'1_sal_mouse165';
            '2_dcz_mouse137';'1_dcz_mouse138';'2_dcz_mouse139';'1_dcz_mouse140';'2_dcz_mouse145';'1_dcz_mouse146';'2_dcz_mouse152';'1_dcz_mouse154';'2_dcz_mouse162';'1_dcz_mouse163';'2_dcz_mouse164';'2_dcz_mouse165'};%'1_dcz_mouse147';%'2_sal_mouse153';'2_sal_mouse147'; %147 had really weird ACa, so hard to confirm dreadds. 153 was anterior to v1. including both does not change result
     groupsBOTHconfirmed=vertcat(ones(12,1),ones(12,1)+1)';
     ordersALL=[1 2 1 2 1 2 1 2 1 2 1 1 ...
                2 1 2 1 2 1 2 1 2 1 2 2];
     mices=(horzcat(1:12,1:12));
   
end
 
for find_times=1;
    times=[-100 0 65 400 700]; %prestim, zero, midstim, endstim
    load(strcat('GLctrl',files{1}),'t_sG1');
    for find_timeindices=1;
        t1=0; t2=0; t3=0; t4=0; t5=0;
        for t=1:size(t_sG1,2);
            if and(t_sG1(t)>times(1),t1==0);
                t1=t;
            end
        end
        for t=1:size(t_sG1,2);
            if and(t_sG1(t)>times(2),t2==0);
                t2=t;
            end
        end
        for t=1:size(t_sG1,2);
            if and(t_sG1(t)>times(3),t3==0);
                t3=t;
            end
        end
        for t=1:size(t_sG1,2);
            if and(t_sG1(t)>times(4),t4==0);
                t4=t;
            end
        end
         for t=1:size(t_sG1,2);
            if and(t_sG1(t)>times(5),t5==0);
                t5=t;
            end
         end
         t5=t; 
    end
end

for GLoddballanalyse=1;
    limtrialsOB=zeros(size(files,1),6); 
     fh = waitbar(0,'Percent completEE...');
    for fila=1:size(files,1);
        %oris=[2 3 5 4 1 8 6 7]
         load(strcat('GLtest',files{fila}));

         load(strcat('GLtest',files{fila},'loco'));
         if exclloco==1; rej1=max(locobintrials(:,1:1000),[],2)'; rej=rej+rej1; end;
         
          mean(rej)
         orders=[]; for zzz=1:200; orders=horzcat(orders,[1:5]); end; orders=orders(1,1:size(trials1,2));
         ersp1=ersp1(:,:,rej==0);     
         tfs1a=tfs1(:,:,rej==0);
         erp1a=erp1(rej==0,:); 
         
         trials1=horzcat(trials1,zeros(1,6));trials=trials1; 
         trials1(1)=0; for t=6:size(trials1,2)-6; if and(trials(t)==4, trials(t)==4); trials1(t)=4; 
         elseif and(trials(t)==3, trials(t)==3); trials(t)=3; 
             elseif and(trials(t)==1, trials(t+6)==3); trials1(t)=1; 
         elseif and(trials(t)==2, trials(t+5)==3); trials1(t)=2; 
         else trials1(t)=0;
         end; end
        trials1=trials1(1,1:end-6); 
        trials1=trials1(rej==0);
        %part 2 limtrials
        for limtrials=limtrials1; %dont need to do this for the 6 and8, but only 4
            if limtrials>0;
                %then, limit to last 6 occurrances
                counts=zeros(1,8);
                for tr=1:size(trials1,2);
                    tr1=tr;
                    aa=trials1(1,tr1);
                    if aa>0
                        counts(1,aa)=counts(1,aa)+1;
                        if counts(1,aa)>limtrials;
                            trials1(1,tr1)=0;
                        end
                    end
                end
            end
        end

        trnums_oddball(fila,1)=sum(trials1==1);
        trnums_oddball(fila,2)=sum(trials1==2);
        trnums_oddball(fila,3)=sum(trials1==3);
        trnums_oddball(fila,4)=sum(trials1==4);
        erspsV1_oddball(:,:,fila,1)=mean(ersp1(:,:,trials1==1),3);
        erspsV1_oddball(:,:,fila,2)=mean(ersp1(:,:,trials1==2),3); 
        erspsV1_oddball(:,:,fila,3)=mean(ersp1(:,:,trials1==3),3);
        erspsV1_oddball(:,:,fila,4)=mean(ersp1(:,:,trials1==4),3);

        limtrialsOB(fila,1)=trnums_oddball(fila,3); 
        limtrialsOB(fila,2)=trnums_oddball(fila,2);
        limtrialsOB(fila,5)=trnums_oddball(fila,4); 
        waitbar((fila/size(files,1)*.5),fh,'doing oddballs');
    end
end
  
for controlanalyse=1;
  
    for fila=1:size(files,1);
        %oris=[2 3 5 4 1 8 6 7]
         load(strcat('GLctrl',files{fila})); 
         load(strcat('GLctrl',files{fila},'loco')); if exclloco==1; rej1=max(locobintrials(:,1:1000),[],2)'; rej=rej+rej1; end;
        
         mean(rej)
         orders=[]; for zzz=1:200; orders=horzcat(orders,[1:5]); end; orders=orders(1,1:size(trials1,2));
         ersp1=ersp1(:,:,rej==0);       
         tfs1a=tfs1(:,:,rej==0);
         erp1a=erp1(rej==0,:); 
         
         trials1(1)=0; for t=2:size(trials1,2); if and(trials1(t)==1, sum(trials1(t-1)==[1])==0); trials1(t)=1; 
         elseif and(trials1(t)==2, sum(trials1(t-1)==[2])==0); trials1(t)=2; 
             elseif and(trials1(t)==5, sum(trials1(t-1)==[5])==0); trials1(t)=5; else trials1(t)=0;
         end; end
         for t=2:size(trials1,2); if orders(1,t)<4; trials1(1,t)=0; end; end;


        trials1=trials1(rej==0);

         if limtrials1>0;
             counts=zeros(1,8);
             for t=1:size(trials1,2);
                 if ordersALL(1,fila)==2; 
                     tt=t;
                 else %for the recordings from the first day, use the latter trials
                     tt=(size(trials1,2)+1)-t;
                 end
                 if trials1(1,tt)>0
                 counts(1,trials1(1,tt))=counts(1,trials1(1,tt))+1;
                 if counts(1,trials1(1,tt))>limtrialsOB(fila,trials1(1,tt))
                     trials1(1,tt)=0;
                 end
                 end
             end
         end

        trnums_control(fila,1)=sum(trials1==1);
        trnums_control(fila,2)=sum(trials1==2);
        trnums_control(fila,3)=sum(trials1==5);
        erspsV1_control(:,:,fila,1)=mean(ersp1(:,:,trials1==1),3);
        erspsV1_control(:,:,fila,2)=mean(ersp1(:,:,trials1==2),3); 
        erspsV1_control(:,:,fila,3)=mean(ersp1(:,:,trials1==5),3);



        itc=tfs1a(:,:,trials1==1)./abs(tfs1a(:,:,trials1==1));nn=trnums_control(fila,1);
        itcsV1_control(:,:,fila,1)=abs(mean(itc,3))-sqrt(-(1/nn)*log(.5));
         itc=tfs1a(:,:,trials1==2)./abs(tfs1a(:,:,trials1==2));nn=trnums_control(fila,2);
        itcsV1_control(:,:,fila,2)=abs(mean(itc,3))-sqrt(-(1/nn)*log(.5));
         itc=tfs1a(:,:,trials1==5)./abs(tfs1a(:,:,trials1==5));nn=trnums_control(fila,3);
        itcsV1_control(:,:,fila,3)=abs(mean(itc,3))-sqrt(-(1/nn)*log(.5));

        waitbar((fila/size(files,1)*.5)+.5,fh,'doing controls');
    end
end
 
%close all force
for visualresponsees=1;
   
    groups=groupsBOTHconfirmed;
  
    figure;
    sc='linear';
    % visual induced oscillations (power
    for grp=1:2;
        subs=groups==grp;
        for power_to_stimV1_control=1;
            tmp1=erspsV1_control(:,:,subs,1);
            tmp2=erspsV1_control(:,:,subs,2);
            tmp3=erspsV1_control(:,:,subs,3);

            subplot(2,4,1); contourf(t_sG1,f_sG1,mean(tmp1,3),50,'linecolor','none'); set(gca,'yscale', sc);title('control 1'); colormap jet;  axis([-150 735 5 40]);cc=[-2 5];ylabel('control context');   caxis(cc);
            subplot(2,4,2); contourf(t_sG1,f_sG1,mean(tmp1,3),50,'linecolor','none'); set(gca,'yscale', sc);title('control 1');  colormap jet;  axis([-150 735 5 40]); caxis(cc);
            subplot(2,4,3);contourf(t_sG1,f_sG1,mean(tmp2,3),50,'linecolor','none'); set(gca,'yscale', sc); title('control 2'); colormap jet;   axis([-150 735 5 40]);caxis(cc);
            subplot(2,4,4);contourf(t_sG1,f_sG1,mean(tmp3,3),50,'linecolor','none'); set(gca,'yscale', sc);title('control 3'); colormap jet;  axis([-150 735 5 40]);caxis(cc);
            make_eps_saveable        
        end

        for powe_to_stimV1_OB=1;
            tmp0=erspsV1_oddball(:,:,subs,1);
            tmp1=erspsV1_oddball(:,:,subs,2);
            tmp2=erspsV1_oddball(:,:,subs,3);
            tmp3=erspsV1_oddball(:,:,subs,4);

            subplot(2,4,5); contourf(t_sG1,f_sG1,mean(tmp0,3),50,'linecolor','none'); set(gca,'yscale', sc); title('redundant'); colormap jet;  axis([-100 735 5 40]);ylabel('oddball');  caxis(cc);
            subplot(2,4,6);contourf(t_sG1,f_sG1,mean(tmp2,3),50,'linecolor','none'); set(gca,'yscale', sc); title('redundant OB'); colormap jet;   axis([-100 735 5 40]);caxis(cc);
            subplot(2,4,7);contourf(t_sG1,f_sG1,mean(tmp1,3),50,'linecolor','none'); set(gca,'yscale', sc); title('predictable OB'); colormap jet;  axis([-100 735 5 40]);caxis(cc);
            subplot(2,4,8);contourf(t_sG1,f_sG1,mean(tmp3,3),50,'linecolor','none'); set(gca,'yscale', sc); title('global OB'); colormap jet; axis([-100 735 5 40]);caxis(cc);
           make_eps_saveable        
        end

    
    freqs=[4:12 20:40];
   
    for powerlineplots=1;
        subs=groups==grp;

        figure;
        for power_to_stimV1_control=1;

            tmp1c=squeeze(mean(erspsV1_control(freqs,:,subs,1),1))';
            tmp2c=squeeze(mean(erspsV1_control(freqs,:,subs,2),1))';
            tmp3c=squeeze(mean(erspsV1_control(freqs,:,subs,3),1))';

            tmp0=squeeze(mean(erspsV1_oddball(freqs,:,subs,1),1))';
            tmp1=squeeze(mean(erspsV1_oddball(freqs,:,subs,2),1))';
            tmp2=squeeze(mean(erspsV1_oddball(freqs,:,subs,3),1))';
            tmp3=squeeze(mean(erspsV1_oddball(freqs,:,subs,4),1))';
       
                    for sub=1:size(tmp0,1);
                    
                        tmp1c(sub,:)=tmp1c(sub,:)-mean(tmp1c(sub,t1:t2),2);
                        tmp2c(sub,:)=tmp2c(sub,:)-mean(tmp2c(sub,t1:t2),2);
                        tmp3c(sub,:)=tmp3c(sub,:)-mean(tmp3c(sub,t1:t2),2);
                        tmp0(sub,:)=tmp0(sub,:)-mean(tmp0(sub,t1:t2),2);
                        tmp1(sub,:)=tmp1(sub,:)-mean(tmp1(sub,t1:t2),2);
                        tmp2(sub,:)=tmp2(sub,:)-mean(tmp2(sub,t1:t2),2);
                        tmp3(sub,:)=tmp3(sub,:)-mean(tmp3(sub,t1:t2),2);
                    end
             

            mn1=mean(tmp1c); st1=std(tmp1c,0,1)./sqrt(size(tmp1c,1));
            mn2=mean(tmp0); st2=std(tmp0,0,1)./sqrt(size(tmp0,1));
            subplot(1,4,1);shadedErrorBar(t_sG1,mn1,st1,'lineProps','k'); hold on
            shadedErrorBar(t_sG1,mn2,st2,'lineProps','b');  title('redundant');axis([-100 650 -1 3.7]);

            mn1=mean(tmp1c); st1=std(tmp1c,0,1)./sqrt(size(tmp1c,1));
            mn2=mean(tmp2); st2=std(tmp2,0,1)./sqrt(size(tmp2,1));
            subplot(1,4,2);shadedErrorBar(t_sG1,mn1,st1,'lineProps','k'); hold on
            shadedErrorBar(t_sG1,mn2,st2,'lineProps','r');  title('redundant OB');axis([-100 650 -1 3.7]);

            mn1=mean(tmp2c); st1=std(tmp2c,0,1)./sqrt(size(tmp2c,1));
            mn2=mean(tmp1); st2=std(tmp1,0,1)./sqrt(size(tmp1,1));
            subplot(1,4,3);shadedErrorBar(t_sG1,mn1,st1,'lineProps','k'); hold on
            shadedErrorBar(t_sG1,mn2,st2,'lineProps','r'); title('predictable OB');axis([-100 650 -1 3.7]);

            mn1=mean(tmp3c); st1=std(tmp3c,0,1)./sqrt(size(tmp3c,1));
            mn2=mean(tmp3); st2=std(tmp3,0,1)./sqrt(size(tmp3,1));
            subplot(1,4,4);shadedErrorBar(t_sG1,mn1,st1,'lineProps','k'); hold on
            shadedErrorBar(t_sG1,mn2,st2,'lineProps','r'); title('global OB'); axis([-100 650 -1 3.7]);


            make_eps_saveable


        end


    end

    tss=t3:t4; %times
    for barplots=1;
        subs=groups==grp;

        figure;
        for power_to_stimV1_control=1;

            tmp1c=squeeze(mean(erspsV1_control(freqs,:,subs,1),1))';
            tmp2c=squeeze(mean(erspsV1_control(freqs,:,subs,2),1))';
            tmp3c=squeeze(mean(erspsV1_control(freqs,:,subs,3),1))';

            tmp0=squeeze(mean(erspsV1_oddball(freqs,:,subs,1),1))';
            tmp1=squeeze(mean(erspsV1_oddball(freqs,:,subs,2),1))';
            tmp2=squeeze(mean(erspsV1_oddball(freqs,:,subs,3),1))';
            tmp3=squeeze(mean(erspsV1_oddball(freqs,:,subs,4),1))';
         
                    for sub=1:size(tmp0,1);
                      
                        tmp1c(sub,:)=tmp1c(sub,:)-mean(tmp1c(sub,t1:t2),2);
                        tmp2c(sub,:)=tmp2c(sub,:)-mean(tmp2c(sub,t1:t2),2);
                        tmp3c(sub,:)=tmp3c(sub,:)-mean(tmp3c(sub,t1:t2),2);
                        tmp0(sub,:)=tmp0(sub,:)-mean(tmp0(sub,t1:t2),2);
                        tmp1(sub,:)=tmp1(sub,:)-mean(tmp1(sub,t1:t2),2);
                        tmp2(sub,:)=tmp2(sub,:)-mean(tmp2(sub,t1:t2),2);
                        tmp3(sub,:)=tmp3(sub,:)-mean(tmp3(sub,t1:t2),2);
                    end
                
            

            dat1=mean(tmp1c(:,tss),2); mn1=mean(dat1);
            dat2=mean(tmp0(:,tss),2); mn2=mean(dat2);
            st1=std(dat1-dat2,0,1)./sqrt(size(dat1,1));
            subplot(1,4,1); errorbar_groups([mn1 mn2]',[st1 st1]', 'bar_colors',[.5 .5 .5; 0 0 1]); hold on;
            x1=(dat1-dat1)+1.05;x2=(dat2-dat2)+1.95;
            subplot(1,4,1);scatter(x1,dat1,'fill','k');
            subplot(1,4,1);scatter(x2,dat2,'fill','k');
            for sub=1:size(dat1,1);
                subplot(1,4,1); line([x1(sub) x2(sub)],[dat1(sub) dat2(sub)],'Color','k');
            end
            ylim([-4 6]);

            dat1=mean(tmp1c(:,tss),2); mn1=mean(dat1);
            dat2=mean(tmp2(:,tss),2); mn2=mean(dat2);
            st1=std(dat1-dat2,0,1)./sqrt(size(dat1,1));
            subplot(1,4,2); errorbar_groups([mn1 mn2]',[st1 st1]', 'bar_colors',[.5 .5 .5; 0 0 1]); hold on;
            x1=(dat1-dat1)+1.05;x2=(dat2-dat2)+1.95;
            subplot(1,4,2);scatter(x1,dat1,'fill','k');
            subplot(1,4,2);scatter(x2,dat2,'fill','k');
            for sub=1:size(dat1,1);
                subplot(1,4,2); line([x1(sub) x2(sub)],[dat1(sub) dat2(sub)],'Color','k');
            end
            ylim([-4 6]);

            dat1=mean(tmp2c(:,tss),2); mn1=mean(dat1);
            dat2=mean(tmp1(:,tss),2); mn2=mean(dat2);
            st1=std(dat1-dat2,0,1)./sqrt(size(dat1,1));
            subplot(1,4,3); errorbar_groups([mn1 mn2]',[st1 st1]', 'bar_colors',[.5 .5 .5; 0 0 1]); hold on;
            x1=(dat1-dat1)+1.05;x2=(dat2-dat2)+1.95;
            subplot(1,4,3);scatter(x1,dat1,'fill','k');
            subplot(1,4,3);scatter(x2,dat2,'fill','k');
            for sub=1:size(dat1,1);
                subplot(1,4,3); line([x1(sub) x2(sub)],[dat1(sub) dat2(sub)],'Color','k');
            end
           ylim([-4 6]);

            dat1=mean(tmp3c(:,tss),2); mn1=mean(dat1);
            dat2=mean(tmp3(:,tss),2); mn2=mean(dat2);
            st1=std(dat1-dat2,0,1)./sqrt(size(dat1,1));
            subplot(1,4,4); errorbar_groups([mn1 mn2]',[st1 st1]', 'bar_colors',[.5 .5 .5; 0 0 1]); hold on;
            x1=(dat1-dat1)+1.05;x2=(dat2-dat2)+1.95;
            subplot(1,4,4);scatter(x1,dat1,'fill','k');
            subplot(1,4,4);scatter(x2,dat2,'fill','k');
            for sub=1:size(dat1,1);
                subplot(1,4,4); line([x1(sub) x2(sub)],[dat1(sub) dat2(sub)],'Color','k');
            end
            ylim([-4 6]);
            make_eps_saveable

            dat1=mean(tmp1c(:,tss),2);
            dat2=mean(tmp0(:,tss),2);

            dat3=mean(tmp1c(:,tss),2);
            dat4=mean(tmp2(:,tss),2);


            dat5=mean(tmp2c(:,tss),2);
            dat6=mean(tmp1(:,tss),2);


            dat7=mean(tmp3c(:,tss),2);
            dat8=mean(tmp3(:,tss),2);

            anova_rm({[dat5 dat6] [dat7 dat8]})
            [h,p,ci,stats]=ttest(dat1,dat2);disp(strcat('adaptation. t=',num2str(stats.tstat),'. p=',num2str(p)))
            [h,p,ci,stats]=ttest(dat3,dat4); disp(strcat('local rdnt global DD. t=',num2str(stats.tstat),'. p=',num2str(p)))
            [h,p,ci,stats]=ttest(dat5,dat6); disp(strcat('local DD global rdnt. t=',num2str(stats.tstat),'. p=',num2str(p)))
            [h,p,ci,stats]=ttest(dat7,dat8); disp(strcat('local DD global DD. t=',num2str(stats.tstat),'. p=',num2str(p)))

        end
        if grp==1;
            cnt1=dat5;
            cnt2=dat7;
            cnt3=dat1; 
            dev1=dat6;
            dev2=dat8; 
        elseif grp==2;
            cnt1=vertcat(cnt1,dat5);
            cnt2=vertcat(cnt2,dat7);
            cnt3=vertcat(cnt3,dat1);
            dev1=vertcat(dev1,dat6);
            dev2=vertcat(dev2,dat8);
        end

 

    end

    end
 
end


drug=groups;
for LMElova=1;
    mask = ismember(drug, 1:2);
    drug_filt = drug(mask);
    orders_filt = ordersALL(mask);
    mice_filt=mices(mask);
    % Number of sessions
    n = numel(cnt1);   % should be 24

    % Stack the dependent variable
    Y = [ ...
        cnt1(:); ...
        cnt2(:); ...
        dev1(:); ...
        dev2(:)];

    % Fixed-effect factors
    Context = categorical([ ...
        repmat("control", n, 1); ...
        repmat("control", n, 1); ...
        repmat("deviant", n, 1); ...
        repmat("deviant", n, 1)]);

    OddballType = categorical([ ...
        repmat("type1", n, 1); ...
        repmat("type2", n, 1); ...
        repmat("type1", n, 1); ...
        repmat("type2", n, 1)]);

    Drug = categorical(repmat(drug_filt(:), 4, 1));
    Mouse = categorical(repmat(mice_filt(:), 4, 1));
    Order = categorical(repmat(orders_filt(:), 4, 1));

    % Create table
    tbl = table(Y, Drug, Context, OddballType, Order, Mouse);

    % Fit linear mixed effects model
    lme = fitlme(tbl, ...
        'Y ~ Drug*Context*OddballType + Order*Context*OddballType + (1|Mouse)');

    % View full model summary
    fe = lme.Coefficients;
    disp(fe)
    anova_results = anova(lme, 'DFMethod', 'satterthwaite');
        disp(anova_results)
end

for plottin_contrasts=1;

%% Prep
cnt1 = cnt1(:); dev1 = dev1(:);
cnt2 = cnt2(:); dev2 = dev2(:);
cnt3=cnt3(:);

drug   = drug_filt(:);    % 1 or 2
mouse  = mice_filt(:);
orders = orders_filt(:);  % 1 or 2

% Deviant minus control
delta_type1 = dev2 - dev1;
delta_type2 = (cnt1 + cnt2+ cnt3)/3;

% Colors
gray   = [0.6 0.6 0.6];
green  = [0.2 0.7 0.2];
ord1_c = [0.85 0.2 0.2];   % noticeable red
ord2_c = [0.2 0.2 0.85];  % noticeable blue

drugs = [1 2];

%% -------- FIGURE 1: TYPE 1 --------
figure; hold on

% Bar means
mean1 = mean(delta_type1(drug == 1));
mean2 = mean(delta_type1(drug == 2));

% SEM
sem1 = std(delta_type1(drug == 1)) / sqrt(sum(drug == 1));
sem2 = std(delta_type1(drug == 2)) / sqrt(sum(drug == 2));

b = bar(drugs, [mean1 mean2], 'FaceColor','flat');
b.CData(1,:) = gray;
b.CData(2,:) = green;

% Error bars (SEM)
errorbar(drugs, [mean1 mean2], [sem1 sem2], ...
    'k', 'LineStyle','none', 'LineWidth',1.5, 'CapSize',8)

% Connected mouse dots
umice = unique(mouse);

for m = 1:numel(umice)
    idx = mouse == umice(m);

    if sum(idx & drug==1)==1 && sum(idx & drug==2)==1

        y1 = delta_type1(idx & drug==1);
        y2 = delta_type1(idx & drug==2);

        o1 = orders(idx & drug==1);
        o2 = orders(idx & drug==2);

        % connect
        plot([1 2], [y1 y2], '-', ...
            'Color',[0.7 0.7 0.7], 'LineWidth',0.75)

        % marker + color by order
        if o1 == 1
            mk1 = '^'; c1 = ord1_c;
        else
            mk1 = 'o'; c1 = ord2_c;
        end

        if o2 == 1
            mk2 = '^'; c2 = ord1_c;
        else
            mk2 = 'o'; c2 = ord2_c;
        end

        scatter(1, y1, 90, c1, mk1, 'filled', ...
            'MarkerEdgeColor','k', 'LineWidth',1.2)
        scatter(2, y2, 90, c2, mk2, 'filled', ...
            'MarkerEdgeColor','k', 'LineWidth',1.2)
    end
end

xlim([0.5 2.5])
xticks([1 2])
xticklabels({'Drug 1','Drug 2'})
ylabel('\Delta deviant2 minus deviant1')
title('Type 1: global minus local deviant')
box off
set(gca,'FontSize',12)
% Legend (dummy handles)
h1 = scatter(nan, nan, 90, ord1_c, '^', 'filled', ...
    'MarkerEdgeColor','k', 'LineWidth',1.2);
h2 = scatter(nan, nan, 90, ord2_c, 'o', 'filled', ...
    'MarkerEdgeColor','k', 'LineWidth',1.2);

legend([h1 h2], {'Saline first','DCZ first'}, ...
    'Location','best', 'Box','off')
%% -------- FIGURE 2: TYPE 2 --------
figure; hold on

% Bar means
mean1 = mean(delta_type2(drug == 1));
mean2 = mean(delta_type2(drug == 2));

% SEM
sem1 = std(delta_type2(drug == 1)) / sqrt(sum(drug == 1));
sem2 = std(delta_type2(drug == 2)) / sqrt(sum(drug == 2));

b = bar(drugs, [mean1 mean2], 'FaceColor','flat');
b.CData(1,:) = gray;
b.CData(2,:) = green;

% Error bars (SEM)
errorbar(drugs, [mean1 mean2], [sem1 sem2], ...
    'k', 'LineStyle','none', 'LineWidth',1.5, 'CapSize',8)

% Connected mouse dots
for m = 1:numel(umice)
    idx = mouse == umice(m);

    if sum(idx & drug==1)==1 && sum(idx & drug==2)==1

        y1 = delta_type2(idx & drug==1);
        y2 = delta_type2(idx & drug==2);

        o1 = orders(idx & drug==1);
        o2 = orders(idx & drug==2);

        plot([1 2], [y1 y2], '-', ...
            'Color',[0.7 0.7 0.7], 'LineWidth',0.75)

        if o1 == 1
            mk1 = '^'; c1 = ord1_c;
        else
            mk1 = 'o'; c1 = ord2_c;
        end

        if o2 == 1
            mk2 = '^'; c2 = ord1_c;
        else
            mk2 = 'o'; c2 = ord2_c;
        end

        scatter(1, y1, 90, c1, mk1, 'filled', ...
            'MarkerEdgeColor','k', 'LineWidth',1.2)
        scatter(2, y2, 90, c2, mk2, 'filled', ...
            'MarkerEdgeColor','k', 'LineWidth',1.2)
    end
end

xlim([0.5 2.5])
xticks([1 2])
xticklabels({'Drug 1','Drug 2'})
ylabel('mean control responses')
title('average control responses')
box off
set(gca,'FontSize',12)

end


for statsforfirstcontrast=1;    
        tmp1=delta_type1(drug == 1);
        tmp2=delta_type1(drug == 2);
        % Combine data into long format
        mouse   = repmat((1:12)', 2, 1);              % mouse ID, repeated for each condition
        cond    = [repmat({'saline'},12,1); repmat({'drug'},12,1)];
        response = [tmp1(:); tmp2(:)];
        orderVar = [ordersALL(1:12)'; ordersALL(13:24)'];              % same order value per mouse, both conditions

        tbl = table(mouse, cond, orderVar, response, ...
            'VariableNames', {'Mouse','Condition','Order','Response'});
        tbl.Mouse = categorical(tbl.Mouse);
        tbl.Condition = categorical(tbl.Condition);
        tbl.Order = categorical(tbl.Order);   % treat order as grouping factor; drop this line if it's continuous

        % Fit LME: fixed effect of condition and order, random of mouse
        lme = fitlme(tbl, 'Response ~ Condition + Order + (1|Mouse) ');
        anova_results = anova(lme, 'DFMethod', 'satterthwaite');
        disp(anova_results)


        disp(lme)
end
for statsforsecondecontrast=1;    
        tmp1=delta_type2(drug == 1);
        tmp2=delta_type2(drug == 2);
        % Combine data into long format
        mouse   = repmat((1:12)', 2, 1);              % mouse ID, repeated for each condition
        cond    = [repmat({'saline'},12,1); repmat({'drug'},12,1)];
        response = [tmp1(:); tmp2(:)];
        orderVar = [ordersALL(1:12)'; ordersALL(13:24)'];              % same order value per mouse, both conditions

        tbl = table(mouse, cond, orderVar, response, ...
            'VariableNames', {'Mouse','Condition','Order','Response'});
        tbl.Mouse = categorical(tbl.Mouse);
        tbl.Condition = categorical(tbl.Condition);
        tbl.Order = categorical(tbl.Order);   % treat order as grouping factor; drop this line if it's continuous

        % Fit LME: fixed effect of condition, random effects of mouse and order
        lme = fitlme(tbl, 'Response ~ Condition + Order + (1|Mouse) ');
anova_results = anova(lme, 'DFMethod', 'satterthwaite');
        disp(anova_results)

end
