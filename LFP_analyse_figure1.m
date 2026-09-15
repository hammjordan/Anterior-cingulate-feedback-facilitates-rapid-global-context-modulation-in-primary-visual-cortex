
clcl
limtrials1=6; %limit control trials to this. if 0, it wont do it;  
for fileload=1;
         files={'pre_mouse34';
           'pre_mouse35'; 
            'pre_mouse37';
            'pre_mouse38';
            'pre_mouse41'; 
            'pre_mouse42';
            'pre_mouse36';'pre_mouse43';'pre_mouse44';
             };
end

for find_times=1;
    times=[-100 0 75 600 700]; %prestim, zero, midstim, endstim
    load(strcat('GLctrl1_',files{1}),'t_sG1');
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
         load(strcat('GLtest1_',files{fila}));
      

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

        itc=tfs1a(:,:,trials1==1)./abs(tfs1a(:,:,trials1==1));nn=trnums_oddball(fila,1);
        itcsV1_oddball(:,:,fila,1)=abs(mean(itc,3))-sqrt(-(1/nn)*log(.5));
         itc=tfs1a(:,:,trials1==2)./abs(tfs1a(:,:,trials1==2));nn=trnums_oddball(fila,2);
        itcsV1_oddball(:,:,fila,2)=abs(mean(itc,3))-sqrt(-(1/nn)*log(.5));
         itc=tfs1a(:,:,trials1==3)./abs(tfs1a(:,:,trials1==3));nn=trnums_oddball(fila,3);
        itcsV1_oddball(:,:,fila,3)=abs(mean(itc,3))-sqrt(-(1/nn)*log(.5));
         itc=tfs1a(:,:,trials1==4)./abs(tfs1a(:,:,trials1==4));nn=trnums_oddball(fila,4);
        itcsV1_oddball(:,:,fila,4)=abs(mean(itc,3))-sqrt(-(1/nn)*log(.5));
        waitbar((fila/size(files,1)*.5),fh,'doing oddballs');
    end
end
  
for controlanalyse=1;
  
    for fila=1:size(files,1);
        %oris=[2 3 5 4 1 8 6 7]
         load(strcat('GLctrl1_',files{fila})); 
       
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
                 tt=(size(trials1,2)+1)-t; 
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
 
close all force
for visualresponsees=1;
    figure;
    sc='linear';
    groups=ones(1,fila);

    figure;
    sc='linear';
    % visual induced oscillations (power
    for grp=1;
        subs=groups==grp;
        for power_to_stimV1_control=1;
            tmp1=erspsV1_control(:,:,subs,1);
            tmp2=erspsV1_control(:,:,subs,2);
            tmp3=erspsV1_control(:,:,subs,3);
            
            subplot(2,4,1); contourf(t_sG1,f_sG1,mean(tmp1,3),50,'linecolor','none'); set(gca,'yscale', sc);title('control 1'); colormap jet;  axis([-150 735 4 50]);cc=[-2 5];ylabel('control context');   caxis(cc);
            subplot(2,4,2); contourf(t_sG1,f_sG1,mean(tmp1,3),50,'linecolor','none'); set(gca,'yscale', sc);title('control 1');  colormap jet;  axis([-150 735 4 50]); caxis(cc);
            subplot(2,4,3);contourf(t_sG1,f_sG1,mean(tmp2,3),50,'linecolor','none'); set(gca,'yscale', sc); title('control 2'); colormap jet;   axis([-150 735 4 50]);caxis(cc);
            subplot(2,4,4);contourf(t_sG1,f_sG1,mean(tmp3,3),50,'linecolor','none'); set(gca,'yscale', sc);title('control 3'); colormap jet;  axis([-150 735 4 50]);caxis(cc);
            make_eps_saveable        
        end

        for powe_to_stimV1_OB=1;
            tmp0=erspsV1_oddball(:,:,subs,1);
            tmp1=erspsV1_oddball(:,:,subs,2);
            tmp2=erspsV1_oddball(:,:,subs,3);
            tmp3=erspsV1_oddball(:,:,subs,4);

            subplot(2,4,5); contourf(t_sG1,f_sG1,mean(tmp0,3),50,'linecolor','none'); set(gca,'yscale', sc); title('redundant'); colormap jet;  axis([-100 735 4 50]);ylabel('oddball');  caxis(cc);
            subplot(2,4,6);contourf(t_sG1,f_sG1,mean(tmp2,3),50,'linecolor','none'); set(gca,'yscale', sc); title('redundant OB'); colormap jet;   axis([-100 735 4 50]);caxis(cc);
            subplot(2,4,7);contourf(t_sG1,f_sG1,mean(tmp1,3),50,'linecolor','none'); set(gca,'yscale', sc); title('predictable OB'); colormap jet;  axis([-100 735 4 50]);caxis(cc);
            subplot(2,4,8);contourf(t_sG1,f_sG1,mean(tmp3,3),50,'linecolor','none'); set(gca,'yscale', sc); title('global OB'); colormap jet; axis([-100 735 4 50]);caxis(cc);
           make_eps_saveable        
        end
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
            shadedErrorBar(t_sG1,mn2,st2,'lineProps','b');  title('redundant');axis([-100 650 -2 5]);

            mn1=mean(tmp1c); st1=std(tmp1c,0,1)./sqrt(size(tmp1c,1));
            mn2=mean(tmp2); st2=std(tmp2,0,1)./sqrt(size(tmp2,1));
            subplot(1,4,2);shadedErrorBar(t_sG1,mn1,st1,'lineProps','k'); hold on
            shadedErrorBar(t_sG1,mn2,st2,'lineProps','r');  title('redundant OB');axis([-100 650 -2 5]);

            mn1=mean(tmp2c); st1=std(tmp2c,0,1)./sqrt(size(tmp2c,1));
            mn2=mean(tmp1); st2=std(tmp1,0,1)./sqrt(size(tmp1,1));
            subplot(1,4,3);shadedErrorBar(t_sG1,mn1,st1,'lineProps','k'); hold on
            shadedErrorBar(t_sG1,mn2,st2,'lineProps','r'); title('predictable OB');axis([-100 650 -2 5]);

            mn1=mean(tmp3c); st1=std(tmp3c,0,1)./sqrt(size(tmp3c,1));
            mn2=mean(tmp3); st2=std(tmp3,0,1)./sqrt(size(tmp3,1));
            subplot(1,4,4);shadedErrorBar(t_sG1,mn1,st1,'lineProps','k'); hold on
            shadedErrorBar(t_sG1,mn2,st2,'lineProps','r'); title('global OB'); axis([-100 650 -2 5]);


           


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
            ylim([-2 5.0]);

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
            ylim([-2 5.0]);

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
            ylim([-2 5.0]);

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
            ylim([-2 5.0]);
            make_eps_saveable

            dat1=mean(tmp1c(:,tss),2);
            dat2=mean(tmp0(:,tss),2);

            dat3=mean(tmp1c(:,tss),2);
            dat4=mean(tmp2(:,tss),2);


            dat5=mean(tmp2c(:,tss),2);
            dat6=mean(tmp1(:,tss),2);


            dat7=mean(tmp3c(:,tss),2);
            dat8=mean(tmp3(:,tss),2);

            anova_rm({[dat1 dat2] [dat3 dat4] [dat5 dat6] [dat7 dat8]})
            [h,p,ci,stats]=ttest(dat1,dat2);disp(strcat('adaptation. t=',num2str(stats.tstat),'. p=',num2str(p)))
            [h,p,ci,stats]=ttest(dat3,dat4); disp(strcat('local rdnt global DD. t=',num2str(stats.tstat),'. p=',num2str(p)))
            [h,p,ci,stats]=ttest(dat5,dat6); disp(strcat('local DD global rdnt. t=',num2str(stats.tstat),'. p=',num2str(p)))
            [h,p,ci,stats]=ttest(dat7,dat8); disp(strcat('local DD global DD. t=',num2str(stats.tstat),'. p=',num2str(p)))

        end


    end

end

    if oldapproach==1;


        for controlanalyse=1;
            fh = waitbar(0,'Percent completEE...');
            for fila=1:size(files,1);
                %oris=[2 3 5 4 1 8 6 7]
                load(strcat('GLctrl1_',files{fila}));
                if doaca==1;
                    ersp1=ersp2; tfs1=tfs2;
                end
                % rej=rej-rej;
                orders=[]; for zzz=1:200; orders=horzcat(orders,[1:5]); end; orders=orders(1,1:size(trials1,2));
                ersp1=ersp1(:,:,rej==0);
                %ersp2=ersp2(:,:,rej==0);
                tfs1a=tfs1(:,:,rej==0);
                %tfs2a=tfs2(:,:,rej==0);
                erp1a=erp1(rej==0,:);
                %erp2a=erp2(rej==0,:);

                trials1(1)=0; for t=2:size(trials1,2); if and(trials1(t)==1, sum(trials1(t-1)==[1])==0); trials1(t)=1;
                elseif and(trials1(t)==2, sum(trials1(t-1)==[2])==0); trials1(t)=2;
                elseif and(trials1(t)==5, sum(trials1(t-1)==[5])==0); trials1(t)=5; else trials1(t)=0;
                end; end
            for t=2:size(trials1,2); if orders(1,t)<4; trials1(1,t)=0; end; end;

            if limtrials>0;
                counts=zeros(1,8);
                for t=1:size(trials1,2);
                    tt=(size(trials1,2)+1)-t;
                    if trials1(1,tt)>0
                        counts(1,trials1(1,tt))=counts(1,trials1(1,tt))+1;
                        if counts(1,trials1(1,tt))>limtrials
                            trials1(1,tt)=0;
                        end
                    end
                end
            end


            trials1=trials1(rej==0);

            trnums_control(fila,1)=sum(trials1==1);
            trnums_control(fila,2)=sum(trials1==2);
            trnums_control(fila,3)=sum(trials1==5);
            erspsV1_control(:,:,fila,1)=mean(ersp1(:,:,trials1==1),3);
            erspsV1_control(:,:,fila,2)=mean(ersp1(:,:,trials1==2),3);
            erspsV1_control(:,:,fila,3)=mean(ersp1(:,:,trials1==5),3);
            %erspsPFC_control(:,:,fila,1)=mean(ersp2(:,:,trials1==1),3);
            %erspsPFC_control(:,:,fila,2)=mean(ersp2(:,:,trials1==2),3);
            % erspsPFC_control(:,:,fila,3)=mean(ersp2(:,:,trials1==5),3);



            itc=tfs1a(:,:,trials1==1)./abs(tfs1a(:,:,trials1==1));nn=trnums_control(fila,1);
            itcsV1_control(:,:,fila,1)=abs(mean(itc,3))-sqrt(-(1/nn)*log(.5));
            itc=tfs1a(:,:,trials1==2)./abs(tfs1a(:,:,trials1==2));nn=trnums_control(fila,2);
            itcsV1_control(:,:,fila,2)=abs(mean(itc,3))-sqrt(-(1/nn)*log(.5));
            itc=tfs1a(:,:,trials1==5)./abs(tfs1a(:,:,trials1==5));nn=trnums_control(fila,3);
            itcsV1_control(:,:,fila,3)=abs(mean(itc,3))-sqrt(-(1/nn)*log(.5));

            for coherencesanalyses=1;
                % clear coherbase
                % for calculatecoherencebetweens=1
                %     for f=1:size(f_sG1,2)
                %         lags=0; cnt=0;
                %         for t=t1:t2
                %             for tr=1:size(tfsALL1,3)
                %                 l1=tfsALL1(f,t,tr)./abs(tfsALL1(f,t,tr));
                %                 l2=tfsALL2(f,t,tr)./abs(tfsALL2(f,t,tr));
                %                 lags=lags+(1*exp(1i*((angle(l1)-angle(l2)))));
                %                 cnt=cnt+1;
                %             end
                %         end
                %         coherbase(f)=abs(lags/cnt);%-(sqrt(-(1/cnt)*log(.5)));
                %     end
                % end
                % cohers_control(:,1,fila)=coherbase;
                %
                % clear coherstim1
                % for calculatecoherencebetweens=1
                %     for f=1:size(f_sG1,2)
                %         lags=0; cnt=0;
                %         for t=t2:t3
                %             for tr=1:size(tfsALL1,3)
                %                 l1=tfsALL1(f,t,tr)./abs(tfsALL1(f,t,tr));
                %                 l2=tfsALL2(f,t,tr)./abs(tfsALL2(f,t,tr));
                %                 lags=lags+(1*exp(1i*((angle(l1)-angle(l2)))));
                %                 cnt=cnt+1;
                %             end
                %         end
                %         coherstim1(f)=abs(lags/cnt);%-(sqrt(-(1/cnt)*log(.5)));
                %     end
                % end
                % cohers_control(:,2,fila)=coherstim1;
                %
                % clear coherstim2
                % for calculatecoherencebetweens=1
                %     for f=1:size(f_sG1,2)
                %         lags=0; cnt=0;
                %         for t=t3:t4
                %             for tr=1:size(tfsALL1,3)
                %                 l1=tfsALL1(f,t,tr)./abs(tfsALL1(f,t,tr));
                %                 l2=tfsALL2(f,t,tr)./abs(tfsALL2(f,t,tr));
                %                 lags=lags+(1*exp(1i*((angle(l1)-angle(l2)))));
                %                 cnt=cnt+1;
                %             end
                %         end
                %         coherstim2(f)=abs(lags/cnt);%-(sqrt(-(1/cnt)*log(.5)));
                %     end
                % end
                % cohers_control(:,3,fila)=coherstim2;
                %
                % clear cohershuff
                % for calculatecoherencebetweens=1
                %     for f=1:size(f_sG1,2)
                %         for z=1:100;
                %             lags=0; cnt=0;
                %             rr=randsample(size(tfsALL2,3),1);
                %             for t=t2:t3
                %                 for tr=1:size(tfsALL1,3)
                %                     if rr==tr; rr=1; if tr==1; rr=5; end; end
                %                     l1=tfsALL1(f,t,tr)./abs(tfsALL1(f,t,tr));
                %                     l2=tfsALL2(f,t,rr)./abs(tfsALL2(f,t,rr));
                %                     lags=lags+(1*exp(1i*((angle(l1)-angle(l2)))));
                %                     cnt=cnt+1;
                %                 end
                %             end
                %             cohershuff(f,z)=abs(lags/cnt);%-(sqrt(-(1/cnt)*log(.5)));
                %         end
                %     end
                % end
                % cohers_control(:,4,fila)=mean(cohershuff,2);
                %
                % clear cohershuff
                % for calculatecoherencebetweens=1
                %     for f=1:size(f_sG1,2)
                %         for z=1:100;
                %             lags=0; cnt=0;
                %             rr=randsample(size(tfsALL2,3),1);
                %             for t=t3:t4
                %                 for tr=1:size(tfsALL1,3)
                %                     if rr==tr; rr=1; if tr==1; rr=5; end; end
                %                     l1=tfsALL1(f,t,tr)./abs(tfsALL1(f,t,tr));
                %                     l2=tfsALL2(f,t,rr)./abs(tfsALL2(f,t,rr));
                %                     lags=lags+(1*exp(1i*((angle(l1)-angle(l2)))));
                %                     cnt=cnt+1;
                %                 end
                %             end
                %             cohershuff(f,z)=abs(lags/cnt);%-(sqrt(-(1/cnt)*log(.5)));
                %         end
                %     end
                % end
                % cohers_control(:,5,fila)=mean(cohershuff,2);
            end
            waitbar((fila/size(files,1)*.5),fh,'doing controls');
            end
        end
        for GLoddballanalyse=1;

            for fila=1:size(files,1);
                %oris=[2 3 5 4 1 8 6 7]
                load(strcat('GLtest1_',files{fila}));
                if doaca==1;
                    ersp1=ersp2; tfs1=tfs2;
                end
                %rej=rej-rej;
                orders=[]; for zzz=1:200; orders=horzcat(orders,[1:5]); end; orders=orders(1,1:size(trials1,2));
                ersp1=ersp1(:,:,rej==0);
                %ersp2=ersp2(:,:,rej==0);
                tfs1a=tfs1(:,:,rej==0);
                %tfs2a=tfs2(:,:,rej==0);
                erp1a=erp1(rej==0,:);
                %erp2a=erp2(rej==0,:);

                trials1=horzcat(trials1,zeros(1,6));
                trials1(1)=0; for t=1:size(trials1,2)-6; if and(trials1(t)==4, trials1(t)==4); trials1(t)=4;
                elseif and(trials1(t)==3, trials1(t)==3); trials1(t)=3;
                elseif and(trials1(t)==1, trials1(t+6)==4); trials1(t)=1;
                elseif and(trials1(t)==2, trials1(t+5)~=2); trials1(t)=2;
                else trials1(t)=0;
                end; end
            trials1=trials1(1,1:end-6);
            trials1=trials1(rej==0);
            %part 2 limtrials
            for limtrials=limtrials; %dont need to do this for the 6 and8, but only 4
                if limtrials>0;
                    %then, limit to last 6 occurrances
                    counts=zeros(1,8);
                    for tr=1:size(trials1,2);
                        tr1=(size(trials1,2)+1)-tr;
                        aa=trials1(1,tr1);
                        if aa>0
                            counts(1,aa)=counts(1,aa)+1;
                            if counts(1,aa)>6;
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
            %erspsPFC_control(:,:,fila,1)=mean(ersp2(:,:,trials1==1),3);
            %erspsPFC_control(:,:,fila,2)=mean(ersp2(:,:,trials1==2),3);
            % erspsPFC_control(:,:,fila,3)=mean(ersp2(:,:,trials1==5),3);



            itc=tfs1a(:,:,trials1==1)./abs(tfs1a(:,:,trials1==1));nn=trnums_oddball(fila,1);
            itcsV1_oddball(:,:,fila,1)=abs(mean(itc,3))-sqrt(-(1/nn)*log(.5));
            itc=tfs1a(:,:,trials1==2)./abs(tfs1a(:,:,trials1==2));nn=trnums_oddball(fila,2);
            itcsV1_oddball(:,:,fila,2)=abs(mean(itc,3))-sqrt(-(1/nn)*log(.5));
            itc=tfs1a(:,:,trials1==3)./abs(tfs1a(:,:,trials1==3));nn=trnums_oddball(fila,3);
            itcsV1_oddball(:,:,fila,3)=abs(mean(itc,3))-sqrt(-(1/nn)*log(.5));
            itc=tfs1a(:,:,trials1==4)./abs(tfs1a(:,:,trials1==4));nn=trnums_oddball(fila,4);
            itcsV1_oddball(:,:,fila,4)=abs(mean(itc,3))-sqrt(-(1/nn)*log(.5));

            waitbar((fila/size(files,1)*.5)+.5,fh,'doing oddballs');
            end
        end
    end
end