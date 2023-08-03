clc
clear
close all
Fun_name='F1'; 
SearchAgents=30; 
Max_iterations=100;
number=30;  
TOTAL=zeros(number,1);  
[lowerbound,upperbound,dimension,fitness]=Hight_Get_Functions_details(Fun_name); 

%% AFT	2021	Ali Baba and the forty thieves
for g=1:number
[Best_score1,Best_pos1,AFT_curve1]=AFT(SearchAgents,Max_iterations,lowerbound,upperbound,dimension,fitness);  
  TOTAL(g,:)=Best_score1;
end
    fprintf('%d\n%d\n%g\n',min(TOTAL),mean(TOTAL),std(TOTAL));
%% DOA	2021	Dingo Optimization Algorithm
for g=1:number
[Best_score2,Best_pos2,DOA_curve2]=DOA(SearchAgents,Max_iterations,lowerbound,upperbound,dimension,fitness);  % 战争策略
  TOTAL(g,:)=Best_score2;
end
    fprintf('%d\n%d\n%g\n',min(TOTAL),mean(TOTAL),std(TOTAL));
%% SCA	2016	Sine Cosine Algorithm
for g=1:number
[Best_score3,Best_pos3,SCA_curve3]=SCA(SearchAgents,Max_iterations,lowerbound,upperbound,dimension,fitness);  
  TOTAL(g,:)=Best_score3;
end
    fprintf('%d\n%d\n%g\n',min(TOTAL),mean(TOTAL),std(TOTAL));
%% SSA	2017	Salp Swarm Algorithm, SSA
for g=1:number
[Best_score4,Best_pos4,SSA_curve4]=SSA(SearchAgents,Max_iterations,lowerbound,upperbound,dimension,fitness);  
  TOTAL(g,:)=Best_score4;
end
    fprintf('%d\n%d\n%g\n',min(TOTAL),mean(TOTAL),std(TOTAL));
%% WSO	2022	White Shark Optimizer
for g=1:number
[Best_score5,Best_pos5,WSO_curve5]=WSO(SearchAgents,Max_iterations,lowerbound,upperbound,dimension,fitness);  
  TOTAL(g,:)=Best_score5;
end
    fprintf('%d\n%d\n%g\n',min(TOTAL),mean(TOTAL),std(TOTAL));
%% WSOA	2022	War Strategy Optimization Algorithm
for g=1:number
[Best_score6,Best_pos6,WSOA_curve6]=WSOA(SearchAgents,Max_iterations,lowerbound,upperbound,dimension,fitness);  % 战争策略
  TOTAL(g,:)=Best_score6;
end
    fprintf('%d\n%d\n%g\n',min(TOTAL),mean(TOTAL),std(TOTAL));
%% YDSE	2023	Young's double-slit experiment optimizer
for g=1:number
[Best_score7,Best_pos7,YDSE_curve7]=YDSE(SearchAgents,Max_iterations,lowerbound,upperbound,dimension,fitness);  % 战争策略
  TOTAL(g,:)=Best_score7;
end
    fprintf('%d\n%d\n%g\n',min(TOTAL),mean(TOTAL),std(TOTAL));
%% TSO	2021
for g=1:number
[Best_score8,Best_pos8,TSO_curve8]=TSO(SearchAgents,Max_iterations,lowerbound,upperbound,dimension,fitness);  % Calculating the solution of the given problem using SABO
  TOTAL(g,:)=Best_score8;
end
    fprintf('%d\n%d\n%g\n',min(TOTAL),mean(TOTAL),std(TOTAL));
%% TSO1	2022
for g=1:number
[Best_score9,Best_pos9,MTSO_curve9]=MTSO(SearchAgents,Max_iterations,lowerbound,upperbound,dimension,fitness);  
  TOTAL(g,:)=Best_score9;
end
    fprintf('%d\n%d\n%g\n',min(TOTAL),mean(TOTAL),std(TOTAL));
%% CTSO	2022
for g=1:number
[Best_score10,Best_pos10,CTSO_curve10]=CTSO(SearchAgents,Max_iterations,lowerbound,upperbound,dimension,fitness);  % 战争策略
  TOTAL(g,:)=Best_score10;
end
    fprintf('%d\n%d\n%g\n',min(TOTAL),mean(TOTAL),std(TOTAL));
%% CETSO
for g=1:number
[Best_score11,Best_pos11,CETSO_curve11]=CETSO(SearchAgents,Max_iterations,lowerbound,upperbound,dimension,fitness);  % 战争策略
  TOTAL(g,:)=Best_score11;
end
    fprintf('%d\t%d\t%g\n',min(TOTAL),mean(TOTAL),std(TOTAL));

%%
color=[1	0.880 	0.473 
1	0.688 	0.543 
1	0.535 	0.669 
0.839	0.465 	0.794 
0.618	0.469 	0.861 
0.273	0.696 	0.924 
0	0.702 	0.912 
0	0.637 	0.837 
0	0.657 	0.708 
0	0.661 	0.578 
0.820 	0.227 	0.157   ];

figure('units','normalized','position', [0.501 0.25 0.38 0.5])
        semilogy(AFT_curve1,'color',color(1,:),'linewidth',5,'linestyle','-.')
            hold on
        semilogy(DOA_curve2,'color',color(2,:),'linewidth',5,'linestyle','-.')
        semilogy(SCA_curve3,'color',color(3,:),'linewidth',5,'linestyle','-.')
        semilogy(SSA_curve4,'color',color(4,:),'linewidth',5,'linestyle','-.')
        semilogy(WSO_curve5,'color',color(5,:),'linewidth',5,'linestyle','-.')
        semilogy(WSOA_curve6,'color',color(6,:),'linewidth',5,'linestyle','-.')
        semilogy(YDSE_curve7,'color',color(7,:),'linewidth',5,'linestyle','-.')
        semilogy(TSO_curve8,'color',color(8,:),'linewidth',6,'linestyle','-')
        semilogy(MTSO_curve9,'color',color(9,:),'linewidth',5,'linestyle','--')
        semilogy(CTSO_curve10,'color',color(10,:),'linewidth',5,'linestyle','--')
        semilogy(CETSO_curve11,'color',color(11,:),'linewidth',6,'linestyle','-')

% legend('AFT','DOA','SCA','SSA','WSO','WSOA','YDSE','TSO','MTSO','CTSO','CETSO')
% grid on
axis tight
% box on
%     xlabel('Iteration');
%     ylabel('Fitness value');
