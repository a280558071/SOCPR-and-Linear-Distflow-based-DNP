%% Distribution Network Planning (DNP) based on Linear Distflow Optimal Power Flow (OPF) Model in 
% M. E. Baran and F. F. Wu, “Optimal capacitor placement on radial distribution systems,” IEEE Trans. Power Delivery, vol. 4, no. 1, pp. 725–734,1989.
% M. E. Baran and F. F. Wu, “Optimal sizing of capacitors placed on a radial distribution system,” IEEE Trans. Power Delivery, vol. 4, no. 1, pp. 735–743, 1989
% A compact form is described in Section III. of 
% Exact convex relaxation of optimal power flow in radial networksL Gan, N Li, U Topcu, SH LowIEEE Transactions on Automatic Control 60 (1), 72-87, 2014

clear all
close all
clc
%% ***********Parameters **********
N=16; % number of load nodes
L=33; % number of dis. lines
Sbase=1e6;  % unit:VA
Ubase=10e3;  % unit:V
Ibase=Sbase/Ubase/1.732;  % unit: A
Zbase=Ubase/Ibase/1.732;  % unit: Ω
LineInf=xlsread('case data\16bus33lines.xlsx','F3:L35');
NodeInf=xlsread('case data\16bus33lines.xlsx','A3:D18');
s=LineInf(:,2);
t=LineInf(:,3);
N_subs=[13,14,15,16];  % Subs nodes
N_loads=1:12; % Load nodes
v_min=0.95^2;
v_max=1.05^2;
S_max=12 ; % max power in any distribution line
M=1e8;
% line investment cost denoted by 
Cost=LineInf(:,7).*LineInf(:,4);
% make load (MVA) become P_load Q_load assuming a power factor n.
S=NodeInf(N_loads,4);
n=0.8;  % power factor
P_load=S*n;
Q_load=S*sqrt(1-n^2);
% make impedance z become r and x assuming a rx rate n_rx
z=LineInf(:,4).*LineInf(:,6)/Zbase; % line impedance unit:p.u. 
n_rx=1;
r=z/sqrt(1+n_rx^2);
x=r/n_rx;
%% ***********Graph of the Distribution System***********
G=graph(s,t);  % G=(V,E), G:graph V:vertex E:edge
idxOut = findedge(G,s,t); % !!!!!!*************** index of edge is not the same with that of in mpc.branch !!!!!!*******
In=myincidence(s,t);  % 节支关联矩阵，node-branch incidence matrix
% Plot the Graph
figure;
p=plot(G,'Layout','force');
p.XData=NodeInf(:,2);
p.YData=NodeInf(:,3);
hold on;
labelnode(p,N_subs,{'13','14','15','16'});
highlight(p,N_loads,'NodeColor','y','Markersize',20,'NodeFontSize',20);
highlight(p,N_subs,'Marker','s','NodeColor','c','Markersize',30,'NodeFontSize',40);
highlight(p,s,t,'EdgeColor','k','LineStyle','-.','LineWidth',2,'EdgeFontSize',8);
text(p.XData, p.YData, p.NodeLabel,'HorizontalAlignment', 'center','FontSize', 15); % put nodes' label in right position.
p.NodeLabel={};
hold off;
%% ***********Variable statement**********
v_i=sdpvar(N,1, 'full'); % v=|U|^2=U·U*, "*" means conjugate
l_ij=sdpvar(L,1, 'full'); % l=|I|^2=I·I*
P_ij=sdpvar(L,1, 'full'); 
Q_ij=sdpvar(L,1, 'full');
P_shed=sdpvar(length(N_loads),1, 'full'); % shedded load
Q_shed=sdpvar(length(N_loads),1, 'full'); % shedded load
g_subs_P=sdpvar(length(N_subs),1,'full'); % generation from substations
g_subs_Q=sdpvar(length(N_subs),1,'full'); % generation from substations
y_ij=binvar(L,1,'full'); % decision variable denoting whether this line should be invested
%% ***********Constraints*************
Cons=[];
%% 1. Power balance
% S_ij=s_i+\sum_(h:h→i) S_hi for any (i,j)∈E, 
% denoted by node-branch incidence matrix I
Cons_S=[];
Cons_S=[In*P_ij==[-(P_load-P_shed);g_subs_P],In*Q_ij==[-(Q_load-Q_shed);g_subs_Q], P_shed>=0,Q_shed==P_shed*sqrt(1-n^2)/n];
Cons=[Cons,Cons_S];
%% 2. Voltage Calculation
% v_i-v_j=2Re(z_ij·S_ij*)=2(r·P_ij+x·Q_ij) 
% → |v_i-v_j-2(r·P_ij+x·Q_ij)|≤(1-y_ij)*M
Cons_V=[v_i(N_subs)==v_max,v_i];
Cons_V=[Cons_V,abs(In'*v_i-2*r.*P_ij-2*x.*Q_ij)<=(1-y_ij)*M];
Cons=[Cons,Cons_V];
%% 3. Voltage limits
% v_min<=v<=v_max
Cons=[Cons,v_min<=v_i<=v_max];
%% 4. P_ij^2+Q_ij^2<=S_ij_max
% Cons_PQ=[abs(P_ij)<=y_ij*S_max,abs(Q_ij)<=y_ij*S_max];
Cons_PQ=[P_ij.^2+Q_ij.^2<=y_ij*S_max.^2];
Cons=[Cons, Cons_PQ];
%% **********Objectives*************
Obj_inv=sum(Cost.*y_ij);
Obj_ope=M*sum(P_shed);
Obj=Obj_inv+Obj_ope;
%% ********* Solve the probelm
ops=sdpsettings('solver','gurobi', 'gurobi.Heuristics',0,'gurobi.Cuts',0); %,'usex0',1,'gurobi.MIPGap',5e-2,
sol=optimize(Cons,Obj,ops)
%% Save the solution with "s_" as start
s_y_ij=value(y_ij);
s_v_i=value(v_i);
s_P_ij=value(P_ij);
s_Q_ij=value(Q_ij);
s_P_shed=value(P_shed);
s_Q_shed=value(Q_shed);
s_g_subs_P=value(g_subs_P);
s_g_subs_Q=value(g_subs_Q);
s_Obj=value(Obj);
s_Obj=value(Obj_inv);
s_Obj=value(Obj_ope);
%% Print the results in the command line
display('————————规划方案如下——————————————');
display(['   建设线路方案： ',num2str(round(s_y_ij',2))]);
display([' (1/0 表示 建设/不建设 该线路, 线路编号参见输入文件)']);
display('————————规划建设成本如下——————————————');
display(['   建设配电线路成本： ',num2str(value(Obj_inv)),' 美元']);
display(['   失负荷成本:  ',num2str(round(value(Obj_ope),2)),'  美元']);
display(['>> 规划建设总成本:  ',num2str(value(Obj)),'  美元']);
%% Plot the results
Gi=graph(s,t);
os=0.1; % offset value to print text in the figure
figure;
pi=plot(Gi,'Layout','force');
pi.LineWidth(idxOut) = 10*abs(s_P_ij)/max(abs(s_P_ij))+0.01;
pi.XData=NodeInf(:,2);
pi.YData=NodeInf(:,3);
labelnode(pi,N_subs,{'13','14','15','16'});
highlight(pi,N_loads,'NodeColor','y','Markersize',20,'NodeFontSize',20);
highlight(pi,N_subs,'Marker','s','NodeColor','c','Markersize',30,'NodeFontSize',40);
text(pi.XData, pi.YData, pi.NodeLabel,'HorizontalAlignment', 'center','FontSize', 15); % put nodes' label in right position.
text(pi.XData+os, pi.YData+os, num2str([round(s_P_shed,2);zeros(length(N_subs),1)]),'HorizontalAlignment', 'center','FontSize', 20, 'Color','r'); % printf the shedded load (if any).
text(pi.XData-os, pi.YData-os, num2str(sqrt(s_v_i)),'HorizontalAlignment', 'center','FontSize', 15, 'Color','g'); % printf the nodes' voltage (sqrt(v_j)).
for l=1:L
    Coor1_PQ=(pi.XData(s(l))+pi.XData(t(l)))/2;
    Coor2_PQ=(pi.YData(s(l))+pi.YData(t(l)))/2;
    text(Coor1_PQ, Coor2_PQ, num2str(round(s_P_ij(l)+s_Q_ij(l)*i,2)),'HorizontalAlignment', 'center','FontSize', 15, 'Color','b'); % printf the complex power in distribution lines(if any).
end 
pi.NodeLabel={};
%% Save the results
warning off
save('LinearDistFlow_planning_results');