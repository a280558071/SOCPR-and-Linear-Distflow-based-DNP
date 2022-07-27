%% This function implement the technique in [1], on exact convex relaxation for distribution system optimal power flow (DOPF)
% It includes the exmaple of SCE 47 bus in [1], and allow comparison
% between with/without s∈S_volt in Sec. IV..
% [1] Gan L.,Li N.,Topcu U.and H. Low S. 2015. Exact Convex Relaxation of Optimal Power Flow in Radial Networks.
% IEEE Transactions on Automatic Control 60(1). 72-87. 10.1109/TAC.2014.2332712.

clear
%% With/without s∈S_volt
WithSvolt=1; % 1 / 0 means with / without s∈S_volt
%% 系统参数
Ubase=12.35e3; % unit:V
Sbase=1e7;   % unit:VA
Ibase=Sbase/sqrt(3)/Ubase; % unit:A
Zbase=Ubase/Ibase/sqrt(3);  % unit: Ω
netpara=xlsread('excel2021','网络参数');
loadpoint=xlsread('excel2021','节点负荷');
L=size(netpara,1);%网络中的支路
r=netpara(:,4)/Zbase;%电阻,unit: p.u.
x=netpara(:,5)/Zbase;%电抗,unit: p.u.
So=loadpoint(:,3)*1e6/Sbase;%节点负荷功率,unit: MVA
judge=loadpoint(:,4); %节点类型
I_max=560.98/Ibase; % Ubase=12.35e3, 假设每条线路最多传输 12MW 功率，则对应电流值I_max为 12e6/Ubase/sqrt(3)= 560.98 A,
v_max=1.1^2;
v_min=0.9^2;

%% 绘制网络情况，便于后续观察
I=netpara(:,2);
J=netpara(:,3);
G=graph(I,J);
h1=figure;
p=plot(G);
highlight(p,find(judge==1),'Marker','s','NodeColor','c','Markersize',10);  % Capacitator
highlight(p,find(judge==0),'Marker','v','NodeColor','y','Markersize',10);  % PV panels
In=myincidence(I,J); % 节点支路关联矩阵
Inn=In;
Inn(Inn>0)=0;    % Inn is the negative part of I，表示从节点出发的支路的编号
%% 节点数
N=47;
%m=4;
%% 决策变量
P_ij=sdpvar(L,1); % Pij=sdpvar(numnodes,numnodes,'full');% 节点1向节点j传送的有功功率
Q_ij=sdpvar(L,1); % Qij=sdpvar(numnodes,numnodes,'full');%节点i向节点j传送的无功功率
u=sdpvar(N,1); % 各节点电压
v=sdpvar(N,1); % u^2
l_ij=sdpvar(L,1) ;% 电流幅值的平方
Pi=sdpvar(N,1);% 节点注入有功功率 ppoint(i)
Qi=sdpvar(N,1);%节点注入wu功功率

%% 阻抗矩阵矩阵节点参数初始化
w=sqrt(-1);
z=r+x*w;
z_c=conj(z); % 对z取共轭

Cons=[];
%% 判断筛选负荷节点类型
Cons_Load=[];
Eta=2.8;  % Scale up the PV penetration rate, 5 means 5 times than Table I in [1]
Pi_max=zeros(N,1);
Qi_max=zeros(N,1);
for i=1:N
    if judge(i)==1     %电容补偿
        Cons_Load=[Cons_Load,Pi(i)==0];
        Qi_max(i)=So(i)*Eta;
        Cons_Load=[Cons_Load,0<=Qi(i)<=Qi_max(i)];   %st=[pl(i)==0,0<=ql(i)<=So(i)];
    elseif judge(i)==2 %普通节点负荷
        Pi_max(i)=-So(i)*0.9;
        Qi_max(i)=-So(i)*sqrt(1-0.9^2);
        Cons_Load=[Cons_Load,Pi(i)==Pi_max];  % cos(0.9)???
        Cons_Load=[Cons_Load,Qi(i)==Qi_max];               %st=[pl(i)==-So(i)*cos(0.9);ql(i)==-So(i)*sin(0.9)];
    elseif judge(i)==0  %PV panel
        Pi_max(i)=So(i)*Eta;
        Cons_Load=[Cons_Load,0<=Pi(i)<=Pi_max(i)];
        Cons_Load=[Cons_Load,0==Qi(i)];
    elseif judge(i)==3  %slack node (首端节点，即变电站)
        Cons_Load=[Cons_Load,-So(i)<=Pi(i)<=So(i)];
        Cons_Load=[Cons_Load,-So(i)*0.5<=Qi(i)<=So(i)*0.5];
    end
end
Cons=[Cons,Cons_Load];
display('Constraints on load type of nodes completed!')
% Cons_Load

%% Check whether C1 holds for the network
De = degree(G); % degreee of each node
LeafNodes=find(De==1);
LeafNodes=LeafNodes(2:end); % delete node 1
Pij_hat=abs(inv(In(2:end,:))*Pi_max(2:end));
Qij_hat=abs(inv(In(2:end,:))*Qi_max(2:end));
IJ=[I J];
A_l=cell(length(LeafNodes),1);
Path_L=cell(length(LeafNodes),1);
Al_uij=zeros(length(LeafNodes),2);
for i=1:length(LeafNodes)
    lt=LeafNodes(i);
    Path = shortestpath(G,lt,1);
    A_l{i}=zeros(2,2,length(Path)-1);
    for j=1:length(Path)-1 % for all node in Path from lt→1, calculate A_l
%         nl=[]; % the set of 1,2,...,nl in C1
        for k=1:L
            if (IJ(k,1)==Path(j) && IJ(k,2)==Path(j+1))||(IJ(k,2)==Path(j) && IJ(k,1)==Path(j+1))
                Path_L{i}=[Path_L{i},Path(j+1)];
                A_l{i}(:,:,j)=diag([1 1])-2/v_min*mtimes([r(k);x(k)],[Pij_hat(k),Qij_hat(k)]);
            end
        end
    end
    Al_uij_temp=diag([1 1]);
    for m=2:length(Path)-1
        Al_uij_temp=Al_uij_temp*A_l{i}(:,:,m);
    end
    Al_uij(i,:)=(Al_uij_temp*[r(find(J==lt));x(find(J==lt))])';
end
if find(Al_uij<0)
    display('Important message: C1 doesn''t hold !')
else
    display('Important message: C1 holds !')
end

%% 电流电压幅值约束
Cons_Limits=[I_max^2>=l_ij>=0, v(2:end)>=v_min];
Cons=[Cons,Cons_Limits];
display('Constraints on voltage and current limits completed!')
% Cons_Limits

%% 节点电压计算约束
Cons_V=[];
Cons_V=[Cons_V, v(1)==v_max];
Cons_V=[Cons_V, In'*v == 2*r.*P_ij+2*x.*Q_ij-(r.^2+x.^2).*l_ij];
Cons=[Cons,Cons_V];
display('Constraints on voltage calculation completed!')
% Cons_V

%% SOC松弛约束
Cons_SOC=[];
for k=1:L
    i=I(k);%首节点
    Cons_SOC=[Cons_SOC,(l_ij(k)+v(i)).^2 >= 4*P_ij(k).^2 + 4*Q_ij(k).^2 + (l_ij(k)-v(i)).^2];
end
Cons=[Cons,Cons_SOC];
display('Constraints on SOCP relaxation completed!')
% Cons_SOC

%%  潮流约束
Cons_PF=[In*P_ij - Inn*(r.*l_ij) == Pi, In*Q_ij - Inn*(x.*l_ij) == Qi];
Cons=[Cons,Cons_PF];
display('Constraints on power flow caculation completed!')
% Cons_PF

%% Constraints on s∈S_{volt}
if WithSvolt==1
    P_ij_s=sdpvar(L,1); % the solution of the Linear DistFlow model, to constitue Svolt
    Q_ij_s=sdpvar(L,1); % the solution of the Linear DistFlow model
    v_s=sdpvar(N,1); % the solution of the Linear DistFlow model
    
    Cons_Svolt=[];
    Cons_Svolt=[Cons_Svolt, v_max>=v_s(2:end)];
    Cons_Svolt=[Cons_Svolt, v_s(1)==v_max];
    Cons_Svolt=[Cons_Svolt, In'*v_s == 2*r.*P_ij_s+2*x.*Q_ij_s];
    Cons_Svolt=[Cons_Svolt, In(2:end,:)*P_ij_s == Pi(2:end), In(2:end,:)*Q_ij_s == Qi(2:end)];
    Cons=[Cons, Cons_Svolt];
    display('Constraints on s∈S_{volt} completed!')
end

%% 目标函数
Pr_sub=0.9;
Pr_pv=[0.8 0.7 0.6 0.7 0.8];
C=Pr_pv*Pi([find(judge==0)])+Pr_sub*Pi(1);
% C=sum(Pi);
%% 求解
ops=sdpsettings('solver', 'gurobi');
result=optimize(Cons,C,ops);
result.info

if result.problem==0
    %% 保存求解结果
    s_P_ij=value(P_ij)*Sbase/1e6;
    s_Q_ij=value(Q_ij)*Sbase/1e6;
    s_v=value(v);
    s_l_ij=value(l_ij);
    s_Pi=value(Pi)*Sbase/1e6;
    s_Qi=value(Qi)*Sbase/1e6;
    display(['光伏发电总出力：',num2str(sum(s_Pi([find(judge==0)])),3),' MW']);
    display(['变电站总出力：',num2str(sum(s_Pi(1)),3),' MW']);
    if WithSvolt==1
        s_P_ij_s=value(P_ij_s)*Sbase/1e6;
        s_Q_ij_s=value(Q_ij_s)*Sbase/1e6;
        s_v_s=value(v_s);
    end
    
    %% 在图上标注一些信息
    for k=1:L
        if norm([s_P_ij(k) s_Q_ij(k)])~=0
            highlight(p,I(k),J(k),'LineWidth',norm([s_P_ij(k) s_Q_ij(k)]));  % 把线路宽度设置为潮流大小
        end
    end
    os=0.1; % offset value to print text in the figure
    text(p.XData-os, p.YData-os, num2str(sqrt(s_v)),'HorizontalAlignment', 'center','FontSize', 10, 'Color','g'); % printf the nodes' voltage (sqrt(v_j)).
    N_cap=find(judge==1);
    N_pv=find(judge==0);
    text(p.XData(N_cap)-os, p.YData(N_cap)+os, num2str(s_Qi(N_cap)*1e6/Sbase,2),'HorizontalAlignment', 'center','FontSize', 10, 'Color','r'); % printf the nodes' Qi if it's a node with Capacitator.
    text(p.XData(N_pv)+os, p.YData(N_pv)-3*os, num2str(s_Pi(N_pv)*1e6/Sbase,2),'HorizontalAlignment', 'center','FontSize', 10, 'Color','k'); % printf the nodes' Pi if it's a node with PVs.
    for k=1:L
        Coor1_PQ=(p.XData(I(k))+p.XData(J(k)))/2;
        Coor2_PQ=(p.YData(I(k))+p.YData(J(k)))/2;
        text(Coor1_PQ, Coor2_PQ, num2str(s_P_ij(k)+s_Q_ij(k)*w,2),'HorizontalAlignment', 'center','FontSize', 10, 'Color','b'); % printf the complex power in distribution lines(if any).
    end
    
    %% 观察SOC约束是否取等号
    for k=1:L
        i=I(k);%首节点
        SOC_gap(k)=(s_l_ij(k)+s_v(i)).^2-(4*(s_P_ij(k)*1e6/Sbase).^2+4*(s_Q_ij(k)*1e6/Sbase).^2+(s_l_ij(k)-s_v(i)).^2);
        if round(SOC_gap(k),3)==0 % always some errors in duality gap
            %display(['线路 ',num2str(I(k)),' ',num2str(J(k)),' SOCP松弛约束取等号！'])
        else
            display(['Line ',num2str(I(k)),' ',num2str(J(k)),' SOCP relaxation is not exact！',' Line impedance z = ',num2str(z(k)*Zbase,3), ' Ω'])
        end
    end
    SOC_gap=SOC_gap';
    
    %% 观察
    
    if WithSvolt==1
        title("Exact SOCP relaxation for DOPF：SCE 47 bus case with s∈S_{volt}")
    else
        title("SOCP relaxation for DOPF：SCE 47 bus case WITHOUT s∈S_{volt}")
    end
end