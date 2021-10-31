function I=myincidence(s,t)
%% my own incidence function to avoid node's order mistakes in MATLAB own function I=incidence(G)
% Inputs: s starting nodes, t: to nodes
% Output: I the incidence matrix I(s(j),j)=1 and I(t(j),j)=-1 if line j is
% starting from s(j) to t(j)
MaxNode=max(max(s),max(t));
I=zeros(MaxNode,length(s));
for j=1:length(s)
    I(s(j),j)=1;
    I(t(j),j)=-1;
end
end