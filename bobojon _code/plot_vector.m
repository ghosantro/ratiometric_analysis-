function [arrowx, arrowy] = plot_vector(x,y,u,v,scale,varargin)

n=size(x,1)*size(x,2);

x=reshape(x,1,n);
y=reshape(y,1,n);
u=reshape(u,1,n);
v=reshape(v,1,n);

rel_head=0.5; % length of arrow head relative to total length
% *****
min_head=0.0; % minimum length of arrow head

alfa=pi/6; % half opening angle of arrow head

ih=ishold;

x1=x;
x2=x+u*scale;
y1=y;
y2=y+v*scale;
r=sqrt((x2-x1).^2 + (y2-y1).^2);

retcos=[cos(alfa), -sin(alfa) ; sin(alfa), cos(alfa)];

% arrow heads
rel_headlength=max(min_head./(r*scale),rel_head*ones(size(r))); 
lvek1=retcos *[(x1-x2) ; (y1-y2)] .* rel_headlength([1;1],:);
lvek2=retcos'*[(x1-x2) ; (y1-y2)] .* rel_headlength([1;1],:);

% arrows
arrowx=[x1; x2; x2+lvek1(1,:); x2; x2+lvek2(1,:)];
arrowy=[y1; y2; y2+lvek1(2,:); y2; y2+lvek2(2,:)];

% if nargin<6 | isempty(linetype)
%    plot(arrowx,arrowy,'blue-')
% else
%    plot(arrowx,arrowy,linetype)
% end

if nargin==7 & ~isempty(mark)
   hold on
   plot(x,y,mark)
end

if nargin==8 & ~isempty(fillcolour)
   hold on
   ha=plot(x,y,mark);
   set(ha,'markerFaceColor',fillcolour)
end

if ih
   hold on
else
   hold off
end
