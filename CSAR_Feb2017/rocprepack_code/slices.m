clear
cla
clf

%close all
for i=1:51
  C(i) = 0.02*(i-1)*2*pi;
end

ncuts=9;
for ii=1:ncuts

if ii==1; load fort.200; end;
if ii==2; load fort.201; end;
if ii==3; load fort.202; end;
if ii==4; load fort.203; end;
if ii==5; load fort.204; end;
if ii==6; load fort.205; end;
if ii==7; load fort.206; end;
if ii==8; load fort.207; end;
if ii==9; load fort.208; end;

n = fort(:,1);
m = fort(:,2);
x = fort(:,3);
y = fort(:,4);
z = fort(:,5);
r = fort(:,6);

tot = size(n);
tot = tot(1,1);

subplot(3,3,ii)
hold on
for i=1:tot
   radius = r(i);
   xshift = x(i);
   yshift = y(i);
    plot(cos(C)*radius+xshift,sin(C)*radius+yshift,'k');
    %if m(i)==1; plot(cos(C)*radius+xshift,sin(C)*radius+yshift,'k'); end;
    %if m(i)==2; plot(cos(C)*radius+xshift,sin(C)*radius+yshift,'r'); end;
    %if m(i)==3; plot(cos(C)*radius+xshift,sin(C)*radius+yshift,'b'); end;
%    plot(xshift,yshift,'x')
end
hold off
axis equal
axis([-1 1 -1 1])
set(gca,'FontSize',[12])
%set(gca,'XGrid','on','YGrid','on')
%set(gca,'YTick',[-1 0 1]);
%set(gca,'XTick',[-1 0 1]); 

if ii==1; title('z=-1'); ylabel('y'); set(gca,'XTick',[]); end;
if ii==2; title('z=-0.75'); set(gca,'XTick',[]); end;
if ii==3; title('z=-0.5'); set(gca,'XTick',[]); end;
if ii==4; title('z=-0.25'); ylabel('y'); set(gca,'XTick',[]); end;
if ii==5; title('z=0'); set(gca,'XTick',[]); end;
if ii==6; title('z=0.25'); set(gca,'XTick',[]); end;
if ii==7; title('z=0.5'); xlabel('x'); ylabel('y'); end;
if ii==8; title('z=0.75'); xlabel('x'); end;
if ii==9; title('z=1'); xlabel('x'); end;

set(gca,'Box','on')

end

print 3d_slices -djpeg90
print 3d_slices -deps
