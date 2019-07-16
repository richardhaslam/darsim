x=[0;6;2;8;4;10;6;12;8];
y=[0;2;4;6;8;10;12;14;16];

t=numel(x);
r=x(t);
nnn=y(t);

s1=subplot(1,1,1);
set(s1,'xlim',[1 numel(x)],'ylim',sort([0 r]),'nextplot','add')

p1=plot(nan,nnn,'r-');
p2=plot(r,nnn,'g*');

axis([0 20 0 20])
set(gca,'ydir','rev')

for t=1:numel(x)
    r=x(t);
    nnn=y(t);

      set(p1,'xdata',[get(p1,'xdata') r],'ydata',[get(p1,'ydata') nnn])
      set(p2,'xdata',r,'ydata',nnn)

      pause(.5)

end
