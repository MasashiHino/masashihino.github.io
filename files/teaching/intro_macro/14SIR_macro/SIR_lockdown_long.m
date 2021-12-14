%-----------------
% code for SIR+Death
% initially provided by Taisuke Nakata
% modified by Masashi Hino
%-----------------

clc;clear all;
close all;
T=150;
beta = 0.5852;

gamma_D = 7/18*0.005;%7/18*0.005;
gamma_R   =7/18*(1-0.005) ;
gamma_S = 7/18*0.0;

S=zeros(1,T);%population of susceptible at all t
I=zeros(1,T);%population of infectious at all t
R=zeros(1,T);% population of recovered at all t
D=zeros(1,T);% number of death at all t

X0    = beta/gamma_R;
Sstar = 1/X0;
Rstar = 1-Sstar;

% Initial Condition
I(1) = 10^(-6);
S(1) = 1 - I(1);
R(1) = 0;
D(1) =0;



% lockdown

for i=1:2
    if i==1      
        L=zeros(1,T);
        T1=51;% timing when lockdown begins
        T2=59;% timing when lockdown ends
        t_L=zeros(1,T);
        t_L(T1:T2)=1;
        L(T1:T2)=0.2;
    elseif i==2
        L=zeros(1,T);
        T1=51;% timing when lockdown begins
        T2=55;% timing when lockdown ends
        t_L=zeros(1,T);
        t_L(T1:T2)=1;
        L(T1:T2)=0.2;
    end
    figure(1)
    hold on
    area(t_L*100)
    
    figure(2)
    subplot(221)
    hold on    
    area(t_L*100)

    subplot(222)
    hold on
    area(t_L*100)
 
    subplot(223)
    hold on
    area(t_L*100)

    subplot(224)
    hold on
    area(t_L*100)    
    
% end
% for i=2:-1:1
%     if i==1     
%         T1=51;% timing when lockdown begins
%         T2=59;% timing when lockdown ends
%         t_L=zeros(1,T);
%         t_L(T1:T2)=1;
%         L(T1:T2)=0.2;
%     elseif i==2
%         T1=51;% timing when lockdown begins
%         T2=55;% timing when lockdown ends
%         t_L=zeros(1,T);
%         t_L(T1:T2)=1;
%         L(T1:T2)=0.2;
%     end

    for t=2:T
        S(t) = S(t-1) - beta*(1-L(t-1))*S(t-1)*(1-L(t-1))*I(t-1)+gamma_S*R(t-1);
        I(t) = I(t-1) + beta*(1-L(t-1))*S(t-1)*(1-L(t-1))*I(t-1)- gamma_R*I(t-1)-gamma_D*I(t-1);
        R(t) = R(t-1) + gamma_R*I(t-1) - gamma_S*R(t-1);
        D(t) = D(t-1) + gamma_D*I(t-1);
    end
    Trans=beta*S.*I;
    
    %%
    figure(1)
    hold on

    plot([1:T],Trans*100,'k','LineWidth',1.5)
    ylim([min(Trans*100),3])
    title('V‹KŠ´õŽÒ”(T)','FontSize',20)
    xlabel('ŽžŠÔ(T)','FontSize',14,'interpreter','latex')
    ylabel('(%)','FontSize',12)
    
    %%
    
    figure(2)
    subplot(2,2,1)
    hold on
    plot([1:T],I*100,'k','LineWidth',1.5);
    ylim([min(I*100),7])
    title('Š´õŽÒ(I)','FontSize',20,'interpreter','latex')
    xlabel('ŽžŠÔ(T)','FontSize',14,'interpreter','latex')
    ylabel('(%)','FontSize',12)
    
    subplot(2,2,2)
    hold on
    plot([1:T],S*100,'k','LineWidth',1.5);
    ylim([min(S*100),100])
    %plot([1:T],ones(1,T)*Sstar,'k--','LineWidth',1.5);
    title('Œ’N‚Èl(S)','FontSize',20,'interpreter','latex')
    
    
    
    
    subplot(2,2,3)
    hold on
    plot([1:T],R*100,'k','LineWidth',1.5);
    ylim([min(R*100),60])
    %plot([1:T],ones(1,T)*Rstar,'k--','LineWidth',1.5);
    title('‰ñ•œ‚µ‚½l(R)','FontSize',20,'interpreter','latex')
    
    subplot(2,2,4)
    hold on
    plot([1:T],D*100,'k','LineWidth',1.5);
    ylim([min(D*100),0.3])
    title('Ž€–SŽÒ”(D)','FontSize',20,'interpreter','latex')
    
end
%%

figure(1)
legend('’·ŠúŠÔ','’·ŠúŠÔ','’ZŠúŠÔ','’ZŠúŠÔ')

figure(2)
subplot(221)
legend('’·ŠúŠÔ','’·ŠúŠÔ','’ZŠúŠÔ','’ZŠúŠÔ')

%print(figure(1),'-dpdf','sir.pdf')
print(figure(1),'-dpng','sir.png')

%%

function applyhatch(h,patterns,colorlist)
%APPLYHATCH Apply hatched patterns to a figure
%  APPLYHATCH(H,PATTERNS) creates a new figure from the figure H by
%  replacing distinct colors in H with the black and white
%  patterns in PATTERNS. The format for PATTERNS can be
%    a string of the characters '/', '\', '|', '-', '+', 'x', '.'
%    a cell array of matrices of zeros (white) and ones (black)
%
%  APPLYHATCH(H,PATTERNS,COLORS) maps the colors in the n by 3
%  matrix COLORS to PATTERNS. Each row of COLORS specifies an RGB
%  color value.
%
%  Note this function makes a bitmap image of H and so is limited
%  to low-resolution, bitmap output.
%
%  Example 1:
%    bar(rand(3,4));
%    applyhatch(gcf,'\-x.');
%
%  Example 2:
%    colormap(cool(6));
%    pie(rand(6,1));
%    legend('Jan','Feb','Mar','Apr','May','Jun');
%    applyhatch(gcf,'|-+.\/',cool(6));
%
%  See also: MAKEHATCH
%  Copyright 2002-2009 The MathWorks, Inc.
  
oldppmode = get(h,'paperpositionmode');
oldunits = get(h,'units');
set(h,'paperpositionmode','auto');
set(h,'units','pixels');
figsize = get(h,'position');
if nargin == 2
  colorlist = [];
end
if verLessThan('matlab','8.4.0')
  bits = hardcopy(h,'-dzbuffer','-r0');
else
  bits = print(h,'-RGBImage','-r0');
end
set(h,'paperpositionmode',oldppmode);
bwidth = size(bits,2);
bheight = size(bits,1);
bsize = bwidth * bheight;
if ~isempty(colorlist)
  colorlist = uint8(255*colorlist);
  [colors,colori] = nextnonbw(0,colorlist,bits);
else
  colors = (bits(:,:,1) ~= bits(:,:,2)) | ...
	   (bits(:,:,1) ~= bits(:,:,3));
end
pati = 1;
colorind = find(colors);
while ~isempty(colorind)
  colorval(1) = bits(colorind(1));
  colorval(2) = bits(colorind(1)+bsize);
  colorval(3) = bits(colorind(1)+2*bsize);
  if iscell(patterns)
    pattern = patterns{pati};
  elseif isa(patterns,'char')
    pattern = makehatch(patterns(pati));
  else
    pattern = patterns;
  end
  pattern = uint8(255*(1-pattern));
  pheight = size(pattern,2);
  pwidth = size(pattern,1);
  ratioh = ceil(bheight/pheight);
  ratiow = ceil(bwidth/pwidth);
  bigpattern = repmat(pattern,[ratioh ratiow]);
  if ratioh*pheight > bheight
    bigpattern(bheight+1:end,:) = [];
  end
  if ratiow*pwidth > bwidth
    bigpattern(:,bwidth+1:end) = [];
  end
  bigpattern = repmat(bigpattern,[1 1 3]);
  color = (bits(:,:,1) == colorval(1)) & ...
	  (bits(:,:,2) == colorval(2)) & ...
	  (bits(:,:,3) == colorval(3));
  color = repmat(color,[1 1 3]);
  bits(color) = bigpattern(color);
  if ~isempty(colorlist)
    [colors,colori] = nextnonbw(colori,colorlist,bits);
  else
    colors = (bits(:,:,1) ~= bits(:,:,2)) | ...
	     (bits(:,:,1) ~= bits(:,:,3));
  end
  colorind = find(colors);
  pati = (pati + 1);
  if pati > length(patterns)
    pati = 1;
  end
end
newfig = figure('units','pixels','visible','off');
imaxes = axes('parent',newfig,'units','pixels');
im = image(bits,'parent',imaxes);
fpos = get(newfig,'position');
set(newfig,'position',[fpos(1:2) figsize(3) figsize(4)+1]);
set(imaxes,'position',[0 0 figsize(3) figsize(4)+1],'visible','off');
set(newfig,'visible','on');

function [colors,out] = nextnonbw(ind,colorlist,bits)
out = ind+1;
colors = [];
while out <= size(colorlist,1)
  if isequal(colorlist(out,:),[255 255 255]) | ...
	isequal(colorlist(out,:),[0 0 0])
    out = out+1;
  else
    colors = (colorlist(out,1) == bits(:,:,1)) & ...
	     (colorlist(out,2) == bits(:,:,2)) & ...
	     (colorlist(out,3) == bits(:,:,3));
    return
  end
end
end
end

function A = makehatch(hatch)
%MAKEHATCH Predefined hatch patterns
%  MAKEHATCH(HATCH) returns a matrix with the hatch pattern for HATCH
%   according to the following table:
%      HATCH        pattern
%     -------      ---------
%        /          right-slanted lines
%        \          left-slanted lines
%        |          vertical lines
%        -          horizontal lines
%        +          crossing vertical and horizontal lines
%        x          criss-crossing lines
%        .          single dots
%
%  See also: APPLYHATCH
%  Copyright 2002-2009 The MathWorks, Inc.
n = 6;
A=zeros(n);
switch (hatch)
 case '/'
  A = fliplr(eye(n));
 case '\'
  A = eye(n);
 case '|'
  A(:,1) = 1;
 case '-'
  A(1,:) = 1;
 case '+'
  A(:,1) = 1;
  A(1,:) = 1;
 case 'x'
  A = eye(n) | fliplr(diag(ones(n-1,1),-1));
 case '.'
  A(1:2,1:2)=1;
 otherwise
  error(['Undefined hatch pattern "' hatch '".']);
end
end

