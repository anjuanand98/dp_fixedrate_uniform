%Xnopt stores quantized regions
%Yn stores quantized points
%x uniform in [a,b]
function main
clc;
close all;
clear all;
%initializing x 
num=50;%number of elements initialized in x range
Nval=[2 4 8 16 32];
numiter=40;%number of iterations
%range of x: [a,b]
a=0;
b=1;
for N=Nval%number of quantization levels
xval=linspace(a,b,num);%initializing samples in x
iter=1;%signifies new value being added to xval from quantization function
Xniter=zeros(N+1,numiter);%quantized region points for each iteration
%arbitrarily initializing stopping critera
endistortion=100;
endistortion1=0;
fe=abs(endistortion-endistortion1)/endistortion;
while fe>10^-3
    [xval Xnopt]=nonstrategic_quantization_hloop(xval,N);
    Xniter(:,iter)=Xnopt;
    %checking fractional difference in cost
    Yn=decoder(Xnopt,b-a);
    endistortion1=endistortion;
    %encoder distortion
    endistortion=encoderdistortion(Xnopt,Yn,b-a);
    fe=abs(endistortion-endistortion1)/endistortion;
    iter=iter+1;
end
%decoder distortion
dedistortion=decoderdistortion(Xnopt,Yn,b-a);
switch N
    case 2
        Xnopt2=Xnopt;
        Yn2=Yn;
        endistortion2=endistortion;
        dedistortion2=dedistortion;
    case 4
        Xnopt4=Xnopt;
        Yn4=Yn;
        endistortion4=endistortion;
        dedistortion4=dedistortion;
    case 8
        Xnopt8=Xnopt;
        Yn8=Yn;
        endistortion8=endistortion;
        dedistortion8=dedistortion;
    case 16
        Xnopt16=Xnopt;
        Yn16=Yn;
        endistortion16=endistortion;
        dedistortion16=dedistortion;
    case 32
        Xnopt32=Xnopt;
        Yn32=Yn;
        endistortion32=endistortion;
        dedistortion32=dedistortion;
end
end
%dataxcubeduniformfixed1.mat holds all the output data
save('dataxcubeduniformfixed.mat','Nval','Xnopt2','Yn2','endistortion2','dedistortion2','Xnopt4','Yn4','endistortion4','dedistortion4','Xnopt8','Yn8','endistortion8','dedistortion8','Xnopt16','Yn16','endistortion16','dedistortion16','Xnopt32','Yn32','endistortion32','dedistortion32');
%plotting encoder and decoder distortions
f=figure;
plot(Nval,[endistortion2, endistortion4, endistortion8, endistortion16, endistortion32],'*-');
hold on;
plot(Nval,[dedistortion2, dedistortion4, dedistortion8, dedistortion16, dedistortion32],'o-');
hold off;
legend({'encoder distortion','decoder distortion'},'FontSize',14);
xlabel('N');
ylabel('distortion');
title('fixed rate xcubed uniform N');
saveas(f,'fixedrate_xcubed_uniform.fig');
saveas(f,'fixedrate_xcubed_uniform.m');

function [y]=decoder(x,w)
N=length(x)-1;
y=zeros(1,N);
fun=@(xv) xv*(1/w);
for i=1:N
    y(i)=integral(fun,x(i),x(i+1),'ArrayValued',true)/((x(i+1)-x(i))*(1/w));%y=integral x fx / integral fx
end

function endistortion=encoderdistortion(x,y,w)
N=length(y);
endistortion=0;
for n=1:N
    fun1=@(xv) ((xv.^3-y(n)).^2).*(1/w);
    endistortion=endistortion+integral(@(xv) fun1(xv),x(n),x(n+1));  
end
    
function dedistortion=decoderdistortion(x,y,w)
N=length(y);
dedistortion=0;
for n=1:N
    fun1=@(xv) ((xv-y(n)).^2).*(1/w);
    dedistortion=dedistortion+integral(@(xv) fun1(xv),x(n),x(n+1));  
end 

function [xval Xnopt]=nonstrategic_quantization_hloop(xval,N)
%initializing x optimal end points
Xnopt=10^7*ones(N+1,1);
Xnopt(1)=xval(1);
Xnopt(end)=xval(end);
w=xval(end)-xval(1);
%initializing [alpha, beta] pairs for x
D1val=zeros((length(xval)-1)*length(xval)/2,2);%stores all possible intervals [alpha,beta]
k=1;%iteration over all possible (alpha,beta) pairs
for i=1:length(xval)-1%for each alpha value possible
    for j=i+1:length(xval)%for each beta value possible, given the alpha value
        D1val(k,1:2)=[xval(i) xval(j)];%(alpha,beta)
        k=k+1;
    end
end
Dnxoxval=D1val(1:length(xval)-1,2);%(xo,x) - possible x values
%initializing D1 for each [alpha, beta] pair in D1val
D1=zeros(size(D1val,1),1);%D1([alpha,beta])=minimum over y integral(alpha to beta (f(x,x-y)p(x)dx))
%calculating D1
for i=1:size(D1val,1)%iterate over possible (alpha,beta) pairs
    q=decoder([D1val(i,1);D1val(i,2)],w);
    D1(i)=encoderdistortion([D1val(i,1);D1val(i,2)],q,w);
end
Dn=10^7*ones(length(Dnxoxval),N-1);%stores Dn values for each n (2:N), for each x in (xo,x)
Xn=10^7*ones(length(Dnxoxval),N-1);%stores X_(n-1) values for each n (1:N-1), for each x in (xo,x)
indmat=10^7*ones(N-1,length(xval),length(xval));%additional x values that give a distortion close to the minimum
xval1=xval;
for n=2:N
    for i=n+1:length(xval)%Dn(xo,x(i)) - there should be sufficient levels between xo and x(i), so i starts from n, e.g., n=2 level quantization implies you need atleast three points in x 
        [minval xboundary ind alpharange]=nonstr_recursion_1(xval,i,n,Xn,Dn,Dnxoxval,D1val,D1);%finding Dn and X_(n-1) values for each (xo,x) for a given n level 
        indmat(n-1,i,alpharange(ind))=1;
        %indexed such that x(i) location in Dnxoxval gives corresponding Dn and Xn values
        Dn(find((Dnxoxval==xval(i))),n-1)=minval;
        Xn(find((Dnxoxval==xval(i))),n-1)=xboundary;
    end
end
xval1=xval;
for n=N:-1:2%backward iteration to find Xnopt
    Xnopt(n)=Xn(find((Dnxoxval==Xnopt(n+1))),n-1);%X_(n-1) optimal=X_(n-1)(xo,X_n opt)
    temp=xval(find(indmat(n-1,find(xval==Xnopt(n+1)),:)==1));
    if length(temp)==1
    Xnopt(n)=temp;
    else
        for it=1:length(temp)
            %sampling more around possible decision levels
            xval1=[xval1 linspace(xval1(find(xval1==temp(it))-1),xval1(find(xval1==temp(it))+1),20)];
            xval1=unique(xval1);
            xval1=sort(xval1);
        end
    end
    %sampling more around decision levels chosen in this iteration
    xvt=linspace(xval1(find(xval1==Xnopt(n))-1),xval1(find(xval1==Xnopt(n))+1),20);
    xval1=[xval1 xvt];
    xval1=unique(xval1);
    xval1=sort(xval1);
end
%removing values far from optimal x region points
t1=1;
for n=2:N
    t3=floor((find(xval1==Xnopt(n))+t1)/2);
    xval1(t1+1:t3)=[];
    t1=floor((find(xval1==Xnopt(n))+find(xval1==Xnopt(n+1)))/2);
end
xval1(t1:end-1)=[];
xval=[xval1 Xnopt'];
xval=sort(xval);
xval=unique(xval); 
%removing values too close to each other
indz=[];
for i=1:length(xval)-1
    if norm(xval(i+1)-xval(i))/norm(xval(i))<10^-4
        indz=[indz i];
    end
end
xval(indz)=[];

function [minval xboundary ind alpharange]=nonstr_recursion_1(xval,xind,n,Xn,Dn,Dnxoxval,D1val,D1)%returns distortion and corresponding value of x and theta
    if n==1%n=1 values are already computed in D1
        minval=D1(find(D1val(:,1)==xval(1) & D1val(:,2)==xval(xind)));
        xboundary=xval(xind);
        return;
    end
    alpharange=[n:xind-1];%indices of alpha values possible in x range (xo<alpha<x)
    arr=zeros(length(alpharange),1);%for each alpha value possible, finding minimum (D_(n)(xo,x)=min over alpha, ao<alpha<x D_(n-1)(xo,alpha)+D1(alpha,x))
    %D_n(xo<alpha<x, thetao<beta<theta)
    for arrx=1:length(alpharange)
        if n-2>=1 && Xn(find((Dnxoxval==xval(alpharange(arrx)))),n-2)~=10^7%referring to a variable for values already computed
            xboundary1=Xn(find((Dnxoxval==xval(alpharange(arrx)))),n-2);
            minval1=Dn(find((Dnxoxval==xval(alpharange(arrx)))),n-2);
        else
        [minval1 xboundary1]=nonstr_recursion_1(xval,alpharange(arrx),n-1,Xn,Dn,Dnxoxval,D1val,D1);%D_(n-1)(xo,alpha) 
        end
        %D_n(xo,x)=min over alpha, xo<alpha<x (D_(n-1)(xo,alpha)+D1(alpha,x))
        arr(arrx)=minval1+D1(find(D1val(:,1)==xval(alpharange(arrx)) & D1val(:,2)==xval(xind)));
    end
    minval=min(arr);%returning minimum distortion value
    in=find(arr==minval);
    %finding possible solutions with distortion negligibly close to (or equal to) the minimum distortion
    ind=find(abs(arr-minval)/minval<10^-4);
    xboundary=xval(alpharange(in(1)));