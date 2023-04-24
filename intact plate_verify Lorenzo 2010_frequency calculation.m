%rectangular plate with uniform elastic supports on boundaries, in-plane vibration
clear all
clc
format long
N=6;									%the number of terms in each direction
a=0.8; 
b=1*a; 									%a, b are the plate dimensions in x and y directions
h=0.008; 								%thickness
nu=0.3; 								%the Poisson's ratio
E=73e9; 								%Young's modulus
Rho=2800; 								%density
A11=E*h/(1-nu^2);
A22=A11;
A12=nu*A11;
A66=E*h/(2*(1+nu));
k_dimensionless=2;					%the dimensionless elastic edge stiffness
k0=k_dimensionless*E*h/(a/2*(1-nu^2));	%don't miss h, very important
k1U=k0;
k1V=0;
k2U=0;
k2V=k0;
k3U=k0;
k3V=0;
k4U=0;
k4V=k0;
%%%%%%%%%%%%%%%%
%admissible functions
%%%%%%%%%%%%%%%%
syms Xi
y(1)=(sin(pi/4*Xi+3*pi/4))^2;
y(2)=sin(pi/4*Xi+3*pi/4)*sin(-pi/2*Xi-3*pi/2);
y(3)=(sin(pi/4*Xi-3*pi/4))^2;
y(4)=sin(pi/4*Xi-3*pi/4)*sin(pi/2*Xi-3*pi/2);
for m=5:N
	y(m)=sin(pi/2*(m-4)*Xi+pi/2*(m-4))*sin(pi/2*Xi+pi/2);
end
%%%%%%%%%%%%%%%%
%Imr00
%%%%%%%%%%%%%%%%
Imr00=zeros(N,N);
for m=1:N
	for n=1:N
		Imr00(m,n)=int(y(m)*y(n),-1,1);
	end
end
%%%%%%%%%%%%%%%%
%Imr11
%%%%%%%%%%%%%%%%
Imr11=zeros(N,N);
for m=1:N
	for n=1:N
		Imr11(m,n)=int(diff(y(m))*diff(y(n)),-1,1);
	end
end
%%%%%%%%%%%%%%%%
%Imr01
%%%%%%%%%%%%%%%%
Imr01=zeros(N,N);
for m=1:N
	for n=1:N
		Imr01(m,n)=int(y(m)*diff(y(n)),-1,1);
	end
end
Ins00=Imr00;
Ins11=Imr11;
Imr10=Imr01';
Ins10=Imr10;
Ins01=Imr01;
%%%%%%%%%%%%%%%%
%Emr & Ens
%%%%%%%%%%%%%%%%
EmrNeg1=zeros(N,N);
EmrPos1=zeros(N,N);
for m=1:N
	for r=1:N
		EmrNeg1(m,r)=subs(y(m),Xi,-1)*subs(y(r),Xi,-1);
		EmrPos1(m,r)=subs(y(m),Xi,1)*subs(y(r),Xi,1);
	end
end
EnsNeg1=EmrNeg1;
EnsPos1=EmrPos1;
%%%%%%%%%%%%%%%%
%Jmr & Jns
%%%%%%%%%%%%%%%%
Jns1U=k1U*Ins00;
Jns3U=k3U*Ins00;
Jmr2U=k2U*Imr00;
Jmr4U=k4U*Imr00;
Jns1V=k1V*Ins00;
Jns3V=k3V*Ins00;
Jmr2V=k2V*Imr00;
Jmr4V=k4V*Imr00;
%%%%%%%%%%%%%%%%%%% stiffness matrix %%%%%%%%%%%%%%%%%%%
K=zeros(2*N^2,2*N^2);
for m=1:N
	for n=1:N
		u=n+(m-1)*N;
		for r=1:N
			for s=1:N
				v=s+(r-1)*N;
				K(u,v)=K(u,v)+A11*(b/a)*Imr11(m,r)*Ins00(n,s)+A66*(a/b)*Imr00(m,r)*Ins11(n,s)+(b/2)*(EmrNeg1(m,r)*Jns1U(n,s)+EmrPos1(m,r)*Jns3U(n,s))+(a/2)*(Jmr2U(m,r)*EnsNeg1(n,s)+Jmr4U(m,r)*EnsPos1(n,s));
				K(u,v+N^2)=K(u,v+N^2)+A12*Imr10(m,r)*Ins01(n,s)+A66*Imr01(m,r)*Ins10(n,s);
				K(u+N^2,v)=K(u+N^2,v)+A12*Imr01(m,r)*Ins10(n,s)+A66*Imr10(m,r)*Ins01(n,s);
				K(u+N^2,v+N^2)=K(u+N^2,v+N^2)+A22*(a/b)*Imr00(m,r)*Ins11(n,s)+A66*(b/a)*Imr11(m,r)*Ins00(n,s)+(b/2)*(EmrNeg1(m,r)*Jns1V(n,s)+EmrPos1(m,r)*Jns3V(n,s))+(a/2)*(Jmr2V(m,r)*EnsNeg1(n,s)+Jmr4V(m,r)*EnsPos1(n,s));
			end
		end
	end
end
%%%%%%%%%%%%%%%%%%% mass matrix %%%%%%%%%%%%%%%%%%%
M=zeros(2*N^2,2*N^2);
for m=1:N
	for n=1:N
		u=n+(m-1)*N;
		for r=1:N
			for s=1:N
				v=s+(r-1)*N;
				M(u,v)=M(u,v)+(Rho*h*a*b/4)*Imr00(m,r)*Ins00(n,s);
				M(u+N^2,v+N^2)=M(u+N^2,v+N^2)+(Rho*h*a*b/4)*Imr00(m,r)*Ins00(n,s);
			end
		end
	end
end
X=zeros(2*N^2,2*N^2);
eigenvalues=zeros(2*N^2,2*N^2);
[X,eigenvalues] = eig(K,M);
Eigenvalues=diag(eigenvalues);
Hz=Eigenvalues.^0.5/(2*pi);
Hz=sort(real(Hz));
Lambda=Hz*2*pi*a/2*(Rho*(1-nu^2)/E)^0.5
