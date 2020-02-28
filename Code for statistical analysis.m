% Number of soil parameters (su, SPT, gama, PI) for all boreholes
m=0;
count_soil=zeros(n_borehole,1);
 for i=1:n_sample
   if depth_soil(i,1)<depth_soil(i+1,1)
     c=i;
     m=m+1;
     count_soil(m,1)=c;
   end
 end
%--------------------------------------------------------------------------
% Extract soil parameters (su, SPT, gama, PI) for each soil layer (SOC, MC,
% SC, FS, VSC, SS) from all boreholes
A = regexp( fileread('soil_type.txt'), '\n', 'split');  % Soil type in Bangkok
for i=1:count_parameter(1,1)
  for j=1:count_soil(1,1)
    if (parameter(i,2)+parameter(i,3))/2>=depth_soil(j,2)&&(parameter(i,2)+parameter(i,3))/2<depth_soil(j,3)
        B{i,:}=A{j};
    end
    if (parameter(i,2)+parameter(i,3))/2>depth_soil(j,3)
        B{i,:}='0';
    end
  end
end
%--------------------------------------------------------------------------
% Simlation of Lognormal distribution for soil parameters
pdf_parameter = fitdist(parameter(:,4),'lognormal')
M_parameter=exp(mu+0.5*sigma^2)
std_parameter=sqrt(M_parameter^2*(exp(sigma^2)-1))
cov_parameter=std_parameter/M_parameter
y = pdf('lognormal',x,mu,sigma);
%--------------------------------------------------------------------------
% Calculation of value of soil parameters in Bangkok by separated into 1000mx2m
% squares
lx=[54000:1000:68000];  % domain size in horizontal direction
ly=[0:2:24];            % domain size in vertical direction
p=0;
q=0;
for i=1:14
 for j=1:12   
   for k=1:n_sample
    if l_p(k,1)>=lx(i,j)&&l_p(k,1)<lx(i,j+1)&&l_p(k,2)>=ly(i,j)&&l_p(k,2)<ly(i,j+1)
    p=p+l_p(i,3);
    q=q+1;
    end  
   end
 end
m_p=p/q;
p(j,i)=m_p;
end
%--------------------------------------------------------------------------
% Simuation of horizontal and vertical correlation length: rho_lx and
% rho_ly
n1=11;      % number of sample in horizontal direction
n2=12;      % number of sample in vertical direction
log_p=log(p);
mu_p=sum(sum(log_p))/(n1*n2);
sigma=sum(sum((log_p-mu_p).^2))/(n1*n2-1);
X=log_p-mu_p;
for i=0:n1-1
 A_lx=X(:,1:n1-i);                                                                          
 B_lx=X(:,1+i:n1);
 C_lx=A_lx.*B_lx;
 D_lx=sum(C_lx);
 E_lx=sum(D_lx);
 F_lx=E_lx/(sigma*(n2*(n1-i)-1));
 rho_lx(i+1,1)=F_lx;
end

for j=0:n2-1
 A_ly=X(1:n2-j,:);
 B_ly=X(1+j:n2,:);
 C_ly=A_ly.*B_ly;
 D_ly=sum(C_ly);
 E_ly=sum(D_ly);
 F_ly=E_ly/(sigma*(n1*(n2-j)-1));
 rho_ly(j+1,1)=F_ly;
end

