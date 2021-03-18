function J=numjacobian(x,measurement,input)

n=length(x); 
[d1,d2,fx] = COST_4SAIL(x,measurement,input); %#ok<ASGLU>
step=1e-6; % 
J = zeros(length(fx),n);

for k=1:n
   xstep = x;
   xstep(k)=x(k)+step;
   [d1,d2,fxstep] = COST_4SAIL(xstep,measurement,input);  %#ok<ASGLU>
   J(:,k)= (fxstep-fx)./step;
end;