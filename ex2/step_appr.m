function out = step_appr(x,subdiv_num,subdiv_length)
%computes the step-function approximation of 'x' by averaging it
%over equally-spaced subintervals of length 'subdiv_length'. 

%Please pre-compute the quantity of subintervals 'subdiv_num', as
%this function may be called multiple times. For the lazy, this
%safeguard is here:
if nargin<3
    subdiv_length = length(x)/subdiv_num;
end
%IF THIS FUNCTION THROWS AN ERROR OR MISBEHAVES:
	%it would probably be because it needs this to be satisfied
	%exactly:
	%length(x) / subdiv_num == subdiv_length.
if (length(x) / subdiv_num) ~= subdiv_length
   fprintf('ERROR: Incompatible size of x and subdiv_num')
end

%Methodology: 
%(1) Reshape into a matrix where each column represents the group
%to be averaged, then 
%(2) Compute averages along the right direction with 'mean', then
%(3) Use a kronecker product to repeat values where they need to
%be.
out = kron(ones(subdiv_length,1),...
      mean(reshape(x,subdiv_length,subdiv_num),1));
out = out(:);
end
