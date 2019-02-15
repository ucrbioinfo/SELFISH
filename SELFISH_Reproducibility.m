function S = SELFISH_Reproducibility(contc1,contc2)
%SELFISH_Reproducibility computes reproducibility score between two contact
%maps
%
%   S =   SELFISH_Reproducibility(contc1,RESOLUTION)
%   Computes reproducibility score between two contact maps given by filenames/(i,j,v) triplets contc1 and
%   contc2 with p-values less than THRESHOLD in the specified INTERVAL. 
%
%   'contc1'           -   contact map [file address] 1 in Rao format
%
%   'contc2'           -   contact map [file address] 2 in Rao format

numOfSteps = 200;
rel_loc  = (0:0.005:1)';
indc = nchoosek(1:numOfSteps,2);
flag = zeros(2,nchoosek(numOfSteps,2));
% read contact maps 1 and 2
if ischar(contc1)  
    disp('Reading contact map 1...');
    fid = fopen(contc1);
    H1 = textscan(fid,'%d \t %d \t %f');

    H1 = [H1{1} H1{2} H1{3}];
    H1 = double(H1);
    fclose(fid);
elseif isnumeric(contc1) & ismatrix(contc1) & (size(contc1,2)==3)
    H1 = contc1; 
end

if ischar(contc2)
    disp('Reading contact map 2...');
    fid = fopen(contc2);
    H2 = textscan(fid,'%d \t %d \t %f');
    H2 = [H2{1} H2{2} H2{3}];
    H2 = double(H2);
    fclose(fid);
elseif isnumeric(contc2) & ismatrix(contc2) & (size(contc2,2)==3)
    H2 = contc2; 
end
            
N = max(max(H1(:,1)),max(H1(:,2)));
vals = zeros(numOfSteps,1);
I = ceil(rel_loc(1:end-1)*N)+1;
J = ceil(I + 2*N/(numOfSteps));
% compute the sum for each block along the diagonal
for j=1:length(I)
    interval = [I(j) J(j)];
    vals(j) = sum(H1(H1(:,1)>=interval(1)&H1(:,1)<=interval(2)&H1(:,2)>=interval(1)&H1(:,2)<=interval(2),3));
end
flag(1,:) = vals(indc(:,1)) > vals(indc(:,2));

vals = zeros(numOfSteps,1);
for j=1:length(I)
    interval = [I(j) J(j)];
    vals(j) = sum(H2(H2(:,1)>=interval(1)&H2(:,1)<=interval(2)&H2(:,2)>=interval(1)&H2(:,2)<=interval(2),3));
end
flag(2,:) = vals(indc(:,1)) > vals(indc(:,2));
 
D = sum((flag(1,:)-flag(2,:)).^2)/size(flag,2);
S = exp(-5*D);
% or --> S = 1 - D;
end