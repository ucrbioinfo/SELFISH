function output = SELFISH_DCI(contc1,contc2,norm1,norm2,THRESHOLD,RESOLUTION,INTERVAL)
%SELFISH_DCI find differential chromatin interactions between two contact
%maps
%
%   output =
%   SELFISH_DCI(contc1,contc2,norm1,norm2,THRESHOLD,RESOLUTION,INTERVAL)
%   finds DCIs between two contact maps given by filenames contc1 and
%   contc2 with p-values less than THRESHOLD in the specified INTERVAL. 
%
%   'contc1'           -   contact map filename 1
%
%   'norm1'            -   normalization vector filename 1
%
%   'norm2'            -   normalization vector filename 2
%
%   'THRESHOLD'        -   THRESHOLD at which DCIs are return
%
%   'RESOLUTION'       -  Data resoultion in bp
%
%   'INTERVAL'         - The interval in bp for which DCIs are detected.

    INTERVAL = ceil(INTERVAL/RESOLUTION);
    % Read contact-map 1 and 2
    disp('Reading contact maps...');
    fid = fopen(contc1);
    H1 = textscan(fid,'%d \t %d \t %f');
    L = min([size(H1{1}(:),1) size(H1{2}(:),1) size(H1{3}(:),1)]);

    H1 = [H1{1}(1:L) H1{2}(1:L) H1{3}(1:L)];
    H1(:,1:2) = ceil(H1(:,1:2)/RESOLUTION)+1;
    H1 = double(H1);
    fid = fopen(contc2);
    H2 = textscan(fid,'%d \t %d \t %f');

    L = min([size(H2{1}(:),1) size(H2{2}(:),1) size(H2{3}(:),1)]);
    H2 = [H2{1}(1:L) H2{2}(1:L) H2{3}(1:L)];
    H2(:,1:2) = ceil(H2(:,1:2)/RESOLUTION)+1;
    H2 = double(H2);
    % % KR Normalize the contact maps
    disp('Normalizing contact maps...');
    KR_norm = load(norm1);
    KR_norm(KR_norm<0.1) = nan;
    H1(:,3) = H1(:,3)./(KR_norm(H1(:,1)).*KR_norm(H1(:,2)));
    H1(isnan( H1(:,3)),3) = 0;

    KR_norm = load(norm2);
    KR_norm(KR_norm<0.1) = nan;
    H2(:,3) = H2(:,3)./(KR_norm(H2(:,1)).*KR_norm(H2(:,2)));
    H2(isnan( H2(:,3)),3) = 0;

    % find distances between interacting loci
    interaction_dist1 = H1(:,2) - H1(:,1);
    interaction_dist2 = H2(:,2) - H2(:,1);

    [a,~,c]  = unique(interaction_dist1);
    avgDiag1 = zeros(max(interaction_dist1)+1,1);
    avgSTD1 = zeros(max(interaction_dist1)+1,1);
    avgDiag1(a+1) = accumarray(c,H1(:,3),[],@(x) mean(x,1));
    avgSTD1(a+1)  = accumarray(c,H1(:,3),[],@(x) std(x,1));

    [a,~,c]  = unique(interaction_dist2);
    avgDiag2 = zeros(max(interaction_dist2)+1,1);
    avgSTD2 = zeros(max(interaction_dist2)+1,1);
    avgDiag2(a+1) = accumarray(c,H2(:,3),[],@(x) mean(x,1));
    avgSTD2(a+1)  = accumarray(c,H2(:,3),[],@(x) std(x,1));

    %% normalize contact maps with respect to their distances up to 2Mb
    normH1 = H1;
    normH1(:,3) = 0;
    normH2 = H2;
    normH2(:,3) = 0;

    [a,~,c]  = unique(interaction_dist1);
    normH1(:,3) = (H1(:,3)-avgDiag1(interaction_dist1+1))./avgSTD1(interaction_dist1+1);
    normH1(isnan( normH1(:,3)),3) = 0;

    [a,~,c]  = unique(interaction_dist2);
    normH2(:,3) = (H2(:,3)-avgDiag2(interaction_dist2+1))./avgSTD2(interaction_dist2+1);
    normH2(isnan( normH2(:,3)),3) = 0;

    %% Find the differences for the specified interval only

    H1 = H1(H1(:,1)>=INTERVAL(1)&H1(:,1)<=INTERVAL(2)&H1(:,2)>=INTERVAL(1)&H1(:,2)<=INTERVAL(2),:);
    normH1 = normH1(normH1(:,1)>=INTERVAL(1)&normH1(:,1)<=INTERVAL(2)&normH1(:,2)>=INTERVAL(1)&normH1(:,2)<=INTERVAL(2),:);
    H2 = H2(H2(:,1)>=INTERVAL(1)&H2(:,1)<=INTERVAL(2)&H2(:,2)>=INTERVAL(1)&H2(:,2)<=INTERVAL(2),:);
    normH2 = normH2(normH2(:,1)>=INTERVAL(1)&normH2(:,1)<=INTERVAL(2)&normH2(:,2)>=INTERVAL(1)&normH2(:,2)<=INTERVAL(2),:);
    H1(:,1:2) = H1(:,1:2) - INTERVAL(1) + 1;
    normH1(:,1:2) = normH1(:,1:2) - INTERVAL(1) + 1;
    H2(:,1:2) = H2(:,1:2) - INTERVAL(1) + 1;
    normH2(:,1:2) = normH2(:,1:2) - INTERVAL(1) + 1;

    %% convert (i,j,v) matrices to full matrices
    intvLen = INTERVAL(2)-INTERVAL(1)+1;

    H1 = sparse(H1(:,1),H1(:,2),H1(:,3),intvLen,intvLen);
    normH1 = sparse(normH1(:,1),normH1(:,2),normH1(:,3),intvLen,intvLen);
    if sum(sum(triu(H1,1)))==0
        H1=tril(H1)+tril(H1,1)';
        normH1=tril(normH1)+tril(normH1,1)';
    else
        H1=triu(H1)+triu(H1,1)';
        normH1=triu(normH1)+triu(normH1,1)';
    end
    normH1 = full(normH1);
    H1 = full(H1);

    H2 = sparse(H2(:,1),H2(:,2),H2(:,3),intvLen,intvLen);
    normH2 = sparse(normH2(:,1),normH2(:,2),normH2(:,3),intvLen,intvLen);
    if sum(sum(triu(H2,1)))==0
        H2=tril(H2)+tril(H2,1)';
        normH2=tril(normH2)+tril(normH2,1)';
    else
        H2=triu(H2)+triu(H2,1)';
        normH2=triu(normH2)+triu(normH2,1)';
    end
    normH2 = full(normH2);
    H2 = full(H2);
    
    disp('Finding differences...');
    %% Gaussian impact radii
    % keep zero-interaction 2D coordinates
    zIndc1 = (H1(:) == 0) ; 
    zIndc2 = (H2(:) == 0) ;
    lidx = tril(true(size(normH1))); % Lower triangular half

    normH1(zIndc1) = 0;
    normH2(zIndc2) = 0;
    normH1(lidx) = 0;
    normH2(lidx) = 0;

    indZ = zIndc1 | zIndc2 | lidx(:);
    indNZ = indZ == 0;

    nNZpts = sum(indNZ); 
    PVAL = ones(nNZpts,1);
    Scales  = zeros(nNZpts,1);

    % Set the method's parameters
    sigma0 = 1.6;
    s = 10;
    scales = zeros(s+2,1);

    % compute the first three filtered images (G1, G2, G3)
    scales(1) = 2*(2*sigma0)+1;
    Gc1 = imgaussfilt(normH1,sigma0*(2^((2-1)/s)));
    scales(2) = 2*(2*(sigma0*(2^((2-1)/s))))+1;
    Gn1 = imgaussfilt(normH1,sigma0*(2^((3-1)/s)));
    scales(3) = 2*(2*(sigma0*(2^((3-1)/s))))+1;

    Gc2 = imgaussfilt(normH2,sigma0*(2^((2-1)/s)));
    Gn2 = imgaussfilt(normH2,sigma0*(2^((3-1)/s)));

    Lc = ((Gc1-Gn1) - (Gc2-Gn2));
    Lc(indZ) = 0;

    % now we iterate for levels of impact radii
    for level= 4:s+2
        Gc1 = Gn1;
        Gn1 = imgaussfilt(normH1,sigma0*(2^((level-1)/s)));
        Gc2 = Gn2;
        Gn2 = imgaussfilt(normH2,sigma0*(2^((level-1)/s)));
        scales(level) = 2*(2*(sigma0*(2^((level-1)/s))))+1;
        Ln = ((Gc1-Gn1) - (Gc2-Gn2));
        Ln(indZ) = 0;

        % compute p-values for the current level
        pd = fitdist(Lc(indNZ),'Normal');
        pval = cdf(pd,Lc(indNZ),'upper');
        pval(pval>0.5) = 1 - pval(pval>0.5);
        pval = pval*2;
        % keep the largest p-value for each paired index as the
        % p-value of the index
        toChangeIndx = pval < PVAL;
        PVAL(toChangeIndx) = pval(toChangeIndx);
        Scales(toChangeIndx) = scales(level); 
        % move to the next impact radius
        %Lp = Lc;
        Lc = Ln;
    end

    correctedPVAL = mafdr(PVAL,'BHFDR',1);
    [ptx,pty] = ind2sub(size(normH1),find(indNZ)); 
    % discard significant DCIs in sparse regions
    %%%
    %
    sigIndx = correctedPVAL < THRESHOLD;
    output = [(ptx(sigIndx)-1)*RESOLUTION+INTERVAL(1) (pty(sigIndx)-1)*RESOLUTION+INTERVAL(1) correctedPVAL(sigIndx)]; 
end    
