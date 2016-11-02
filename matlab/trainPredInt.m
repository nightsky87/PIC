function trainPredInt(N)
    % Reset the random number generator for repeatability
    rng(0);
    
    % Define the number of training patches
    numPatches = 100000;

    % Collect traning patches from the Kodak database
    P = collectPatches(numPatches);
    
    % Reduce the patches to their respective sums
    S = reducePatch(P);
    
    % Initialize the interpolants by clustering
    I = initInterpolants(P,S,N);
    
    % Refine the interpolants iteratively
    I = refineInterpolants(P,S,I,100);
    
    % Quantize the interpolants
    for i = 1:3
        for j = 1:N
            I{i,j} = round(64 * I{i,j});
        end
    end
    
    % Generate the C lookup table
    generateLookup(I);
end

function P = collectPatches(numPatches)
    % Divide the number of patches evenly between the 24 Kodak images
    npImage = ones(1,24) * round(numPatches / 24);
    npImage(end) = numPatches - sum(npImage(1:23));

    % Allocate memory for all training patches
    P = zeros(36,numPatches,'single');

    % Process each Kodak image
    a = 1;
    for i = 1:24
        % Load the image to memory
        X = imread(sprintf('../images/Kodak/kodim%02d.bmp',i));

        % Map to the YCbCr and drop the chroma components
        X = rgb2ycbcr(X);
        X = double(X(:,:,1));

        % Collect training patches
        P(:,a:a+npImage(i)-1) = single(im2colrand(X,[6 6],npImage(i)));
        a = a + npImage(i);
    end
end

function S = reducePatch(P)
    % Define a summing matrix for 2x2 sub-patches
    a = 1;
    M = zeros(9,36);
    for j = 1:2:6
        for i = 1:2:6
            T = zeros(6);
            T(i:i+1,j:j+1) = 1;
            M(a,:) = T(:);
            a = a + 1;
        end
    end
    
    % Remove the top-left average
    M = M(2:end,:);

    % Calculate the sums
    S = M * P;
end

function I = initInterpolants(P,S,N)
    % Find the central average
    A = S(4,:) / 4;
    
    % Copy the neighborhood sums for each position
    S1 = S(3:8,:);
    S2 = S([1 2 4 5 7 8],:);
    S3 = S([4 5 7 8],:);
    
    % Copy the target pixels for each position
    P1 = P(21,:);
    P2 = P(16,:);
    P3 = P(22,:);
    
    % Find the pixel differences from the central average
    D = [P1; P2; P3] - repmat(A,[3 1]);
    
    % Generate an random set of cluster centers
    C = randn(3,N);
    C = C ./ repmat(sqrt(sum(C .^ 2)),[3 1]);
    
    % Cluster the differences
    for n = 1:10
        % Find the nearest cluster
        [proj,clust] = max(C' * D,[],1);
        
        % Update each cluster center
        regen = false(1,N);
        for i = 1:N
            % Locate all relevant data
            mask = (clust == i);
            
            % Regenerate the cluster center if unused
            if nnz(mask) == 0
                regen(i) = true;
                continue;
            end
            
            % Update the cluster center
            c = D(:,mask) * proj(mask)';
            C(:,i) = c / sqrt(c' * c);
        end
        
        % Regenerate unused clusters
        T = randn(3,nnz(regen));
        T = T ./ repmat(sqrt(sum(T .^ 2)),[3 1]);
        C(:,regen) = T;
    end
    
    % Determine the final cluster membership
    [~,clust] = max(C' * D,[],1);
    
    % Generate the initial interpolants using the final clusters
    I = cell(3,N);
    for i = 1:N
        % Locate all relevant data
        mask = (clust == i);
        
        % Generate the interpolants using least-squares regression
        I{1,i} = P1(:,mask) / S1(:,mask);
        I{2,i} = P2(:,mask) / S2(:,mask);
        I{3,i} = P3(:,mask) / S3(:,mask);
    end
end

function I = refineInterpolants(P,S,I,numIter)
    % Determine the number of interpolants
    N = size(I,2);
    
    % Copy the neighborhood sums for each position
    S1 = S(3:8,:);
    S2 = S([1 2 4 5 7 8],:);
    S3 = S([4 5 7 8],:);
    
    % Copy the target pixels for each position
    P1 = repmat(P(21,:),[N 1]);
    P2 = repmat(P(16,:),[N 1]);
    P3 = repmat(P(22,:),[N 1]);

    % Perform refinement iteratively
    curErr = Inf;
    for n = 1:numIter
        % Stack the interpolants
        I1 = cell2mat(I(1,:)');
        I2 = cell2mat(I(2,:)');
        I3 = cell2mat(I(3,:)');
    
        % Copy the previous error
        prevErr = curErr;
        
        % Find the best interpolant for each patch
        E1 = (P1 - I1 * S1) .^ 2;
        E2 = (P2 - I2 * S2) .^ 2;
        E3 = (P3 - I3 * S3) .^ 2;
        [E,int] = min(E1 + E2 + E3,[],1);
        
        % Check for termination
        curErr = sum(E);
        if curErr > prevErr
            I = Iprev;
            break;
        end
        fprintf('Iteration %d: %.0f\n',n,prevErr);
        
        % Update each interpolant
        Iprev = I;
        for i = 1:N
            % Locate all relevant data
            mask = (int == i);

            % Generate the interpolants using least-squares regression
            I{1,i} = P1(1,mask) / S1(:,mask);
            I{2,i} = P2(1,mask) / S2(:,mask);
            I{3,i} = P3(1,mask) / S3(:,mask);
        end
    end
end

function generateLookup(I)
    % Determine the number of interpolants
    N = size(I,2);

    % Reorder the coefficients for row-major formating
    for i = 1:N
        I{1,i} = I{1,i}([1 4 2 5 3 6]);
        I{2,i} = I{2,i}([1 3 5 2 4 6]);
        I{3,i} = I{3,i}([1 3 2 4]);
    end
    
    % Create the lookup header
    fh = fopen('../source/PICCom/PICComPredLookup.h','wt');
    
    fprintf(fh,'#ifndef PIC_COM_PRED_LOOKUP_H\n');
    fprintf(fh,'#define PIC_COM_PRED_LOOKUP_H\n\n');
    fprintf(fh,'#include "ComDef.h"\n\n');
    
    fprintf(fh,'static const s16 pred_w1[%d] = {\\\n',6*N);
    for i = 1:N-1
        fprintf(fh,'\t% 3d, % 3d, % 3d, % 3d, % 3d, % 3d,\\\n',I{1,i});
    end
    fprintf(fh,'\t% 3d, % 3d, % 3d, % 3d, % 3d, % 3d };\n\n',I{1,N});
    
    fprintf(fh,'static const s16 pred_w2[%d] = {\\\n',6*N);
    for i = 1:N-1
        fprintf(fh,'\t% 3d, % 3d, % 3d, % 3d, % 3d, % 3d,\\\n',I{2,i});
    end
    fprintf(fh,'\t% 3d, % 3d, % 3d, % 3d, % 3d, % 3d };\n\n',I{2,N});

    fprintf(fh,'static const s16 pred_w3[%d] = {\\\n',4*N);
    for i = 1:N-1
        fprintf(fh,'\t% 3d, % 3d, % 3d, % 3d,\\\n',I{3,i});
    end
    fprintf(fh,'\t% 3d, % 3d, % 3d, % 3d };\n\n',I{3,N});

    fprintf(fh,'#endif\n');
    
    fclose(fh);
end