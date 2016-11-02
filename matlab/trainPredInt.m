function trainPredInt(N)
    % Reset the random number generator for repeatability
    rng(0);
    
    % Define the number of training patches
    numPatches = 100000;

    % Collect traning patches from the Kodak database
    P = collectPatches(numPatches);
    
    % Reduce the patches to their respective sums
    S = reducePatchVert(P);
    
    % Initialize the interpolants by clustering
    I = initInterpolants(P,S,N);
    
    % Refine the interpolants iteratively
    I = refineInterpolants(P,S,I,100);
    
    % Quantize the interpolants
    for i = 1:N
        I{i} = round(64 * I{i});
    end
    
    % Generate the C lookup table
    generateLookup(I);
end

function P = collectPatches(numPatches)
    % Divide the number of patches evenly between the 24 Kodak images
    npImage = ones(1,24) * round(numPatches / 24);
    npImage(end) = numPatches - sum(npImage(1:23));

    % Allocate memory for all training patches
    P = zeros(12,numPatches,'single');

    % Process each Kodak image
    a = 1;
    for i = 1:24
        % Load the image to memory
        X = imread(sprintf('../images/Kodak/kodim%02d.bmp',i));

        % Map to the YCbCr and drop the chroma components
        X = rgb2ycbcr(X);
        X = double(X(:,:,1));

        % Collect training patches
        P(:,a:a+npImage(i)-1) = single(im2colrand(X,[4 3],npImage(i)));
        a = a + npImage(i);
    end
end

function S = reducePatchVert(P)
    % Define the summing matrix for 2x1 sub-patches
    a = 1;
    M = kron(eye(6),ones(1,2));

    % Calculate the vertical sums
    S = M * P;
end

function I = initInterpolants(P,S,N)
    % Copy the target pixels
    P = P(6,:);
    
    % Find the difference of the pixel from the neighboring averages
    D = repmat(P,[6 1]) - S / 2;
    
    % Generate an random set of cluster centers
    C = randn(6,N);
    C = C ./ repmat(sqrt(sum(C .^ 2)),[6 1]);
    
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
        T = randn(6,nnz(regen));
        T = T ./ repmat(sqrt(sum(T .^ 2)),[6 1]);
        C(:,regen) = T;
    end
    
    % Determine the final cluster membership
    [~,clust] = max(C' * D,[],1);
    
    % Generate the initial interpolants using the final clusters
    I = cell(N,1);
    for i = 1:N
        % Locate all relevant data
        mask = (clust == i);
        
        % Generate the interpolants using least-squares regression
        I{i} = P(:,mask) / S(:,mask);
    end
end

function I = refineInterpolants(P,S,I,numIter)
    % Determine the number of interpolants
    N = size(I,1);
    
    % Copy the target pixels for each position
    P = repmat(P(6,:),[N 1]);

    % Perform refinement iteratively
    curErr = Inf;
    for n = 1:numIter
        % Stack the interpolants
        Is = cell2mat(I);
    
        % Copy the previous error
        prevErr = curErr;
        
        % Find the best interpolant for each patch
        E = (P - Is * S) .^ 2;
        [E,int] = min(E,[],1);
        
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
            I{i} = P(1,mask) / S(:,mask);
        end
    end
end

function generateLookup(I)
    % Determine the number of interpolants
    N = size(I,1);

    % Reorder the coefficients for row-major formating
    for i = 1:N
        I{i} = I{i}([1 3 5 2 4 6]);
    end
    
    % Create the lookup header
    fh = fopen('../source/PICCom/PICComPredLookup.h','wt');
    
    fprintf(fh,'#ifndef PIC_COM_PRED_LOOKUP_H\n');
    fprintf(fh,'#define PIC_COM_PRED_LOOKUP_H\n\n');
    fprintf(fh,'#include "ComDef.h"\n\n');
    
    fprintf(fh,'static const s16 pred_weights[%d] = {\\\n',6*N);
    for i = 1:N-1
        fprintf(fh,'\t% 3d, % 3d, % 3d, % 3d, % 3d, % 3d,\\\n',I{i});
    end
    fprintf(fh,'\t% 3d, % 3d, % 3d, % 3d, % 3d, % 3d };\n\n',I{N});
    fprintf(fh,'#endif\n');
    
    fclose(fh);
end