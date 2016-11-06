function trainSmallBasis(N)
    % Reset the random number generator for repeatability
    rng(0);
    
    % Define the number of training patches
    numPatches = 100000;

    % Collect traning patches from the Kodak database
    P = collectPatches(numPatches);
    
    % Cluster the patches
    B = trainBasis(P,N);
    
    % Quantize the basis vectors and reorder for row-major format
    for i = 1:N
        B{i} = round(64 * B{i}([1 3 2 4],:));
    end
    
    % Generate the C lookup table
    generateLookup(B);
end

function P = collectPatches(numPatches)
    % Divide the number of patches evenly between the 24 Kodak images
    npImage = ones(1,24) * round(numPatches / 24);
    npImage(end) = numPatches - sum(npImage(1:23));

    % Allocate memory for all training patches
    P = zeros(4,numPatches,'single');

    % Process each Kodak image
    a = 1;
    for i = 1:24
        % Load the image to memory
        X = imread(sprintf('../images/Kodak/kodim%02d.bmp',i));

        % Map to the YCbCr and drop the chroma components
        X = rgb2ycbcr(X);
        X = double(X(:,:,1));

        % Collect training patches
        P(:,a:a+npImage(i)-1) = single(im2colrand(X,[2 2],npImage(i)));
        a = a + npImage(i);
    end
end

function B = trainBasis(P,N)
    % Remove the DC component of the patches
    P = P - repmat(mean(P),[4 1]);

    % Generate an random set of cluster centers
    C = randn(4,N);
    C = C ./ repmat(sqrt(sum(C .^ 2)),[4 1]);
    
    % Cluster the patches
    for n = 1:10
        % Find the nearest cluster
        proj = C' * P;
        [~,clust] = max(abs(proj),[],1);
        
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
            c = P(:,mask) * proj(mask)';
            C(:,i) = c / sqrt(c' * c);
        end
        
        % Regenerate unused clusters
        T = randn(4,nnz(regen));
        T = T ./ repmat(sqrt(sum(T .^ 2)),[4 1]);
        C(:,regen) = T;
    end
    
    % Determine the final cluster membership
    [~,clust] = max(abs(C' * P),[],1);
    
    % Generate an orthogonal basis for each cluster
    B = cell(N,1);
    for i = 1:N
        % Locate all relevant data
        mask = (clust == i);
        
        % Generate the interpolants using least-squares regression
        [U,~,~] = svd(P(:,mask),'econ');
        U = U(:,1:3);
        
        % Enforce that the largest coefficient in each basis is positive
        for j = 1:3
            [~,ind] = max(abs(U(:,j)));
            U(:,j) = U(:,j) / sign(U(ind,j));
        end
        B{i} = U;
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

    % Create the lookup header
    fh = fopen('../source/PICCom/PICComBasisLookup.h','wt');
    
    fprintf(fh,'#ifndef PIC_COM_BASIS_LOOKUP_H\n');
    fprintf(fh,'#define PIC_COM_BASIS_LOOKUP_H\n\n');
    fprintf(fh,'#include "ComDef.h"\n\n');
    fprintf(fh,'#define NUM_BASIS %d\n\n',N);
    
    fprintf(fh,'static const s32 basis_weights[%d] = {\\\n',12*N);
    for i = 1:N-1
        fprintf(fh,'\t% 3d, % 3d, % 3d, % 3d, % 3d, % 3d, % 3d, % 3d, % 3d, % 3d, % 3d, % 3d,\\\n',I{i}(:));
    end
    fprintf(fh,'\t% 3d, % 3d, % 3d, % 3d, % 3d, % 3d, % 3d, % 3d, % 3d, % 3d, % 3d, % 3d};\n\n',I{N}(:));
    fprintf(fh,'#endif\n');
    
    fclose(fh);
end