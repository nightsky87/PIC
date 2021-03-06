function trainSmallBasis(N)
    % Reset the random number generator for repeatability
    rng(0);
    
    % Define the number of training patches
    numPatches = 100000;

    % Collect traning patches from the Kodak database
    P = collectPatches(numPatches);
    
    % Initialize the interpolants
    I = initInterpolants(P,N);
    
    % Iteratively refine the interpolants
    [I,B] = refineInterpolants(P,I,100);
    
    % Generate the C lookup table
    generateLookup(I,B);
end

function P = collectPatches(numPatches)
    % Divide the number of patches evenly between the 24 Kodak images
    npImage = ones(1,24) * round(numPatches / 24);
    npImage(end) = numPatches - sum(npImage(1:23));

    % Allocate memory for all training patches
    P = zeros(16,numPatches,'single');

    % Process each Kodak image
    a = 1;
    for i = 1:24
        % Load the image to memory
        X = imread(sprintf('../images/Kodak/kodim%02d.bmp',i));

        % Map to the YCbCr and drop the chroma components
        X = rgb2ycbcr(X);
        X = double(X(:,:,1));

        % Collect training patches
        P(:,a:a+npImage(i)-1) = single(im2colrand(X,[4 4],npImage(i)));
        a = a + npImage(i);
    end
end

function I = initInterpolants(P,N)
    % Generate the sums of the four sub-patches
    S(1,:) = sum(P([1 2 5 6],:));
    S(2,:) = sum(P([9 10 13 14],:));
    S(3,:) = sum(P([3 4 7 8],:));
    S(4,:) = sum(P([11 12 15 16],:));
    
    % Generate a subset of patches and remove the DC component of the subset
    P = P([1 5 2 6],:);
    Z = P - repmat(mean(P),[4 1]);

    % Generate an random set of cluster centers
    C = randn(4,N);
    C = C ./ repmat(sqrt(sum(C .^ 2)),[4 1]);
    
    % Cluster the patches
    for n = 1:10
        % Find the nearest cluster
        proj = C' * Z;
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
            c = Z(:,mask) * proj(mask)';
            C(:,i) = c / sqrt(c' * c);
        end
        
        % Regenerate unused clusters
        T = randn(4,nnz(regen));
        T = T ./ repmat(sqrt(sum(T .^ 2)),[4 1]);
        C(:,regen) = T;
    end
    
    % Determine the final cluster membership
    [~,clust] = max(abs(C' * Z),[],1);
    
    % Generate an interpolant for each cluster
    I = cell(N,1);
    for i = 1:N
        % Locate all relevant data
        mask = (clust == i);
        
        % Generate the interpolants using least-squares regression
        I{i} = P(:,mask) / S(:,mask);
    end
end

function [I,B] = refineInterpolants(P,I,numIter)
    % Define the quantization scale
    qScale = 4096;
    
    % Determine the number of interpolants
    N = size(I,1);
    
    % Generate the sums of the four sub-patches
    S(1,:) = sum(P([1 2 5 6],:));
    S(2,:) = sum(P([9 10 13 14],:));
    S(3,:) = sum(P([3 4 7 8],:));
    S(4,:) = sum(P([11 12 15 16],:));
    
    % Generate a subset of patches
    P = repmat(P([1 5 2 6],:),[N 1]);

    % Perform refinement iteratively
    curErr = Inf;
    for n = 1:numIter
        % Stack the interpolants
        Is = cell2mat(I);
    
        % Copy the previous error
        prevErr = curErr;
        
        % Find the best interpolant for each patch
        E = reshape(sum(reshape((P - Is * S) .^ 2,4,[])),N,[]);
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
            I{i} = P(1:4,mask) / S(:,mask);
        end
    end
    
    % Quantize the interpolants
    for i = 1:N
        I{i} = round(qScale * I{i});
    end
    
    % Stack the interpolants
    Is = cell2mat(I);
    
    % Find the best interpolant for each patch
    E = reshape(sum(reshape((P - Is * S / qScale) .^ 2,4,[])),N,[]);
    [~,int] = min(E,[],1);
    
    % Find the error basis for each interpolant
    B = cell(N,1);
    for i = 1:N
        % Locate all relevant data
        mask = (int == i);

        % Calculate the interpolation errors
        E = P(1:4,mask) - round(I{i} * S(:,mask) / qScale);
        [U,~,~] = svd(E,'econ');
        
        % Standardize the orientations of the basis vectors
        for j = 1:3
            [~,ind] = max(abs(U(:,j)));
            U(:,j) = U(:,j) / sign(U(ind,j));
        end
        B{i} = round(qScale * U(:,1:3));
    end
end

function generateLookup(I,B)
    % Determine the number of interpolants
    N = size(I,1);

    % Create the lookup header
    fh = fopen('../source/PICCom/PICComBasisLookup.h','wt');
    
    fprintf(fh,'#ifndef PIC_COM_BASIS_LOOKUP_H\n');
    fprintf(fh,'#define PIC_COM_BASIS_LOOKUP_H\n\n');
    fprintf(fh,'#include "ComDef.h"\n\n');
    fprintf(fh,'#define NUM_BASIS %d\n\n',N);
    
    fprintf(fh,'static const s32 interp_weights[%d] = {\\\n',16*N);
    for i = 1:N-1
        fprintf(fh,'\t% 4d, % 4d, % 4d, % 4d, % 4d, % 4d, % 4d, % 4d, % 4d, % 4d, % 4d, % 4d, % 4d, % 4d, % 4d, % 4d,\\\n',reshape(I{i}',[],1));
    end
    fprintf(fh,'\t% 4d, % 4d, % 4d, % 4d, % 4d, % 4d, % 4d, % 4d, % 4d, % 4d, % 4d, % 4d, % 4d, % 4d, % 4d, % 4d};\n\n',reshape(I{i}',[],1));

    fprintf(fh,'static const s32 basis_weights[%d] = {\\\n',12*N);
    for i = 1:N-1
        fprintf(fh,'\t% 4d, % 4d, % 4d, % 4d, % 4d, % 4d, % 4d, % 4d, % 4d, % 4d, % 4d, % 4d,\\\n',B{i}(:));
    end
    fprintf(fh,'\t% 4d, % 4d, % 4d, % 4d, % 4d, % 4d, % 4d, % 4d, % 4d, % 4d, % 4d, % 4d};\n\n',B{N}(:));
    fprintf(fh,'#endif\n');
    
    fclose(fh);
end