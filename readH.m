function code = readH(path)
    %   read .alist    
    %   Note: For irregular LDPC the node degrees (used later in DE for 
    %   calculating threshold)  are ceiled mean values
    %
    
    [h, q] = alist2sparse(path);
    [code.H, code.G] = ldpc_h2g(h,q);
    code.H = full(code.H);
    code.H = code.H';
    [code.VN,code.CN] = size(code.H);
    code.K = size(code.G,1);%code.VN - code.CN;
    code.R = size(code.G,1)/size(code.H,1);

    code.dv = ceil(mean(sum(code.H,2)));
    code.dc = ceil(mean(sum(code.H,1)));
    code.CN2VN = zeros(code.VN,code.dv);
    
    for i = 1:code.VN
        code.CN2VN(i,[1:length(find(code.H(i,:)))]) =find(code.H(i,:));
    end
    for i = 1:code.CN
        code.VN2CN(i,[1:length(find(code.H(:,i)))]) = find(code.H(:,i));
    end

end