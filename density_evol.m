function [its, mi] = density_evol( channel, dv, dc, NIter, K)
    %   Calculating for codes (dv, dc) in bi-awgn channel with (sigma^2) 
    %   the MI between input and output
    %
    %   its - number of iterations
    %   mi - mutual inf.
    %
    %
    
    [PChannel,QAwgn,bounds]= biAwgn2Dmc(channel.var,channel.M); %PChannel - a* - channel transition probs 2 x M
    [jPChannel,ind] = channelSort(double(PChannel));
    [QChannel,MIinit] = quantBiDmc(jPChannel,channel.Z); %opt. Quantizer (Z x M)
    [~,reInd] = sort(ind);
    jPChannel = jPChannel(:,reInd);
    QChannel = QChannel(:,reInd);
    rInit = PChannel*QChannel';

    mi = 0;
    its = NIter;

    rIn = rInit;
    tIn = rInit;

    for i = 1:NIter
        for j = 1:dc - 2
            [tIn] = CnDecomp(tIn, rIn, K);
        end
        lIn = tIn;
        tIn = rInit;

        for j = 1:dv - 2
            [tIn,tTild] = VnDecomp(tIn, lIn, K);
        end

        [tOut, tTild,mi] = VnDecomp(tIn, lIn, K);

        if (mi) >= 1
            its = i;
            return;
        end

        rIn = tOut;
        tIn = tOut;
    end
end

function [tOut] = CnDecomp( tIn, rIn, K)
    %       Decomposed Check
    %               Node
    % V0: tIn --> O --> tOut->tIn --> O --> tOut->tIn -->... --> Lij
    %             ^                   ^
    %    V1: rIn  |          V2: rIn  |
    %               (vn to cn msgs)

    %   rIn - Input pr distr for quantized VN to CN msg at current iteration
    %   (for initial iteration: transition probs for DMC channel produced from BI-AWGN ch)
    %   tIn - Input pr distr from output of prev decomp. node (2-input 1-output, tIn=rIn for
    %   the first 2-input node of decomposed CN)
    %   Output: tOut - output pr distr for CN to VN msgs
    %          
    %

    m = [[1/2 0 0 1/2];
        [0 1/2 1/2 0]];

    tTild = m*kron(tIn, rIn);
    [tTild,ind] = channelSort(tTild);

    Q = quantBiDmc(tTild, K); %quantizing decomposition-part node into K levels
    [~,reInd] = sort(ind);
    Q = Q(:,reInd);
    %output prob distr

    tTild = tTild(:,reInd);
    tOut = tTild*Q';
end

function [tOut, tTild,mi] = VnDecomp(tIn,lIn, K)
    tTild = [[kron(tIn(1,:),lIn(1,:))];[kron(tIn(2,:),lIn(2,:))]];
    [tTild,ind] = channelSort(tTild);
     [Q,mi] = quantBiDmc(tTild,K);
    [~,reInd] = sort(ind);
    Q = Q(:,reInd);
    tTild = tTild(:,reInd);

    tOut = tTild*Q';
end

