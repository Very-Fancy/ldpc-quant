function thr = threshold_sigma(quant, code, NIter)
    %   Decoding threshold  sigma (std) for code ensemble in awgn channel 
    %
    %   quant = [M Z K], M - uniform Q levels, Z - DMC Q levels, K - decoder Q
    %   levels ( Set Z = K )
    %   code: result of readH()
    %
    %   Output:
    %   thr.channel: M, Z - levels; sigma - decoding threshold (std); var -
    %   threshold noise variance (sigma^2); QBounds - Bounds for
    %   non-uniform Q
    %
    %
    
    QBounds=[];
    channel.M = quant(1,1);
    channel.Z = quant(1,2);
    thr.K = quant(1,3);
    
    sigmaMax = 10;
    sigmaMin = 0;
    d_sigma0 =100;

    while (sigmaMax - sigmaMin) < d_sigma0
        d_sigma0 = (sigmaMax - sigmaMin);
        channel.sigma = ((sigmaMax + sigmaMin)/2);
        channel.var = (channel.sigma^2);

        [ Iters, conv] = density_evol( channel, code.dv, code.dc, NIter, thr.K);

        if conv >= 1.0
            sigmaMin = channel.sigma;
        else
            sigmaMax = channel.sigma;
        end
    end

    thr.channel = channel;
    thr.conv = conv;
    thr.Niter = Iters;

    [PChannel,QAwgn,bounds] = biAwgn2Dmc(channel.var, channel.M);
    [QChannel,MIinit] = quantBiDmc(PChannel, channel.Z);
    
    for i = 2:channel.Z
        lvl = (find(QChannel(i,:)==1));%max(QChannel(i,:));
        QBounds = [QBounds bounds(1,min(lvl))];
    end
    thr.QBounds = QBounds;

end



