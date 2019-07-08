function simulate(M, thrf, ldpc_filename, iter, snr_array, q_uniform, q_nonuniform, fe)
    %   Estimation of unquantized, uniformly and non-uniformly
    %   quantized min-sum LDPC decoders performance.
    %
    %   M = 2 for bpsk
    %   thrf - set 0 to calculte new decoding treshold and quantization bounds, 
    %   or put the name of the file with saved decoding threshold
    %   ldpc_filename - alist file
    %   snr_array - snr values to simulate
    %   q_uniform - number of bits for uniform Q
    %   q_nonuniform - bits for non-uniform Q (< q_uniform)
    %   fe - number of frames with error (the same for unquantized, unif. Q and non
    %   unif. Q)
    %
    %   simulate(2, 0,'codes/LDPC_5G_K120_R075.alist', 20,[1 1.5 2 2.5 3 3.5 4 4.5 5], 6, 4, 10)
    %
    %
    
    addpath('quantDMC','codes','binary_soft_decoder');
    %   q_uniform = 6;
    %   q_nonuniform = 4;

    init_ldpc = @(x) decode_soft(-2, x);
    free_ldpc = @(x) decode_soft(-1, x);

    code = readH(ldpc_filename);

    if thrf~=0
        thr = load(thrf);
        thr = thr.thr;
    else
        thr = threshold_sigma([2^q_uniform-1 2^q_nonuniform-1 2^q_nonuniform-1], code, iter);
    end

    disp(sprintf('\nNoise decoding threshold sigma^2 = %f\n', thr.channel.var));
    % initialize LDPC
    [h, q] = alist2sparse(ldpc_filename);
    [H, G] = ldpc_h2g(h,q);
    [K, N] = size(G);
    R = K/(N);

    [ldpc, ~, ~] = init_ldpc(ldpc_filename);
    Es = 1;

    in_ber = zeros(1, length(snr_array));
    ber_unq = zeros(1, length(snr_array));
    ber_uniq= zeros(1, length(snr_array));
    ber_q = zeros(1, length(snr_array));

    fer_unq = zeros(1, length(snr_array));
    fer_q = zeros(1, length(snr_array));
    fer_uniq = zeros(1, length(snr_array));

    hDemod = comm.BPSKDemodulator('DecisionMethod', 'Log-likelihood ratio');

    for ii = 1:length(snr_array)
        %   NEW SNR
        snr = snr_array(ii);

        noise_var = 0.5*10^(-(snr)/10)*(1/R);
        sigma = sqrt(noise_var);

        %   Uniform Q bounds
        [~,~,B] = biAwgn2Dmc(sigma^2,2^q_uniform - 1);

        hDemod.Variance = sigma^2;
        hDemod.release();

        disp(sprintf('\n\n========== SNR = %f, EbN0 = %f, sigma = %1.2e ===========', snr,  snr, sigma/sqrt(Es)));

        tests = 0;
        in_errors = 0;

        errors_unq = 0;
        errors_uniq = 0;
        errors_q = 0;
        wrong_dec_unq = 0;
        wrong_dec_uniq = 0;
        wrong_dec_q=0;

        while  wrong_dec_unq < fe ||  wrong_dec_uniq < fe || wrong_dec_q < fe%
            tests = tests + 1;
            iwd = zeros(1,K); %randi(q, 1, K) - 1 %%zeros(1,K); %

            cwd = zeros(1, N); %ldpc_encode(iwd, G, 2);%
            modulated = 1-cwd.*2;

            rx = modulated + sigma*(randn(1,N));

            in_llrs = hDemod(rx');
            
            %   Clipping the values to be in [-2, 2]
            in_llrs_unq=[];
            for i = 1:N
                in_llrs_unq=[in_llrs_unq sign(rx(i))*(min([abs(rx(i)) 2]))];
            end

            est_cwd = 1*(in_llrs < 0);

            in_error_indexes = find(cwd ~= est_cwd');
            in_errors = in_errors + length(in_error_indexes);

            scale_array =ones(1, N);
            offset_array =0*ones(1,N);

            %   Uniform Quantization
            [~,in_llrs_uniq] = quantiz(in_llrs_unq, B(2:end-1), -(2^(q_uniform-1)-1):2^(q_uniform-1)-1);

            %   Non-Uniform Quantization
            [~,in_llrs_q] = quantiz(in_llrs_unq, thr.QBounds, -(2^(q_nonuniform-1)-1):2^(q_nonuniform-1)-1);
            
            %   decode_soft(0, ...) - unquantized min-sum decoder
            %   decode_soft(Q, ...) - min-sum decoder quantized to Q bits
            %   Unquantized:
            [result, number_of_iter, ecwd_unq, out_llrs] = decode_soft(0, ldpc,  in_llrs,  iter, scale_array, offset_array);
           
            %   Uniform Quantization:
            [result, number_of_iter, ecwd_uniq, out_llrs] = decode_soft(q_uniform, ldpc, in_llrs_uniq, iter,  scale_array, offset_array);
           
            %   Non-Uniform Quantization:
            [result, number_of_iter, ecwd_q, out_llrs] = decode_soft(q_nonuniform, ldpc, in_llrs_q,  iter,  scale_array, offset_array);

            eiwd_unq = ecwd_unq(N-K+1:end);
            eiwd_uniq = ecwd_uniq(N-K+1:end); % information bits are stored in the last positions of cwd
            eiwd_q = ecwd_q(N-K+1:end);

            error_indexes_unq = find(iwd ~= eiwd_unq);%find(cwd ~= ecwd_unq);%
            error_indexes_uniq = find(iwd ~= eiwd_uniq);%find(cwd ~= ecwd_uniq);%
            error_indexes_q = find(iwd ~= eiwd_q);%find(cwd ~= ecwd_q);%

            if ~isempty(error_indexes_unq) || ~isempty(error_indexes_uniq) || ~isempty(error_indexes_q)
                if ~isempty(error_indexes_unq)
                    wrong_dec_unq = wrong_dec_unq + 1;
                    errors_unq = errors_unq + length(error_indexes_unq);
                    fer_unq(ii) = wrong_dec_unq/tests;
                    ber_unq(ii) = errors_unq/N/tests;
                end

                if ~isempty(error_indexes_uniq)
                    wrong_dec_uniq = wrong_dec_uniq + 1;
                    errors_uniq = errors_uniq + length(error_indexes_uniq);
                    fer_uniq(ii) = wrong_dec_uniq/tests;
                    ber_uniq(ii) = errors_uniq/N/tests;
                end

                if ~isempty(error_indexes_q)
                    wrong_dec_q = wrong_dec_q + 1;
                    errors_q = errors_q + length(error_indexes_q);
                    fer_q(ii) = wrong_dec_q/tests;
                    ber_q(ii) = errors_q/N/tests;
                end

                in_ber(ii) = in_errors/N/tests;
                disp(sprintf('\tin_ber = %f, FER: non-Q = %f, Q-%d = %f, Q-%d = %f, BER: non-Q = %f, Q-%d = %f, Q-%d = %f' , ...
                    in_ber(ii), fer_unq(ii), q_uniform, fer_uniq(ii), q_nonuniform, ...
                    fer_q(ii), ber_unq(ii), q_uniform, ber_uniq(ii), q_nonuniform, ber_q(ii)));
            end
        end

    end
    free_ldpc(ldpc);
end