function [results] = simulate_noncoherent_joint_const(params)
%SIMULATE_COHERENT function that simulates the non-coherent EEP constellation scheme for two
%users or for a single user. 

%% Parameters
N_users = params.N_users; 
M = params.M; % Number of Rx antennas (BS)
L = params.L; % Tx length (in bits)
bps = params.bps; % 2 bits/symbol in QPSK
L_sym = params.L_sym; % Tx length in syms
N_subcarriers = params.N_subcarriers; % Number of dft points
L_ofdm_syms = params.L_ofdm_syms; % Length in ofdm symbols
pwr = params.user_pwr;
width = params.width;
phase_dist = params.phase_dist;
K = params.K;
N_taps = params.N_taps;
angles = params.user_angles;
n_ch_uses = params.n_channel_uses;

%% Create bit strings
bits = round(rand(L, N_users));

%% Constellation modulation (QPSK)
syms = QPSK_modulation(bits); 
if N_users == 2
    syms(:,2) = syms(:,2) * exp(j* pi/4); 
    pwr = [1 1];
    joint_syms = sum(syms, 2);
end

%% Differential OFDM encoding/modulation/carrier alocation
switch params.diff_decoding_dimension
    case 'time'
        [ofdm_signal] = OFDM_diff_modulation_time(syms, N_subcarriers);
    case 'freq'
        [ofdm_signal] = OFDM_diff_modulation_freq(syms, N_subcarriers);
end
%% Rician channel
H = rician_channel(angles, N_subcarriers, M, N_taps, K, phase_dist);

%% SNR sweep loop
SNR_sweep = params.SNR_sweep;
SER_total_mtx = zeros(length(SNR_sweep), n_ch_uses);
BER_total_mtx = zeros(length(SNR_sweep), n_ch_uses);
SINR_total_mtx = zeros(length(SNR_sweep), n_ch_uses);

for SNR_idx = 1:length(SNR_sweep)
    for ch_use = 1:n_ch_uses
        
        SNR_dB = SNR_sweep(SNR_idx);
        N0 = (10.^(-SNR_dB/10)); 
        
        %% Transmission 
        y = tx_ofdm_signal(ofdm_signal, H, N0);
        
        %% Differential OFDM decoding/demodulation/carrier dealocation
        switch params.diff_decoding_dimension
            case 'time'
                rx_syms = OFDM_diff_demodulation_time(y); 
            case 'freq'
                rx_syms = OFDM_diff_demodulation_freq(y); 
        end
        
        rx_syms = rx_syms(1:L_sym, :);% Neglect zero padded symbols due to fixed N_subcarriers
        if isequal(params.diff_decoding_dimension,'freq')
            rx_syms_nm = rx_syms./mean(abs(rx_syms),1); 
        else
            rx_syms_nm = rx_syms;
        end

        if N_users == 1
            det_syms = QPSK_detector(rx_syms_nm); % Min distance QPSK detection 
            det_bits = QPSK_demodulator(det_syms); % Map symbols to bits
            
            % SER
            error_sym = (det_syms - syms);
            error_sym_flags = (error_sym~=0);
            SER_total = sum(error_sym_flags, 'all')/(L * N_users);
            
            % SINR (from EVM) 
            evm = sqrt(sum(abs(rx_syms_nm - syms).^2, 'all')/(L_sym*N_users));
            SINR_dB = -10*log10(evm);

        elseif N_users == 2
            [det_bits, det_joint_syms] = joint_const_detection(rx_syms_nm); % Joint constellation
            % Metrics
            % SER Joint
            error_sym = (transpose(det_joint_syms) - joint_syms);
            error_sym_flags = (abs(error_sym)>0.1);
            SER_user = sum(error_sym_flags, 1)/L;
            SER_total = sum(error_sym_flags, 'all')/(length(joint_syms));
    
            % SINR Joint (from EVM) 
            evm = sqrt(sum(abs(rx_syms_nm - joint_syms).^2, 'all')/(L_sym*N_users));
            SINR_dB = -10*log10(evm);

        else
            assert(N_users < 3, 'Number of users must be 1 or 2 for the legacy NCH.')
        end

            
        % BER
        error_bits = abs(bits - det_bits);
        BER_total = sum(error_bits, 'all')/(L * N_users);
            
        SER_total_mtx(SNR_idx, ch_use) = SER_total;
        BER_total_mtx(SNR_idx, ch_use) = BER_total;
        SINR_total_mtx(SNR_idx, ch_use) = SINR_dB;

        text = strcat("Ch use number: ", int2str(ch_use) ,"; SNR :", int2str(SNR_dB));
        disp(text)
    end
end

results.SER_total_mtx = mean(SER_total_mtx, 2);
results.BER_total_mtx = mean(BER_total_mtx, 2);
results.SINR_total_mtx = mean(SINR_total_mtx, 2);

end

