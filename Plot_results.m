clear, clc, close all;

load("Results\CH_err_wspace.mat", 'BER_total_mtx', 'SINR_total_mtx', 'SNR_sweep')
CH_err_BER  = BER_total_mtx;
CH_err_SINR = SINR_total_mtx;

load("Results\CH_perfect_wspace.mat", 'BER_total_mtx', 'SINR_total_mtx', 'SNR_sweep')
CH_perf_BER  = BER_total_mtx;
CH_perf_SINR = SINR_total_mtx;

load("Results\NC_wspace.mat", 'BER_total_mtx', 'SINR_total_mtx')
NC_BER  = BER_total_mtx;
NC_SINR = SINR_total_mtx;

%%
figure(1)
    hold on, grid on
    title('Bit Error Rate')
    xlabel('SNR (dB)')
    ylabel('BER')
    plot(SNR_sweep, CH_perf_BER,'--^', 'DisplayName','Coherenet MRC: Perfect Channel Knowledge')
    plot(SNR_sweep, CH_err_BER, '--o', 'DisplayName','Coherenet MRC: Channel error var = N0')
    plot(SNR_sweep, NC_BER, '-.square', 'DisplayName','Non-coherent')
    yscale log
    legend()


figure(2)
    hold on, grid on
    title('Signal to Interference and Noise Ratio')
    xlabel('SNR (dB)')
    ylabel('SINR (dB)')
    plot(SNR_sweep, CH_perf_SINR, '--^', 'DisplayName','Coherenet MRC: Perfect Channel Knowledge')
    plot(SNR_sweep, CH_err_SINR, '--o', 'DisplayName','Coherenet MRC: Channel error var = N0')
    plot(SNR_sweep, NC_SINR, '-.square', 'DisplayName','Non-coherent')
    legend()