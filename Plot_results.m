clear, clc, close all;

folder = "Results\14-Apr-2025 12_43_21/";

load(folder + "params")

load(folder + "results_nch_eep")
NC_SER_eep  = results_nch_eep.SER_total_mtx;
NC_SINR_eep = results_nch_eep.SINR_total_mtx;

load(folder + "results_nch_freq")
NC_SER_freq  = results_nch_freq.SER_total_mtx;
NC_SINR_freq = results_nch_freq.SINR_total_mtx;

load(folder + "results_mrc_perf.mat")
CH_SER_mrc_perf  = results_ch_perfect_mrc.SER_total_mtx;
CH_SINR_mrc_perf = results_ch_perfect_mrc.SINR_total_mtx;

load(folder + "results_zf_perf.mat")
CH_SER_zf_perf  = results_ch_perfect_zf.SER_total_mtx;
CH_SINR_zf_perf = results_ch_perfect_zf.SINR_total_mtx;

load(folder + "results_mmse_perf.mat")
CH_SER_mmse_perf  = results_ch_perfect_mmse.SER_total_mtx;
CH_SINR_mmse_perf = results_ch_perfect_mmse.SINR_total_mtx;

load(folder + "results_mrc_imperf.mat")
CH_SER_mrc_imperf  = results_ch_imperfect_mrc.SER_total_mtx;
CH_SINR_mrc_imperf = results_ch_imperfect_mrc.SINR_total_mtx;

load(folder + "results_zf_imperf.mat")
CH_SER_zf_imperf  = results_ch_imperfect_zf.SER_total_mtx;
CH_SINR_zf_imperf = results_ch_imperfect_zf.SINR_total_mtx;

load(folder + "results_mmse_imperf.mat")
CH_SER_mmse_imperf  = results_ch_imperfect_mmse.SER_total_mtx;
CH_SINR_mmse_imperf = results_ch_imperfect_mmse.SINR_total_mtx;

SNR_sweep = params.SNR_sweep;

close all
figure(1)
    hold on, grid on

    xlabel('\textbf{SNR (dB)}','Interpreter','latex');
    ylabel('\textbf{SER}','Interpreter','latex');
    % plot(SNR_sweep, CH_BER_mrc_perf,'--^', 'DisplayName','Coherenet MRC: Perfect Channel Knowledge')
    % plot(SNR_sweep, CH_SER_mrc_imperf, '--o', 'DisplayName','Coherenet MRC: Channel error var = N0')
    % plot(SNR_sweep, CH_SER_zf_perf,'--^', 'DisplayName','Coherenet ZF: Perfect Channel Knowledge')
    % plot(SNR_sweep, CH_SER_zf_imperf, '--o', 'DisplayName','Coherenet ZF: Channel error var = N0')
    % plot(SNR_sweep, CH_SER_mmse_perf,'--^', 'DisplayName','MMSE with Perfect Channel Knowledge')
    plot(SNR_sweep, CH_SER_mmse_imperf, '-.r^', 'DisplayName','\textbf{MMSE with Channel errors}')
    % plot(SNR_sweep, NC_SER_time, '-.square', 'DisplayName','Non-coherent time')
    plot(SNR_sweep, NC_SER_freq, '-ksquare', 'DisplayName','\textbf{Proposed non-coherent}')
    % plot(SNR_sweep, NC_SER_leg, '-.^', 'DisplayName','\textbf{Non-coherent single user from [7]}')
    plot(SNR_sweep, NC_SER_eep, '-.x', 'DisplayName','\textbf{Non-coherent EEP from [11]}')
    yscale log
    xlim([-12, 1])
    legend('Interpreter','latex', 'location','southwest', 'FontWeight','bold');
    set(gca,'TickLabelInterpreter','latex')
    set(findall(gcf,'-property','FontSize'),'FontSize',12.5)
    box on;

figure(2)
    hold on, grid on
    % title('Signal to Interference and Noise Ratio')
    xlabel('\textbf{SNR (dB)}','Interpreter','latex');
    ylabel('\textbf{SINR (dB)}','Interpreter','latex');
    % plot(SNR_sweep, CH_SINR_mrc_perf, '--^', 'DisplayName','Coherenet MRC: Perfect Channel Knowledge')
    % plot(SNR_sweep, CH_SINR_mrc_imperf, '--o', 'DisplayName','Coherenet MRC: Channel error var = N0')
    % plot(SNR_sweep, CH_SINR_zf_perf, '--^', 'DisplayName','Coherenet ZF: Perfect Channel Knowledge')
    % plot(SNR_sweep, CH_SINR_zf_imperf, '--o', 'DisplayName','Coherenet ZF: Channel error var = N0')
    % plot(SNR_sweep, -CH_SINR_mmse_perf, '--^', 'DisplayName','MMSE with Perfect Channel Knowledge')
    plot(SNR_sweep, CH_SINR_mmse_imperf, '-.r^', 'DisplayName','\textbf{MMSE with Channel errors}')
    % plot(SNR_sweep, NC_SINR_time, '-.square', 'DisplayName','Non-coherent time')
    plot(SNR_sweep, NC_SINR_freq, '-ksquare', 'DisplayName','\textbf{Proposed non-coherent}')
    % plot(SNR_sweep, NC_SINR_leg, '-^', 'DisplayName','\textbf{Non-coherent single user from [7]}')
    plot(SNR_sweep, NC_SINR_eep, '-.x', 'DisplayName','\textbf{Non-coherent EEP from [11]}')
    legend('Interpreter','latex', 'location','northwest');
    set(gca,'TickLabelInterpreter','latex')
    set(findall(gcf,'-property','FontSize'),'FontSize',12.5)
    box on;

    %% Plot antenna scaling
    clear, close all;
    results_folder = "Results\antenna_scaling/";
    load(results_folder + "workspace")

    SNR_sweep = params.SNR_sweep;

    figure(1)
    hold on, grid on
    xlabel('\textbf{SNR (dB)}','Interpreter','latex');
    ylabel('\textbf{SER}','Interpreter','latex');
    for i = 1:length(N_ant_array)-1

        plot(SNR_sweep, SER_ant_scaling_nch_freq(i,:), '-ksquare')
        % plot(SNR_sweep, SER_ant_scaling_perfect_mmse(i,:), '-.bx')
        plot(SNR_sweep, SER_ant_scaling_imperfect_mmse(i,:), '-.r^')
  
    end
    legend_nch = ['\textbf{Non-coherent scheme}'];
    legend_mmse = ['MMSE with channel errors'];
    legend_labels = {'\textbf{Non-coherent scheme}', '\textbf{MMSE with channel errors}'};
    l= legend(legend_nch, legend_mmse, 'Interpreter','latex', 'location','southwest');

    labels = {'\textbf{M=50}' '\textbf{M=100}' '\textbf{M=200}' '\textbf{M=400}'};
    x_nch = [-2.8 -5.5 -8.2 -12];
    y_nch = [1e-1 1e-2 10^(-3.2) 1e-4];
    x_ch = [-1.15 -2.6 -5 -8];
    y_ch = [1e-4 1e-5 10^(-4.3) 10^(-4.6)];
    text(x_nch, y_nch, labels, 'Color', 'black', 'Interpreter','latex')
    text(x_ch, y_ch, labels, 'Color', 'red', 'Interpreter','latex')
    box on;
    yscale log;
    set(gca,'TickLabelInterpreter','latex')
    set(findall(gcf,'-property','FontSize'),'FontSize',12.5)
    xlim([-20 2])

    figure(2)
    hold on, grid on
    xlabel('SNR (dB)','Interpreter','latex');
    ylabel('SINR (dB)','Interpreter','latex');
    for i = 1:length(N_ant_array)-1
        plot(SNR_sweep, SINR_ant_scaling_nch_freq(i, :), '-.ksquare')
        plot(SNR_sweep, SINR_ant_scaling_perfect_mmse(i, :), '-.bx')
        plot(SNR_sweep, SINR_ant_scaling_imperfect_mmse(i, :), '-.r^')
    end
    l= legend(legend_nch, legend_mmse, 'Interpreter','latex', 'location','southwest');

    set(gca,'TickLabelInterpreter','latex')
    set(findall(gcf,'-property','FontSize'),'FontSize',12.5)
    box on;