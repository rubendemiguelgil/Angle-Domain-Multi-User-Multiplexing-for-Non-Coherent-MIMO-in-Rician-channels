clear, clc, close all;

folder = "Results\subfolder_name_here\";

load(folder + "params")

load(folder + "results_nch_eep")
NC_SER_eep  = results_nch_eep.SER_total_mtx;
NC_SINR_eep = results_nch_eep.SINR_total_mtx;

load(folder + "results_nch_freq")
NC_SER_freq  = results_nch_freq.SER_total_mtx;
NC_SINR_freq = results_nch_freq.SINR_total_mtx;

load(folder + "results_mmse_imperf.mat")
CH_SER_mmse_imperf  = results_ch_imperfect_mmse.SER_total_mtx;
CH_SINR_mmse_imperf = results_ch_imperfect_mmse.SINR_total_mtx;

SNR_sweep = params.SNR_sweep;
M = params.M;
K = params.K;
U = params.N_users;
W = params.width;
sigma = sqrt(10.^(-SNR_sweep/10));

SINR_analytical = ((M * K + W)/(1+K)) ./ (((W/M*(U^2-U) + 2*(K*U -K))/(1+K)^2) + ((W*U - W + 2 * sigma.^2 * W/M * U + 2* K * sigma.^2)/(1+K)) + W/M * sigma.^4);
SINR_dB_analytical = 10*log10(SINR_analytical);

close all
figure(1)
    hold on, grid on
    xlabel('\textbf{SNR (dB)}','Interpreter','latex');
    ylabel('\textbf{SER}','Interpreter','latex');
    plot(SNR_sweep, CH_SER_mmse_imperf, '-.r^', 'DisplayName','\textbf{MMSE with Channel errors}')
    plot(SNR_sweep, NC_SER_freq, '-ksquare', 'DisplayName','\textbf{Proposed non-coherent}')
    plot(SNR_sweep, NC_SER_eep, '-.x', 'DisplayName','\textbf{Non-coherent EEP from [11]}')
    yscale log
    xlim([-12, 0])
    legend('Interpreter','latex', 'location','southwest', 'FontWeight','bold');
    set(gca,'TickLabelInterpreter','latex')
    set(findall(gcf,'-property','FontSize'),'FontSize',12.5)
    box on;

figure(2)
    hold on, grid on
    xlabel('\textbf{SNR (dB)}','Interpreter','latex');
    ylabel('\textbf{SINR (dB)}','Interpreter','latex');
    plot(SNR_sweep, CH_SINR_mmse_imperf, '-.r^', 'DisplayName','\textbf{MMSE with Channel errors}')
    plot(SNR_sweep, NC_SINR_freq, '-ksquare', 'DisplayName','\textbf{Proposed non-coherent}')
    plot(SNR_sweep, NC_SINR_eep, '-.x', 'DisplayName','\textbf{Non-coherent EEP from [11]}')
    plot(SNR_sweep, SINR_dB_analytical, '-b+', 'DisplayName','\textbf{Analytical SINR from (18)}')
    legend('Interpreter','latex', 'location','southeast');
    set(gca,'TickLabelInterpreter','latex')
    set(findall(gcf,'-property','FontSize'),'FontSize',12.5)
    box on;

    %% Plot antenna scaling
    clear, clf;
    results_folder = "Results\antenna_scaling\";
    load(results_folder + "workspace")

    SNR_sweep = params.SNR_sweep;

    figure(1)
    hold on, grid on
    xlabel('\textbf{SNR (dB)}','Interpreter','latex');
    ylabel('\textbf{SER}','Interpreter','latex');
    for i = 1:length(N_ant_array)

        plot(SNR_sweep, SER_ant_scaling_nch_freq(i,:), '-ksquare')
        plot(SNR_sweep, SER_ant_scaling_imperfect_mmse(i,:), '-.r^')
  
    end
    legend_nch = ['\textbf{Non-coherent scheme}'];
    legend_mmse = ['\textbf{MMSE with channel errors}'];
    legend_labels = {'\textbf{Non-coherent scheme}', '\textbf{MMSE with channel errors}'};
    l= legend(legend_nch, legend_mmse, 'Interpreter','latex', 'location','southwest');

    labels = {'\textbf{M=50}' '\textbf{M=100}' '\textbf{M=200}' '\textbf{M=400}'};
    x_nch = [-2.8 -6.2 -7.6 -12];
    y_nch = [1e-1 10^(-2.2) 10^(-3) 1e-4];
    x_ch = [-1 -2 -2.7 -5.7];
    y_ch = [10^(-1.5) 10^(-3.3) 10^(-5) 10^(-4.3)];
    text(x_nch, y_nch, labels, 'Color', 'black', 'Interpreter','latex')
    text(x_ch, y_ch, labels, 'Color', 'red', 'Interpreter','latex')
    box on;
    yscale log;
    set(gca,'TickLabelInterpreter','latex')
    set(findall(gcf,'-property','FontSize'),'FontSize',12.5)
    xlim([-20 3])
    ylim([3*10^(-7) 1])

