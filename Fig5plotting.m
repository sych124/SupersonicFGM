SiNSUS_MTS = [-412	-273	-198	-149	-117	-81.2	-52	-24.9	-10.7	-6.3	-4.4	-2.8	-1.36	-1.3	-1.1	-1	-0.9	-0.8	-0.6	-0.38	-0.3	-0.18	-0.16];
SiNSUS_ROM= [-394	-272	-200	-152	-120	-83.6	-53.2	-25	-10.5	-6.1	-4.3	-2.72	-1.68	-1.3	-1.1	-1	-0.9	-0.8	-0.6	-0.38	-0.3	-0.18	-0.16];
k = [ 0.1	 0.2	 0.3	 0.4	 0.5	 0.75	 1	 2	 3	 4	 5	 7.5	 10	 11	 12	 13	 14	 15	 20	 25	 30	 35	 40];

figure();

    hold on;
    plot(k, SiNSUS_MTS,'-k','LineWidth',2);


    hold on;
    plot(k, SiNSUS_ROM,'-r','LineWidth',2);

    set(gcf, 'Position', [10, 10, 900+160, 900+180]);
    set(gca, 'GridAlpha',0.4,'FontSize',18);
    %set(gca, 'XScale','log');
    %set(gca, 'YScale','log');
    %title('Frequency derivation of AlO/Al');
    ylabel('S_{T_{cr}}','FontSize',20);
    xlabel('Volum fraction (k)','FontSize',20);
    xlim([0 40]);
    grid on;
    
legend('MTS','ROM','FontSize',26,'Location','northwest');