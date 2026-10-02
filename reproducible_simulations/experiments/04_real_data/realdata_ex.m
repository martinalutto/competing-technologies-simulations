years = [2006 2007 2008 2009];

blackberry = [7.0 9.6 16.6 20.8];
windows    = [14.0 12.2 11.8 8.8];
ios        = [0 2.7 8.2 15.1];

figure;
hold on;

plot(years, blackberry, '-o', 'LineWidth', 2.3, ...
    'MarkerSize', 10, 'DisplayName', 'BlackBerry');

plot(years, windows, '--s', 'LineWidth', 2.3, ...
    'MarkerSize', 10, 'DisplayName', 'Windows Mobile');

plot(years, ios, '-.^', 'LineWidth', 2.3, ...
    'MarkerSize', 10, 'DisplayName', 'iOS');

xlabel('Year');
ylabel('Global market share');

xticks(years);
xlim([2005.8 2009.2]);
ylim([0 23]);

legend('Location','best','Box','off');
box on;
grid on;

set(gca,'FontSize',16);

exportgraphics(gcf,'smartphone_market_share.pdf', ...
    'ContentType','vector');