case_list = {'case01.mat', 'case02.mat', 'case03.mat', 'case04.mat'};
xrange = [0.15 0.3];

for i = 1 : 4
	subplot(2, 2, i);
	casename1 = case_list{i};
	load(casename1);
	plot(t_axis, d2 * 4, 'b', 'linewidth', 2); hold on;
	plot(t_axis, d, 'r--', 'linewidth', 2);
	xlim(xrange);
end

