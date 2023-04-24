function fig = animate_system_seq(xhist,thist,VIDEO,video_title,t_inc)

if VIDEO, video=video_writer(video_title, t_inc); end
fig = plot_system_seq(xhist);
set(gcf,'Visible','on')

t_range = 1:10:size(thist);

mav_view = uav_viewer();
h = animatedline('LineWidth', 1.2, 'Color', 'r');

for j = t_range
        path_index = ~isnan(xhist(:,9));
        if path_index(j)
            mav_view.update(thist(j), xhist(j,:));
            hold on
            addpoints(h, xhist(j, 2), xhist(j, 1),-xhist(j, 3))
            hold on
        end
%     pause(0.002)
    if VIDEO, video.update(thist(j));  end
end

if VIDEO, video.close(); end
