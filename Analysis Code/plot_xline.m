
function plot_xline(line_x)
for l = 1:length(line_x)
handle = xline (line_x(l),'--',num2str(line_x(l)));
handle.HandleVisibility = "Off";
end
end

