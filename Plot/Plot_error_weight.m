function Plot_error_weight(w_lms,w_rls,e_lms,e_rls) 

%% 画图
    %抽头对比图
    figure
    plot(w_lms)
    hold on
    plot(w_rls)
    legend("FFE-LMS","FFE-RLS")
    title("均衡器抽头")
    
    %训练的对比图
    figure
    plot(abs(e_lms))
    hold on
    plot(abs(e_rls))
    legend("FFE-LMS","FFE-RLS")
    xlabel("迭代次数")
    ylabel("误差")

end
