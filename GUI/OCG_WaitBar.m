classdef OCG_WaitBar

    properties (SetAccess=private)

        Fig
        Txt
        maxIter
        Line1
        Line2
        tstart

    end

    methods

        function obj = OCG_WaitBar(maxIter)

            obj.tstart = tic;
            obj.maxIter = maxIter;

            % 界面
            WS = get(0,'ScreenSize');
            wt = 500;
            ht = 100;
            obj.Fig = figure('Position',[(WS(3)-wt)/2,(WS(4)-ht)/2,wt,ht],...
                'Name','GUI','NumberTitle','off','menu','none',...
                'Color','white','Resize','off');

            % 头像
            barH = 35;
            AxesIcon = axes(obj.Fig,'Units','pixels','Position',[1,1,ht,ht]);
            try
                icon = imread('今天不飞了.png');
            catch
                icon = ind2rgb(round(255*rescale(peaks(100))+1),[1-hot(128); hot(128)]);
            end
            imshow(icon,'Parent',AxesIcon)

            % 提示区
            PnlInfo = uipanel(obj.Fig,'Units','pixels','Position', [1+ht,barH,wt-ht,ht-barH]);
            strlist = {'当前次数',0,'总次数',obj.maxIter,'当前耗时','0'};
            obj.Txt = cell(6,1);
            % 设置控件的列宽和行高
            colWidth = 1/4; % 因为要分为四列，所以每列宽度为1/4
            rowHeight = 1/2; % 两行控件，所以每行高度为1/2
            for i = 1:2 % 两行
                for j = 1:4 % 第一行四个控件，第二行两个控件
                    if i == 2 && j > 2
                        break; % 第二行只需创建两个控件，所以这里退出循环
                    end
                    % 计算索引
                    n = (i - 1) * 4 + j; % 根据行和列计算索引
                    obj.Txt{n} = uicontrol(PnlInfo,'style','edit','Enable','off',...
                        'Units','normalized','Position',[(j-1)*colWidth, (i-1)*rowHeight , colWidth, rowHeight],...
                        'String',strlist{n},'Fontsize',16);
                end
            end

            % 进度区
            AxesBar = axes(obj.Fig,'Units','pixels','Position',[1+ht,1,wt-ht,barH]);
            axis(AxesBar,[-0.05,1.05,-0.2,0.2])
            axis(AxesBar,'off')
            hold(AxesBar,'on')
            obj.Line1 = plot(AxesBar,[0,1],[0,0],'-','LineWidth',15,'Color',[0.9,0.9,0.9]);
            obj.Line2 = plot(AxesBar,[0,0],[0,0],'-','LineWidth',15,'Color',[0.1,0.9,0.1]);
            hold(AxesBar,'off')
            drawnow

        end

        function updata(obj,iter)

            % 迭代次数
            set(obj.Txt{2},'string',iter)

            % 时间
            tnow = toc(obj.tstart);
            set(obj.Txt{6},'string',...
                [int2str(floor(tnow/60)),':',int2str(floor(rem(tnow,60)))])



            % 进度条
            obj.Line1.XData = [iter/obj.maxIter,1];
            obj.Line2.XData = [0,iter/obj.maxIter];
            drawnow

        end

        function closeWaitBar(obj)
            close(obj.Fig)
            clear("obj")
        end

    end




end