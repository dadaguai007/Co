classdef DataSaver < handle
    properties
        % 存储数据的变量
        data;
        % 存储文件路径
        filePath;
        % 存储文件性能
        name;
    end
    
    methods
        % 构造函数，用于初始化数据和文件路径
        function obj = DataSaver(data, filePath,name)
            obj.data         = data;
            obj.filePath     = filePath;
            obj.name         = name;
        end

        % 创建文件夹的方法
        function createFolder(obj)
            if ~exist(obj.filePath, 'dir')
                mkdir(obj.filePath);
                fprintf('Folder %s created successfully.\n', obj.filePath);
            else
                fprintf('Folder %s already exists.\n', obj.filePath);
            end
        end


        % 保存数据到.mat 文件的方法
        function saveToMat(obj)
            % 组成生成路径
            txt=sprintf('%s\\%s.mat',obj.filePath,obj.name);
            try
                % 检查文件路径是否以.mat 结尾
                if ~endsWith(txt, '.mat')
                    txt = [txt, '.mat'];
                end
                save_data=obj.data;
                % 保存数据到.mat 文件
                save(txt, 'save_data');
                fprintf('Data saved to %s successfully.\n', obj.filePath);
            catch ME
                fprintf('Error saving data to .mat file: %s\n', ME.message);
            end
        end
        
    end
end