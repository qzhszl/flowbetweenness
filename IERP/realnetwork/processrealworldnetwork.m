clear; clc;

% === 修改为你的 GML 文件路径 ===
gmlFile = 'D:\\data\\flow betweenness\\IERP\\realnetwork\\NewSpain_18c_travelmap.gml';

txt = fileread(gmlFile);
lines = regexp(txt, '\r\n|\r|\n', 'split');

edges_src = [];
edges_tgt = [];
edges_w   = [];

i = 1;
nline = numel(lines);

% --- 解析 edges（source,target,weight/value） ---
while i <= nline
    s = strtrim(lines{i});
    if startsWith(s, 'edge')
        i = i + 1;
        src = NaN;
        tgt = NaN;
        w   = NaN;

        while i <= nline
            t = strtrim(lines{i});

            % source
            tok = regexp(t, 'source\s+(-?\d+)', 'tokens','once');
            if ~isempty(tok)
                src = str2double(tok{1});
            end

            % target
            tok = regexp(t, 'target\s+(-?\d+)', 'tokens','once');
            if ~isempty(tok)
                tgt = str2double(tok{1});
            end

            % weight / value
            tok = regexp(t, '(weight|value)\s+([+-]?\d+(\.\d+)?([eE][+-]?\d+)?)', 'tokens','once');
            if ~isempty(tok)
                w = str2double(tok{2});
            end

            if startsWith(t, ']')
                break;
            end

            i = i + 1;
        end

        % 默认权重为 1
        if isnan(w), w = 1; end

        edges_src(end+1,1) = src;
        edges_tgt(end+1,1) = tgt;
        edges_w(end+1,1)   = w;
    end
    i = i + 1;
end

% --- 无自环 ---
valid = edges_src ~= edges_tgt;
edges_src = edges_src(valid);
edges_tgt = edges_tgt(valid);
edges_w   = edges_w(valid);

% --- 无向 + 合并重边 (u < v) ---
pairMap = containers.Map();

for k = 1:numel(edges_src)
    u = edges_src(k);
    v = edges_tgt(k);
    w = edges_w(k);

    % 排序无向边
    a = min(u,v);
    b = max(u,v);
    key = sprintf('%d_%d', a, b);

    % 合并重边（权重相加）
    if isKey(pairMap, key)
        pairMap(key) = pairMap(key) + w;
    else
        pairMap(key) = w;
    end
end

% --- 输出 txt 文件: node1  node2  weight ---
keysList = keys(pairMap);
fid = fopen('edges_final.txt','w');

for k = 1:numel(keysList)
    key = keysList{k};
    nums = sscanf(key, '%d_%d');
    u = nums(1);
    v = nums(2);
    w = pairMap(key);

    fprintf(fid, '%d\t%d\t%g\n', u, v, w);
end

fclose(fid);

fprintf('已生成 edges_final.txt（三列：node1 node2 weight）\n');
