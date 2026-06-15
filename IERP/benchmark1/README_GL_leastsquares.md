# GL_leastsquares

调用接口：

```matlab
[output_Atilde, output_Omega] = GL_leastsquares(D);
```

`D` 被解释为论文中的目标 effective resistance matrix：

- `D` 必须为实数、有限、非负、对称方阵；
- 对角线必须为 `0`；
- 与作者原始实现一致，非对角线上的 `0` 表示该节点对没有约束；
- 输出 `output_Atilde` 是无向加权图的邻接矩阵（权重是 conductance）；
- 输出 `output_Omega` 是由 `output_Atilde` 重新计算的 effective resistance matrix。

实现对应论文 Problem 2：

```text
min_H sum_(i,j in S) (r_H(i,j) - D(i,j))^2
```

流程：

1. 使用最短路补全缺失的 resistance constraints；
2. 使用论文 Theorem 1 的闭式公式初始化，并投影掉负边权；
3. 如果初始化已达到零误差，则它已经是 least-squares 全局最优解；
4. 否则调用作者原始 `effResGDSmall.m` 或 `effResGD.m`，执行投影梯度下降/随机坐标下降和线搜索；
5. 从最终邻接矩阵计算 `output_Omega`。

作者原始代码完整保存在 `upstream_graph_similarity_learning`，适配器没有修改其中任何文件。原代码运行时生成的 `iter*.mat` 被隔离在临时目录并自动清理。

作者的大图函数调用 Statistics Toolbox 的 `randsample`。若当前 MATLAB
没有该函数，适配器会临时使用 `compat/randsample.m` 中等价的无放回采样实现；
该兼容文件不在作者原始代码目录内。

运行测试：

```matlab
cd('C:\Users\86748\Documents\MATLAB\flowbetweenness\IERP\benchmark1')
test_GL_leastsquares
```

论文：

J. G. Hoskins, C. Musco, C. Musco, and C. E. Tsourakakis,
"Inferring Networks From Random Walk-Based Node Similarities,"
NeurIPS 2018.

作者代码：

https://github.com/cnmusco/graph-similarity-learning
