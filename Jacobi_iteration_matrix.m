%Jacobi迭代法的矩阵形式
function x = Jacobi_iteration_matrix(A, b, x0, tol, maxIter)
% A: 系数矩阵
% b: 右侧向量
% x0: 初始解
% tol: 容忍度
% maxIter: 最大迭代次数

L = -tril(A,-1); %取矩阵A的严格下三角矩阵的相反数
U = -triu(A,1); %取矩阵A的严格上三角矩阵的相反数
D = diag(diag(A)); %取矩阵A的对角矩阵
B = D\(L+U);
g = D\b;

x = x0;%初始化解向量
for k = 1:maxIter
    x_old = x;  % 保存上一次的解
    x = B*x+g;
    if norm(x - x_old, inf) < tol% 以无穷范数为标准，检查是否达到容忍度
        fprintf('矩阵形式Jacobi迭代法经过%d次迭代后找到解向量\n', k);
        break;
    end
end
