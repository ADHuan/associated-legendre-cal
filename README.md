# 缔合勒让德函数计算（Associated Legendre Functions Calculation）

## 简介

编写了计算缔合勒让德函数的程序，供参考，共有四种方法：

- 标准前向列推法[<sup>1</sup>](#refer-anchor-1)
- 标准前向行推法[<sup>1</sup>](#refer-anchor-1)
- Belikov递推法[<sup>2</sup>](#refer-anchor-2)
- 跨阶次递推法（Swarztrauber 法）[<sup>2</sup>](#refer-anchor-2)

其中标准前向行、列推法使用MATLAB编写，其余使用C++编写

## 注意
高阶项由于超过了double类型的最小值范围$(P_{nm}<-1.7E-308)$，计算不精准，需注意

## 参考

<div id="refer-anchor-1"></div>

- [1] [吴星,刘雁雨. 多种超高阶次缔合勒让德函数计算方法的比较[J]. 测绘科学技术学报,2006,23(3):188-191. DOI:10.3969/j.issn.1673-6338.2006.03.010.](https://d.wanfangdata.com.cn/periodical/chxyxb200603010)

<div id="refer-anchor-2"></div>

- [2] [欧阳明达,张敏利,于亮. 采用Belikov列推和跨阶次递推方法计算超高阶缔合勒让德函数[J]. 测绘工程,2017,26(7):12-15,21. DOI:10.19349/j.cnki.issn1006-7949.2017.07.003.](https://d.wanfangdata.com.cn/periodical/chgc201707003)

