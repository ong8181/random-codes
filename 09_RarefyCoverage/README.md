# `rarefy_even_coverage` 関数
`rarefy_even_coverage` は `phyloseq` オブジェクトをインプットとして、coverage-based rarefaction を行うための関数です。内部で `iNEXT` package と `phyloseq:::rarefaction_subsample` を利用して coverage-based rarefaction を実行します. `rarefy_even_coverage` 関数は `metagMisc` package の `phyloseq_coverage_raref` 関数と同様の結果を返します (https://github.com/vmikk/metagMisc). ただし、`phyloseq` オブジェクト内の分類群が行でも列でも (`taxa_are_rows = TRUE` or `taxa_are_rows = FALSE`) 使用可能です. さらに、`rarefy_even_coverage` 関数は rarefaction カーブの可視化のためのオブジェクトも返します.

詳しくは `demo_rarefy.R` を御覧ください.

## 関数の説明
``` r
rarefy_even_coverage <-  function(ps_obj,
                                  coverage = 0.97,
                                  remove_not_rarefied = FALSE,
                                  include_iNEXT_results = FALSE,
                                  nboot = 40,
                                  knots = 50,
                                  n_rarefy_iter = 1,
                                  rarefy_average_method = "round",
                                  sample_method = "phyloseq",
                                  ran_seed = 1234)
```

### 重要パラメータ
- `ps_obj`: `phyloseq` オブジェクト.
- `coverage`:  カバレッジ (デフォルト = 0.97 [= 97%])
- `include_iNEXT_results`: rarefaction カーブを図示するための `iNEXT` の結果を出力するかどうか. `TRUE` とすると返り値はリストとなり、１つ目の要素が `phyloseq` オブジェクト、２つ目の要素が `iNEXT` の結果となります. また、`TRUE` とすると計算時間が増大します. `FALSE` の場合は `phyloseq` オブジェクトのみが出力されます.

### その他のパラメータ
- `remove_not_rarefied`: 指定されたカバレッジに届かないサンプルを出力から除くかどうか.
- `nboot`: `iNEXT` 関数の `nboot` を指定 (`include_iNEXT_results = TRUE` の場合のみ有効).
- `knots`: `iNEXT` 関数の `knots` を指定 (`include_iNEXT_results = TRUE` の場合のみ有効). rarefaction カーブを描くための点 (knot) の数. `knots` の数が少ないと rarefaction カーブがカクカクします.
- `n_rarefy_iter`: OTU テーブルを何回 rarefy するか. 数を増やすとランダムサンプリングに起因する誤差を軽減できます (デフォルト = 1).
- `rarefy_average_method`: `n_rarefy_iter >= 2` のとき、複数回の rarefaction の結果をどのように要約するか. `rarefy_average_method = "round"` のときは `round()` を使って平均値を丸めます. `rarefy_average_method = "floor"` のときは `floor()` を使って丸めるので、複数回の rarefaction で平均して 1 回以上その OTU が選択されない場合はその OTU は 0 となります. `rarefy_average_method = "ceiling"` のときは `ceiling()` を使って丸めるので、複数回の rarefaction のうち、1 回でもその OTU が選択された場合はその OTU は 1 となります. 
- `sample_method`: rarefaction の際に `vegan::rrarefy()` を使うか (`sample_method = "vegan"`) 、`phyloseq:::rarefaction_subsample()` を使うか (`sample_method = "phyloseq"`).
- `ran_seed`: ランダムシード値.


## 実行例
```{r}
# OTU テーブルの rarefaction と可視化のためのオブジェクトを同時に返す
## include_iNEXT_results = TRUE とすると可視化のための計算を実行します
## OTU テーブルの rarefaction のみが必要であればここを FALSE とします
ps_rare_raw <- rarefy_even_coverage(ps_sample,
                                    coverage = 0.97,
                                    include_iNEXT_results = TRUE)
# phyloseq オブジェクトのみを抽出
ps_rare <- ps_rare_raw[[1]]                      
# 図示
plot_rarefy(ps_rare_raw)
```

<img src="img/rarefy_plot.png" width="800px">


## 実行例




# `rarefy_even_coverage` function
This repository includes convenient functions to perform coverage-based rarefaction. `phyloseq` object can be easily rarefied based on a user-specified coverage by the following command. Functions implemented in `iNEXT` package is used. `rarefy_even_coverage` returns almost identical results with `phyloseq_coverage_raref` function in `metagMisc` (https://github.com/vmikk/metagMisc), but a phyloseq object with `taxa_are_rows = TRUE` or `taxa_are_rows = FALSE` is accepted. In addition, `rarefy_even_coverage` returns a rarefaction curve for visualization.

For detail, please run `demo_rarefy.R`.

## Description
``` r
rarefy_even_coverage <-  function(ps_obj,
                                  coverage = 0.97,
                                  remove_not_rarefied = FALSE,
                                  include_iNEXT_results = FALSE,
                                  nboot = 40,
                                  knots = 50,
                                  n_rarefy_iter = 1,
                                  rarefy_average_method = "round",
                                  sample_method = "phyloseq",
                                  ran_seed = 1234)
```

### Important parameters
- `ps_obj`: `phyloseq` object.
- `coverage`:  User-specified parameter (default = 0.97 [= 97%])
- `include_iNEXT_results`: Include `iNEXT` results or not. If `TRUE`, the function returns a list which contains two elements. The first object is a rarefied phyloseq object, and the second object is an iNEXT result. Also, if `TRUE`, computation time will increase. If`FALSE`, it returns a rarefied`phyloseq` object only.

### Other parameters
- `remove_not_rarefied`: Remove samples of which coverage is lower than `coverage`.
- `nboot`: Specifi `nboot` of `iNEXT` function (valid if `include_iNEXT_results = TRUE`).
- `knots`: Specifi `knots` of `iNEXT` function (valid if `include_iNEXT_results = TRUE`).
- `n_rarefy_iter`: The number of iterations of rarefactions (default = 1).
- `rarefy_average_method`: If `n_rarefy_iter >= 2`, this argument determines how the multiple rarefactions are summarized. `rarefy_average_method = "round"` uses `round()`. `rarefy_average_method = "floor"` uses `floor()`. `rarefy_average_method = "ceiling"` uses `ceiling()`.
- `sample_method`: Specify which function is used for rarefaction. `sample_method = "vegan"` uses `vegan::rrarefy()`, while `sample_method = "phyloseq"` uses `phyloseq:::rarefaction_subsample()`.
- `ran_seed`: Random seed.


## Example
```{r}
# Calculate rarefied matrix and rarefaction curve simultaneously
ps_rare_raw <- rarefy_even_coverage(ps_sample,
                                    coverage = 0.97,
                                    include_iNEXT_results = TRUE)
# Extract phyloseq object
ps_rare <- ps_rare_raw[[1]]                      
```

In addition, results of the coverage-based rarefaction can be checked by visualizing the rarefaction curves.

```{r}
plot_rarefy(ps_rare_raw)
```

<img src="img/rarefy_plot.png" width="800px">
