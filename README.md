Replication materials for IJF ILC 2014 paper
================

**[Irregular Leadership Changes in 2014: Forecasts using ensemble, split-population duration models](http://www.sciencedirect.com/science/article/pii/S0169207015000485)**

For questions contact the corresponding author [Michael Ward](mailto:michael.don.ward@gmail.com) or [Andreas Beger](mailto:adbeger@gmail.com).

This article is a summary of a longer technical report for the [PITF](http://en.wikipedia.org/wiki/Political_Instability_Task_Force). The complete [original report](http://arxiv.org/abs/1409.7105) is available on arXiv.org, and contains a large amount of additional information on the method we used for forecasting, accuracy assessments, etc.

**Citation:**

``` bibtex
@article{beger2016irregular,
  title={Irregular Leadership Changes in 2014: Forecasts using ensemble, split-population duration models},
  author={Beger, Andreas and Dorff, Cassy L. and Ward, Michael D.},
  journal={International Journal of Forecasting},
  year={2016},
  volume={32},
  issue={1},
  pages={98--111}
  }    
```

Getting the code and data
-------------------------

The easiest way to get the replication code is to [download a zip](https://github.com/andybega/ijf-ilc2014/archive/master.zip). Alternatively, you can clone the repository through the Github GUI client ([OS X](https://mac.github.com/), [Windows](https://windows.github.com/)).

The data are available on dataverse: <http://dx.doi.org/10.7910/DVN/28942>. Several smaller intermediate results are included in the git data folder, but replicating the full analysis will require the larger raw data from dataverse.

Running the replication
-----------------------

1.  [Download](https://github.com/andybega/ijf-ilc2014/archive/master.zip) or [clone](github-mac://openRepo/https://github.com/andybega/ijf-ilc2014) this repository.

2.  Download the data sets on [Dataverse](http://dx.doi.org/10.7910/DVN/28942), at least the 2 beginning with `irc-data` and place them in `replication/data`.

3.  In `runme.R`, change the working directory path on line 33.

4.  Source or run the code in `runme.R`. We recommend running through the code block by block rather than sourcing. The original analysis was run on OS X using R 3.0.2 and 3.1.1.

The script relies on two packages, `EBMAforecastbeta` and `spduration` that are not available on CRAN. They are included in `replication/R/packages` with both OS X and Windows versions. The replication script will attempt to install them if they are not already present, but you may have to do so manually if this fails.

Files and scripts
-----------------

`data`:

-   `all_preds.rda` - contains all theme/ensemble predictions from 2001 to 2014-09; used throughout `runme.r` to replicate figures in the same order as in the article, even though the models needed to create it are estimated in the same script
-   `ensemble_data.rda` - calibration/test data to estimate ensemble
-   `irc_data_mod.rda` - imputed data
-   `ensemble.rda` - saved ensemble model object
-   `irc-data-v3.rda` - raw, unimputed source data
-   `model_estimates.rda` - saved estimates for the 7 theme models

`graphics`:

-   Contains the graphics used in the article.

`R/packages`:

-   `EBMAforecastbeta_0.44.tar.gz` – OS X source package
-   `EBMAforecastbeta_0.44.zip` – Windows source package
-   `spduration_0.12.tar.gz` – OS X source package
-   `spduration_0.12.zip` – Windows source package

`R/utilities`:

-   `ensemble_forecast.r` - helper functions to calculate ensemble forecast
-   `gather_preds.r` - gathers all theme/ensemble predictions from 2001 to 2014-09 in one data frame, `all_preds.rda`
-   `theme_models.r` - helper functions for theme model fit
-   `varDecomp.r` - helpfer functions for variable variance decomposition
-   `worldMap.r` - function for choropleth worldmap

2019-04-05 Update
-----------------

Checked replication and updated several issues. See `runme.R` for more details in the notes at the top.

To replicate the exact results, use the saved fitted models and predictions.

``` r
sessionInfo()
```

    ## R version 3.5.2 (2018-12-20)
    ## Platform: x86_64-apple-darwin15.6.0 (64-bit)
    ## Running under: macOS Mojave 10.14.4
    ## 
    ## Matrix products: default
    ## BLAS: /Library/Frameworks/R.framework/Versions/3.5/Resources/lib/libRblas.0.dylib
    ## LAPACK: /Library/Frameworks/R.framework/Versions/3.5/Resources/lib/libRlapack.dylib
    ## 
    ## locale:
    ## [1] en_US.UTF-8/en_US.UTF-8/en_US.UTF-8/C/en_US.UTF-8/en_US.UTF-8
    ## 
    ## attached base packages:
    ## [1] stats     graphics  grDevices utils     datasets  methods   base     
    ## 
    ## loaded via a namespace (and not attached):
    ##  [1] compiler_3.5.2  magrittr_1.5    tools_3.5.2     htmltools_0.3.6
    ##  [5] yaml_2.2.0      Rcpp_1.0.1      stringi_1.4.3   rmarkdown_1.11 
    ##  [9] knitr_1.22      stringr_1.4.0   xfun_0.5        digest_0.6.18  
    ## [13] evaluate_0.13
