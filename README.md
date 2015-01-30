**Replication materials for:<br/>
Irregular Leadership Changes in 2014: Forecasts using ensemble, split-population duration models**
***

For questions contact the corresponding author [Michael Ward](mailto:michael.d.ward@duke.edu) or [Andreas Beger](andreas.beger@duke.edu).

This article is a summary of a longer technical report for the [PITF](http://en.wikipedia.org/wiki/Political_Instability_Task_Force). The complete [original report](http://arxiv.org/abs/1409.7105) is available on arXiv.org, and contains a large amount of additional information on the method we used for forecasting, accuracy assessments, etc.


**Citation:**

Beger, Andreas, Cassy L. Dorff, and Michael D. Ward, 2015, "Irregular Leadership Changes in 2014: Forecasts using ensemble, split-population duration models," *International Journal of Forecasting*. 

```bibtex
@article{beger2015irregular,
  title={Irregular Leadership Changes in 2014: Forecasts using ensemble, split-population duration models},
  author={Beger, Andreas and Dorff, Cassy L. and Ward, Michael D.},
  journal={International Journal of Forecasting},
  year={2015},
  volume={},
  issue={},
  pages={}
  }    
```

Getting the code and data
-----

The easiest way to get the replication code is to [download a zip](https://github.com/andybega/ijf-ilc2014/archive/master.zip). Alternatively, you can clone the repository through the Github GUI client ([OS X](https://mac.github.com/), [Windows](https://windows.github.com/)).

The data, including several intermediate results, are available on dataverse: [http://dx.doi.org/10.7910/DVN/28942](http://dx.doi.org/10.7910/DVN/28942).


Running the replication
-----

1. [Download](https://github.com/andybega/ijf-ilc2014/archive/master.zip) or [clone](github-mac://openRepo/https://github.com/andybega/ijf-ilc2014) this repository. 

2. Download the data sets on [Dataverse](http://dx.doi.org/10.7910/DVN/28942), at least the 2 beginning with `irc-data` and place them in `replication/data`.

3. In `runme.R`, change the working directory path on line 33.

4. Source or run the code in `runme.R`. We recommend running through the code block by block rather than sourcing. The original analysis was run on OS X using R 3.0.2 and 3.1.1.

The script relies on two packages, `EBMAforecastbeta` and `spduration` that are not available on CRAN. They are included in `replication/R/packages` with both OS X and Windows versions. The replication script will attempt to install them if they are not already present, but you may have to do so manually if this fails.

Files and scripts
------

`data`
* `all_preds.rda` - contains all theme/ensemble predictions from 2001 to 2014-09; used throughout `runme.r` to replicate figures in the same order as in the article, even though the models needed to create it are estimated in the same script   
* `ensemble_data.rda` - calibration/test data to estimate ensemble
* `irc_data_mod.rda` - imputed data
* `ensemble.rda` - saved ensemble model object
* `irc-data-v3.rda` - raw, unimputed source data
* `model_estimates.rda` - saved estimates for the 7 theme models

`graphics`
* Contains the graphics used in the article.

`R/packages`
* `EBMAforecastbeta_0.44.tar.gz` – OS X source package
* `EBMAforecastbeta_0.44.zip` – Windows source package
* `spduration_0.12.tar.gz` – OS X source package
* `spduration_0.12.zip` – Windows source package

`R/utilities`
* `ensemble_forecast.r` - helper functions to calculate ensemble forecast
* `gather_preds.r` - gathers all theme/ensemble predictions from 2001 to 2014-09 in one data frame, `all_preds.rda`
* `theme_models.r` - helper functions for theme model fit
* `varDecomp.r` - helpfer functions for variable variance decomposition
* `worldMap.r` - function for choropleth worldmap



