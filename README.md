CoolBox
=======

[![Binder](https://mybinder.org/badge_logo.svg)](https://mybinder.org/v2/gh/GangCaoLab/CoolBox/master?filepath=demo%2FTestRegion.ipynb)
[![Python Package using Conda](https://github.com/GangCaoLab/CoolBox/workflows/Python%20Package%20using%20Conda/badge.svg)](https://github.com/GangCaoLab/CoolBox/actions?query=workflow%3A%22Python+Package+using+Conda%22)

WIP

Flexible, user-friendly genomic data visualization toolkit. 

![](imgs/title.png)

Highlights:

* Multi-omics data interactively visualization
* Show within Jupyter notebook
* User-friendly API (ggplot2-like Python EDSL)

Documents:
* [Wiki](https://github.com/Nanguage/CoolBox/wiki)
* Jupyter notebook [walkthrough](demo/coolbox_guide.ipynb)

TODO List:

+ Using in CLI and show in Qt Window.
    + fire CLI to compose tracks.
    + Panel in Qt.
        - ref [qtvoila](https://luiztauffer.github.io/guacamole-data-science/posts/2020-04-20-qtvoila/)
+ Install with conda
    + Upload to [Bioconda](https://bioconda.github.io/)
+ Better documents.
+ Interactively show in independent web-page.
    + ref this [example](https://github.com/jupyter-widgets/ipywidgets/tree/master/examples/web3)
    + [voila](https://github.com/voila-dashboards/voila)
+ Speed up loading.
    + speed up Arcs/BEDGraph loading with make index.
+ Optimizing plot
    + other ploting system(JS based)
+ Support more data type
    + Image type
    + BEDPE/Pairs (For Arcs plot)


Thanks

+ [pyGenomeTracks](https://github.com/deeptools/pyGenomeTracks),
CoolBox's plot system is fork from it.


