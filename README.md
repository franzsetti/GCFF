# Graph-Cuts for F-Formation (GCFF)

by Francesco Setti, Chris Russell, Chiara Bassetti and Marco Cristani


### Introduction

*An F-formation arises whenever two or more people sustain a spatial and orientational relationship in which the space between them is one to which they have equal, direct, and exclusive access.*  
[A. Kendon, 1990]

GCFF is a method to detect F-formations in static images, leveraging on the power of graph-cuts algorithms in clustering graphs.
In our formulation the nodes are represented by people in the scene and the candidate o-space centres, while edges are defined between each pair of nodes of different type (i.e. between a person and a candidate o-space centre).
We model the probability of each individual to belong to a specific o-space as: Pr(C,i) = O(Gi | μ ,i) ∝ exp (- ||xGi-μi||^2/σ^2) .
We then build the cost function by adding a Minimum Description Length (MDL) prior and considering the log function of the probability obtained.
Moreover, we introduce an additive term which acts as the visibility constraint on the individual i regardless of the group person j is assigned to.

You can find more details about the method on the original PloS One [paper](http://journals.plos.org/plosone/article?id=10.1371/journal.pone.0123783) or on the dedicated [webpage](http://vips.sci.univr.it/research/fformation/).


### License

GCFF is released under the MIT License (refer to the LICENSE file for details).


### Citing GCFF

If you find GCFF useful in your research, please consider citing:

    @article{setti2015f,
        title={F-formation detection: Individuating free-standing conversational groups in images},
        author={Setti, Francesco and Russell, Chris and Bassetti, Chiara and Cristani, Marco},
        journal={PloS One},
        volume={10},
        number={5},
        pages={e0123783},
        year={2015},
        publisher={Public Library of Science}
    }


### Data

In order to run the example provided and to reproduce the experiments of the original paper, you may want to download the datasets [here](http://vips.sci.univr.it/research/fformation/download/data.zip).
