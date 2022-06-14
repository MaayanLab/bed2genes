# bed2genes
Small python package to map bed peaks to genes. The input is a bed file loaded as a pandas dataframe. It will filter bed entries to only show peaks that have a gene start site close by as defined by `upstream_distance` and `downstream_distance`. Currently the package supports `human` and `mouse`.


# Installation

bed2genes is currently only available as a Python package in this GitHub repository. You can install the bed2genes Python package and its dependencies through pip by using the following command:

```
$ pip install git+https://github.com/MaayanLab/bed2genes.git
```

### Python example

```python
import bed2genes
import pandas as pd
import urllib.request

# download example gene expression signature
url = "https://raw.githubusercontent.com/MaayanLab/bed2genes/main/bed2genes/data/example.bed"
urllib.request.urlretrieve(url, "example.bed")

bed = pd.read_csv("example.bed", sep="\t")

bed_genes = bed2genes.map_genes(bed, species="human", upstream_distance=-2000, downstream_distance=500, reload=False)
```

### Optional Parameters

The main function of `bed2genes.map_genes()` supports several optional parameters. The default parameters should work well for most use cases.

| parameter name | type | default | description |
|:-----|:---------|:-------------|:------|
| `species`	| string | human | Reference genome. (human/mouse)|
| `upstream_distance` | int | -2000 | max bp upstream of gene start. |
| `downstream_distance` | int | 500 | max bp downstream (in gene) of gene start. |
| `reload` | bool | False | Reload gene information or use cache. |
